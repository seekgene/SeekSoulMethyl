import sys
from multiprocessing import Pipe, Process, Queue
import multiprocessing.connection
import io
import dnaio
from xopen import xopen

class Reader(Process):
    def __init__(self, file1:str, file2:str, connections:multiprocessing.connection, queue:Queue, buffer_size:int):
        """

        Args:
            file1: read1 fastq file
            file2: read2 fastq file
            connections: 
            queue: 
            buffer_size: buffer
        
        Returns:
            Reader

        """
        super().__init__()
        self.file1 = file1
        self.file2 = file2
        self.connections = connections
        self.queue = queue
        self.buffer_size = buffer_size

    def run(self):
        try:
            chunk_index = 0
            for file1, file2 in zip(self.file1, self.file2):
                with xopen(file1, 'rb') as f1:
                    with xopen(file2, 'rb') as f2:
                        for (chunk1, chunk2) in dnaio.read_paired_chunks(f1, f2, self.buffer_size):
                            worker_index = self.queue.get()
                            pipe = self.connections[worker_index]
                            pipe.send(chunk_index)
                            pipe.send_bytes(chunk1)
                            pipe.send_bytes(chunk2)
                            chunk_index += 1
            for _ in range(len(self.connections)):
                worker_index = self.queue.get()
                self.connections[worker_index].send(-1)
        except Exception as e:
            for worker_index in range(len(self.connections)):
                self.connections[worker_index].send(-2)
            raise e

class Writer:
    def __init__(self, file:str, file_multi:str, paired_out:bool=False):
        """
        Args:
            file: read2 fastq file
            file_muti: unknow barcode
        """
        self.paired_out = paired_out

        if paired_out:
            self.file1 = f"{file}_1.fq.gz"
            self.file2 = f"{file}_2.fq.gz"
            self._fh1 = xopen(self.file1, mode='wb')
            self._fh2 = xopen(self.file2, mode='wb')

            self.file1_multi = f"{file_multi}_1.fq.gz"
            self.file2_multi = f"{file_multi}_2.fq.gz"
            self._fh_multi1 = xopen(self.file1_multi, mode='wb')
            self._fh_multi2 = xopen(self.file2_multi, mode='wb')
        else:
            self.file = f"{file}.fq.gz"
            self._fh1 = xopen(self.file, mode='wb')

            self.file_multi = f"{file_multi}.fq.gz"
            self._fh_multi1 = xopen(self.file_multi, mode='wb')

        self._chunks = dict()
        self._current_index = 0

    def write(self, data, index):
        self._chunks[index] = data
        if self.paired_out:
            while self._current_index in self._chunks:
                (_r1, _r2), (_multi_r1, _multi_r2) = self._chunks[self._current_index]
                self._fh1.write(_r1)
                self._fh2.write(_r2)
                self._fh_multi1.write(_multi_r1)
                self._fh_multi2.write(_multi_r2)
                del self._chunks[self._current_index]
                self._current_index += 1
        else:
            while self._current_index in self._chunks:
                (_r1,), (_multi_r1,) = self._chunks[self._current_index]
                self._fh1.write(_r1)
                self._fh_multi1.write(_multi_r1)
                del self._chunks[self._current_index]
                self._current_index += 1            

    def wrote_everything(self):
        return not self._chunks

    def close(self):
        self._fh1.close()
        self._fh_multi1.close()
        if self.paired_out:
            self._fh2.close()
            self._fh_multi2.close()

class Worker(Process):
    """
    """
    def __init__(self, id_, read_pipe, write_pipe, need_work_queue, func, paired_out=False):
        """
        """
        super().__init__()
        self._id = id_
        self.read_pipe = read_pipe
        self.write_pipe = write_pipe
        self.need_work_queue = need_work_queue
        self.func = func
        self.paired_out = paired_out

    def run(self):
        try:
            while True:
                self.need_work_queue.put(self._id)
                chunk_index = self.read_pipe.recv()
                if chunk_index == -1:
                    break
                elif chunk_index == -2:
                    e, tb_str = self.read_pipe.recv()
                    raise e
                try:
                    data = self.read_pipe.recv_bytes()
                    input = io.BytesIO(data)
                    data = self.read_pipe.recv_bytes()
                    input2 = io.BytesIO(data)
                    if self.paired_out:
                        tmp = (io.BytesIO(), io.BytesIO())
                        tmp_multi = (io.BytesIO(), io.BytesIO())
                        _ = self.func(fq1=input, fq2=input2, fq_out=tmp, fqout_multi=tmp_multi)
                        self.write_pipe.send(chunk_index)
                        self.write_pipe.send_bytes(tmp[0].getvalue())
                        self.write_pipe.send_bytes(tmp[1].getvalue())
                        self.write_pipe.send_bytes(tmp_multi[0].getvalue())
                        self.write_pipe.send_bytes(tmp_multi[1].getvalue())
                        self.write_pipe.send(_)
                    else:
                        tmp = (io.BytesIO(),)
                        tmp_multi = (io.BytesIO(),)
                        _ = self.func(fq1=input, fq2=input2, fq_out=tmp, fqout_multi=tmp_multi)
                        self.write_pipe.send(chunk_index)
                        self.write_pipe.send_bytes(tmp[0].getvalue())
                        self.write_pipe.send_bytes(tmp_multi[0].getvalue())
                        self.write_pipe.send(_)
                except Exception as e:
                    self.write_pipe.send(-2)
                    self.write_pipe.send((e, traceback.format_exc()))
                    raise e  
            self.write_pipe.send(-1)
        except Exception as e:
            self.write_pipe.send(-2)
            self.write_pipe.send((e, traceback.format_exc()))
            raise e
        ##except Exception as e:
        ##    self.write_pipe.send(-2)
        ##    raise e
            

class Pipeline:
    def __init__(self, func, fq1, fq2, fqout, fqout_multi, core, stat=None, paired_out=False, buffer_size=16*1024**2):
        self.n_workers = core
        self.fq1 = fq1
        self.fq2 = fq2
        self.fqout = fqout
        self.fqout_multi = fqout_multi
        self.buffer_size = buffer_size
        self.need_work_queue = Queue()
        self.func = func
        self.paired_out = paired_out
        self.stat = stat

    def run(self):
        try:
            # start reader process
            reader_connections = [Pipe(duplex=False) for _ in range(self.n_workers)]
            _pipes, _conn = zip(*reader_connections)
            _reader_process = Reader(self.fq1, self.fq2, _conn, self.need_work_queue, self.buffer_size)
            _reader_process.daemon = True
            _reader_process.start()

            self.workers = []
            self.connections = []
            self.writer = Writer(self.fqout, self.fqout_multi, self.paired_out)
            for index in range(self.n_workers):
                conn_r, conn_w = Pipe(duplex=False)
                self.connections.append(conn_r)
                worker = Worker(index, _pipes[index], conn_w, self.need_work_queue,
                                self.func, self.paired_out)
                worker.daemon = True
                worker.start()
                self.workers.append(worker)

            # write output
            while self.connections:
                ready_connections = multiprocessing.connection.wait(self.connections)
                for connection in ready_connections:
                    try:
                        chunk_index = connection.recv()
                        if chunk_index == -1:
                            self.connections.remove(connection)
                            continue
                        elif chunk_index == -2:
                            ##sys.stderr.write('err!!!\n')
                            error, traceback_str = connection.recv()
                            print(f"error :\n{traceback_str}", file=sys.stderr)
                            self._cleanup(_reader_process)
                            raise error 
                        # if single?
                        data1 = connection.recv_bytes()
                        data2 = connection.recv_bytes()
                        data_multi1 = connection.recv_bytes()
                        data_multi2 = connection.recv_bytes()
                        self.writer.write([(data1, data2), (data_multi1, data_multi2)], chunk_index)
                        _stat = connection.recv()
                        self.stat.update(**_stat)
                    except EOFError:
                        self.connections.remove(connection)
            assert self.writer.wrote_everything()
            for w in self.workers:
                w.join()
            _reader_process.join()
            self.writer.close()

        except Exception as e:
            for w in self.workers:
                try:
                    w.terminate()
                except:
                    pass
            try:
                _reader_process.terminate()
            except:
                pass
            try:
                self.writer.close()
            except:
                pass
            raise  
