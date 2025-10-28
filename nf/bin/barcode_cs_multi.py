#!/usr/bin/env python

import os
import sys
import json
import dnaio
import itertools
from functools import partial
from multiprocessing import Pipe, Process, Queue
import multiprocessing.connection
import io
from xopen import xopen
from collections import defaultdict, Counter
import re
import gzip
from loguru import logger
from cutadapt.adapters import FrontAdapter, BackAdapter, Sequence
import click
import random
import string

# Pre-compiled regular expressions for CH mode filtering
# Forward chain direction patterns
R1_FORWARD_CH_PATTERN = re.compile(r'C[ATC]')  # R1: C[ATC]
R2_FORWARD_CH_PATTERN = re.compile(r'[ATG]G')  # R2: [ATG]G

# Reverse chain direction patterns  
R1_REVERSE_CH_PATTERN = re.compile(r'[ATG]G')  # R1: [ATG]G
R2_REVERSE_CH_PATTERN = re.compile(r'C[ATC]')  # R2: C[ATC]


R1_MINLEN = 20
R2_MINLEN = 60
ATAC_R1_MINLEN = 20
ATAC_R2_MINLEN = 60


# Key position arrays for CT conversion rate calculation
# Default positions (non-ME5 mode)
DEFAULT_POSITIONS = [24, 27, 28, 31, 36, 38, 57]
# ME5 mode positions
ME5_POSITIONS = [42, 45, 46, 49, 54, 56, 75]

# Non-insert sequence length configuration
# Default mode: 52
DEFAULT_NON_INSERT_LEN = 52
# ME5 mode: 70
ME5_NON_INSERT_LEN = 70

# Sequence extraction parameter configuration
# Default mode parameters
DEFAULT_SEQ_7F_START = 17
DEFAULT_SEQ_7F_LEN = 7
DEFAULT_L17ME_START = 39
DEFAULT_L17ME_LEN = 18

# ME5 mode parameters
ME5_SEQ_7F_START = 0
ME5_SEQ_7F_LEN = 0
ME5_L17ME_START = 57
ME5_L17ME_LEN = 18

# Maximum sequence length configuration
DEFAULT_MAX_READ1_LEN = 81
ME5_MAX_READ1_LEN = 63


CHEMISTRY = {
    "DD-M":{
        'shift': False,
        'structure': 'B17',
        'adapter1': [["AGATGTGTATAAGAGAYAG", "5", 0.1, 9],
                    ["CTGTCTCTTATACACATCT","3",0.1, 9]], # ME,ME-rev Y:C/T
        'adapter2': [["AGATGTGTATAAGAGACAG", "5", 0.1, 9],
                     ["CTRTCTCTTATACACATCT","3",0.1, 9]], # ME,ME-rev R:G/A
        'match_type': (1,),
    },
    "ME5":{
        'shift': False,
        'structure': 'B17U12',
        'adapter1': [["AGATGTGTATAAGAGAYAG", "5", 0.1, 9],
                    ["CTGTCTCTTATACACATCT","3",0.1, 9]], # ME,ME-rev Y:C/T
        'adapter2': [["AGATGTGTATAAGAGACAG", "5", 0.1, 9],
                     ["CTRTCTCTTATACACATCT","3",0.1, 9]], # ME,ME-rev R:G/A
        'match_type': (1,),
    }    
}

class Reader(Process):
    """
    Read paired fastq files
    """
    def __init__(self, file1:str, file2:str, connections:multiprocessing.connection, queue:Queue, buffer_size:int):
        """Initialize reader process

        Args:
            file1: read1 fastq file
            file2: read2 fastq file
            connections: 
            queue: Queue storing idle worker process indices
            buffer_size: Buffer size for reading fastq files
        
        Returns:
            Reader object

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
                with xopen(file1, 'rb', threads=4) as f1:
                    with xopen(file2, 'rb', threads=4) as f2:
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
    """
    Handle output content
    """
    def __init__(self, file:str, file_multi:str, paired_out:bool=False, split_fastq:int=0, samplename:str=""):
        """Output processed sequences
        Args:
            file: read2 fastq file
            file_muti: Fastq filename for undetermined barcodes
            paired_out: Whether to output paired reads
            split_fastq: Split files by first n bases of barcode (0=no split)
            samplename: Sample name for file naming
        """
        self.paired_out = paired_out
        self.split_fastq = split_fastq
        self.samplename = samplename
        self.file_handles = {}  # Store file handles for different barcode prefixes

        if paired_out:
            # Create output files for both chain directions
            self.file1_forward = f"{file}_forward_1.fq.gz"
            self.file2_forward = f"{file}_forward_2.fq.gz"
            self._fh1_forward = xopen(self.file1_forward, mode='wb', threads=4, compresslevel=1)
            self._fh2_forward = xopen(self.file2_forward, mode='wb', threads=4, compresslevel=1)
            
            self.file1_reverse = f"{file}_reverse_1.fq.gz"
            self.file2_reverse = f"{file}_reverse_2.fq.gz"
            self._fh1_reverse = xopen(self.file1_reverse, mode='wb', threads=4, compresslevel=1)
            self._fh2_reverse = xopen(self.file2_reverse, mode='wb', threads=4, compresslevel=1)

            self.file1_multi = f"{file_multi}_1.fq.gz"
            self.file2_multi = f"{file_multi}_2.fq.gz"
            self._fh_multi1 = xopen(self.file1_multi, mode='wb', threads=4, compresslevel=1)
            self._fh_multi2 = xopen(self.file2_multi, mode='wb', threads=4, compresslevel=1)
        else:
            # Single-end output
            self.file_forward = f"{file}_forward.fq.gz"
            self._fh_forward = xopen(self.file_forward, mode='wb', threads=4, compresslevel=1)
            
            self.file_reverse = f"{file}_reverse.fq.gz"
            self._fh_reverse = xopen(self.file_reverse, mode='wb', threads=4, compresslevel=1)

            self.file_multi = f"{file_multi}.fq.gz"
            self._fh_multi1 = xopen(self.file_multi, mode='wb', threads=4, compresslevel=1)

        self._chunks = dict()
        self._current_index = 0

    def _get_barcode_prefix(self, read_name):
        """Extract barcode prefix from read name"""
        # read name format: barcode_umi_chain_direction_alt_original_name
        parts = read_name.split('_')
        if len(parts) > 0:
            full_barcode = parts[0]  # Complete barcode
            return full_barcode[:self.split_fastq] if self.split_fastq > 0 else ""
        return ""

    def _extract_barcode_prefix_from_bytes(self, fastq_bytes):
        """Extract barcode prefix from fastq byte data"""
        try:
            # Convert byte data to string
            fastq_str = fastq_bytes.decode('utf-8')
            lines = fastq_str.strip().split('\n')
            if len(lines) >= 1:
                # First line is read name
                read_name = lines[0]
                if read_name.startswith('@'):
                    read_name = read_name[1:]  # Remove @ symbol
                return self._get_barcode_prefix(read_name)
        except Exception as e:
            # If parsing fails, return empty string
            pass
        return ""

    def _get_file_handle(self, prefix, direction):
        """Get or create file handle for corresponding prefix and direction"""
        if self.split_fastq == 0:
            # No splitting, use original file handles
            if direction == 'forward':
                return (self._fh1_forward, self._fh2_forward) if self.paired_out else (self._fh_forward,)
            elif direction == 'reverse':
                return (self._fh1_reverse, self._fh2_reverse) if self.paired_out else (self._fh_reverse,)
        
        # Split mode, use dynamic file handles
        key = f"{prefix}_{direction}"
        if key not in self.file_handles:
            # Create new file handles with correct file paths
            if self.paired_out:
                # Extract directory and base name from original file path
                base_path = self.file1_forward.replace('_forward_1.fq.gz', '')
                file1 = f"{base_path}_{direction}_{prefix}_1.fq.gz"
                file2 = f"{base_path}_{direction}_{prefix}_2.fq.gz"
                self.file_handles[key] = (xopen(file1, mode='wb', threads=4, compresslevel=1), xopen(file2, mode='wb', threads=4, compresslevel=1))
            else:
                base_path = self.file_forward.replace('_forward.fq.gz', '')
                file1 = f"{base_path}_{direction}_{prefix}.fq.gz"
                self.file_handles[key] = (xopen(file1, mode='wb', threads=4, compresslevel=1),)
        return self.file_handles[key]

    def _process_split_data(self, r1_data, r2_data, direction):
        """Process split data, parse each read and write to corresponding files"""
        try:
            # Convert byte data to string
            r1_str = r1_data.decode('utf-8')
            r2_str = r2_data.decode('utf-8')
            
            # Split by lines
            r1_lines = r1_str.strip().split('\n')
            r2_lines = r2_str.strip().split('\n')
            
            # Process each read (4 lines per read)
            for i in range(0, len(r1_lines), 4):
                if i + 3 < len(r1_lines):
                    try:
                        # Extract read name
                        read_name = r1_lines[i]
                        if read_name.startswith('@'):
                            read_name = read_name[1:]  # Remove @ symbol
                        
                        # Extract barcode prefix
                        barcode_prefix = self._get_barcode_prefix(read_name)
                        
                        # If unable to extract barcode prefix, write to default file
                        if not barcode_prefix:
                            if direction == 'forward':
                                read_data = '\n'.join(r1_lines[i:i+4]) + '\n'
                                self._fh1_forward.write(read_data.encode())
                                if i + 3 < len(r2_lines):
                                    read_data2 = '\n'.join(r2_lines[i:i+4]) + '\n'
                                    self._fh2_forward.write(read_data2.encode())
                            elif direction == 'reverse':
                                read_data = '\n'.join(r1_lines[i:i+4]) + '\n'
                                self._fh1_reverse.write(read_data.encode())
                                if i + 3 < len(r2_lines):
                                    read_data2 = '\n'.join(r2_lines[i:i+4]) + '\n'
                                    self._fh2_reverse.write(read_data2.encode())
                            continue
                        
                        # Get corresponding file handle
                        fh1, fh2 = self._get_file_handle(barcode_prefix, direction)
                        
                        # Write read data
                        read_data = '\n'.join(r1_lines[i:i+4]) + '\n'
                        fh1.write(read_data.encode())
                        
                        if i + 3 < len(r2_lines):
                            read_data2 = '\n'.join(r2_lines[i:i+4]) + '\n'
                            fh2.write(read_data2.encode())
                            
                    except Exception as e:
                        # Single read parsing failed, write to default file
                        print(f"Warning: Failed to process read at line {i+1}: {e}")
                        if direction == 'forward':
                            read_data = '\n'.join(r1_lines[i:i+4]) + '\n'
                            self._fh1_forward.write(read_data.encode())
                            if i + 3 < len(r2_lines):
                                read_data2 = '\n'.join(r2_lines[i:i+4]) + '\n'
                                self._fh2_forward.write(read_data2.encode())
                        elif direction == 'reverse':
                            read_data = '\n'.join(r1_lines[i:i+4]) + '\n'
                            self._fh1_reverse.write(read_data.encode())
                            if i + 3 < len(r2_lines):
                                read_data2 = '\n'.join(r2_lines[i:i+4]) + '\n'
                                self._fh2_reverse.write(read_data2.encode())
        except Exception as e:
            # Entire chunk parsing failed, write to default file
            print(f"Error: Failed to process chunk for {direction}: {e}")
            if direction == 'forward':
                self._fh1_forward.write(r1_data)
                self._fh2_forward.write(r2_data)
            elif direction == 'reverse':
                self._fh1_reverse.write(r1_data)
                self._fh2_reverse.write(r2_data)

    def write(self, data, index):
        self._chunks[index] = data
        if self.paired_out:
            while self._current_index in self._chunks:
                (_r1_forward, _r2_forward, _r1_reverse, _r2_reverse), (_multi_r1, _multi_r2) = self._chunks[self._current_index]
                
                # If splitting is enabled, need to parse each read and write to corresponding files separately
                if self.split_fastq > 0:
                    # Process forward reads
                    if _r1_forward:
                        self._process_split_data(_r1_forward, _r2_forward, 'forward')
                    
                    # Process reverse reads
                    if _r1_reverse:
                        self._process_split_data(_r1_reverse, _r2_reverse, 'reverse')
                else:
                    # No splitting, use original file handles
                    self._fh1_forward.write(_r1_forward)
                    self._fh2_forward.write(_r2_forward)
                    self._fh1_reverse.write(_r1_reverse)
                    self._fh2_reverse.write(_r2_reverse)
                
                self._fh_multi1.write(_multi_r1)
                self._fh_multi2.write(_multi_r2)
                del self._chunks[self._current_index]
                self._current_index += 1
        else:
            while self._current_index in self._chunks:
                (_r1_forward, _r1_reverse), (_multi_r1,) = self._chunks[self._current_index]
                
                # If splitting is enabled, need to parse read name to determine file handle
                if self.split_fastq > 0:
                    # Parse forward reads' read name
                    if _r1_forward:
                        forward_prefix = self._extract_barcode_prefix_from_bytes(_r1_forward)
                        if forward_prefix:
                            fh, = self._get_file_handle(forward_prefix, 'forward')
                            fh.write(_r1_forward)
                    
                    # Parse reverse reads' read name
                    if _r1_reverse:
                        reverse_prefix = self._extract_barcode_prefix_from_bytes(_r1_reverse)
                        if reverse_prefix:
                            fh, = self._get_file_handle(reverse_prefix, 'reverse')
                            fh.write(_r1_reverse)
                else:
                    # No splitting, use original file handles
                    self._fh_forward.write(_r1_forward)
                    self._fh_reverse.write(_r1_reverse)
                
                self._fh_multi1.write(_multi_r1)
                del self._chunks[self._current_index]
                self._current_index += 1            

    def wrote_everything(self):
        return not self._chunks

    def close(self):
        if self.paired_out:
            self._fh1_forward.close()
            self._fh2_forward.close()
            self._fh1_reverse.close()
            self._fh2_reverse.close()
            self._fh_multi1.close()
            self._fh_multi2.close()
        else:
            self._fh_forward.close()
            self._fh_reverse.close()
            self._fh_multi1.close()
        
        # Close all dynamically created file handles
        for file_handles in self.file_handles.values():
            for fh in file_handles:
                fh.close()

class Worker(Process):
    """Worker process class
    """
    def __init__(self, id_, read_pipe, write_pipe, need_work_queue, func, paired_out=False, split_fastq=0):
        """
        """
        super().__init__()
        self._id = id_
        self.read_pipe = read_pipe
        self.write_pipe = write_pipe
        self.need_work_queue = need_work_queue
        self.func = func
        self.paired_out = paired_out
        self.split_fastq = split_fastq

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
                data = self.read_pipe.recv_bytes()
                input = io.BytesIO(data)
                data = self.read_pipe.recv_bytes()
                input2 = io.BytesIO(data)
                if self.paired_out:
                    tmp_forward = (io.BytesIO(), io.BytesIO())
                    tmp_reverse = (io.BytesIO(), io.BytesIO())
                    tmp_multi = (io.BytesIO(), io.BytesIO())
                    _ = self.func(fq1=input, fq2=input2, fq_out_forward=tmp_forward, fq_out_reverse=tmp_reverse, fqout_multi=tmp_multi)
                    self.write_pipe.send(chunk_index)
                    self.write_pipe.send_bytes(tmp_forward[0].getvalue())
                    self.write_pipe.send_bytes(tmp_forward[1].getvalue())
                    self.write_pipe.send_bytes(tmp_reverse[0].getvalue())
                    self.write_pipe.send_bytes(tmp_reverse[1].getvalue())
                    self.write_pipe.send_bytes(tmp_multi[0].getvalue())
                    self.write_pipe.send_bytes(tmp_multi[1].getvalue())
                    self.write_pipe.send(_)
                else:
                    tmp_forward = (io.BytesIO(),)
                    tmp_reverse = (io.BytesIO(),)
                    tmp_multi = (io.BytesIO(),)
                    _ = self.func(fq1=input, fq2=input2, fq_out_forward=tmp_forward, fq_out_reverse=tmp_reverse, fqout_multi=tmp_multi)
                    self.write_pipe.send(chunk_index)
                    self.write_pipe.send_bytes(tmp_forward[0].getvalue())
                    self.write_pipe.send_bytes(tmp_reverse[0].getvalue())
                    self.write_pipe.send_bytes(tmp_multi[0].getvalue())
                    self.write_pipe.send(_)
            self.write_pipe.send(-1)
        except Exception as e:
            self.write_pipe.send(-2)
            raise e

class Pipeline:
    def __init__(self, func, fq1, fq2, fqout, fqout_multi, core, stat=None, paired_out=False, buffer_size=16*1024**2, split_fastq=0, samplename=""):
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
        self.split_fastq = split_fastq
        self.samplename = samplename

    def run(self):
        # start reader process
        reader_connections = [Pipe(duplex=False) for _ in range(self.n_workers)]
        _pipes, _conn = zip(*reader_connections)
        _reader_process = Reader(self.fq1, self.fq2, _conn, self.need_work_queue, self.buffer_size)
        _reader_process.daemon = True
        _reader_process.start()

        # start worker processes
        self.workers = []
        self.connections = []
        self.writer = Writer(self.fqout, self.fqout_multi, self.paired_out, self.split_fastq, self.samplename)
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
                chunk_index = connection.recv()
                if chunk_index == -1:
                    self.connections.remove(connection)
                    continue
                elif chunk_index == -2:
                    sys.stderr.write('err!!!\n')
                # if single?
                data1_forward = connection.recv_bytes()
                data2_forward = connection.recv_bytes()
                data1_reverse = connection.recv_bytes()
                data2_reverse = connection.recv_bytes()
                data_multi1 = connection.recv_bytes()
                data_multi2 = connection.recv_bytes()
                self.writer.write([(data1_forward, data2_forward, data1_reverse, data2_reverse), (data_multi1, data_multi2)], chunk_index)
                _stat = connection.recv()
                self.stat.update(**_stat)
        assert self.writer.wrote_everything()
        for w in self.workers:
            w.join()
        _reader_process.join()
        self.writer.close()

class AdapterFilter:
    """Filter adapters"""
    def __init__(self, adapter1:list=[], adapter2:list=[], non_insert_len:int=52, chemistry:str="DD-M"):
        self.adapter1 = []
        self.non_insert_len = non_insert_len
        self.chemistry = chemistry
        
        # Support error rate and overlap parameters
        # sequence: str,
        # max_errors: float = 0.1,
        # min_overlap: int = 16,
        for p in adapter1:
            if len(p) >= 4:  # If there are additional parameters (error rate and overlap)
                if p[1] == "3":
                    self.adapter1.append(BackAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3]))
                elif p[1] == "5":
                    self.adapter1.append(FrontAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3]))
            else:  # Default parameters
                if p[1] == "3":
                    self.adapter1.append(BackAdapter(sequence=p[0], min_overlap=10))
                elif p[1] == "5":
                    self.adapter1.append(FrontAdapter(sequence=p[0], min_overlap=7))
        
        self.adapter2 = []
        for p in adapter2:
            if len(p) >= 4:  # If there are additional parameters (error rate and overlap)
                if p[1] == "3":
                    self.adapter2.append(BackAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3]))
                elif p[1] == "5":
                    self.adapter2.append(FrontAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3]))
            else:  # Default parameters
                if p[1] == "3":
                    self.adapter2.append(BackAdapter(sequence=p[0], min_overlap=10))
                elif p[1] == "5":
                    self.adapter2.append(FrontAdapter(sequence=p[0], min_overlap=10))
    
    def filter(self, r1=None, r2=None) -> tuple:
        flag = False
        r1_me_left = False
        r1_me_right = False
        if r1 and self.adapter1:
            sp2_target = False  # Flag for whether target adapter is detected
            for _ in self.adapter1:
                # print(_)
                m = _.match_to(r1.sequence)
                # print(f'm: {m}')
                if m:
                    if _.sequence == "AGATGTGTATAAGAGAYAG": 
                        r1_me_left = True                    
                    if _.sequence == "CTGTCTCTTATACACATCT":
                        r1_me_right = True
                    flag = True
                    r1 =  m.trimmed(r1)
                    # print(r1.sequence)
            #r1_start = 9
            r1_start = 0 # qiuxia 
            r1_end = len(r1.sequence)
            if len(r1.sequence) > 18:
                if r1_me_left:
                    r1_start = 9 # cutadapter rm me left, rm additial 9bp for uncorret methylation of Transposase
                else:
                    r1_start = self.non_insert_len # cutadapter not rm me left, start from non_insert_len
                if r1_me_right:
                    r1_end = len(r1.sequence) - 9 # cutadapter cut me right, rm additial 9bp for uncorret methylation of Transposase
                max_length = ME5_MAX_READ1_LEN if self.chemistry == "ME5" else DEFAULT_MAX_READ1_LEN
                if r1_end - r1_start > max_length:
                    r1_end = r1_start + max_length
                r1.sequence = r1.sequence[r1_start:r1_end]
                r1.qualities = r1.qualities[r1_start:r1_end]
        r2_me_left = False
        r2_me_right = False
        if r2 and self.adapter2:
            me_target = False
            for _ in self.adapter2:
                m = _.match_to(r2.sequence)
                if m:
                    if _.sequence == "AGATGTGTATAAGAGACAG":  # 'adapter2': ["CTGTCTCTTATACACATCT", "3"],  me-rev
                        #print(f'CTRTCTCTTATACACATCT')
                        r2_me_left = True
                    if _.sequence == "CTRTCTCTTATACACATCT":  # 'adapter2': ["AGATGTGTATAAGAGACAG", "5"],  me
                        r2_me_right = True
                    flag = True
                    r2 =  m.trimmed(r2)
            r2_start = 9 # trim 9bp, if left adapter detected or not  
            r2_end = len(r2.sequence)
            if len(r2.sequence) > 18:
                if r2_me_right:
                    r2_end = len(r2.sequence) - 9  # cutadapter cut me right, rm additial 9bp for uncorret methylation of Transposase
            r2.sequence = r2.sequence[r2_start:r2_end]
            r2.qualities = r2.qualities[r2_start:r2_end] 
        return flag, r1, r2

class QcStat:
    """Summary statistics"""
    def __init__(self):
        self.data = { }

    def update(self, **d):
        if not self.data:
            self.data = d
        else:
            for k, v in d.items():
                self.data[k] += v

    @staticmethod
    def _sort_gc(d):
        if not d:  # If dictionary is empty, return empty result structure
            return {b: [] for b in 'ATCGN'}
        idx_max = max([k[0] for k in d])
        return {
            b: [d.get((i, b), 0) for i in range(idx_max+1)] for b in 'ATCGN'
        }

    @staticmethod
    def _sort_q(d, phred=33):
        if not d:  # If dictionary is empty, return empty result structure
            return {}
        idx_max = max([k[0] for k in d])
        q_max = max([ord(k[1])-phred for k in d])
        return {
            i: [d.get((i, chr(q+phred)), 0) for q in range(q_max+1)] for i in range(idx_max+1)
        }

    def save(self, path='summary.json'):
        tmp = {'__version__': 'v1.0.0'}
        for k in self.data:
            if k.endswith('_gc'):
                tmp[k] = self._sort_gc(self.data[k])
            elif k.endswith('_q'):
                tmp[k] = self._sort_q(self.data[k])
            else:
                tmp[k] = dict(self.data[k])
        with open(path, 'w') as fh:
            json.dump(tmp, fh, indent=4)

def parse_structure(string:str) -> tuple:
    """Parse adapter structure

    Use letters B, L, U, X and T with numbers to represent reads structure.
    B represents barcode bases;
    L represents linker bases;
    U represents UMI bases;
    X represents any bases, used as placeholder;
    T represents T bases;
    Numbers after letters indicate base length.

    Args:
        string: Adapter structure description

    Returns:
        Two-dimensional tuple containing ordered structure parts and lengths.
        For example:
        When string is B8L8B8L10B8U8, returns:
            (('B', 8), ('L', 8), ('B', 8), ('L', 10), ('B', 8), ('U', 8))
    """
    regex = re.compile(r'([BLUXT])(\d+)')
    groups = regex.findall(string)
    return tuple([(_[0], int(_[1])) for _ in groups])

def read_file(file_list: list) -> dict:
    """Prepare whitelist set
        Args:
            file_list: Path to each segment whitelist file
    """
    wl_dict = dict()
    for i, wl_file in enumerate(file_list):
        white_list = set()
        with xopen(wl_file, "r", threads=4) as fh:
            for l in fh:
                if l.startswith("#"):
                    continue
                la = l.strip()
                if not la:
                    continue
                white_list.add(la)
        wl_dict[i] = white_list
    return wl_dict

def get_new_bc(bc:str, white_list:set, distance:int)->set:
    """Return the intersection of the set after mismatching at each position of the original barcode with the whitelist set"""

    if distance == 1:
        BASE_LIST = ["T", "C", "G", "A"]
        mm_dict = dict()
        for i, c in enumerate(bc):
            if c == "N":
                mm_dict = { bc[:i] + base + bc[i+1:]:f"{i}{base}" for base in BASE_LIST }
                break  
            else:
                mm_dict.update({ bc[:i] + base + bc[i+1:]:f"{i}{base}" for base in BASE_LIST if base!=c })
                
        bc_set = set(mm_dict.keys()).intersection(white_list)
        # return {k: mm_dict[k] for k in bc_set}
    else:
        bc_dict = defaultdict(set)
        for bc_true in white_list:
            hmm = sum(ch1 != ch2 for ch1,ch2 in zip(bc_true,bc))
            if hmm <= distance:
                bc_dict[hmm].add(bc_true)
                
        bc_set = set()
        if len(bc_dict) != 0:
            sorted_items = sorted(bc_dict.items(), key=lambda x: x[0])
            bc_set = sorted_items[0][1]
    
    return bc_set

def hamming_distance(s1, s2):
    return len([(i, j) for i, j in zip(s1, s2) if i != j])

def calculate_average_quality(quality_string):
    # Assuming Phred+33 encoding
    quality_scores = [ord(char) - 33 for char in quality_string]
    average_quality = sum(quality_scores) / len(quality_scores)
    return average_quality

def calculate_ct_conversion_rate(r1, stat_Dict, positions=None, chemistry="DD-M"):
    """
    Calculate CT conversion rate for methylation analysis
    
    This function analyzes specific positions in the R1 sequence to detect methylation conversion patterns:
    - Check 7f sequence pattern (TTGCTGT or TTGTTGT)
    - Check 17bp linker methylase sequence (GTAGATGTGTATAAGAGA)
    - Extract bases and quality values at specified positions
    - Calculate C and T conversion rates based on insert sequence type (CC or TT)
    
    Args:
        r1: R1 sequence object containing sequence and qualities attributes
        stat_Dict: Statistics dictionary for recording various counts
        positions: Key position array, uses default positions if None
        chemistry: Chemistry type for selecting different parameters
        
    Returns:
        dict: Contains updated statistical information
    """
    # If no positions provided, use default positions
    if positions is None:
        positions = DEFAULT_POSITIONS
    
    # Set sequence extraction parameters based on chemistry parameter
    if chemistry == "ME5":
        seq_7f_start = ME5_SEQ_7F_START
        seq_7f_len = ME5_SEQ_7F_LEN
        l17me_start = ME5_L17ME_START
        l17me_len = ME5_L17ME_LEN
    else:
        seq_7f_start = DEFAULT_SEQ_7F_START
        seq_7f_len = DEFAULT_SEQ_7F_LEN
        l17me_start = DEFAULT_L17ME_START
        l17me_len = DEFAULT_L17ME_LEN
    
    # Check if sequence length is sufficient, use maximum value in positions for dynamic checking
    if positions and len(positions) >= 7:
        max_position = max(positions)
        if len(r1.sequence) <= max_position:
            return stat_Dict
    else:
        # If no positions information, use conservative length checking
        if len(r1.sequence) < 80:
            return stat_Dict
    
    # Extract key sequence fragments
    seq_7f = r1.sequence[seq_7f_start:(seq_7f_start+seq_7f_len)]
    l17me = r1.sequence[l17me_start:(l17me_start+l17me_len)]
    
    # Check 17bp linker methylase sequence
    if l17me == 'GTAGATGTGTATAAGAGA':
        stat_Dict["num_17lme"] += 1
    
    # Check 7f sequence pattern
    if chemistry != "ME5" and (seq_7f != 'TTGCTGT' and seq_7f != 'TTGTTGT'):
        # In non-ME5 mode, if 7f sequence doesn't match, return directly
        return stat_Dict
    
    # 7f sequence matches, increase count
    if chemistry != "ME5":
        stat_Dict["num_7f"] += 1

    
    # If 17lme condition is met, perform detailed analysis
    if l17me == 'GTAGATGTGTATAAGAGA':
        # If both 7f and 17lme conditions are met, increase count
        stat_Dict["num_7f17lme"] += 1
        
        # Extract bases and quality values at key positions
        bases = []
        qualities = []
        
        for pos in positions:
            bases.append(r1.sequence[pos])
            qualities.append(r1.qualities[pos])
        
        # Combine sequence and quality values
        allseq = ''.join(bases)
        # Now both modes have 7 positions, handle uniformly
        countseq = allseq[:5]  # First 5 positions for counting
        insert = bases[-2] + bases[-1]  # Second-to-last and last positions as insert sequence
        allqua = ''.join(qualities)
        
        # Calculate average quality value
        average_quality = calculate_average_quality(allqua)
        
        # Only perform conversion rate calculation when average quality >= 30
        if average_quality >= 30:
            if insert == 'CC': 
                # CC insert type: calculate C and T conversion rates
                B_count_C = countseq.count('C')
                B_count_T = countseq.count('T')
                if B_count_C + B_count_T == len(countseq):
                    stat_Dict["line_B"] += 1
                    stat_Dict["B_all_C"] += B_count_C
                    
            elif insert == 'TT':
                # TT insert type: calculate T and C conversion rates
                A_count_T = countseq.count('T')
                A_count_C = countseq.count('C')
                if A_count_T + A_count_C == len(countseq):  # Must be C or T, other bases indicate sequencing errors
                    stat_Dict["line_A"] += 1
                    stat_Dict["A_all_T"] += A_count_T
    
    return stat_Dict

def generate_random_string(length=12):
    return ''.join(random.choices(string.ascii_uppercase, k=length))

def should_filter_read_ch_pattern(r1_sequence, r2_sequence, chain_direction, threshold=0):
    """
    Check if reads should be filtered due to CH pattern
    
    Args:
        r1_sequence: R1 sequence string
        r2_sequence: R2 sequence string
        chain_direction: Chain direction ('forward' or 'reverse')
        threshold: Filter threshold, 0 means no filtering, >0 means filter when pattern appears more than this number of times
        
    Returns:
        bool: True means the read should be filtered, False means keep it
    """
    # If threshold is 0, no filtering is performed
    if threshold == 0:
        return False
    
    if chain_direction == 'forward':
        # Forward chain direction filter conditions
        # R1: C[ATC] exceeds threshold times
        # R2: [ATG]G exceeds threshold times
        r1_count = len(R1_FORWARD_CH_PATTERN.findall(r1_sequence))
        r2_count = len(R2_FORWARD_CH_PATTERN.findall(r2_sequence))
    elif chain_direction == 'reverse':
        # Reverse chain direction filter conditions
        # R1: [ATG]G exceeds threshold times
        # R2: C[ATC] exceeds threshold times
        r1_count = len(R1_REVERSE_CH_PATTERN.findall(r1_sequence))
        r2_count = len(R2_REVERSE_CH_PATTERN.findall(r2_sequence))
    else:
        # Unknown chain direction, no filtering
        return False
    
    # If any pattern appears more than threshold times, filter the read
    return r1_count > threshold or r2_count > threshold

def determine_chain_direction(r1_sequence, chemistry="DD-M", positions=None):
    """
    Determine chain direction based on R1 sequence
    
    Args:
        r1_sequence: R1 sequence string
        chemistry: Chemistry type for determining chain direction
        positions: Key position array for extracting insert sequence
        
        Returns:
        str: chain_direction
            - chain_direction: 'forward', 'reverse', or 'unknown'
    """
    chain_direction = 'unknown'
    
    # Use maximum value in positions for dynamic length checking
    if positions and len(positions) >= 7:
        max_position = max(positions)
        if len(r1_sequence) > max_position:
            # Extract positions 1, 6, 7 from positions (indices 0, 5, 6)
            insert_positions = [positions[0], positions[5], positions[6]]
            insert = ''.join([r1_sequence[pos] for pos in insert_positions])
            
            # Determine chain direction based on C count in insert
            C_count = insert.count('C')
            if C_count == 3:  # Contains 3 Cs
                chain_direction = 'reverse'
            else:  # Other cases
                chain_direction = 'forward'
    
    return chain_direction


def summary(seq, seq_q, seq_dict, qua_dict):
    
    for i, (base, q) in enumerate(zip(seq, seq_q)):
        seq_dict[(i,base)] += 1
        qua_dict[(i,q)] += 1
        
    return seq_dict, qua_dict

def process_barcode(fq1, fq2, fq_out_forward, fq_out_reverse, fqout_multi, r1_structure, shift, shift_pattern,
                    barcode_wl_dict, linker_wl_dict, match_type_dict, adapter1=[["AAAAAAAAAAAA", "3"],],
                    adapter2=[["AAAAAAAAAAAA", "3"],], do_B_correction=True, do_L_correction=True,
                    use_multi=True, use_short_read=False, paired_out=True, positions=None, non_insert_len=52, chemistry="DD-M", split_fastq=0, filter_ch=0):
    
    barcode_list_flag = False
    linker_list_flag = False
    if len(barcode_wl_dict)>0:
        barcode_list_flag = True

    if len(linker_wl_dict)>0:
        linker_list_flag = True

    stat_Dict = defaultdict(int)
    Barcode_Counter = Counter()
    
    Barcode_GC_Counter = Counter()
    UMI_GC_Counter = Counter()
    R2_GC_Counter = Counter()
    Barcode_Q_Counter = Counter()
    UMI_Q_Counter = Counter()
    R2_Q_Counter = Counter()
    
    adapter_filter = AdapterFilter(adapter1=adapter1, adapter2=adapter2, non_insert_len=non_insert_len, chemistry=chemistry)
    
    fh = dnaio.open(fq1, fq2, fileformat="fastq", mode="r", open_threads=4)
    if paired_out:
        outfh_forward = dnaio.open(fq_out_forward[0], fq_out_forward[1], fileformat="fastq", mode="w", open_threads=4)
        outfh_reverse = dnaio.open(fq_out_reverse[0], fq_out_reverse[1], fileformat="fastq", mode="w", open_threads=4)
    else:
        outfh_forward = dnaio.open(fq_out_forward[0], fileformat="fastq", mode="w", open_threads=4)
        outfh_reverse = dnaio.open(fq_out_reverse[0], fileformat="fastq", mode="w", open_threads=4)

    if use_multi:
        if paired_out:
            outfh_multi = dnaio.open(fqout_multi[0], fqout_multi[1], fileformat="fastq", mode="w", open_threads=4)
        else:
            outfh_multi = dnaio.open(fqout_multi[0], fileformat="fastq", mode="w", open_threads=4)
    
    for r1, r2 in fh:
        stat_Dict["total"] += 1
        
        start_pos = 0
        end_pos = 0
        sequence = r1.sequence
        qualities = r1.qualities
        
        # deal with shift
        if shift:
            shift_pos = sequence[:7].find(shift_pattern)
            if shift_pos < 0:
                stat_Dict["no_anchor"] += 1
                # logger.debug(f"{r1.name},{sequence},{sequence[:7]} no anchor!")
                continue
            else:
                start_pos = shift_pos + 1
        
        # get barcode/umi/quality sequence          
        old_seqs = defaultdict(list)
        new_seqs = defaultdict(list)
        seq_quals = defaultdict(list)
        B = 0
        L = 0
        is_valid = True
        is_multi = False
        is_correct = False
        is_B_no_correction = False
        is_L_no_correction = False

        for _, (code, n) in enumerate(r1_structure):
            end_pos = start_pos + n
            seq = sequence[start_pos:end_pos]
            quals = qualities[start_pos:end_pos]

            if code == "B":
                old_seqs["B"].append(seq)
                seq_quals["B"].append(quals)

                if barcode_list_flag: # match barcode in whitelist
                    if seq in barcode_wl_dict.get(B, barcode_wl_dict[0]):
                        new_seqs["B"].append({seq})
                    else:
                        if do_B_correction:
                            bc_set = get_new_bc(seq, barcode_wl_dict.get(B, barcode_wl_dict[0]), match_type_dict.get(B, match_type_dict[0]))

                            if len(bc_set) == 0:
                                is_valid = False
                                is_B_no_correction = True
                                # logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} no barcode!")
                                break
                            elif len(bc_set) == 1:
                                new_seqs["B"].append(bc_set)
                                is_correct = True
                                # logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(bc_set)} do_B_correction!")
                            else:
                                new_seqs["B"].append(bc_set)
                                # logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(bc_set)} do_B_correction!")
                                is_multi = True
                        else:
                            is_valid = False
                            break
                else:
                    new_seqs["B"].append({seq})
                B += 1

            elif code == "L":   
                if linker_list_flag: # linker correction step
                    if seq in linker_wl_dict.get(L, linker_wl_dict[0]):
                        pass
                    else:
                        if do_L_correction:
                            lk_set = get_new_bc(seq, linker_wl_dict.get(L, linker_wl_dict[0]))
                            if len(lk_set) == 0:
                                is_valid = False
                                is_L_no_correction = True
                                # logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(lk_set)} no linker!")
                                break
                        else:
                            is_valid = False
                            break
                L += 1
                
            elif code == "U":
                new_seqs["U"].append(seq)
                seq_quals["U"].append(quals)
                
            start_pos = start_pos + n

        # check double instances and calculate CT conversion rate
        if is_valid:
            # Use new CT conversion rate calculation function
            stat_Dict = calculate_ct_conversion_rate(r1, stat_Dict, positions, chemistry)
            
            # Use function to determine chain direction
            chain_direction = determine_chain_direction(r1.sequence, chemistry, positions)
            barcode_old = "".join(old_seqs["B"])
            Barcode_Counter[barcode_old] += 1

            #get base summary for umi/r2
            umi = "".join(new_seqs["U"])
            umi_q = "".join(seq_quals["U"])
            barcode_q = "".join(seq_quals["B"])
            
            UMI_GC_Counter, UMI_Q_Counter = summary(umi, umi_q, UMI_GC_Counter, UMI_Q_Counter)
            R2_GC_Counter, R2_Q_Counter = summary(r2.sequence, r2.qualities, R2_GC_Counter, R2_Q_Counter)
            
            r1.sequence = sequence[start_pos:]
            r1.qualities = qualities[start_pos:]
                        
            if is_multi: #write r2 multi files
                if use_multi:         
                    #update barcode quality
                    Barcode_Q_Counter.update(enumerate(barcode_q))
                    bc_new_lst = []
                    for element in itertools.product(*new_seqs["B"]):          
                        barcode_new = "".join(element)
                        bc_new_lst.append(barcode_new)
                        
                    # Simplify multi file read name format, excluding chain direction information
                    # Store candidate barcode list and original read name separately
                    bc_new_all = ":".join(bc_new_lst)
                    # Format: barcode_old_candidate_list_umi_original_read_name
                    r2.name = "_".join([barcode_old, bc_new_all, umi, r2.name])
                    r1.name = "_".join([barcode_old, bc_new_all, umi, r1.name])
                    r1.sequence = sequence[start_pos:]
                    r1.qualities = qualities[start_pos:]
                    
                    outfh_multi.write(r1, r2)
            else:  #write r2 files based on chain direction
                stat_Dict["valid"] += 1
                flag, r1, r2 = adapter_filter.filter(r1, r2)
                if flag:
                    if (not use_short_read) or len(r1) == 0 or len(r2) == 0:
                        if len(r1) < R1_MINLEN or len(r2) < R2_MINLEN:
                            stat_Dict["too_short"] += 1
                            continue
                    else:
                        stat_Dict["trimmed"] += 1

                barcode_new = "".join([_.pop() for _ in new_seqs["B"]])
                Barcode_GC_Counter, Barcode_Q_Counter = summary(barcode_old, barcode_q, Barcode_GC_Counter, Barcode_Q_Counter)
                
                #find alterations
                if is_correct:
                    _alt = "".join([str(i)+o for i, (o,n) in enumerate(zip(barcode_old, barcode_new)) if o != n])
                else:
                    _alt = "M"

                # Determine UMI format based on chemistry type
                if chemistry == "ME5":
                    final_umi = umi  # ME5 mode uses original UMI directly
                else:
                    # DD-M mode uses first 12bp of final output r1 sequence as umi
                    final_umi = r1.sequence[:12] if len(r1.sequence) >= 12 else r1.sequence

                # Add barcode to end of read name
                new_r2_read_name = f'{r2.name}:{barcode_new}'
                new_r1_read_name = f'{r1.name}:{barcode_new}'
                
                r2.name = "_".join([barcode_new, final_umi, chain_direction, _alt, new_r2_read_name])
                r1.name = "_".join([barcode_new, final_umi, chain_direction, _alt, new_r1_read_name])

                # First count all reads determined as forward and reverse
                if chain_direction == 'forward':
                    stat_Dict["forward"] += 1
                elif chain_direction == 'reverse':
                    stat_Dict["reverse"] += 1
                
                # CH mode filtering: perform final filtering before output
                if should_filter_read_ch_pattern(r1.sequence, r2.sequence, chain_direction, filter_ch):
                    # Count filtered reads by chain direction separately
                    if chain_direction == 'forward':
                        stat_Dict["forward_chimeric"] += 1
                    elif chain_direction == 'reverse':
                        stat_Dict["reverse_chimeric"] += 1
                else:
                    # Output to different files by chain direction (only output filtered reads)
                    if chain_direction == 'forward':
                        outfh_forward.write(r1, r2)
                    elif chain_direction == 'reverse':
                        outfh_reverse.write(r1, r2)
            if is_correct:
                stat_Dict["B_corrected"] += 1
                # Count B_corrected by chain direction
                if chain_direction == 'forward':
                    stat_Dict["B_corrected_forward"] += 1
                elif chain_direction == 'reverse':
                    stat_Dict["B_corrected_reverse"] += 1
        else:
            if is_B_no_correction:
                stat_Dict["B_no_correction"] += 1

            if is_L_no_correction:
                stat_Dict["L_no_correction"] += 1
    
    # Close file handles
    fh.close()
    if use_multi:
        outfh_multi.close()
    outfh_forward.close()
    outfh_reverse.close()

    return {
            "stat": Counter(stat_Dict),
            "barcode_count": Barcode_Counter,
            "barcode_gc": Barcode_GC_Counter,
            "umi_gc": UMI_GC_Counter,
            "r2_gc": R2_GC_Counter,
            "barcode_q": Barcode_Q_Counter,
            "umi_q": UMI_Q_Counter,
            "r2_q": R2_Q_Counter
        }

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("--fq1", "fq1", required=True, type=click.Path(), multiple=True, help="Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(), multiple=True, help="Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--barcode", multiple=True, required=True, help="Barcode white list file, can specify multiple times.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--chemistry", required=True, type=click.Choice(["DD-M","ME5"]), help="chemistry")
@click.option("--split_fastq", default=0, type=int, show_default=True, help="Split output by first n bases of barcode (0=no split)")
@click.option("--filter_ch", default=0, type=int, show_default=True, help="CH pattern filtering threshold (0=no filtering, >0=filter when pattern appears more than this number of times)")

def barcode_main(chemistry, fq1:list, fq2:list, samplename: str, outdir:str,
                 barcode:list=[], match_type:list=[], shift:str=True, shift_pattern:str="A",
                 structure:str="B17U12X7", linker: list=[],
                 core:int=4, do_B_correction=True, do_L_correction=True,
                 use_multi=True, use_short_read=False, adapter1=[["TTTTTTTTTTTT", "5"], ],
                 adapter2=[["AAAAAAAAAAAA", "3"], ], paired_out=True, split_fastq:int=0, filter_ch:int=0):
    shift = CHEMISTRY[chemistry]['shift']
    structure = CHEMISTRY[chemistry]['structure']
    adapter1 = CHEMISTRY[chemistry]['adapter1']
    adapter2 = CHEMISTRY[chemistry]['adapter2']
    match_type = CHEMISTRY[chemistry]['match_type']

    logger.info("extract barcode start!")
    #parse r1 structure
    r1_structure = parse_structure(structure)
  
    #get wl dict for bc/linker
    # When barcode or linker is empty, return empty dictionary
    barcode_wl_dict = read_file(barcode)
    linker_wl_dict = read_file(linker)
    print(f'linker_wl_dict : {linker_wl_dict}')
    match_type_dict = {ind: val for ind, val in enumerate(match_type)}
    print(f'match_type_dict : {match_type_dict}')

    if len(barcode_wl_dict)>0 and do_B_correction:
        logger.info("barcode one base mismatch allowed.")
    else:
        logger.info("barcode mismatch NOT allowed.")

    if "L" in structure:
        if len(linker_wl_dict)>0 and do_L_correction:
            logger.info("linker one base mismatch allowed.")
        else:
            logger.info("linker mismatch NOT allowed.")

    if use_multi:
        logger.info("rescue barcode match multi barcode in whitelist.")
    else:
        logger.info("ignore barcode match multi barcode in whitelist.")
    
    # Select positions and non_insert_len based on chemistry parameter
    if chemistry == "ME5":
        positions = ME5_POSITIONS
        non_insert_len = ME5_NON_INSERT_LEN
        logger.info(f"Chemistry '{chemistry}' detected - using ME5 mode positions for CT conversion rate calculation.")
        logger.info(f"Using ME5 mode non_insert_len: {non_insert_len}")
    else:
        positions = DEFAULT_POSITIONS
        non_insert_len = DEFAULT_NON_INSERT_LEN
        logger.info(f"Chemistry '{chemistry}' detected - using default positions for CT conversion rate calculation.")
        logger.info(f"Using default non_insert_len: {non_insert_len}")
    
    #worker function
    worker_func = partial(
        process_barcode,
        r1_structure=r1_structure,
        shift=shift,
        shift_pattern=shift_pattern,
        barcode_wl_dict=barcode_wl_dict,
        linker_wl_dict=linker_wl_dict,
        match_type_dict=match_type_dict,
        do_B_correction=do_B_correction,
        do_L_correction=do_L_correction,
        use_multi=use_multi,
        use_short_read=use_short_read,
        adapter1=adapter1,
        adapter2=adapter2,
        positions=positions,
        non_insert_len=non_insert_len,
        chemistry=chemistry,
        split_fastq=split_fastq,
        filter_ch=filter_ch,
    )
    
    stat = QcStat()
    
    os.makedirs(f"{outdir}/step1", exist_ok=True)
    fqout = os.path.join(outdir, f"step1/{samplename}")
    fqout_multi = os.path.join(outdir, f"step1/{samplename}_multi")
    json_multi = os.path.join(outdir, f"step1/{samplename}_multi.json")
    
    pipeline = Pipeline(
        func=worker_func,
        fq1=fq1,
        fq2=fq2,
        fqout=fqout,
        fqout_multi=fqout_multi,
        stat=stat,
        core=core,
        paired_out=paired_out,
        split_fastq=split_fastq,
        samplename=samplename
    )
    pipeline.run()

    fqout1 = f"{fqout}_1.fq.gz"
    fqout2 = f"{fqout}_2.fq.gz"
    if use_multi:
        # find the multiple barcodes
        logger.info("deal multi start!")
        fqout_multi1 = f"{fqout_multi}_1.fq.gz"
        fqout_multi2 = f"{fqout_multi}_2.fq.gz"
        adapter_filter = AdapterFilter(adapter1=adapter1, adapter2=adapter2, non_insert_len=non_insert_len, chemistry=chemistry)
        multi_stat = defaultdict(int)
        # Create output files for both chain directions for multi rescue
        fqout_forward = f"{fqout}_forward"
        fqout_reverse = f"{fqout}_reverse"
        
        # If splitting is enabled, need to create dynamic file handles
        if split_fastq > 0:
            forward_files = {}  # {prefix: (file1, file2)}
            reverse_files = {}  # {prefix: (file1, file2)}
            
            def get_multi_output_handle(prefix, direction):
                """Get or create output file handle for multi rescue"""
                if direction == 'forward':
                    if prefix not in forward_files:
                        file1 = f"{fqout_forward}_{prefix}_1.fq.gz"
                        file2 = f"{fqout_forward}_{prefix}_2.fq.gz"
                        forward_files[prefix] = dnaio.open(file1, file2, mode="a", open_threads=4)
                    return forward_files[prefix]
                elif direction == 'reverse':
                    if prefix not in reverse_files:
                        file1 = f"{fqout_reverse}_{prefix}_1.fq.gz"
                        file2 = f"{fqout_reverse}_{prefix}_2.fq.gz"
                        reverse_files[prefix] = dnaio.open(file1, file2, mode="a", open_threads=4)
                    return reverse_files[prefix]
                return None
        
        with dnaio.open(fqout_forward + "_1.fq.gz", fqout_forward + "_2.fq.gz", mode="a", open_threads=4) as f_forward:
            with dnaio.open(fqout_reverse + "_1.fq.gz", fqout_reverse + "_2.fq.gz", mode="a", open_threads=4) as f_reverse:
                fh = dnaio.open(fqout_multi1, fqout_multi2, fileformat="fastq", mode="r", open_threads=4)
                for r1, r2 in fh:
                    multi_stat["total"] += 1
                    final_barcode = None
                    
                    # Parse multi file read name format: barcode_old_candidate_list_umi_original_read_name
                    bc_old, r2_candidate, umi, r2_name = r2.name.split("_", 3)
                    r2_candidate = r2_candidate.split(":")
                    # r2_name now only contains original read name, not candidate barcode list
                    
                    # Similarly parse r1.name to ensure format consistency
                    _, _, _, r1_name = r1.name.split("_", 3)
                    
                    read_num = 0
                    for _ in sorted(r2_candidate):
                        v = stat.data["barcode_count"].get(_, 0)
                        if v > read_num:
                            read_num = v
                            final_barcode = _

                    if not final_barcode:
                        multi_stat["B_no_correction"] += 1
                        stat.data["stat"]["B_no_correction"] += 1
                        continue
                        
                    multi_stat["valid"] += 1
                    stat.data["stat"]["valid"] += 1

                    flag, r1, r2 = adapter_filter.filter(r1, r2)
                    if flag:
                        if (not use_short_read) or len(r1) == 0 or len(r2) == 0:
                            if len(r1) < R1_MINLEN or len(r2) < R2_MINLEN:
                                multi_stat["too_short"] += 1
                                stat.data["stat"]["too_short"] += 1
                                continue
                        else:
                            multi_stat["trimmed"] += 1
                            stat.data["stat"]["trimmed"] += 1
                    else: # input r1 or r2 witout adapter, but was short after fastp
                        if len(r1) < R1_MINLEN or len(r2) < R2_MINLEN:
                            multi_stat["too_short"] += 1
                            stat.data["stat"]["too_short"] += 1
                            continue                   
                    alt_l = [str(i)+o for i, (o,n) in enumerate(zip(bc_old, final_barcode)) if o != n]
                    _alt = "".join([alt for alt in alt_l])

                    # Re-determine chain direction (based on finally determined barcode)
                    chain_direction = determine_chain_direction(r1.sequence, chemistry, positions)
                    
                    # Determine UMI format based on chemistry type
                    if chemistry == "ME5":
                        final_umi = umi  # ME5 mode uses original UMI directly
                    else:
                        # DD-M mode uses first 12bp of final output r1 sequence as umi
                        final_umi = r1.sequence[:12] if len(r1.sequence) >= 12 else r1.sequence
                    
                    # Add barcode to read name in specified format, consistent with non-multi format
                    # Note: r2_name and r1_name do not contain candidate barcode list, only original read name
                    new_r2_read_name = f'{r2_name}:{final_barcode}'
                    new_r1_read_name = f'{r1_name}:{final_barcode}'
                    
                    r2.name = "_".join([final_barcode, final_umi, chain_direction, _alt, new_r2_read_name])
                    r1.name = "_".join([final_barcode, final_umi, chain_direction, _alt, new_r1_read_name])
                    
                    # First count all reads determined as forward and reverse
                    if chain_direction == 'forward':
                        multi_stat["forward"] += 1
                        stat.data["stat"]["forward"] += 1
                    elif chain_direction == 'reverse':
                        multi_stat["reverse"] += 1
                        stat.data["stat"]["reverse"] += 1
                    
                    # CH mode filtering: perform final filtering before output
                    if should_filter_read_ch_pattern(r1.sequence, r2.sequence, chain_direction, filter_ch):
                        # Count filtered reads by chain direction separately
                        if chain_direction == 'forward':
                            multi_stat["forward_chimeric"] += 1
                            stat.data["stat"]["forward_chimeric"] += 1
                        elif chain_direction == 'reverse':
                            multi_stat["reverse_chimeric"] += 1
                            stat.data["stat"]["reverse_chimeric"] += 1
                    else:
                        # Output to different files based on chain direction (only output filtered reads)
                        if split_fastq > 0:
                            # Split mode: write to corresponding split file based on barcode prefix
                            barcode_prefix = final_barcode[:split_fastq] if len(final_barcode) >= split_fastq else final_barcode
                            if chain_direction == 'forward':
                                outfh = get_multi_output_handle(barcode_prefix, 'forward')
                                if outfh:
                                    outfh.write(r1, r2)
                                else:
                                    f_forward.write(r1, r2)  # Fallback option
                            elif chain_direction == 'reverse':
                                outfh = get_multi_output_handle(barcode_prefix, 'reverse')
                                if outfh:
                                    outfh.write(r1, r2)
                                else:
                                    f_reverse.write(r1, r2)  # Fallback option
                        else:
                            # No split mode: write to default files
                            if chain_direction == 'forward':
                                f_forward.write(r1, r2)
                            elif chain_direction == 'reverse':
                                f_reverse.write(r1, r2)

        # Close all dynamically created file handles
        if split_fastq > 0:
            for outfh in forward_files.values():
                outfh.close()
            for outfh in reverse_files.values():
                outfh.close()

        with open(json_multi, "w") as fh:
            json.dump(multi_stat, fp=fh, indent=4)

    del stat.data["barcode_count"]
    logger.info("deal multi done!")
    stat.data["stat"]["chemistry"] = chemistry
    stat.data["stat"]["gexname"] = samplename
    # ct_mean = A_all_T / (5 * line_A)
    stat.data["stat"]["ct_mean"] = stat.data["stat"]["A_all_T"] / (5 * stat.data["stat"]["line_A"])
    # cc_mean = B_all_C / (5 * line_B)
    stat.data["stat"]["cc_mean"] = stat.data["stat"]["B_all_C"] / (5 * stat.data["stat"]["line_B"])
    # rate_7f = num_7f / total
    stat.data["stat"]["rate_7f"] = stat.data["stat"]["num_7f"] / stat.data["stat"]["total"]
    # rate_17lme = num_17lme / total
    stat.data["stat"]["rate_17lme"] = stat.data["stat"]["num_17lme"] / stat.data["stat"]["total"]
    # rate_7f17lme = num_7f17lme / total
    stat.data["stat"]["rate_7f17lme"] = stat.data["stat"]["num_7f17lme"] / stat.data["stat"]["total"]
    stat.save(os.path.join(outdir, f"{samplename}_summary.json"))
    logger.info("extract barcode done!")
    # Return all output file paths
    return (fqout_forward + "_1.fq.gz", fqout_forward + "_2.fq.gz", 
            fqout_reverse + "_1.fq.gz", fqout_reverse + "_2.fq.gz")

if __name__ == '__main__':
    barcode_main()
