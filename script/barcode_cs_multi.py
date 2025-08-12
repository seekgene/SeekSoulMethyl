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


R1_MINLEN = 20
R2_MINLEN = 60
ATAC_R1_MINLEN = 20
ATAC_R2_MINLEN = 60


CHEMISTRY = {
    "DD-M":{
        'shift': False,
        'structure': 'B17',
        'adapter1': [["AGATGTGTATAAGAGAYAG", "5", 0.1, 9],
                    ["CTGTCTCTTATACACATCT","3",0.1, 9]], # ME,ME-rev Y:C/T
        'adapter2': [["AGATGTGTATAAGAGACAG", "5", 0.1, 9],
                     ["CTRTCTCTTATACACATCT","3",0.1, 9]], # ME,ME-rev R:G/A
        'match_type': (1,),
    }    
}

class Reader(Process):
    """
    读取paired fastq
    """
    def __init__(self, file1:str, file2:str, connections:multiprocessing.connection, queue:Queue, buffer_size:int):
        """初始化读入进程

        Args:
            file1: read1 fastq file
            file2: read2 fastq file
            connections: 
            queue: 存放空闲工作进程index的队列
            buffer_size: 读取fastq时的buffer大小
        
        Returns:
            Reader对象

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
    """
    处理输出内容
    """
    def __init__(self, file:str, file_multi:str, paired_out:bool=False):
        """输出处理好的序列
        Args:
            file: read2 fastq file
            file_muti: 不能确定barcode的fastq文件名称
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
    """工作进程类
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
            self.write_pipe.send(-1)
        except Exception as e:
            self.write_pipe.send(-2)
            raise e

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
        # start reader process
        reader_connections = [Pipe(duplex=False) for _ in range(self.n_workers)]
        _pipes, _conn = zip(*reader_connections)
        _reader_process = Reader(self.fq1, self.fq2, _conn, self.need_work_queue, self.buffer_size)
        _reader_process.daemon = True
        _reader_process.start()

        # start worker processes
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
                chunk_index = connection.recv()
                if chunk_index == -1:
                    self.connections.remove(connection)
                    continue
                elif chunk_index == -2:
                    sys.stderr.write('err!!!\n')
                # if single?
                data1 = connection.recv_bytes()
                data2 = connection.recv_bytes()
                data_multi1 = connection.recv_bytes()
                data_multi2 = connection.recv_bytes()
                self.writer.write([(data1, data2), (data_multi1, data_multi2)], chunk_index)
                _stat = connection.recv()
                self.stat.update(**_stat)
        assert self.writer.wrote_everything()
        for w in self.workers:
            w.join()
        _reader_process.join()
        self.writer.close()

class AdapterFilter:
    """过滤接头"""
    def __init__(self, adapter1:list=[], adapter2:list=[],):
        self.adapter1 = []
        
        # 支持容错率和overlap参数
        # sequence: str,
        # max_errors: float = 0.1,
        # min_overlap: int = 16,
        for p in adapter1:
            if len(p) >= 4:  # 如果有额外的参数（容错率和overlap）
                if p[1] == "3":
                    self.adapter1.append(BackAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3]))
                elif p[1] == "5":
                    self.adapter1.append(FrontAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3]))
            else:  # 默认参数
                if p[1] == "3":
                    self.adapter1.append(BackAdapter(sequence=p[0], min_overlap=10))
                elif p[1] == "5":
                    self.adapter1.append(FrontAdapter(sequence=p[0], min_overlap=7))
        
        self.adapter2 = []
        for p in adapter2:
            if len(p) >= 4:  # 如果有额外的参数（容错率和overlap）
                if p[1] == "3":
                    self.adapter2.append(BackAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3]))
                elif p[1] == "5":
                    self.adapter2.append(FrontAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3]))
            else:  # 默认参数
                if p[1] == "3":
                    self.adapter2.append(BackAdapter(sequence=p[0], min_overlap=10))
                elif p[1] == "5":
                    self.adapter2.append(FrontAdapter(sequence=p[0], min_overlap=10))
    
    def filter(self, r1=None, r2=None) -> tuple:
        flag = False
        r1_me_left = False
        r1_me_right = False
        if r1 and self.adapter1:
            sp2_target = False  # 是否检测到目标接头标志
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
                    r1_start = 52 # cutadapter not rm me left, start from 52
                if r1_me_right:
                    r1_end = len(r1.sequence) - 9 # cutadapter cut me right, rm additial 9bp for uncorret methylation of Transposase
                if r1_end - r1_start > 81:
                    r1_end = r1_start + 81
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
    """汇总统计"""
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
        if not d:  # 如果字典为空，返回空的结果结构
            return {b: [] for b in 'ATCGN'}
        idx_max = max([k[0] for k in d])
        return {
            b: [d.get((i, b), 0) for i in range(idx_max+1)] for b in 'ATCGN'
        }

    @staticmethod
    def _sort_q(d, phred=33):
        if not d:  # 如果字典为空，返回空的结果结构
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
    """解析接头结构

    使用字母B、L、U、X和T以及数字表示reads结构。
    B表示barcode部分碱基；
    L表示linker部分碱基；
    U表示umi部分碱基；
    X表示任意碱基，用于占位；
    T表示T碱基；
    字母后数字表示碱基长度。

    Args:
        string: 接头结构描述

    Returns:
        返回二维tuple,内容为按顺序的各部分结构和长度。
        例如：
            当string是B8L8B8L10B8U8,返回:
            (('B', 8), ('L', 8), ('B', 8), ('L', 10), ('B', 8), ('U', 8))
    """
    regex = re.compile(r'([BLUXT])(\d+)')
    groups = regex.findall(string)
    return tuple([(_[0], int(_[1])) for _ in groups])

def read_file(file_list: list) -> dict:
    """准备白名单set
        Args:
            file_list: 每段白名单文件的路径
    """
    wl_dict = dict()
    for i, wl_file in enumerate(file_list):
        white_list = set()
        with xopen(wl_file, "r") as fh:
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
    """返回原始barcode各位置错配后的set与白名单set的交集"""

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
    # 假设使用 Phred+33 编码
    quality_scores = [ord(char) - 33 for char in quality_string]
    average_quality = sum(quality_scores) / len(quality_scores)
    return average_quality

def generate_random_string(length=12):
    return ''.join(random.choices(string.ascii_uppercase, k=length))

def summary(seq, seq_q, seq_dict, qua_dict):
    
    for i, (base, q) in enumerate(zip(seq, seq_q)):
        seq_dict[(i,base)] += 1
        qua_dict[(i,q)] += 1
        
    return seq_dict, qua_dict

def process_barcode(fq1, fq2, fq_out, fqout_multi, r1_structure, shift, shift_pattern,
                    barcode_wl_dict, linker_wl_dict, match_type_dict, adapter1=[["AAAAAAAAAAAA", "3"],],
                    adapter2=[["AAAAAAAAAAAA", "3"],], do_B_correction=True, do_L_correction=True,
                    use_multi=True, use_short_read=False, paired_out=True):
    
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
    
    adapter_filter = AdapterFilter(adapter1=adapter1, adapter2=adapter2)
    
    fh = dnaio.open(fq1, fq2, fileformat="fastq", mode="r")
    if paired_out:
        outfh = dnaio.open(fq_out[0], fq_out[1], fileformat="fastq", mode="w")
    else:
        outfh = dnaio.open(fq_out[0], fileformat="fastq", mode="w")

    if use_multi:
        if paired_out:
            outfh_multi = dnaio.open(fqout_multi[0], fqout_multi[1], fileformat="fastq", mode="w")
        else:
            outfh_multi = dnaio.open(fqout_multi[0], fileformat="fastq", mode="w")
    
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

        # check double instances
        if is_valid:
            if len(r1.sequence) >= 70:
                seq_7f = r1.sequence[17:24]
                l17me = r1.sequence[39:57]
                if l17me == 'GTAGATGTGTATAAGAGA':
                    stat_Dict["num_17lme"] += 1
                
                if seq_7f == 'TTGCTGT' or seq_7f == 'TTGTTGT':
                    stat_Dict["num_7f"] += 1
                    if l17me == 'GTAGATGTGTATAAGAGA':
                        stat_Dict["num_7f17lme"] += 1

                        pos1 = r1.sequence[24]
                        qua1 = r1.qualities[24]
                        pos2 = r1.sequence[27]
                        qua2 = r1.qualities[27]
                        pos3 = r1.sequence[28]
                        qua3 = r1.qualities[28]
                        pos4 = r1.sequence[31]
                        qua4 = r1.qualities[31]
                        pos5 = r1.sequence[36]
                        qua5 = r1.qualities[36]
                        pos6 = r1.sequence[38]
                        qua6 = r1.qualities[38]
                        pos7 = r1.sequence[57]
                        qua7 = r1.qualities[57]
                        allseq = pos1 + pos2 + pos3 + pos4 + pos5 + pos6 + pos7
                        countseq = allseq[0:5]
                        insert = pos6 + pos7
                        allqua = qua1 + qua2 + qua3 + qua4 + qua5 + qua6 + qua7
                        average_quality = calculate_average_quality(allqua)
                        if average_quality >= 30:
                            if insert == 'CC': 
                                B_count_C = countseq.count('C')
                                B_count_T = countseq.count('T')
                                if B_count_C + B_count_T == len(countseq):
                                    stat_Dict["line_B"] += 1
                                    stat_Dict["B_all_C"] += B_count_C
                            elif insert == 'TT':
                                A_count_T = countseq.count('T')
                                A_count_C = countseq.count('C')
                                if A_count_T + A_count_C == len(countseq): # 必需是C或者T，如果是其他碱基说明是测错了
                                    stat_Dict["line_A"] += 1
                                    stat_Dict["A_all_T"] += A_count_T
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
                        
                    bc_new_all = ":".join(bc_new_lst)
                    r2.name = "_".join([barcode_old, bc_new_all, umi, r2.name])
                    r1.name = "_".join([barcode_old, bc_new_all, umi, r1.name])
                    r1.sequence = sequence[start_pos:]
                    r1.qualities = qualities[start_pos:]
                    outfh_multi.write(r1, r2)
            else:  #write r2 files
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

                umi12 = generate_random_string()
                umi19 = umi+umi12

                r2.name = "_".join([barcode_new, umi19, _alt, r2.name])
                r1.name = "_".join([barcode_new, umi19, _alt, r1.name])

                # r2.name = "_".join([barcode_new, umi, _alt, r2.name])
                # r1.name = "_".join([barcode_new, umi, _alt, r1.name])
                outfh.write(r1, r2)
            if is_correct:
                stat_Dict["B_corrected"] += 1
        else:
            if is_B_no_correction:
                stat_Dict["B_no_correction"] += 1

            if is_L_no_correction:
                stat_Dict["L_no_correction"] += 1
    if use_multi:
        outfh_multi.close()
    outfh.close()

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
@click.option("--chemistry", required=True, type=click.Choice(["DD-Q","DD_AG","DD-M"]), help="chemistry")

def barcode_main(chemistry, fq1:list, fq2:list, samplename: str, outdir:str,
                 barcode:list=[], match_type:list=[], shift:str=True, shift_pattern:str="A",
                 structure:str="B17U12X7", linker: list=[],
                 core:int=4, do_B_correction=True, do_L_correction=True,
                 use_multi=True, use_short_read=False, adapter1=[["TTTTTTTTTTTT", "5"], ],
                 adapter2=[["AAAAAAAAAAAA", "3"], ], paired_out=True):
    shift = CHEMISTRY[chemistry]['shift']
    structure = CHEMISTRY[chemistry]['structure']
    adapter1 = CHEMISTRY[chemistry]['adapter1']
    adapter2 = CHEMISTRY[chemistry]['adapter2']
    match_type = CHEMISTRY[chemistry]['match_type']

    logger.info("extract barcode start!")
    #parse r1 structure
    r1_structure = parse_structure(structure)
  
    #get wl dict for bc/linker
    # 当barcode或linker为空时，返回空字典
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
        paired_out=paired_out
    )
    pipeline.run()

    fqout1 = f"{fqout}_1.fq.gz"
    fqout2 = f"{fqout}_2.fq.gz"
    if use_multi:
        # find the multiple barcodes
        logger.info("deal multi start!")
        fqout_multi1 = f"{fqout_multi}_1.fq.gz"
        fqout_multi2 = f"{fqout_multi}_2.fq.gz"
        adapter_filter = AdapterFilter(adapter1=adapter1, adapter2=adapter2)
        multi_stat = defaultdict(int)
        with dnaio.open(fqout1, fqout2, mode="a") as f:
            fh = dnaio.open(fqout_multi1, fqout_multi2, fileformat="fastq", mode="r")
            for r1, r2 in fh:
                multi_stat["total"] += 1
                final_barcode = None
                
                bc_old, r2_candidate, umi, r2_name = r2.name.split("_", 3)
                r2_candidate = r2_candidate.split(":")
                
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
                # stat.data["barcode_count"][_] += 1

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

                umi12 = generate_random_string()
                umi19 = umi+umi12

                r2.name = "_".join([final_barcode, umi19, _alt, r2_name])
                r1.name = r2.name.replace(' 2:',' 1:')
                f.write(r1, r2)

        with open(json_multi, "w") as fh:
            json.dump(multi_stat, fp=fh, indent=4)

    del stat.data["barcode_count"]
    logger.info("deal multi done!")
    stat.data["stat"]["chemistry"] = chemistry
    stat.data["stat"]["gexname"] = samplename
    # ct_mean = (A_all_T + B_all_C) / (7 * line_A + 7 * line_B)
    #stat.data["stat"]["ct_mean"] = (stat.data["stat"]["A_all_T"] + stat.data["stat"]["B_all_C"]) / (7 * stat.data["stat"]["line_A"] + 7 * stat.data["stat"]["line_B"])
    stat.data["stat"]["ct_mean"] = stat.data["stat"]["A_all_T"]  / (5 * stat.data["stat"]["line_A"])
    stat.data["stat"]["cc_mean"] = stat.data["stat"]["B_all_C"]  / (5 * stat.data["stat"]["line_B"])
    # rate_7f = num_7f / total
    stat.data["stat"]["rate_7f"] = stat.data["stat"]["num_7f"] / stat.data["stat"]["total"]
    # rate_17lme = num_17lme / total
    stat.data["stat"]["rate_17lme"] = stat.data["stat"]["num_17lme"] / stat.data["stat"]["total"]
    # rate_7f17lme = num_7f17lme / total
    stat.data["stat"]["rate_7f17lme"] = stat.data["stat"]["num_7f17lme"] / stat.data["stat"]["total"]
    stat.save(os.path.join(outdir, f"{samplename}_summary.json"))
    logger.info("extract barcode done!")
    return fqout1, fqout2

if __name__ == '__main__':
    barcode_main()
