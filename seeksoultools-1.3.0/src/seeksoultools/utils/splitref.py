import os
from collections import defaultdict
import click
from intervaltree import IntervalTree
from loguru import logger

def fa_reader(fa):
    """read fa"""
    chr = ""
    seq = ""
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    yield chr, seq
                chr = line.strip().split(">")[1].split(" ")[0]
                seq = ""
            else:
                seq += line.strip()
        if seq:
            yield chr, seq

def filter_long_chr(fa, max_len=2**29):
    """scan fa"""
    chrs_too_long_dict = {}
    reader = fa_reader(fa)
    for chr, seq in reader:
        chr_len = len(seq)
        if chr_len > max_len:
            chrs_too_long_dict[chr] = chr_len
    return chrs_too_long_dict

def find_split_point(chr, tree_dict, start=0, chunk=2**29, step=1000):
    """find split point"""
    n = 0
    while tree_dict[chr].at(chunk+start):
        chunk -= step
        n += 1
        if n>10000:
            raise Exception(f"no good split point for {chr} in 10Mb base.")
    return chunk+start

def loop_whole_chr(chr, chr_len, tree_dict, chunk=2**29, step=1000):
    """loop whole chr"""
    split_point = []
    start = 0
    while chr_len > chunk:
        point = find_split_point(chr, tree_dict, start, chunk, step)
        split_point.append(point)
        chr_len -= point
        start = start + point
    return split_point

def scan_gtf(gtf, chrs_too_long_dict, chunk=2**29, step=1000):
    """scan gtf and pick split point"""
    tree_dict = defaultdict(IntervalTree)
    split_point_dict = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.strip() == "":
                continue
            tmp = line.strip().split("\t")
            if tmp[0] in chrs_too_long_dict:
                tree_dict[tmp[0]].addi(int(tmp[3]), int(tmp[4]), 1)

    for chr, chr_len in chrs_too_long_dict.items():
        split_point = loop_whole_chr(chr, chr_len, tree_dict, chunk, step)
        split_point_dict[chr] = split_point
    return split_point_dict

def new_fa(fa, split_point_dict, outdir):
    filename = os.path.join(outdir, os.path.basename(fa))
    assert os.path.abspath(filename)!=os.path.abspath(fa), f"output file should not be same as input file: {os.path.abspath(filename)}"
    with open(filename, "w") as fh:
        reader = fa_reader(fa)
        for chr, seq in reader:
            if chr in split_point_dict:
                start = 0
                for n, p in enumerate(split_point_dict[chr]):
                    tmp_seq = seq[start:p]
                    fh.write(f">{chr}_s{n} {len(tmp_seq)}\n")
                    fh.write(tmp_seq)
                    fh.write("\n")
                    start = p
                tmp_seq = seq[p:]
                if tmp_seq:
                    fh.write(f">{chr}_s{n+1} {len(tmp_seq)}\n")
                    fh.write(tmp_seq)
                    fh.write("\n")
            else:
                fh.write(f">{chr}\n")
                fh.write("\n")
    logger.info(f"new fa: {filename}")

def new_gtf(gtf, split_point_dict, outdir):
    filename = os.path.join(outdir, os.path.basename(gtf))
    with open(filename, "w") as fh, \
            open(gtf, "r") as fh_gtf:
        for line in fh_gtf:
            if line.startswith("#"):
                fh.write(line)
                continue
            if line.strip() == "":
                fh.write(line)
                continue
            tmp = line.strip().split("\t")
            tmp[3] = int(tmp[3])
            tmp[4] = int(tmp[4])
            if tmp[0] in split_point_dict:
                tmp_end = [tmp[4] - p for p in split_point_dict[tmp[0]] if tmp[4]-p>=0]
                tmp[0] = f"{tmp[0]}_s{len(tmp_end)}"
                if tmp_end:
                    tmp_len = tmp[4] - tmp[3]
                    tmp[4] = tmp_end[-1]
                    tmp[3] = tmp[4] - tmp_len
                new_line = "\t".join([str(_) for _ in tmp])
                fh.write(f"{new_line}\n")
            else:
                fh.write(line)
    logger.info(f"new gtf: {filename}")

# @click.command()
# @click.option("-f", "--fa", required=True)
# @click.option("-g", "--gtf", required=True)
# @click.option("-o", "--outdir", required=True)
def split_ref(fa, gtf, outdir):
    assert os.path.abspath(outdir)!=os.path.dirname(os.path.abspath(fa)), f"output dir should not be same as input file for fa: {os.path.abspath(outdir)}"
    assert os.path.abspath(outdir)!=os.path.dirname(os.path.abspath(gtf)), f"output dir should not be same as input file for gtf: {os.path.abspath(outdir)}"

    chrs_too_long_dict = filter_long_chr(fa)
    logger.info(f"chr with long length: {chrs_too_long_dict}")

    split_point_dict = scan_gtf(gtf, chrs_too_long_dict)
    logger.info(f"split point: {split_point_dict}")

    logger.info(f"split gtf.")
    new_gtf(gtf, split_point_dict, outdir)

    logger.info(f"split fa.")
    new_fa(fa, split_point_dict, outdir)
