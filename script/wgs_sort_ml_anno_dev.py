import pandas as pd
import numpy as np
import gzip
import threading
from subprocess import run
from collections import defaultdict
import click
import json
import os


def mklist(cbfile):
    cb_list = []
    if cbfile.endswith(".gz"):
        _open = gzip.open
    else:
        _open = open
    with _open(cbfile, "rt") as f:
        for line in f:
            cb = line.strip()
            cb_list.append(cb)
    return cb_list

def countcov(gcovfile):
    d = {}
    with open(gcovfile, "r") as fh:
        for line in fh:
            chrname = line.strip().split('\t')[0]
            depth = int(line.strip().split('\t')[1])
            cov = int(line.strip().split('\t')[2])
            total = int(line.strip().split('\t')[3])
            d[chrname] = d.get(chrname, {'cov_1x':0, 'cov_2x':0, 'cov_4x':0, 'cov_10x':0, 'total':0})
            if depth >= 1:
                d[chrname]['cov_1x'] += cov
            if depth >= 2:
                d[chrname]['cov_2x'] += cov
            if depth >= 4:
                d[chrname]['cov_4x'] += cov
            if depth >= 10:
                d[chrname]['cov_10x'] += cov
            d[chrname]['total'] = total

    genome_cov_1x = 0
    genome_cov_2x = 0
    genome_cov_4x = 0
    genome_cov_10x = 0
    genome_all = 0
    for item in d:
        genome_cov_1x += d[item]['cov_1x']
        genome_cov_2x += d[item]['cov_2x']
        genome_cov_4x += d[item]['cov_4x']
        genome_cov_10x += d[item]['cov_10x']
        genome_all += d[item]['total']

    gcov_1x = genome_cov_1x / genome_all
    gcov_2x = genome_cov_2x / genome_all
    gcov_4x = genome_cov_4x / genome_all
    gcov_10x = genome_cov_10x / genome_all
    
    return gcov_1x, gcov_2x, gcov_4x, gcov_10x

def print_barcodes(cbfile, n, outdir, sample, genomefa, allcpath):
    cblist = mklist(cbfile=cbfile)
    df = pd.DataFrame(columns=['Barcode', 'genomecov_1x', 'genomecov_2x', 'genomecov_4x', 'genomecov_10x'])

    # 将列表分成n个子列表
    chunk_size = len(cblist) // n
    sublists = [cblist[i:i + chunk_size] for i in range(0, len(cblist), chunk_size)]

    # 定义一个线程任务，用于输出子列表中的元素
    def task(sublist, outdir, sample, genomefa):
        bedtools_path = 'bedtools'
        for barcode in sublist:
            # allc position 1-based
            cmd = (
                "mkdir -p {outdir}/bed;"
                "mkdir -p {outdir}/cb_cov;"
                "gunzip -dc {allcpath}/{barcode}_allc.gz | awk -v OFS='\t' '($4 ~/^CG/ ){{print $1,$2-1,$2}}' > {outdir}/bed/{barcode}_CpG.bed; "
                "{bedtools_path} genomecov -ibam {allcpath}/{barcode}_dedup_sort.bam > {outdir}/cb_cov/{sample}_{barcode}_genomecov.txt").format(
                outdir=outdir, sample=sample, barcode=barcode, genomefa=genomefa, bedtools_path=bedtools_path, allcpath = allcpath)
            run(cmd, shell=True)
    # 为每个子列表创建一个线程并启动它
    threads = []
    for sublist in sublists:
        print(len(sublist))
        t = threading.Thread(target=task, args=(sublist, outdir, sample, genomefa))
        t.start()
        threads.append(t)
    print('所有线程已启动。。。')

    # 等待所有线程完成
    for t in threads:
        t.join()
    print('所有线程已完成！')

    for barcode in cblist:
        gcovfile = f'{outdir}/cb_cov/{sample}_{barcode}_genomecov.txt'
        gcov_1x, gcov_2x, gcov_4x, gcov_10x = countcov(gcovfile=gcovfile)
        df = pd.concat([df, pd.DataFrame({'Barcode': [barcode], 'genomecov_1x': [gcov_1x], 'genomecov_2x': [gcov_2x], 'genomecov_4x': [gcov_4x], 'genomecov_10x': [gcov_10x]})], ignore_index=True)
    
    outcsv = os.path.join(outdir, sample+'_cb_genomecov.xls')
    df.to_csv(outcsv, sep='\t', index=False)
    print('所有细胞覆盖度统计已完成！')
    # return df

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--filtercb', required=True,help='step3/filtercb')
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option('--outdir', required=True,help='outdir')
@click.option('--samplename', required=True,help='samplename')
@click.option('--genomefa', required=True,help='/PROJ/k8s/DATA/nas-9b08bbd0-5dcd-4f3a-bdc8-6f0c22252535/database/refdata-gex-GRCh38-UCSC-hg38/GRCh38.d1.vd1.fa')
@click.option('--allcpath', required=True,help='allcools result path')
def main(filtercb, core, outdir, samplename, genomefa, allcpath):
    print_barcodes(cbfile=filtercb, n=core, outdir=outdir, sample=samplename, genomefa=genomefa, allcpath = allcpath)

if __name__ == '__main__':
    main()