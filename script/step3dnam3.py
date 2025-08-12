import os
import re
import gzip
from collections import defaultdict
from itertools import groupby
import numpy as np
import pandas as pd
import pysam
from loguru import logger
import click


def hamming_distance(s1, s2):
    return len([(i, j) for i, j in zip(s1, s2) if i != j])

def umi_correct(geneid_umi_dict, bc, umi_correct_detail_fh):
    counts_dict = defaultdict(lambda: [0, 0])
    for gene_id, umi_dict in geneid_umi_dict.items():
        # 直接计数，不进行UMI矫正
        counts_dict[gene_id][0] = len(umi_dict.keys())  # UMI数量
        counts_dict[gene_id][1] = sum(umi_dict.values())  # 总reads数
    return counts_dict, geneid_umi_dict

def umi_count(samfile, reads_group, umi_correct_detail_fh):
    assigned_dict = defaultdict(lambda: defaultdict(int))
    reads_idlist = []
    for barcode_umi, g in groupby(reads_group, key=lambda x: x.qname):
        _, umi = barcode_umi.split("_")[:2]
        # quchu poly
        if umi == umi[0]*len(umi):
            continue
        tmp_dict = defaultdict(int)
        n = 0
        reads_idlist = []
        # print(reads_idlist)
        for r in g:
            chromosome = samfile.get_reference_name(r.reference_id)
            start_position = (r.reference_start + 1) // 100
            r_id = r.qname
            # print(f'--{n}--r_id:{r_id}')
            n += 1
            if r.has_tag("XT"):
                XT =r.get_tag("XT").split("XT:Z:")[-1]
                if "," in XT:
                    continue
                else:
                    gene_id = XT.split("XT:Z:")[-1]
            elif r.get_tag("XS") == "Unassigned_NoFeatures" or r.get_tag("XS") == "Unassigned_Overlapping_Length":
                gene_id = chromosome + "_" + str(start_position)
            else:
                print(r)
            if r_id in reads_idlist:
                continue
            else:
                reads_idlist.append(r_id)
                # 使用比对位置作为gene_id和umi
                umi =  chromosome + "_" + str(r.reference_start + 1)
                assigned_dict[gene_id][umi] += 1
            # assigned_dict[gene_id][umi] += 1
        
    counts_dict, geneid_umi_dict = umi_correct(assigned_dict, _, umi_correct_detail_fh)
    return counts_dict, geneid_umi_dict

def bam2table(bam, detail_file, counts_file, umi_correct_detail):
    sam_file = pysam.AlignmentFile(bam, "rb")

    umi_correct_detail_fh = open(umi_correct_detail, "w")
    with open(detail_file, "w") as fh1, open(counts_file, "w") as fh2:
        fh1.write("\t".join(["cellID", "geneID", "UMI", "Num"]) + "\n")
        fh2.write("\t".join(["cellID", "geneID", "UMINum", "ReadsNum"]) + "\n")
        for barcode, g in groupby(sam_file, key=lambda x: x.qname.split("_", 1)[0]):
            # print('-------1----------------------------------------')
            counts_dict, geneid_umi_dict = umi_count(sam_file, g, umi_correct_detail_fh)
            for gene_id in geneid_umi_dict:
                for umi in geneid_umi_dict[gene_id]:
                    raw_umi_count = geneid_umi_dict[gene_id][umi]
                    fh1.write(f"{barcode}\t{gene_id}\t{umi}\t{raw_umi_count}\n")
            for gene_id in counts_dict:
                umi_num, reads_num = counts_dict[gene_id]
                fh2.write(f"{barcode}\t{gene_id}\t{umi_num}\t{reads_num}\n")
    umi_correct_detail_fh.close()

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--outdir', required=True,help='outdir')
@click.option('--bam', required=True,help='input file. step2/target/sample_SortByName.bam.')

def count(bam, outdir):
    #basedir = os.path.join(outdir, "step3")
    #os.makedirs(basedir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    detail_file = os.path.join(outdir, "detail.xls")
    counts_file = os.path.join(outdir, "counts.xls")
    umi_file = os.path.join(outdir, "umi.xls")

    bam2table(bam=bam, detail_file=detail_file, counts_file=counts_file, umi_correct_detail=umi_file)
    print('输出完成')

if __name__ == '__main__':
    count()