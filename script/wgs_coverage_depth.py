import re
import os
import click
import json
import numpy as np
from collections import defaultdict


def mapping_pairs(bismarkfile):
    with open(bismarkfile, "r") as bsmap_f:
        for line in bsmap_f:
            if "Sequence pairs analysed in total: " in line:
                totalpairs = int(re.findall("Sequence pairs analysed in total:\t([0-9]+) ",line)[0])
                # print(totalreads)
            if "Number of paired-end alignments with a unique best hit: " in line:
                uniquepairs = int(re.findall("Number of paired-end alignments with a unique best hit:\t([0-9]+) ",line)[0])
                # print(uniquereads)
            if "Sequence pairs did not map uniquely: " in line:
                not_uniq_alignedpairs = int(re.findall("Sequence pairs did not map uniquely:\t([0-9]+) ",line)[0])
                # print(alignedreads)
            alignedpairs = uniquepairs + not_uniq_alignedpairs
    return totalpairs, alignedpairs, uniquepairs

def mapping(bismarkfile):
    totalpairs = 0
    uniquepairs = 0
    not_uniq_alignedpairs = 0
    with open(bismarkfile, "r") as bsmap_f:
        for line in bsmap_f:
            if "Sequence pairs analysed in total:" in line:
                totalpairs = int(line.strip().split('\t')[-1])
            if "Number of paired-end alignments with a unique best hit:" in line:
                uniquepairs = int(line.strip().split('\t')[-1])
            if "Sequence pairs did not map uniquely:" in line:
                not_uniq_alignedpairs = int(line.strip().split('\t')[-1])
                
    alignedpairs = uniquepairs + not_uniq_alignedpairs
    totalreads = totalpairs * 2
    alignedreads = alignedpairs * 2
    uniquereads = uniquepairs * 2
    
    maprato = alignedreads / totalreads
    uniquerato = uniquereads / totalreads
    print(f'map : {maprato:.2%}')
    print(f'unique : {uniquerato:.2%}')

    return totalreads, alignedreads, uniquereads


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--outdir', required=True,help='outdir')
@click.option('--align_summary', required=True,help='bismark align summary')
@click.option('--samplename', required=True,help='sample name')
@click.option('--summary_json', required=True,help='summary json')
@click.option('--coveragefile', required=True,help='coveragefile')
def count(align_summary, outdir, samplename,summary_json,coveragefile):
    total_reads, aligned_reads, unique_reads = mapping(bismarkfile = align_summary)

    #coveragefile = os.path.join(outdir, "wgs", "genomecov", "bedtools_genomecoverage.txt")

    d = {}
    with open(coveragefile, "r") as fh:
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
    # print(d['chr1'])

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
    
    with open(summary_json, "r") as fh:
        summary = json.load(fh)

    with open(summary_json, "w") as fhout:
        summary_map = defaultdict()
        summary_map["Reads Mapped to Genome"] = aligned_reads / total_reads
        summary_map["Reads Mapped Confidently to Genome"] = unique_reads / total_reads
        summary_cov = defaultdict()
        summary_cov["Genome Coverage rate"] = genome_cov_1x / genome_all
        summary_cov["Genome 2X coverage rate"] = genome_cov_2x / genome_all
        summary_cov["Genome 4X coverage rate"] = genome_cov_4x / genome_all
        summary_cov["Genome 10X coverage rate"] = genome_cov_10x / genome_all
        # print(summary_map)
        # print(summary_cov)
        summary["mapping"] = summary_map
        summary["coverage"] = summary_cov
        json.dump(
            summary,
            fhout,
            indent=4,
            default=lambda o: int(o) if isinstance(o, np.int64) else o
        )

if __name__ == '__main__':
    count()