import os
import sys
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import re
import dnaio
import pysam
import pandas as pd
import click
from loguru import logger
from collections import defaultdict
from itertools import groupby


def count_reads(bam:str, samplename:str, outdir:str, max_cells:int = 12000) -> list:
    '''
    1. Count reads for each barcode
    2. Select reads number top 12000 barcodes to analysis
    '''
    bam_file = pysam.AlignmentFile(bam, "rb")
    barcode_counts = {}
    for read in bam_file:
        barcode = read.query_name.split('_')[0]
        if read.is_read1:
            continue
        if barcode in barcode_counts:
            barcode_counts[barcode] += 1
        else:
            barcode_counts[barcode] = 1
    bam_file.close()
    
    count_df = pd.DataFrame.from_dict(barcode_counts, orient = 'index')
    count_df.columns = ['reads_counts']
    count_df['barcode'] = count_df.index
    count_df = count_df.sort_values(by = 'reads_counts', ascending = False)
    count_df.to_csv(f'{outdir}/{samplename}_reads_counts.txt', sep = '\t')
    top_barcodes = count_df.head(n = max_cells).index
    return top_barcodes

def split_sortbyname_bam(
    bam:str,
    outdir:str,
    keep_barcodes:list):
    '''
    Split bam according to keep barcodes
    input bam must be sorted by read name!
    bam: path to bam file that sorted by name
    outdir: path to save split bams
    keep_barcodes: a list include all keep barcodes.
    '''
    os.makedirs(outdir, exist_ok = True)
    barcode_set = set(keep_barcodes)
    
    barcode_done = 0
    input_bam = pysam.AlignmentFile(bam, "rb")
    pbar = tqdm(total = len(barcode_set))
    for barcode, g in groupby(input_bam, key = lambda x: x.qname.split("_", 1)[0]):
        if barcode in barcode_set:
            with pysam.AlignmentFile(
                f'{outdir}/{barcode}.bam', 
                'wb', 
                template = input_bam) as outfh:
                for r in g:
                    outfh.write(r)
            pbar.update(1)
    input_bam.close()
    pbar.close()
    
@click.command()
@click.option('--bam', help = 'BAM file path.', required=True)
@click.option('--samplename', help = 'Samplename.', required=True)
@click.option('--outdir', default = '.', help = 'Output directory.')
@click.option('--max_cells', 
              default = 12000, 
              type = int, 
              help = 'The number of cells to extract.',
              show_default = True)
@click.option('--filtered_barcode', 
              help = 'Filtered barcode file path. Each line is a barcode.')
def main(bam:str, samplename: str, outdir:str, max_cells = 12000, filtered_barcode = None):
    if not filtered_barcode:
        top_barcodes = count_reads(bam, samplename, outdir, max_cells)
    else:
        with open(filtered_barcode, 'r') as f:
            top_barcodes = [line.strip() for line in f.readlines()]
    os.makedirs(f'{outdir}/split_bams', exist_ok = True)
    split_sortbyname_bam(bam, f'{outdir}/split_bams', top_barcodes)
    
if __name__ == '__main__':
    main()
    
    
    
