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

@click.command()
@click.option('--fq1', help = 'R1 fastq file path.', required=True)
@click.option('--fq2', help = 'R2 fastq file path.', required=True)
@click.option('--samplename', help = 'samplename.', required=True)
@click.option('--outdir', default = '.', help = 'Output directory.')
def rename_fastq(fq1:str, fq2:str, samplename:str, outdir:str) -> str:
    '''
    Add raw UMI to the end of read name
    fq1: path to r1 fastq
    fq2: path to r2 fastq
    samplename: prefix for output fastq
    outdir: path to save edited fastq files.
    '''
    logger.info('Add raw UMI to the end of read name of fastq')
    outfq1 = os.path.join(
        outdir, 
        f'{samplename}_1_rename.fq.gz'
    )
    outfq2 = os.path.join(
        outdir, 
        f'{samplename}_2_rename.fq.gz'
    )

    with dnaio.open(file1 = outfq1, file2 = outfq2, mode = 'w') as writer:
        with dnaio.open(file1 = fq1, file2 = fq2, mode = 'r') as reader:
            for r1, r2 in tqdm(reader):
                read_names_1 = r1.name.split(' ')
                #barcode, umi, map_status, platform_info = read_names_1[0].split('_')
                barcode = read_names_1[0].split('_')[0]
                new_r1_read_name = f'{read_names_1[0]}:{barcode} {read_names_1[1]}:{barcode}' # use barcode
                
                read_names_2 = r2.name.split(' ')
                new_r2_read_name = f'{read_names_2[0]}:{barcode} {read_names_2[1]}:{barcode}'
                
                r1.name = new_r1_read_name
                r2.name = new_r2_read_name
                writer.write(r1, r2)
    logger.info('Add raw UMI finish')

if __name__ == '__main__':
    rename_fastq()