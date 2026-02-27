#!/usr/bin/env python

import os
# Set OpenBLAS and OpenMP thread limits before importing numpy/pandas
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import sys
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import re
import logging
import click
import shutil
import subprocess

logger = logging.getLogger(__name__)
def run_allcools(args) -> str:
    '''
    Run ALLCools to extract CpG context, 
    and generate methylation reads counts and coverage counts for each CpG.
    
    bam: path to single cell barcode bam file
    genome_fa: path to reference fasta file
    barcode: barcode name
    outdir: path to save results
    chrom_size_path: path to a file records chromosome size
    align_method: alignment tools. If use bismark, will use --convert_bam_strandness parameter.
    tag: BAM tag to correction.
    '''
    bam, genome_fa, barcode, outdir, chrom_size_path, align_method, tag = args
    other_args = ''
    if align_method == 'bismark':
        other_args = '--convert_bam_strandness'
    if tag:
        other_args += f' --tag {tag}'
    samtools_sort_cmd = (f'samtools sort -o {outdir}/{barcode}_sort.bam {bam};')
    samtools_index_cmd = (f'samtools index {outdir}/{barcode}_sort.bam;')
    bam_to_allc_cmd = (
        f'allcools bam-to-allc '
        f'--bam_path {outdir}/{barcode}_sort.bam '
        f'--reference_fasta {genome_fa} '
        f'--output_path {outdir}/{barcode}_allc {other_args} '
        f'--save_count_df;'
    )
    try:
        subprocess.run(samtools_sort_cmd, check=True, shell = True)
        subprocess.run(samtools_index_cmd, check=True, shell = True)
        subprocess.run(bam_to_allc_cmd, check=True, shell = True)
        # Remove intermediate BAM files
        if os.path.exists(f'{outdir}/{barcode}_sort.bam'):
            os.remove(f'{outdir}/{barcode}_sort.bam')
        if os.path.exists(f'{outdir}/{barcode}_sort.bam.bai'):
            os.remove(f'{outdir}/{barcode}_sort.bam.bai')
        return f'{barcode} done'
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed processing {barcode}: {e}")
        raise e
   
@click.command()
@click.option('--indir', 
              help = 'Directory containing all single-cell barcode bam files', 
              required = True)
@click.option('--samplename', 
              help = 'Sample name', 
              required = True)
@click.option('--outdir', help = 'Output directory', default = './', show_default = True)
@click.option('--genomefa', help = 'Path to reference genome fasta file', required = True)
@click.option('--chrom_size_path', help = 'Path to chromosome size file', required = True)
@click.option('--align_method', 
              help = 'Alignment tool type, default is bismark', 
              default = 'bismark',
              show_default = True)
@click.option('--core', 
              help = 'core', 
              show_default=True, 
              default=1)
@click.option('--tag', 
              help = 'BAM tag to correction', 
              show_default=True, 
              default=None)
@click.option('--filtered_barcode', 
              help = 'filtered barcode list', 
              show_default=True, 
              default=None)
def main(
    indir: str, 
    samplename: str,
    outdir: str, 
    genomefa: str, 
    chrom_size_path: str, 
    align_method: str,
    core: int,
    tag: str = None,
    filtered_barcode:str = None):
    """
    Main function for methylation data analysis pipeline
    """
    samplename = re.sub(' ', '', samplename)
    if not core:
        n_cores = min(os.cpu_count() - 1, 30)
    else:
        n_cores = core
    if n_cores <= 0:
        ncores = 1
    # bam list
    bams = [f for f in os.listdir(indir) if f.endswith('.bam')]
    if filtered_barcode:
        filtered_barcode_list = []
        with open(filtered_barcode, 'r') as fh:
            for line in fh:
                filtered_barcode_list.append(line.strip())
        filtered_barcode_set = set(filtered_barcode_list)
        bams_new = []
        for i in bams:
            bam_name = re.sub('.bam', '', i)
            if bam_name in filtered_barcode_set: 
                bams_new.append(i)
        bams = bams_new
                
    logger.info('Extract CpG context')
    os.makedirs(f'{outdir}/{os.path.basename(indir)}_allcools', exist_ok = True)
    process_args = [
        (os.path.join(indir, bam),
        genomefa,
        re.sub('.bam', '', bam),
        f'{outdir}/{os.path.basename(indir)}_allcools',
        chrom_size_path,
        align_method, tag) 
        for bam in bams]
    logger.info('ALLCools running')
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        for col_idx, (info) in enumerate(
            tqdm(executor.map(run_allcools, process_args),
                total = len(bams),
                desc = "ALLCools processing barcodes")):
            pass
    logger.info('Extract CpG context finish')

if __name__ == '__main__':
    main()                
        