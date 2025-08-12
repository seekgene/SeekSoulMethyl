import os
# Set OpenBLAS and OpenMP thread limits before importing numpy/pandas
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
# 设置环境变量以解决scikit-learn和scipy兼容性问题
os.environ['SCIPY_ARRAY_API'] = '1'

import sys
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import re
import pandas as pd
import logging
import click
import matplotlib.pyplot as plt
from ALLCools.mcds import MCDS
from ALLCools.plot import *
import seaborn as sns
import anndata
from scipy.io import mmwrite
import scipy.sparse
import shutil
import subprocess
import gzip

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
    '''
    bam, genome_fa, barcode, outdir, chrom_size_path, align_method = args
    other_args = ''
    if align_method == 'bismark':
        other_args = '--convert_bam_strandness'
    samtools_sort_cmd = (f'samtools sort -o {outdir}/{barcode}_dedup_sort.bam {bam};')
    bam_to_allc_cmd = (
        f'allcools bam-to-allc '
        f'--bam_path {outdir}/{barcode}_dedup_sort.bam '
        f'--reference_fasta {genome_fa} '
        f'--output_path {outdir}/{barcode}_allc {other_args}; '
    )
    try:
        subprocess.run(samtools_sort_cmd, check=True, shell = True)
        subprocess.run(bam_to_allc_cmd, check=True, shell = True)
        return f'{barcode} done'
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed processing {barcode}: {e}")
        raise e

def run_methyldackel(args) -> str:
    '''
    Run MethylDackel to extract CpG context, 
    and generate methylation reads counts and unmethylation reads counts for each CpG.
    
    bam: path to single cell barcode bam file
    genome_fa: path to reference fasta file
    barcode: barcode name
    outdir: path to save results
    '''
    bam, genome_fa, barcode, outdir = args
    
    methyldackel_cmd = (f'MethylDackel extract '
           f'--keepSingleton '
           f'-o {outdir}/{barcode} '
           f'-@ 1 '
           f'{genome_fa} {bam}'
           )
    subprocess.run(methyldackel_cmd, check = True)
    return f'{barcode} done'

def allc_generate_datasets(
    indir: str, 
    outdir: str, 
    chrom_size_path: str, 
    samplename: str, 
    ncores:int = None, 
    bed:list = None, 
    region_name:list = None,
    out_sparse:bool = False):
    '''
    Run ALLCools to generate methylation matrix
    
    indir: directory contains all single cell allcools bam-to-allc results.
    outdir: path to save results
    chrom_size_path: path to a file records chromosome size
    samplename: sample name
    ncores: cores
    bed: bed file contain region
    region_name: specify a name to bed region
    out_sparse: wether output sparse matrix
    '''
    os.makedirs(outdir, exist_ok = True)
    with open(f'{outdir}/allc_file_paths.tsv', 'w') as outfh:
        for i in os.listdir(indir):
            if i.endswith('_allc.gz'):
                outfh.write(f'{re.sub("_allc.gz", "", i)}\t{indir}/{i}\n')
    
    if ncores == None:
        ncores = max(os.cpu_count() - 1, 1)
    other_args = ' '
    region_list = ['chrom1M', 'chrom100k', 'chrom500k']
    if bed and len(bed) > 0 and region_name and len(region_name) > 0:
        for bed_path, name in zip(bed, region_name):
            region_name_tmp = re.sub(' ', '_', name)
            other_args = other_args + (
                f'--regions {region_name_tmp} {bed_path} '
                f'--quantifiers {region_name_tmp} hypo-score CGN cutoff=0.9 '
                f'--quantifiers {region_name_tmp} count CGN '
            ) 
            region_list.append(region_name_tmp)
    generate_cmd = (
        f'allcools generate-dataset --allc_table {outdir}/allc_file_paths.tsv '
        f'--output_path {outdir}/{samplename}.mcds '
        f'--chrom_size_path {chrom_size_path} '
        f'--obs_dim cell '
        f'--cpu {ncores} '
        f'--chunk_size 50 '
        f'--regions chrom1M 1000000 '
        f'--regions chrom100k 100000 '
        f'--regions chrom500k 500000 '
        f'--quantifiers chrom1M count CGN '
        f'--quantifiers chrom100k count CGN '
        f'--quantifiers chrom500k count CGN '
        f'--quantifiers chrom1M hypo-score CGN cutoff=0.9 '
        f'--quantifiers chrom100k hypo-score CGN cutoff=0.9 '
        f'--quantifiers chrom500k hypo-score CGN cutoff=0.9 '
        f'{other_args};'
    )
    logger.info(generate_cmd)
    subprocess.run(generate_cmd, check = True, shell = True)
    '''
    if out_sparse:
        for var_dim in region_list:
            out_sparse_matrix(
                f'{outdir}/{samplename}.mcds', 
                f'{outdir}/sparse_matrix/{var_dim}',
                var_dim
            )
    '''

def out_sparse_matrix(mcds_path:str, outdir:str, var_dim:str = 'chrom100k'):
    obs_dim = 'cell'  # observation
    # load to memory or not
    load = True
    mcds = MCDS.open(
    mcds_path, 
    obs_dim = obs_dim, 
    var_dim=var_dim
    )
    total_feature = mcds.get_index(var_dim).size
    mcds.add_feature_cov_mean(var_dim=var_dim)
    # get mc reads count and cov reads count
    var_dim_mc = mcds[f'{var_dim}_da'].sel({  
        'count_type': 'mc'
    })
    var_dim_cov = mcds[f'{var_dim}_da'].sel({
        'count_type': 'cov'
    })
    if load and (mcds.get_index('cell').size <= 20000):
        var_dim_mc.load()
        var_dim_cov.load()
        
    barcodes_df = pd.DataFrame(mcds.cell.values)
    features_df = pd.DataFrame(mcds[f'{var_dim}'].values)
    os.makedirs(f'{outdir}/coverage_sparse_mat/', exist_ok = True)
    with gzip.open(f'{outdir}/coverage_sparse_mat/matrix.mtx.gz', 'w') as mat_out:
        mmwrite(
            mat_out, 
            scipy.sparse.coo_matrix(var_dim_cov.values), 
            field = 'integer')  
    barcodes_df.to_csv(
        f'{outdir}/coverage_sparse_mat/barcodes.tsv.gz', 
        sep = '\t', 
        index = False, 
        header = False,
        compression='gzip'
    )
    features_df_name = []
    for n in range(len(mcds[f'{var_dim}'].values)):
        _chr = mcds[f'{var_dim}'][f'{var_dim}_chrom'].values[n]
        _start = mcds[f'{var_dim}'][f'{var_dim}_start'].values[n]
        _end = mcds[f'{var_dim}'][f'{var_dim}_end'].values[n]
        _name = f'{_chr}_{_start}_{_end}'
        features_df_name.append(_name)
    features_df['name'] = features_df_name
    features_df['Type'] = 'Methy Expression'
    features_df.to_csv(
        f'{outdir}/coverage_sparse_mat/features.tsv.gz', 
        sep = '\t', 
        index = False, 
        header = False,
        compression='gzip'
    )
    
    os.makedirs(f'{outdir}/methy_sparse_mat/', exist_ok = True)
    with gzip.open(f'{outdir}/methy_sparse_mat/matrix.mtx.gz', 'w') as mat_out:
        mmwrite(
            mat_out, 
            scipy.sparse.coo_matrix(var_dim_mc.values), 
            field = 'integer'
        ) 
    barcodes_df.to_csv(
        f'{outdir}/methy_sparse_mat/barcodes.tsv.gz', 
        sep = '\t', 
        index = False, 
        header = False,
        compression='gzip'
    )
    features_df['name'] = features_df_name
    features_df['Type'] = 'Methy Expression'
    features_df.to_csv(
        f'{outdir}/methy_sparse_mat/features.tsv.gz', 
        sep = '\t', 
        index = False, 
        header = False,
        compression='gzip'
    )
    
@click.command()
@click.option('--indir', 
              help = 'Directory containing all single-cell barcode bam files', 
              required = True)
@click.option('--samplename', 
              help = 'Sample name', 
              required = True)
@click.option('--outdir', help = 'Output directory', default = './', show_default = True)
@click.option('--bed', 
              nargs=2, 
              type=click.Tuple([str, str]),
              multiple = True,
              help = '''
              First item is the bed name to save results.
              Second item is the path to bed file containing regions of interest. 
              ''', 
              default = None,
              show_default = True)
@click.option('--genomefa', help = 'Path to reference genome fasta file', required = True)
@click.option('--chrom_size_path', help = 'Path to chromosome size file', required = True)
@click.option('--run_type', 
              help = 'Run type: methydackel, allcools, or both. Select tool to extract CpG context.', 
              type = click.Choice(['methydackel', 'allcools', 'both']), 
              default = 'allcools',
              show_default = True)
@click.option('--align_method', 
              help = 'Alignment tool type, default is bismark', 
              default = 'bismark',
              show_default = True)
@click.option('--out_sparse', 
              help = 'Wether save sparse matrix', 
              is_flag = True,
              show_default=True, 
              default=False)
@click.option('--filtered_barcode', 
              help = 'filtered barcode list', 
              show_default=True, 
              default=None)
def main(
    indir: str, 
    samplename: str,
    outdir: str, 
    bed: list, 
    genomefa: str, 
    chrom_size_path: str, 
    run_type: str, 
    align_method: str,
    out_sparse: bool,
    filtered_barcode:str = None):
    """
    Main function for methylation data analysis pipeline
    """
    samplename = re.sub(' ', '', samplename)
    
    n_cores = min(os.cpu_count() - 1, 30)
    if n_cores <= 0:
        ncores = 1
    
    # parse bed parameter
    bed_paths = []
    region_names = []
    if bed and len(bed) > 0:
        for name, path in bed:
            region_names.append(name)
            bed_paths.append(path)
    # bam list
    bams = [f for f in os.listdir(indir) if f.endswith('deduplicated.bam')]
    if len(bams) == 0:
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
    if run_type == 'methyldackel':
        os.makedirs(f'{outdir}/methyldackel', exist_ok = True)
        process_args = [
            (os.path.join(indir, bam),
            genomefa,
            re.sub('.bam', '', bam),
            f'{outdir}/methyldackel') 
            for bam in bams]
        logger.info('MethylDackel running')
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            for col_idx, (info) in enumerate(
                tqdm(executor.map(run_methyldackel, process_args),
                    total = len(bams),
                    desc = "MethylDackel processing barcodes")):
                pass
    elif run_type == 'allcools':
        os.makedirs(f'{outdir}/allcools', exist_ok = True)
        process_args = [
            (os.path.join(indir, bam),
            genomefa,
            re.sub('.bam', '', bam),
            f'{outdir}/allcools',
            chrom_size_path,
            align_method) 
            for bam in bams]
        logger.info('ALLCools running')
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            for col_idx, (info) in enumerate(
                tqdm(executor.map(run_allcools, process_args),
                    total = len(bams),
                    desc = "ALLCools processing barcodes")):
                pass
        allc_generate_datasets(
            indir = f'{outdir}/allcools',
            outdir = f'{outdir}/allcools_generate_datasets',
            samplename = samplename,
            chrom_size_path = chrom_size_path,
            ncores = n_cores,
            bed = bed_paths, 
            region_name = region_names
        )
    else:
        os.makedirs(f'{outdir}/methydackel', exist_ok = True)
        process_args = [
            (os.path.join(indir, bam),
            genomefa,
            re.sub('.bam', '', bam),
            f'{outdir}/methydackel'
            ) 
            for bam in bams]
        logger.info('MethylDackel running')
        with ProcessPoolExecutor(max_workers = n_cores) as executor:
            for col_idx, (info) in enumerate(
                tqdm(executor.map(run_methyldackel, process_args),
                    total = len(bams),
                    desc="MethylDackel processing barcodes")):
                pass
        
        os.makedirs(f'{outdir}/allcools', exist_ok = True)
        process_args = [
            (os.path.join(indir, bam),
            genomefa,
            re.sub('.bam', '', bam),
            f'{outdir}/allcools',
            chrom_size_path,
            align_method) 
            for bam in bams]
        logger.info('ALLCools running')
        with ProcessPoolExecutor(max_workers = n_cores) as executor:
            for col_idx, (info) in enumerate(
                tqdm(executor.map(run_allcools, process_args),
                    total = len(bams),
                    desc = "Processing sample")):
                pass
    logger.info('Extract CpG context finish')
    
    

if __name__ == '__main__':
    main()                
        
        
        
    
    