#!/usr/bin/env python3

import os
import sys
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
import pysam
import pandas as pd
import click
from loguru import logger
from itertools import groupby
import numpy as np


def count_reads(bam: str, samplename: str, outdir: str, max_cells: int = 12000) -> list:
    """
    Count reads for each barcode and select top barcodes for analysis.
    
    Args:
        bam: Path to BAM file
        samplename: Sample name for output files
        outdir: Output directory
        max_cells: Maximum number of cells to extract (default: 12000)
    
    Returns:
        List of top barcodes based on read counts
    """
    try:
        bam_file = pysam.AlignmentFile(bam, "rb")
        barcode_counts = {}
        
        for read in bam_file:
            # Extract barcode from read name (format: barcode_other_info)
            barcode = read.query_name.split('_')[0]
            # Only count read2 to avoid double counting paired reads
            #if read.is_read1:
            #    continue
            barcode_counts[barcode] = barcode_counts.get(barcode, 0) + 1
        
        bam_file.close()
    except Exception as e:
        logger.error(f"Error reading BAM file {bam}: {e}")
        raise
    
    # Create DataFrame from barcode counts
    count_df = pd.DataFrame.from_dict(barcode_counts, orient='index')
    count_df.columns = ['reads_counts']
    count_df['barcode'] = count_df.index
    count_df = count_df.sort_values(by='reads_counts', ascending=False)
    
    # Save read counts to file
    os.makedirs(outdir, exist_ok=True)
    count_df.to_csv(f'{outdir}/{samplename}_reads_counts.txt', sep='\t')
    
    # Return top barcodes
    top_barcodes = count_df.head(n=max_cells).index.tolist()
    logger.info(f"Selected top {len(top_barcodes)} barcodes from {len(barcode_counts)} total barcodes")
    return top_barcodes

def get_barcodes_from_gexcb_and_cbcsv(gexcb: str, cbcsv: str = None) -> list:
    gexcb = pd.read_csv(gexcb, header = None, names = ['barcode'], sep = '\t')
    if cbcsv:
        cbcsv_map = pd.read_csv(cbcsv, header = 0, sep = ',')
        return cbcsv_map[cbcsv_map['gex_cb'].isin(gexcb['barcode'])]['m_cb'].tolist()
    else:
        return gexcb['barcode'].tolist()

def split_sortbyname_bam(
    bam: str,
    outdir: str,
    keep_barcodes: list,
    batch_id: int):
    """
    Split BAM file according to specified barcodes and count reads.
    
    Note: Input BAM file must be sorted by read name!
    
    Args:
        bam: Path to BAM file sorted by read name
        outdir: Output directory to save split BAM files
        keep_barcodes: List of barcodes to extract
        batch_id: Batch identifier for logging
        
    Returns:
        Dictionary with barcode read counts for this batch
    """
    os.makedirs(outdir, exist_ok=True)
    barcode_set = set(keep_barcodes)
    processed_barcodes = set()
    barcode_read_counts = {barcode: 0 for barcode in keep_barcodes}
    
    try:
        input_bam = pysam.AlignmentFile(bam, "rb")
        for barcode, reads_group in groupby(input_bam, key=lambda x: x.qname.split("_", 1)[0]):
            # Stop if we've processed all target barcodes
            if len(processed_barcodes) >= len(barcode_set):
                break
                
            if barcode in barcode_set and barcode not in processed_barcodes:
                output_path = f'{outdir}/{barcode}.bam'
                read_count = 0
                with pysam.AlignmentFile(output_path, 'wb', template=input_bam) as outfh:
                    for read in reads_group:
                        outfh.write(read)
                        # Only count read2 to avoid double counting paired reads
                        #if not read.is_read1:
                        read_count += 1
                
                barcode_read_counts[barcode] = read_count
                processed_barcodes.add(barcode)
                #logger.debug(f"Processed barcode {barcode} in batch {batch_id}")
        input_bam.close()
        #logger.info(f'Finished splitting BAM file for batch {batch_id}, processed {len(processed_barcodes)} barcodes, total mapped reads: {total_mapped_reads}')
        
        return barcode_read_counts
        
    except Exception as e:
        logger.error(f"Error splitting BAM file in batch {batch_id}: {e}")
        raise

@click.command()
@click.option('--bam', help='BAM file path (must be sorted by read name)', required=True)
@click.option('--samplename', help='Sample name for output files', required=True)
@click.option('--outdir', default='.', help='Output directory')
@click.option('--max_cells', 
              default=20000, 
              type=int, 
              help='Maximum number of cells to extract',
              show_default=True)
@click.option('--gexcb', 
              help='Path to RNA filtered barcode file (one barcode per line)',
              show_default=True)
@click.option('--cbcsv', 
              help='Path to bUCB3_whitelist.csv',
              show_default=True)
@click.option('--core', 
              default=1,
              type=int,
              help='Number of CPU cores to use for parallel processing')
def main(bam: str, samplename: str, outdir: str, max_cells: int = 20000, core: int = 1, gexcb: str = None, cbcsv: str = None):
    # Get barcodes either from counting or from provided file
    if not gexcb:
        logger.info("Counting reads and selecting top barcodes...")
        top_barcodes = count_reads(bam, samplename, outdir, max_cells)
    else:
        logger.info(f"Loading barcodes through file: {gexcb} and {cbcsv}")
        try:
            top_barcodes = get_barcodes_from_gexcb_and_cbcsv(gexcb, cbcsv)
        except FileNotFoundError:
            logger.error(f"Barcode file not found: {gexcb} or {cbcsv}")
            raise
        except Exception as e:
            logger.error(f"Error reading barcode file: {e}")
            raise
    
    # Sort barcodes and split into batches for parallel processing
    sorted_barcodes = sorted(top_barcodes)
    batches = np.array_split(sorted_barcodes, core)
    logger.info(f'Split {len(top_barcodes)} barcodes into {core} batches')
    
    # Create output directory
    split_bams_dirname = re.sub('_bismark_.*', '', os.path.basename(bam))
    split_output_dir = f'{outdir}/{split_bams_dirname}'
    os.makedirs(split_output_dir, exist_ok=True)
    
    # Process batches in parallel and collect read counts
    logger.info("Starting parallel BAM splitting...")
    all_barcode_counts = {}
    
    with ProcessPoolExecutor(max_workers=core) as executor:
        futures = []
        for i, barcode_batch in enumerate(batches):
            if len(barcode_batch) > 0:  # Skip empty batches
                futures.append(
                    executor.submit(split_sortbyname_bam, bam, split_output_dir, barcode_batch.tolist(), i)
                )
        
        # Wait for all tasks to complete with progress bar and collect results
        for future in as_completed(futures):
            try:
                batch_counts = future.result()  # Get the barcode counts from this batch
                all_barcode_counts.update(batch_counts)
            except Exception as e:
                logger.error(f"Error in batch processing: {e}")
                raise
    
    # Save filtered barcode read counts to table if using filtered_barcode
    if gexcb:
        count_df = pd.DataFrame.from_dict(all_barcode_counts, orient='index')
        count_df.columns = ['reads_counts']
        count_df['barcode'] = count_df.index
        count_df = count_df.sort_values(by='reads_counts', ascending=False)
        count_df = count_df[count_df['reads_counts'] > 0]
        
        # Save to file
        count_df.to_csv(f'{split_output_dir}/{split_bams_dirname}_filtered_barcode_reads_counts.csv', sep=',', index=False)
        count_df[['barcode']].to_csv(f'{split_output_dir}/{split_bams_dirname}_filtered_barcode', sep='\t', index=False, header=False)
        
        logger.info(f"Total filtered barcodes: {count_df.shape[0]}")
        logger.info(f"Total reads for filtered barcodes: {count_df['reads_counts'].sum()}")
    
    logger.info(f'Successfully finished splitting BAM file: {bam}')
if __name__ == '__main__':
    main()
    
    
    
