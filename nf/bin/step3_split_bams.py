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

def split_bam_single_pass(
    bam: str,
    outdir: str,
    keep_barcodes: list,
    max_open_files: int = 256):
    """
    Split BAM file in a single sequential pass using groupby.

    Since the input BAM is sorted by read name (samtools sort -n),
    reads with the same barcode are contiguous. We can traverse
    the BAM once, opening output files only when needed, and
    closing them when that barcode's reads are exhausted.

    This avoids the previous approach where N workers each traversed
    the entire BAM from start (O(N) redundant I/O).

    Args:
        bam: Path to BAM file sorted by read name
        outdir: Output directory to save split BAM files
        keep_barcodes: List of barcodes to extract
        max_open_files: Maximum number of simultaneously open output files

    Returns:
        Dictionary with barcode read counts
    """
    os.makedirs(outdir, exist_ok=True)
    barcode_set = set(keep_barcodes)
    barcode_read_counts = {barcode: 0 for barcode in keep_barcodes}

    # Track open file handles, close old ones when pool exceeds max_open_files
    open_handles = {}  # barcode -> AlignmentFile
    handle_order = []  # FIFO order of opened barcodes for LRU closing

    try:
        input_bam = pysam.AlignmentFile(bam, "rb")
        template = input_bam

        processed_barcodes = set()

        for barcode, reads_group in groupby(input_bam, key=lambda x: x.qname.split("_", 1)[0]):
            # Skip barcodes not in our target set
            if barcode not in barcode_set:
                continue

            # Stop if we've processed all target barcodes
            if len(processed_barcodes) >= len(barcode_set):
                break

            # Get or create output file handle
            if barcode not in open_handles:
                output_path = f'{outdir}/{barcode}.bam'
                open_handles[barcode] = pysam.AlignmentFile(output_path, 'wb', template=template)
                handle_order.append(barcode)

                # LRU close: if too many open files, close the oldest ones
                while len(open_handles) > max_open_files:
                    oldest_bc = handle_order.pop(0)
                    if oldest_bc in open_handles and oldest_bc != barcode:
                        open_handles[oldest_bc].close()
                        del open_handles[oldest_bc]

            outfh = open_handles[barcode]
            read_count = 0

            for read in reads_group:
                outfh.write(read)
                read_count += 1

            barcode_read_counts[barcode] = read_count
            processed_barcodes.add(barcode)

            # Close this barcode's handle immediately since all its reads are done
            # (sorted BAM guarantees contiguous reads for same barcode)
            if barcode in open_handles:
                open_handles[barcode].close()
                del open_handles[barcode]
                if barcode in handle_order:
                    handle_order.remove(barcode)

        input_bam.close()

        return barcode_read_counts

    except Exception as e:
        logger.error(f"Error splitting BAM file: {e}")
        # Close any remaining open handles
        for fh in open_handles.values():
            try:
                fh.close()
            except:
                pass
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
              help='Number of CPU cores (used for read counting when no gexcb provided)')
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

    # Create output directory
    split_bams_dirname = re.sub('_bismark_.*', '', os.path.basename(bam))
    split_output_dir = f'{outdir}/{split_bams_dirname}'
    os.makedirs(split_output_dir, exist_ok=True)

    logger.info(f"Starting single-pass BAM splitting for {len(top_barcodes)} barcodes...")

    # Use single-pass strategy for BAM splitting
    # This reads the sorted BAM once instead of N times (one per worker)
    all_barcode_counts = split_bam_single_pass(
        bam, split_output_dir, top_barcodes)

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
