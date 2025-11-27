#!/usr/bin/env python
import os
import sys
import re
import pysam
import pandas as pd
import click
from loguru import logger
import numpy as np
from collections import defaultdict

def counts_reads(bam:str, outdir:str):
    """
    Count the number of aligned reads per barcode in a paired-end BAM file.

    Args:
        bam: Path to the input BAM file
        outdir: Output directory
    """
    input_bam = pysam.AlignmentFile(bam, "rb")
    barcode_counts = defaultdict(int)
    
    for read in input_bam.fetch(until_eof=True):
        if read.has_tag("CB"):
            barcode = read.get_tag("CB")
            barcode_counts[barcode] += 1
    
    input_bam.close()
    
    os.makedirs(outdir, exist_ok=True)
    
    df = pd.DataFrame({
        'barcode': list(barcode_counts.keys()),
        'aligned_reads': list(barcode_counts.values())
    })
    
    df = df.sort_values('aligned_reads', ascending=False)
    
    # Write results to file
    output_file = os.path.join(outdir, os.path.basename(bam).replace('.bam', '_cb_aligned_reads_counts.csv'))
    df.to_csv(output_file, index=False)
    
    logger.info(f"Counting finished: {len(barcode_counts)} barcodes; results saved to {output_file}")
    
    return output_file

@click.command()
@click.option('-b', '--bam', required=True, help='Path to the input BAM file')
@click.option('-o', '--outdir', required=True, help='Output directory')
def main(bam, outdir):
    """Count aligned reads per barcode in a paired-end BAM file"""
    logger.info(f"Start processing BAM file: {bam}")
    output_file = counts_reads(bam, outdir)
    logger.info(f"Completed; results saved to: {output_file}")

if __name__ == '__main__':
    main()
