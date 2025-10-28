#!/usr/bin/env python3
"""
FASTQ File Merging Tool
Efficiently merge multiple FASTQ files using the dnaio library

Features:
- Support merging any number of FASTQ files
- Automatic detection of compression formats (.gz, .bz2, etc.)
- Support for single-end and paired-end sequencing data
- Memory-efficient streaming processing
- Detailed progress monitoring and statistics
- Support for quality control and filtering options

Author: AI Assistant
Version: 1.0
"""

import os
import sys
import time
import argparse
from pathlib import Path
from typing import List, Optional, Tuple
import dnaio
from loguru import logger
from tqdm import tqdm

# Configure logging
logger.remove()
logger.add(sys.stderr, format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}", level="INFO")

def validate_fastq_files(input_files: List[str]) -> List[str]:
    """
    Validate whether input FASTQ files exist and are readable
    
    Args:
        input_files: List of input file paths
    
    Returns:
        List of validated file paths
    
    Raises:
        FileNotFoundError: If file does not exist
        ValueError: If file format is incorrect
    """
    valid_files = []
    
    for file_path in input_files:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File does not exist: {file_path}")
        
        if not os.path.isfile(file_path):
            raise ValueError(f"Path is not a file: {file_path}")
        
        # Check file extensions
        valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fastq.bz2', '.fq.bz2']
        if not any(file_path.lower().endswith(ext) for ext in valid_extensions):
            logger.warning(f"File extension may be incorrect: {file_path}")
        
        valid_files.append(file_path)
        logger.info(f"Validation passed: {file_path}")
    
    return valid_files

def detect_paired_end(input_files: List[str]) -> Tuple[bool, List[str], List[str]]:
    """
    Detect whether data is paired-end sequencing
    
    Args:
        input_files: List of input file paths
    
    Returns:
        (is_paired, r1_files, r2_files)
    """
    r1_files = []
    r2_files = []
    
    for file_path in input_files:
        filename = os.path.basename(file_path)
        
        # Detect R1/R2 patterns
        if '_R1_' in filename or '_1.f' in filename or '.R1.' in filename:
            r1_files.append(file_path)
        elif '_R2_' in filename or '_2.f' in filename or '.R2.' in filename:
            r2_files.append(file_path)
        else:
            # If unable to determine, assume single-end data
            r1_files.append(file_path)
    
    is_paired = len(r1_files) > 0 and len(r2_files) > 0 and len(r1_files) == len(r2_files)
    
    if is_paired:
        logger.info(f"Detected paired-end sequencing data: {len(r1_files)} file pairs")
        # Sort to ensure correct pairing
        r1_files.sort()
        r2_files.sort()
    else:
        logger.info(f"Detected single-end sequencing data: {len(r1_files)} files")
        r2_files = []
    
    return is_paired, r1_files, r2_files

def merge_single_end_fastq(input_files: List[str], output_file: str, 
                          min_length: int = 0, max_length: int = 0) -> dict:
    """
    Merge single-end FASTQ files
    
    Args:
        input_files: List of input file paths
        output_file: Output file path
        min_length: Minimum sequence length filter
        max_length: Maximum sequence length filter (0 means no limit)
    
    Returns:
        Statistics dictionary
    """
    logger.info(f"Starting to merge {len(input_files)} single-end FASTQ files")
    
    stats = {
        'total_reads': 0,
        'filtered_reads': 0,
        'output_reads': 0,
        'total_bases': 0
    }
    
    # Ensure output directory exists
    #os.makedirs(os.path.dirname(output_file), exist_ok=True)
    print(input_files)
    with dnaio.open(file1 = output_file, mode='w') as writer:
        for input_file in tqdm(input_files, desc="Processing files"):
            logger.info(f"Processing: {input_file}")
            
            with dnaio.open(file1 = input_file, mode = 'r') as reader:
                for record in tqdm(reader, desc=f"Reading {os.path.basename(input_file)}", leave=False):
                    stats['total_reads'] += 1
                    
                    # Length filtering
                    seq_length = len(record.sequence)
                    if min_length > 0 and seq_length < min_length:
                        stats['filtered_reads'] += 1
                        continue
                    
                    if max_length > 0 and seq_length > max_length:
                        stats['filtered_reads'] += 1
                        continue
                    
                    # Write to output file
                    writer.write(record)
                    stats['output_reads'] += 1
                    stats['total_bases'] += seq_length
    
    return stats

def merge_paired_end_fastq(r1_files: List[str], r2_files: List[str], 
                          output_r1: str, output_r2: str,
                          min_length: int = 0, max_length: int = 0) -> dict:
    """
    Merge paired-end FASTQ files
    
    Args:
        r1_files: List of R1 file paths
        r2_files: List of R2 file paths
        output_r1: R1 output file path
        output_r2: R2 output file path
        min_length: Minimum sequence length filter
        max_length: Maximum sequence length filter
    
    Returns:
        Statistics dictionary
    """
    logger.info(f"Starting to merge {len(r1_files)} pairs of paired-end FASTQ files")
    
    stats = {
        'total_pairs': 0,
        'filtered_pairs': 0,
        'output_pairs': 0,
        'total_bases_r1': 0,
        'total_bases_r2': 0
    }
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(output_r1), exist_ok=True)
    os.makedirs(os.path.dirname(output_r2), exist_ok=True)
    
    with dnaio.open(output_r1, mode='w') as writer_r1, \
         dnaio.open(output_r2, mode='w') as writer_r2:
        
        for r1_file, r2_file in tqdm(zip(r1_files, r2_files), 
                                    total=len(r1_files), desc="Processing file pairs"):
            logger.info(f"Processing: {r1_file} and {r2_file}")
            
            with dnaio.open(r1_file) as reader_r1, \
                 dnaio.open(r2_file) as reader_r2:
                
                for record_r1, record_r2 in tqdm(zip(reader_r1, reader_r2), 
                                                desc=f"Reading {os.path.basename(r1_file)}", 
                                                leave=False):
                    stats['total_pairs'] += 1
                    
                    # Length filtering
                    r1_length = len(record_r1.sequence)
                    r2_length = len(record_r2.sequence)
                    
                    filter_pair = False
                    if min_length > 0 and (r1_length < min_length or r2_length < min_length):
                        filter_pair = True
                    
                    if max_length > 0 and (r1_length > max_length or r2_length > max_length):
                        filter_pair = True
                    
                    if filter_pair:
                        stats['filtered_pairs'] += 1
                        continue
                    
                    # Write to output files
                    writer_r1.write(record_r1)
                    writer_r2.write(record_r2)
                    stats['output_pairs'] += 1
                    stats['total_bases_r1'] += r1_length
                    stats['total_bases_r2'] += r2_length
    
    return stats

def print_statistics(stats: dict, is_paired: bool, processing_time: float):
    """
    Print statistics information
    
    Args:
        stats: Statistics dictionary
        is_paired: Whether data is paired-end
        processing_time: Processing time
    """
    logger.info("=== Merge Completion Statistics ===")
    logger.info(f"Processing time: {processing_time:.2f} seconds")
    
    if is_paired:
        logger.info(f"Total input sequence pairs: {stats['total_pairs']:,}")
        logger.info(f"Filtered sequence pairs: {stats['filtered_pairs']:,}")
        logger.info(f"Output sequence pairs: {stats['output_pairs']:,}")
        logger.info(f"R1 total bases: {stats['total_bases_r1']:,}")
        logger.info(f"R2 total bases: {stats['total_bases_r2']:,}")
        logger.info(f"Total bases: {stats['total_bases_r1'] + stats['total_bases_r2']:,}")
        
        if stats['total_pairs'] > 0:
            filter_rate = (stats['filtered_pairs'] / stats['total_pairs']) * 100
            logger.info(f"Filter rate: {filter_rate:.2f}%")
        
        if processing_time > 0:
            speed = stats['output_pairs'] / processing_time
            logger.info(f"Processing speed: {speed:.0f} sequence pairs/second")
    else:
        logger.info(f"Total input sequences: {stats['total_reads']:,}")
        logger.info(f"Filtered sequences: {stats['filtered_reads']:,}")
        logger.info(f"Output sequences: {stats['output_reads']:,}")
        logger.info(f"Total bases: {stats['total_bases']:,}")
        
        if stats['total_reads'] > 0:
            filter_rate = (stats['filtered_reads'] / stats['total_reads']) * 100
            logger.info(f"Filter rate: {filter_rate:.2f}%")
        
        if processing_time > 0:
            speed = stats['output_reads'] / processing_time
            logger.info(f"Processing speed: {speed:.0f} sequences/second")

def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
        description="Merge multiple FASTQ files using dnaio",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # Merge single-end data
  python merge_fastq_files.py -i file1.fastq file2.fastq file3.fastq -o merged.fastq
  
  # Merge paired-end data
  python merge_fastq_files.py -i sample1_R1.fastq sample1_R2.fastq sample2_R1.fastq sample2_R2.fastq -o merged
  
  # Add length filtering
  python merge_fastq_files.py -i *.fastq -o merged.fastq --min-length 50 --max-length 300
        """
    )
    
    parser.add_argument('-i', '--input', nargs='+', required=True,
                       help='Input FASTQ file paths (supports wildcards)')
    
    parser.add_argument('-o', '--output', required=True,
                       help='Output file path (single-end data) or prefix (paired-end data)')
    
    parser.add_argument('--min-length', type=int, default=0,
                       help='Minimum sequence length filter (default: 0, no filtering)')
    
    parser.add_argument('--max-length', type=int, default=0,
                       help='Maximum sequence length filter (default: 0, no limit)')
    
    parser.add_argument('--force-single-end', action='store_true',
                       help='Force processing as single-end data')
    
    parser.add_argument('--compression-level', type=int, default=6,
                       help='Compression level (1-9, default: 6)')
    
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output mode')
    
    args = parser.parse_args()
    
    # Set log level
    if args.verbose:
        logger.remove()
        logger.add(sys.stderr, format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}", level="DEBUG")
    
    start_time = time.time()
    
    try:
        # Validate input files
        logger.info(f"Validating {len(args.input)} input files")
        valid_files = validate_fastq_files(args.input)
        
        # Detect data type
        if args.force_single_end:
            is_paired = False
            r1_files = valid_files
            r2_files = []
        else:
            is_paired, r1_files, r2_files = detect_paired_end(valid_files)
        
        # Execute merge
        if is_paired:
            # Paired-end data merge
            output_r1 = f"{args.output}_R1.fastq.gz"
            output_r2 = f"{args.output}_R2.fastq.gz"
            
            logger.info(f"Output files: {output_r1}, {output_r2}")
            
            stats = merge_paired_end_fastq(
                r1_files, r2_files, output_r1, output_r2,
                args.min_length, args.max_length
            )
        else:
            # Single-end data merge
            output_file = args.output
            if not output_file.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                output_file += '.fastq.gz'
            
            logger.info(f"Output file: {output_file}")
            
            stats = merge_single_end_fastq(
                r1_files, output_file,
                args.min_length, args.max_length
            )
        
        # Calculate processing time and print statistics
        end_time = time.time()
        processing_time = end_time - start_time
        
        print_statistics(stats, is_paired, processing_time)
        
        logger.info("✅ FASTQ file merge completed!")
        
    except Exception as e:
        logger.error(f"Error occurred during merge process: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()