#!/usr/bin/env python
import pandas as pd
import click
import os
import glob
from pathlib import Path

@click.command()
@click.option('--raw_aligned_reads_counts', '-r', required=True, help='Path to the raw aligned reads counts file or directory')
@click.option('--expected_cell_num', '-e', default = 3000, type=int, help='Expected number of cells')
@click.option('--outdir', '-o', required=True, help='Output directory')
def estimated_cells(raw_aligned_reads_counts:str, expected_cell_num:int, outdir:str):
    """
    Estimate the number of cells and filter barcodes.
    If `raw_aligned_reads_counts` is a directory, process all files ending with
    `_cb_aligned_reads_counts.csv` in that directory and aggregate by the `barcode` column.
    """
    os.makedirs(outdir, exist_ok=True)
    
    barcode_total_reads = {}
    
    # Check whether the input is a file or a directory
    if os.path.isdir(raw_aligned_reads_counts):
        csv_files = glob.glob(os.path.join(raw_aligned_reads_counts, "*_cb_aligned_reads_counts.csv"))
        
        if not csv_files:
            print(f"No *_cb_aligned_reads_counts.csv files found in {raw_aligned_reads_counts}")
            return
            
        print(f"Found {len(csv_files)} files to process")
        
        # Read all files and aggregate reads per barcode
        for csv_file in csv_files:
            try:
                df = pd.read_csv(csv_file)
                
                if 'barcode' not in df.columns or 'aligned_reads' not in df.columns:
                    print(f"Warning: {csv_file} missing required columns (barcode, aligned_reads), skipping.")
                    continue
                
                for _, row in df.iterrows():
                    barcode = row['barcode']
                    reads = row['aligned_reads']
                    if barcode in barcode_total_reads:
                        barcode_total_reads[barcode] += reads
                    else:
                        barcode_total_reads[barcode] = reads
                        
            except Exception as e:
                print(f"Error processing file {csv_file}: {str(e)}")
                continue
    else:
        try:
            df = pd.read_csv(raw_aligned_reads_counts)
            if 'barcode' not in df.columns or 'aligned_reads' not in df.columns:
                print(f"Warning: {raw_aligned_reads_counts} missing required columns (barcode, aligned_reads).")
                return
                
            for _, row in df.iterrows():
                barcode = row['barcode']
                reads = row['aligned_reads']
                barcode_total_reads[barcode] = reads
                
        except Exception as e:
            print(f"Error reading file {raw_aligned_reads_counts}: {str(e)}")
            return
    
    if not barcode_total_reads:
        print("No valid data found. Exiting.")
        return
        
    merged_df = pd.DataFrame({
        'barcode': list(barcode_total_reads.keys()),
        'aligned_reads': list(barcode_total_reads.values())
    })
    
    merged_file = os.path.join(outdir, 'merged_barcode_counts.csv')
    merged_df.to_csv(merged_file, index=False)
    
    merged_df = merged_df.sort_values(by='aligned_reads', ascending=False)
    
    # Compute threshold
    percentile = 99
    threshold_index = int(expected_cell_num * (1 - percentile / 100))
    if threshold_index >= len(merged_df):
        threshold_index = len(merged_df) - 1
        
    readscut = merged_df.iloc[threshold_index]['aligned_reads'] * 0.1
    
    filtered_df = merged_df[merged_df['aligned_reads'] > readscut]
    
    filtered_file = os.path.join(outdir, 'filtered_barcode_read_counts.csv')
    filtered_df.to_csv(filtered_file, sep=',', index=False, header=False)
    filtered_df["barcode"].to_csv(os.path.join(outdir, 'filtered_barcode'), sep=',', index=False, header=False)
    
if __name__ == '__main__':
    estimated_cells()
