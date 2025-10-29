#!/usr/bin/env python3

import sys
import os
import pandas as pd
import json
import glob
import argparse
from pathlib import Path

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Merge single cell methylation metrics from multiple CSV files')
    parser.add_argument('input_dir', help='Input directory containing count.csv files')
    parser.add_argument('--cbcsv', default=None, help='Gex cb and meth cb map csv.')
    parser.add_argument('-o', '--output', default='merged_metrics', help='Output prefix for CSV and JSON files')
    parser.add_argument('--delete-source', action='store_true', 
                       help='Delete source _all.gz.count.csv files after successful aggregation')
    return parser.parse_args()

def generate_cb_map_dcit(cbcsv):
    """Generate mapping dictionary for gex cb and meth cb"""
    if cbcsv is None:
        return None
    df = pd.read_csv(cbcsv, header = 0, sep = ',')
    return dict(zip(df['m_cb'], df['gex_cb']))

def extract_cell_barcode(filename):
    """Extract cell barcode from filename"""
    # Filename format: AAGAGTTAGAAGATTTG_all.gz.count.csv
    basename = os.path.basename(filename)
    if basename.endswith('_all.gz.count.csv'):
        return basename.replace('_all.gz.count.csv', '')
    elif basename.endswith('.count.csv'):
        return basename.replace('.count.csv', '')
    else:
        return basename

def read_single_cell_metrics(csv_file):
    """Read metrics file for a single cell"""
    try:
        df = pd.read_csv(csv_file, index_col=0, header = 0)
        return df
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return None

def aggregate_metrics_to_csv(input_dir, output_csv, gex_cb_map):
    """Aggregate metrics from all cells to CSV file"""
    csv_files = glob.iglob(f"{input_dir}/**/*allc.gz.count.csv", recursive=True)
    
    #if not csv_files:
    #    print(f"No count.csv files found in {input_dir}")
    #    return None
    
    #print(f"Found {len(csv_files)} CSV files to process")
    
    # Store aggregated data for all cells
    all_cells_data = []
    
    for csv_file in csv_files:
        cell_barcode = extract_cell_barcode(csv_file)
        df = read_single_cell_metrics(csv_file)
        
        if df is None:
            continue
            
        # Calculate total count of CpG-related contexts
        cpg_contexts = ['CGT', 'CGG', 'CGC', 'CGA']
        total_cpg_number = 0
        for context in cpg_contexts:
            if context in df.index:
                total_cpg_number += df.loc[context, 'number']
        try:
            if gex_cb_map:
                gex_cb = gex_cb_map.get(cell_barcode, None)
            else:
                gex_cb = cell_barcode
            # Calculate summary metrics for each cell
                cell_summary = {
                    'cell_barcode': cell_barcode,
                    'total_mc': df['mc'].sum(),
                    'total_cov': df['cov'].sum(),
                    'total_number': df['number'].sum(),
                    'total_cpg_number': total_cpg_number,
                    'weighted_mc_rate': df['mc'].sum() / df['cov'].sum() if df['cov'].sum() > 0 else 0,
                    'genome_cov': df['genome_cov'].iloc[0] if len(df) > 0 else 0,
                    'genome_cov_raw_umi': df['genome_cov_raw_umi'].iloc[0] if len(df) > 0 else 0,
                    'cell_saturation': df['cell_saturation'].iloc[0] if len(df) > 0 else 0,
                    'genome_cov_new_umi': df['genome_cov_new_umi'].iloc[0] if len(df) > 0 else 0,
                    'gex_cb': gex_cb
                }
        except:
            print(f"skip {cell_barcode}", flush = True)
            continue
        
        # Add detailed metrics for each methylation context
        for context in df.index:
            cell_summary[f'{context}_mc'] = df.loc[context, 'mc']
            cell_summary[f'{context}_cov'] = df.loc[context, 'cov']
            cell_summary[f'{context}_number'] = df.loc[context, 'number']
            cell_summary[f'{context}_mc_rate'] = df.loc[context, 'mc_rate']
        
        all_cells_data.append(cell_summary)
    
    # Convert to DataFrame and save
    if all_cells_data:
        summary_df = pd.DataFrame(all_cells_data)
        summary_df.to_csv(output_csv, index=False)
        print(f"CSV summary saved to {output_csv}")
        return summary_df
    else:
        print("No valid data to aggregate")
        return None

def aggregate_metrics_to_json(input_dir, output_json, gex_cb_map):
    """Aggregate metrics from all cells to JSON file"""
    csv_files = glob.iglob(f"{input_dir}/**/*allc.gz.count.csv", recursive=True)
    
    if not csv_files:
        print(f"No reads_counts.csv files found in {input_dir}")
        return None
    
    # Store detailed data for all cells
    all_cells_json = {}
    
    # Global statistics
    global_stats = {
        'total_cells': 0,
        'total_mc_all_cells': 0,
        'total_cov_all_cells': 0,
        'contexts_summary': {}
    }
    
    for csv_file in csv_files:
        cell_barcode = extract_cell_barcode(csv_file)
        df = read_single_cell_metrics(csv_file)
        
        if df is None:
            continue
            
        # Calculate total count of CpG-related contexts
        cpg_contexts = ['CGT', 'CGG', 'CGC', 'CGA']
        total_cpg_number = 0
        for context in cpg_contexts:
            if context in df.index:
                total_cpg_number += df.loc[context, 'number']
        try:    
            if gex_cb_map:
                gex_cb = gex_cb_map.get(cell_barcode, None)
            else:
                gex_cb = cell_barcode
            # Store complete data for individual cell
            cell_data = {
                'file_path': csv_file,
                'total_mc': int(df['mc'].sum()),
                'total_cov': int(df['cov'].sum()),
                'total_number': int(df['number'].sum()),
                'total_cpg_number': int(total_cpg_number),
                'weighted_mc_rate': float(df['mc'].sum() / df['cov'].sum()) if df['cov'].sum() > 0 else 0.0,
                'genome_cov': float(df['genome_cov'].iloc[0]) if len(df) > 0 else 0.0,
                'genome_cov_raw_umi': int(df['genome_cov_raw_umi'].iloc[0]) if len(df) > 0 else 0,
                'cell_saturation': float(df['cell_saturation'].iloc[0]) if len(df) > 0 else 0.0,
                'genome_cov_new_umi': int(df['genome_cov_new_umi'].iloc[0]) if len(df) > 0 else 0,
                'gex_cb': gex_cb,
                'contexts': {}
            }
        except:
            print(f"skip {cell_barcode}", flush = True)
            continue
        
        # Add data for each methylation context
        for context in df.index:
            cell_data['contexts'][context] = {
                'mc': int(df.loc[context, 'mc']),
                'cov': int(df.loc[context, 'cov']),
                'number': int(df.loc[context, 'number']),
                'mc_rate': float(df.loc[context, 'mc_rate'])
            }
            
            # Update global statistics
            if context not in global_stats['contexts_summary']:
                global_stats['contexts_summary'][context] = {
                    'total_mc': 0,
                    'total_cov': 0,
                    'total_number': 0,
                    'cell_count': 0
                }
            
            global_stats['contexts_summary'][context]['total_mc'] += int(df.loc[context, 'mc'])
            global_stats['contexts_summary'][context]['total_cov'] += int(df.loc[context, 'cov'])
            global_stats['contexts_summary'][context]['total_number'] += int(df.loc[context, 'number'])
            global_stats['contexts_summary'][context]['cell_count'] += 1
        
        all_cells_json[cell_barcode] = cell_data
        global_stats['total_cells'] += 1
        global_stats['total_mc_all_cells'] += int(df['mc'].sum())
        global_stats['total_cov_all_cells'] += int(df['cov'].sum())
    
    # Calculate global average methylation rate
    if global_stats['total_cov_all_cells'] > 0:
        global_stats['global_mc_rate'] = global_stats['total_mc_all_cells'] / global_stats['total_cov_all_cells']
    else:
        global_stats['global_mc_rate'] = 0.0
    
    # Calculate average methylation rate for each context
    for context in global_stats['contexts_summary']:
        context_data = global_stats['contexts_summary'][context]
        if context_data['total_cov'] > 0:
            context_data['average_mc_rate'] = context_data['total_mc'] / context_data['total_cov']
        else:
            context_data['average_mc_rate'] = 0.0
    
    # Build final JSON structure
    final_json = {
        'summary': global_stats,
        'cells': all_cells_json
    }
    
    # Save JSON file
    with open(output_json, 'w') as f:
        json.dump(final_json, f, indent=2)
    
    print(f"JSON summary saved to {output_json}")
    return final_json

def main():
    args = parse_arguments()
    
    # Check input directory
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist")
        sys.exit(1)
    
    # Set output file paths
    output_csv = f"{args.output}.csv"
    output_json = f"{args.output}.json"
    
    print(f"Processing files in: {args.input_dir}")
    print(f"Output CSV: {output_csv}")
    print(f"Output JSON: {output_json}")
    if args.delete_source:
        print("Warning: Source CSV files will be deleted after successful aggregation")
    if not args.cbcsv:
        gex_cb_map = None
    else:
        gex_cb_map = generate_cb_map_dcit(args.cbcsv)
    # Generate CSV summary
    csv_result = aggregate_metrics_to_csv(args.input_dir, output_csv, gex_cb_map)
    
    # Generate JSON summary
    json_result = aggregate_metrics_to_json(args.input_dir, output_json, gex_cb_map)
    
    if csv_result is not None and json_result is not None:
        print(f"\nSuccessfully processed {len(csv_result)} cells")
        print(f"CSV file: {output_csv}")
        print(f"JSON file: {output_json}")
        
        # If delete source files option is specified, delete original CSV files
        if args.delete_source:
            csv_files = glob.glob(os.path.join(args.input_dir, "*count.csv"))
            deleted_count = 0
            for csv_file in csv_files:
                try:
                    os.remove(csv_file)
                    deleted_count += 1
                    #print(f"Deleted: {csv_file}")
                except Exception as e:
                    print(f"Warning: Failed to delete {csv_file}: {e}")
            print(f"Successfully deleted {deleted_count} source CSV files")
    else:
        print("Failed to process the data")
        sys.exit(1)

if __name__ == "__main__":
    main()

