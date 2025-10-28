#!/usr/bin/env python

import json
import sys
import os
import gzip
import bz2
from Bio import SeqIO

def open_compressed_file(filepath):
    """
    Intelligently open files, supporting compressed formats (.gz, .bz2) and regular files
    """
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    elif filepath.endswith('.bz2'):
        return bz2.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def read_chrom_sizes(chrom_sizes_file):
    """
    Read chrom.sizes file and return list of contained chromosomes
    """
    chromosomes = set()
    try:
        with open(chrom_sizes_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    chrom = line.split('\t')[0]
                    chromosomes.add(chrom)
    except FileNotFoundError:
        print(f"Error: chrom.sizes file {chrom_sizes_file} not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading chrom.sizes file: {str(e)}", file=sys.stderr)
        sys.exit(1)
    
    return chromosomes

def count_cpg_by_chromosome(fasta_file, chrom_sizes_file=None):
    """
    Count CpG sites for each chromosome and return JSON format result
    Supports compressed FASTA files (.gz, .bz2)
    If chrom_sizes_file is provided, only count chromosomes contained within
    """
    # Get filename (without path)
    filename = os.path.basename(fasta_file)
    
    # Read list of chromosomes to count
    target_chromosomes = None
    if chrom_sizes_file:
        target_chromosomes = read_chrom_sizes(chrom_sizes_file)
        print(f"Found {len(target_chromosomes)} chromosomes in chrom.sizes file", file=sys.stderr)
    
    # Initialize result dictionary
    result = {
        filename: {
            "chromosomes": {},
            "total_cg_sites": 0,
            "filtered_by_chrom_sizes": chrom_sizes_file is not None
        }
    }
    
    total_cpg = 0
    processed_chroms = 0
    skipped_chroms = 0
    
    try:
        # Use intelligent file opening function
        with open_compressed_file(fasta_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                chromosome = record.id
                
                # If chrom.sizes file is specified, only process chromosomes contained within
                if target_chromosomes and chromosome not in target_chromosomes:
                    skipped_chroms += 1
                    continue
                
                sequence = str(record.seq).upper()  # Convert to uppercase for uniform processing
                
                # Count CG sites in current chromosome
                cg_count = sequence.count("CG")
                
                # Add to results
                result[filename]["chromosomes"][chromosome] = cg_count
                total_cpg += cg_count
                processed_chroms += 1
        
        # Set total count
        result[filename]["total_cg_sites"] = total_cpg
        
        # Add statistics information
        if chrom_sizes_file:
            result[filename]["stats"] = {
                "processed_chromosomes": processed_chroms,
                "skipped_chromosomes": skipped_chroms,
                "chrom_sizes_file": os.path.basename(chrom_sizes_file)
            }
            print(f"Processed {processed_chroms} chromosomes, skipped {skipped_chroms} chromosomes", file=sys.stderr)
        
    except FileNotFoundError:
        result[filename] = {"error": f"File {fasta_file} not found"}
    except Exception as e:
        result[filename] = {"error": f"Error processing file: {str(e)}"}
    
    return result

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python count_cg_sites.py <fasta_file> [chrom_sizes_file]")
        print("  fasta_file: Input FASTA file (supports .gz and .bz2 compression)")
        print("  chrom_sizes_file: Optional chrom.sizes file to filter chromosomes")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    chrom_sizes_file = sys.argv[2] if len(sys.argv) == 3 else None
    
    # Count CpG sites
    result = count_cpg_by_chromosome(fasta_file, chrom_sizes_file)
    
    # Output JSON format result
    print(json.dumps(result, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()