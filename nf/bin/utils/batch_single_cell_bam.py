#!/usr/bin/env python3

import os
import sys
import argparse
import glob
import subprocess
import shutil
import logging
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def run_command(cmd, shell=False):
    """Run a shell command and check for errors."""
    try:
        # logger.info(f"Running command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
        if shell:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        else:
            subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e}")
        raise

def sort_bam(bam_file, outdir, threads):
    """
    Step 1: Sort BAM by name (samtools sort -n) with safety checks.
    """
    bam_path = Path(bam_file)
    basename = bam_path.name
    if basename.endswith('.bam'):
        base_core = basename[:-4]
    else:
        base_core = basename
    
    sorted_bam_name = f"{base_core}_sortn.bam"
    sorted_bam_path = Path(outdir) / sorted_bam_name
    
    # Check existing
    if sorted_bam_path.exists():
        try:
            subprocess.run(["samtools", "quickcheck", str(sorted_bam_path)], check=True, stderr=subprocess.DEVNULL)
            logger.info(f"Skipping existing valid sorted BAM: {sorted_bam_path}")
            return str(sorted_bam_path)
        except subprocess.CalledProcessError:
            logger.warning(f"Found incomplete/invalid sorted BAM: {sorted_bam_path}. Re-sorting.")
            try:
                sorted_bam_path.unlink()
            except OSError:
                pass
    
    temp_bam = sorted_bam_path.with_suffix(".tmp.bam")
    
    cmd = [
        "samtools", "sort", "-n",
        "-@", str(threads),
        "-o", str(temp_bam),
        str(bam_path)
    ]
    
    try:
        run_command(cmd)
        temp_bam.rename(sorted_bam_path)
    except Exception as e:
        logger.error(f"Sorting failed for {bam_file}: {e}")
        if temp_bam.exists():
            try:
                temp_bam.unlink()
            except OSError:
                pass
        raise
        
    return str(sorted_bam_path)

def split_bam(sorted_bam, outdir, threads, assay_type, gex_barcodes, cbcsv, script_path):
    """
    Step 2: Split BAM into single cells using step3_split_bams.py
    """
    sorted_bam_path = Path(sorted_bam)
    samplename = sorted_bam_path.stem # removes .bam
    
    # Calculate expected output marker to check for completion
    # Replicating logic from step3_split_bams.py:
    # split_bams_dirname = re.sub('_bismark_.*', '', os.path.basename(bam))
    split_bams_dirname = re.sub('_bismark_.*', '', sorted_bam_path.name)
    marker_file = Path(outdir) / split_bams_dirname / f"{split_bams_dirname}_filtered_barcode_reads_counts.csv"
    
    if marker_file.exists():
        logger.info(f"Split results already exist and appear complete: {marker_file}. Skipping.")
        return sorted_bam
    
    cmd = [
        "python3", "-u", str(script_path),
        "--bam", str(sorted_bam_path),
        "--outdir", str(outdir),
        "--samplename", samplename,
        "--core", str(threads),
        "--gexcb", gex_barcodes
    ]
    
    if assay_type == "DD-MET3":
        if not cbcsv:
            raise ValueError("DD-MET3 requires --cbcsv")
        cmd.extend(["--cbcsv", cbcsv])
        
    run_command(cmd)
    return sorted_bam

def safe_merge_bam(f_bam, r_bam, out_bam, threads=1):
    """
    Merge forward and reverse BAMs for a single cell with safety checks.
    1. Check if output exists and is valid (samtools quickcheck).
    2. Write to temp file.
    3. Rename temp to final on success.
    """
    out_bam = Path(out_bam)
    
    # 1. Check existing
    if out_bam.exists():
        try:
            # Quick check for BAM validity
            subprocess.run(["samtools", "quickcheck", str(out_bam)], check=True, stderr=subprocess.DEVNULL)
            # logger.info(f"Skipping existing valid BAM: {out_bam}")
            return
        except subprocess.CalledProcessError:
            logger.warning(f"Found incomplete/invalid BAM: {out_bam}. Re-processing.")
            try:
                out_bam.unlink()
            except OSError:
                pass

    # 2. Prepare temp output
    temp_bam = out_bam.with_suffix(".tmp.bam")
    
    try:
        if f_bam and r_bam:
            cmd = ["samtools", "merge", "-n", "-f", "-@", str(threads), str(temp_bam), str(f_bam), str(r_bam)]
            run_command(cmd)
        elif f_bam:
            shutil.copy2(f_bam, temp_bam)
        elif r_bam:
            shutil.copy2(r_bam, temp_bam)
        
        # 3. Rename to final
        temp_bam.rename(out_bam)
        
    except Exception as e:
        logger.error(f"Failed to merge {out_bam}: {e}")
        if temp_bam.exists():
            try:
                temp_bam.unlink()
            except OSError:
                pass
        raise

def merge_single_cells(outdir, merged_root_dir, parallel_jobs=1):
    """
    Step 3: Merge forward and reverse BAMs for each single cell.
    Parallelized version.
    """
    subdirs = [d for d in Path(outdir).iterdir() if d.is_dir()]
    forward_dirs = {}
    reverse_dirs = {}
    
    for d in subdirs:
        name = d.name
        if "_forward_" in name:
            key = name.replace("_forward_", "_PLACEHOLDER_")
            forward_dirs[key] = d
        elif "_reverse_" in name:
            key = name.replace("_reverse_", "_PLACEHOLDER_")
            reverse_dirs[key] = d
            
    all_keys = set(forward_dirs.keys()) | set(reverse_dirs.keys())
    
    # Collect all merge tasks
    merge_tasks = []
    # Collect post-processing tasks (barcodes/counts merging) - keep these sequential per sample or fast enough to do in main loop
    
    # We need to setup directories first
    sample_infos = []

    for key in all_keys:
        f_dir = forward_dirs.get(key)
        r_dir = reverse_dirs.get(key)
        
        sample_base_name = key.replace("_PLACEHOLDER_", "_merged_")
        sample_base_name = sample_base_name.split("_bismark")[0]
        
        sample_merged_dir = Path(merged_root_dir) / sample_base_name
        sample_merged_dir.mkdir(parents=True, exist_ok=True)
        
        f_bams = {f.stem: f for f in f_dir.glob("*.bam")} if f_dir else {}
        r_bams = {f.stem: f for f in r_dir.glob("*.bam")} if r_dir else {}
        
        all_bcs = set(f_bams.keys()) | set(r_bams.keys())
        
        for bc in all_bcs:
            f_bam = f_bams.get(bc)
            r_bam = r_bams.get(bc)
            out_bam = sample_merged_dir / f"{bc}.bam"
            
            merge_tasks.append((f_bam, r_bam, out_bam))
            
        sample_infos.append({
            'sample_base_name': sample_base_name,
            'sample_merged_dir': sample_merged_dir,
            'f_dir': f_dir,
            'r_dir': r_dir
        })

    logger.info(f"Found {len(merge_tasks)} barcode merge tasks across {len(sample_infos)} samples.")
    
    # Run merge tasks in parallel
    # Use max_workers=parallel_jobs
    # Since these are independent, we can run them.
    
    with ProcessPoolExecutor(max_workers=parallel_jobs) as executor:
        futures = []
        for f_bam, r_bam, out_bam in merge_tasks:
            # We use 1 thread per samtools merge since we are parallelizing by file
            futures.append(executor.submit(safe_merge_bam, f_bam, r_bam, out_bam, threads=1))
            
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logger.error(f"Merge task failed: {e}")

    # After BAMs are merged, process auxiliary files (barcodes, counts)
    logger.info("Merging auxiliary files (barcodes, counts)...")
    for info in sample_infos:
        sample_base_name = info['sample_base_name']
        sample_merged_dir = info['sample_merged_dir']
        f_dir = info['f_dir']
        r_dir = info['r_dir']
        
        # 3.2 Merge filtered_barcode
        f_barcodes = list(f_dir.glob("*_filtered_barcode")) if f_dir else []
        r_barcodes = list(r_dir.glob("*_filtered_barcode")) if r_dir else []
        
        if f_barcodes or r_barcodes:
            merged_barcode_file = sample_merged_dir / f"{sample_base_name}_filtered_barcode"
            # Check if exists (simple check, or force overwrite? Auxiliary files are fast to regenerate)
            # Let's regenerate to ensure consistency with BAMs
            with open(merged_barcode_file, 'w') as out_f:
                files_str = " ".join([str(p) for p in f_barcodes + r_barcodes])
                cmd = f"cat {files_str} | sort | uniq"
                subprocess.run(cmd, shell=True, stdout=out_f, check=True)
        
        # 3.3 Merge reads counts
        f_counts = list(f_dir.glob("*_filtered_barcode_reads_counts.csv")) if f_dir else []
        r_counts = list(r_dir.glob("*_filtered_barcode_reads_counts.csv")) if r_dir else []
        
        if f_counts or r_counts:
            merged_counts_file = sample_merged_dir / f"{sample_base_name}_filtered_barcode_reads_counts.csv"
            barcode_counts = {}
            
            for csv_file in f_counts + r_counts:
                try:
                    with open(csv_file, 'r') as f:
                        header = next(f) 
                        for line in f:
                            parts = line.strip().split(',')
                            if len(parts) >= 2:
                                count = int(parts[0])
                                bc = parts[1]
                                barcode_counts[bc] = barcode_counts.get(bc, 0) + count
                except Exception as e:
                    logger.warning(f"Error reading {csv_file}: {e}")
            
            with open(merged_counts_file, 'w') as f:
                f.write("reads_counts,barcode\n")
                for bc, count in barcode_counts.items():
                    f.write(f"{count},{bc}\n")


def main():
    parser = argparse.ArgumentParser(description="Batch process single cell BAMs (Sort -> Split -> Merge)")
    parser.add_argument("--assay_type", required=True, choices=["DD-MET5", "DD-MET3"], help="Assay type")
    parser.add_argument("--bam_dir", required=True, help="Directory containing Bismark BAM files")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--gex_barcodes", required=True, help="Path to RNA barcodes.tsv.gz")
    parser.add_argument("--cbcsv", help="Path to whitelist CSV (required for DD-MET3)")
    parser.add_argument("--threads", type=int, default=4, help="Threads per job")
    parser.add_argument("--parallel_jobs", type=int, default=8, help="Number of parallel jobs")
    parser.add_argument("--keep_temp", action="store_true", help="Keep intermediate files (sorted BAMs)")
    
    args = parser.parse_args()
    
    # Locate step3 script
    script_dir = Path(__file__).parent.parent
    step3_script = script_dir / "step3_split_bams.py"
    if not step3_script.exists():
        script_dir = Path(__file__).parent
        step3_script = script_dir / "step3_split_bams.py"
        if not step3_script.exists():
            logger.error(f"step3_split_bams.py not found at {step3_script}")
            sys.exit(1)
        
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    
    # Define subdirectories
    sortn_dir = Path(args.outdir) / "sortn_bam"
    split_dir = Path(args.outdir) / "split_bam"
    merge_dir = Path(args.outdir) / "merge_fr_bam"
    
    sortn_dir.mkdir(parents=True, exist_ok=True)
    split_dir.mkdir(parents=True, exist_ok=True)
    merge_dir.mkdir(parents=True, exist_ok=True)
    
    # Gather BAMs
    bams = glob.glob(os.path.join(args.bam_dir, "*_bismark_bt2_pe.bam"))
    if not bams:
        logger.error(f"No *_bismark_bt2_pe.bam files found in {args.bam_dir}")
        sys.exit(1)
    
    logger.info(f"Found {len(bams)} BAM files to process.")
    
    # Step 1 & 2: Sort and Split (Pipeline per BAM)
    # We can chain them or do all sorts then all splits.
    # To maximize resource usage and match user request for "Parallel", 
    # we can run the pipeline for each BAM in parallel.
    # However, Step 3 requires BOTH forward and reverse to be done.
    
    # Let's run Sort+Split in parallel tasks.
    
    with ProcessPoolExecutor(max_workers=args.parallel_jobs) as executor:
        # 1. Sort
        logger.info("Starting Step 1: Sort BAMs...")
        future_to_bam = {
            executor.submit(sort_bam, bam, sortn_dir, args.threads): bam 
            for bam in bams
        }
        
        sorted_bams = []
        for future in as_completed(future_to_bam):
            bam = future_to_bam[future]
            try:
                sorted_bam = future.result()
                sorted_bams.append(sorted_bam)
                logger.info(f"Sorted: {sorted_bam}")
            except Exception as e:
                logger.error(f"Sorting failed for {bam}: {e}")
        
        # 2. Split
        logger.info("Starting Step 2: Split BAMs...")
        future_to_split = {
            executor.submit(split_bam, sb, split_dir, args.threads, args.assay_type, args.gex_barcodes, args.cbcsv, step3_script): sb
            for sb in sorted_bams
        }
        
        for future in as_completed(future_to_split):
            sb = future_to_split[future]
            try:
                future.result()
                logger.info(f"Split completed for: {sb}")
                
                # Cleanup intermediate sorted BAM
                if not args.keep_temp:
                    try:
                        os.remove(sb)
                        logger.info(f"Deleted intermediate sorted BAM: {sb}")
                    except OSError as e:
                        logger.warning(f"Failed to delete {sb}: {e}")
                        
            except Exception as e:
                logger.error(f"Splitting failed for {sb}: {e}")

    # Step 3: Merge
    logger.info("Starting Step 3: Merge Forward/Reverse...")
    merge_single_cells(split_dir, merge_dir, args.parallel_jobs)
    logger.info("All processing completed.")

if __name__ == "__main__":
    main()
