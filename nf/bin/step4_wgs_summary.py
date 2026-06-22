#!/usr/bin/env python

import json
import os
import re
import click
import gzip
import pandas as pd

def extract_percentage(text, pattern):
    """Extract percentage from text"""
    match = re.search(pattern, text)
    return float(match.group(1)) if match else 0.0


def extract_number(text, pattern):
    """Extract number from text"""
    match = re.search(pattern, text)
    return int(match.group(1)) if match else 0

def parse_bismark_report(report_file):
    """
    Parse bismark report file and extract methylation context metrics.

    Bismark PE reports are line-oriented. Parsing line by line avoids
    repeatedly scanning the full report with many regular expressions.
    """
    metrics = {
        'totalreads': 0,
        'uniquereads': 0,
        'not_uniq_alignedpairs': 0,
        'alignedreads': 0,
        'cpg_methylation_rate': 0.0,
        'chg_methylation_rate': 0.0,
        'chh_methylation_rate': 0.0,
        'unknown_methylation_rate': 0.0,
        'methylated_cpg': 0,
        'unmethylated_cpg': 0,
        'methylated_chg': 0,
        'unmethylated_chg': 0,
        'methylated_chh': 0,
        'unmethylated_chh': 0,
        'methylated_unknown': 0,
        'unmethylated_unknown': 0,
    }

    def _value_after_colon(line):
        return line.split(':', 1)[1].strip() if ':' in line else ''

    def _parse_int(line):
        value = _value_after_colon(line).replace(',', '')
        match = re.search(r'\d+', value)
        return int(match.group(0)) if match else 0

    def _parse_pct(line):
        value = _value_after_colon(line)
        match = re.search(r'([\d.]+)\s*%', value)
        return float(match.group(1)) if match else 0.0

    try:
        with open(report_file, 'r') as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith('Sequence pairs analysed in total:'):
                    metrics['totalreads'] = _parse_int(line) * 2
                elif line.startswith('Number of paired-end alignments with a unique best hit:'):
                    metrics['uniquereads'] = _parse_int(line) * 2
                elif line.startswith('Sequence pairs did not map uniquely:'):
                    metrics['not_uniq_alignedpairs'] = _parse_int(line) * 2
                elif line.startswith('C methylated in CpG context:'):
                    metrics['cpg_methylation_rate'] = _parse_pct(line)
                elif line.startswith('C methylated in CHG context:'):
                    metrics['chg_methylation_rate'] = _parse_pct(line)
                elif line.startswith('C methylated in CHH context:'):
                    metrics['chh_methylation_rate'] = _parse_pct(line)
                elif line.startswith('C methylated in Unknown context'):
                    metrics['unknown_methylation_rate'] = _parse_pct(line)
                elif line.startswith("Total methylated C's in CpG context:"):
                    metrics['methylated_cpg'] = _parse_int(line)
                elif line.startswith("Total unmethylated C's in CpG context:"):
                    metrics['unmethylated_cpg'] = _parse_int(line)
                elif line.startswith("Total methylated C's in CHG context:"):
                    metrics['methylated_chg'] = _parse_int(line)
                elif line.startswith("Total unmethylated C's in CHG context:"):
                    metrics['unmethylated_chg'] = _parse_int(line)
                elif line.startswith("Total methylated C's in CHH context:"):
                    metrics['methylated_chh'] = _parse_int(line)
                elif line.startswith("Total unmethylated C's in CHH context:"):
                    metrics['unmethylated_chh'] = _parse_int(line)
                elif line.startswith("Total methylated C's in Unknown context:"):
                    metrics['methylated_unknown'] = _parse_int(line)
                elif line.startswith("Total unmethylated C's in Unknown context:"):
                    metrics['unmethylated_unknown'] = _parse_int(line)

        metrics['alignedreads'] = metrics['uniquereads'] + metrics['not_uniq_alignedpairs']
    except Exception as e:
        print(f"Error parsing {report_file}: {e}")
        return None

    return metrics


def get_methylation_metrics(samplename, outdir):
    """
    Get methylation context metrics
    
    Args:
        samplename: Sample name
        outdir: Output directory (containing step2/bismark subdirectory)
        
    Returns:
        dict: Dictionary containing methylation metrics, returns default values if unable to obtain
    """
    bismark_dir = os.path.join(outdir, 'step2', 'bismark')
    if not os.path.exists(bismark_dir):
        bismark_dir = outdir
    # Case 1: Check if single report file exists
    single_report = os.path.join(bismark_dir, f'{samplename}_1_bismark_bt2_PE_report.txt')
    if os.path.exists(single_report):
        print(f"Found single bismark report: {single_report}")
        metrics = parse_bismark_report(single_report)
        if metrics:
            return {
                'cpg_methylation_rate': metrics['cpg_methylation_rate'],
                'chg_methylation_rate': metrics['chg_methylation_rate'],
                'chh_methylation_rate': metrics['chh_methylation_rate'],
                'unknown_methylation_rate': metrics['unknown_methylation_rate'],
                'mapping': {
                    'totalreads': metrics['totalreads'],
                    'alignedreads': metrics['alignedreads'],
                    'uniquereads': metrics['uniquereads'],
                    'Reads Mapped to Genome': metrics['alignedreads'] / metrics['totalreads'] if metrics['totalreads'] > 0 else 0.0,
                    'Reads Mapped Confidently to Genome': metrics['uniquereads'] / metrics['totalreads'] if metrics['totalreads'] > 0 else 0.0 
                }
            }
        else:
            print(f"Warning: Single bismark report file is empty for {samplename}")
            return None
    
    # Case 2: Check if forward and reverse report files exist
    forward_report = os.path.join(bismark_dir, f'{samplename}_forward_1_bismark_bt2_PE_report.txt')
    reverse_report = os.path.join(bismark_dir, f'{samplename}_reverse_1_bismark_bt2_PE_report.txt')
        
    if os.path.exists(forward_report) and os.path.exists(reverse_report):
        print(f"Found forward and reverse bismark reports: {forward_report}, {reverse_report}")
        forward_metrics = parse_bismark_report(forward_report)
        reverse_metrics = parse_bismark_report(reverse_report)
            
        if forward_metrics and reverse_metrics:
            # Calculate average methylation rate
            # For CpG: (forward_methylated + reverse_methylated) / (forward_total + reverse_total) * 100
            cpg_methylated_total = forward_metrics['methylated_cpg'] + reverse_metrics['methylated_cpg']
            cpg_unmethylated_total = forward_metrics['unmethylated_cpg'] + reverse_metrics['unmethylated_cpg']
            cpg_total = cpg_methylated_total + cpg_unmethylated_total
            cpg_rate = (cpg_methylated_total / cpg_total * 100) if cpg_total > 0 else 0.0
                
            chg_methylated_total = forward_metrics['methylated_chg'] + reverse_metrics['methylated_chg']
            chg_unmethylated_total = forward_metrics['unmethylated_chg'] + reverse_metrics['unmethylated_chg']
            chg_total = chg_methylated_total + chg_unmethylated_total
            chg_rate = (chg_methylated_total / chg_total * 100) if chg_total > 0 else 0.0
                
            chh_methylated_total = forward_metrics['methylated_chh'] + reverse_metrics['methylated_chh']
            chh_unmethylated_total = forward_metrics['unmethylated_chh'] + reverse_metrics['unmethylated_chh']
            chh_total = chh_methylated_total + chh_unmethylated_total
            chh_rate = (chh_methylated_total / chh_total * 100) if chh_total > 0 else 0.0
                
            unknown_methylated_total = forward_metrics['methylated_unknown'] + reverse_metrics['methylated_unknown']
            unknown_unmethylated_total = forward_metrics['unmethylated_unknown'] + reverse_metrics['unmethylated_unknown']
            unknown_total = unknown_methylated_total + unknown_unmethylated_total
            unknown_rate = (unknown_methylated_total / unknown_total * 100) if unknown_total > 0 else 0.0
                
            return {
                'cpg_methylation_rate': cpg_rate,
                'chg_methylation_rate': chg_rate,
                'chh_methylation_rate': chh_rate,
                'unknown_methylation_rate': unknown_rate,
                'mapping': {
                    'totalreads': forward_metrics['totalreads'] + reverse_metrics['totalreads'],
                    'alignedreads': forward_metrics['alignedreads'] + reverse_metrics['alignedreads'],
                    'uniquereads': forward_metrics['uniquereads'] + reverse_metrics['uniquereads'],
                    'Reads Mapped to Genome': (forward_metrics['alignedreads'] + reverse_metrics['alignedreads']) / (forward_metrics['totalreads'] + reverse_metrics['totalreads']) if (forward_metrics['totalreads'] + reverse_metrics['totalreads']) > 0 else 0.0,
                    'Reads Mapped Confidently to Genome': (forward_metrics['uniquereads'] + reverse_metrics['uniquereads']) / (forward_metrics['totalreads'] + reverse_metrics['totalreads']) if (forward_metrics['totalreads'] + reverse_metrics['totalreads']) > 0 else 0.0 
                }
            }
        else:
            print(f"Warning: Forward or reverse bismark report file is empty for {samplename}")
            return None
    # Case 3: Check for multiple report files
    report_files = [f for f in os.listdir(bismark_dir) if f.startswith(f'{samplename}_') and f.endswith('_bismark_bt2_PE_report.txt')]
    if report_files:
        print(f"Found multiple bismark report files for {samplename}: {report_files}")
        # Merge metrics from all report files
        total_metrics = {
            'totalreads': 0,
            'alignedreads': 0,
            'uniquereads': 0,
            'methylated_cpg': 0,
            'unmethylated_cpg': 0,
            'methylated_chg': 0,
            'unmethylated_chg': 0,
            'methylated_chh': 0,
            'unmethylated_chh': 0,
            'methylated_unknown': 0,
            'unmethylated_unknown': 0
        }
        for report_file in report_files:
            report_path = os.path.join(bismark_dir, report_file)
            metrics = parse_bismark_report(report_path)
            if metrics:
                for key, value in metrics.items():
                    try: 
                        total_metrics[key] += value
                    except KeyError:
                        pass
        
        # Calculate average methylation rates after processing all files
        cpg_methylated_total = total_metrics['methylated_cpg']
        cpg_unmethylated_total = total_metrics['unmethylated_cpg']
        cpg_total = cpg_methylated_total + cpg_unmethylated_total
        cpg_rate = (cpg_methylated_total / cpg_total * 100) if cpg_total > 0 else 0.0
        
        chg_methylated_total = total_metrics['methylated_chg']
        chg_unmethylated_total = total_metrics['unmethylated_chg']
        chg_total = chg_methylated_total + chg_unmethylated_total
        chg_rate = (chg_methylated_total / chg_total * 100) if chg_total > 0 else 0.0
        
        chh_methylated_total = total_metrics['methylated_chh']
        chh_unmethylated_total = total_metrics['unmethylated_chh']
        chh_total = chh_methylated_total + chh_unmethylated_total
        chh_rate = (chh_methylated_total / chh_total * 100) if chh_total > 0 else 0.0
        
        unknown_methylated_total = total_metrics['methylated_unknown']
        unknown_unmethylated_total = total_metrics['unmethylated_unknown']
        unknown_total = unknown_methylated_total + unknown_unmethylated_total
        unknown_rate = (unknown_methylated_total / unknown_total * 100) if unknown_total > 0 else 0.0
        
        return {
            'cpg_methylation_rate': cpg_rate,
            'chg_methylation_rate': chg_rate,
            'chh_methylation_rate': chh_rate,
            'unknown_methylation_rate': unknown_rate,
            'mapping': {
                'totalreads': total_metrics['totalreads'],
                'alignedreads': total_metrics['alignedreads'],
                'uniquereads': total_metrics['uniquereads'],
                'Reads Mapped to Genome': total_metrics['alignedreads'] / total_metrics['totalreads'] if total_metrics['totalreads'] > 0 else 0.0,
                'Reads Mapped Confidently to Genome': total_metrics['uniquereads'] / total_metrics['totalreads'] if total_metrics['totalreads'] > 0 else 0.0
            }
        }
    
    else:
        # If none are found, return default values
        print(f"Warning: No bismark report files found for {samplename}")
        return {
            'cpg_methylation_rate': 0.0,
            'chg_methylation_rate': 0.0,
            'chh_methylation_rate': 0.0,
            'unknown_methylation_rate': 0.0,
            'mapping': {
                'totalreads': 0,
                'alignedreads': 0,
                'uniquereads': 0,
                'Reads Mapped to Genome': 0.0,
                'Reads Mapped Confidently to Genome': 0.0
            }
        }

def get_total_cpg(allcfile: str, genome_total_cpg_file: str = None) -> dict:
    """
    Get total CpG count from allc file
    
    Args:
        allcfile: Path to allc file
        
    Returns:
        int: Total CpG count
    """
    total_cpg = 0
    with gzip.open(allcfile, 'rt') as f:
        for line in f:
            total_cpg += 1
    # Get total CpG count from genome_total_cpg file
    genome_total_cpg = 0
    
    if genome_total_cpg_file and os.path.exists(genome_total_cpg_file):
        with open(genome_total_cpg_file, 'r') as f:
            genome_info = json.load(f)
            genome_total_cpg = next(iter(genome_info.values()))["total_cg_sites"]
    # Convert total_cpg / genome_total_cpg to percentage format %.2f
    cpg_methylation_rate = total_cpg / genome_total_cpg if genome_total_cpg > 0 else 0.0
    return {"Total CPGs Detected": total_cpg, "CpG Coverage rate": cpg_methylation_rate}

def parse_cell_info(cells_reads_csv: str, cells_allc_metric_csv: str) -> dict:
    """
    Summarize cell-level quality control information
    
    Args:
        cells_reads_csv: Path to filtered_barcode_reads_counts.csv file
        cells_allc_metric_csv: Path to sample_cells.csv file
        
    Returns:
        dict: Dictionary containing cell information
    """
    # cells_reads has two columns: barcode, reads_counts
    print(cells_reads_csv, flush = True)
    # sort value by reads counts, ascending = True
    cells_reads = pd.read_csv(cells_reads_csv, header = 0, sep = ",").sort_values(by = "reads_counts", ascending = False)
    cells_reads.index = cells_reads["barcode"]
    print(cells_allc_metric_csv, flush = True)
    cells_allc_metrics = pd.read_csv(cells_allc_metric_csv, index_col=0, header = 0, sep = ",").sort_values(by = "genome_cov", ascending = False)
    cells_allc_metrics["barcode"] = cells_allc_metrics.index.to_series().astype(str).str.replace("_allc.gz", "", regex=False)
    aligned_cells = cells_allc_metrics[cells_allc_metrics["barcode"].isin(cells_reads.index)].copy()
    if aligned_cells.shape[0] == 0:
        aligned_cells = cells_allc_metrics.copy()
    aligned_cells = aligned_cells.sort_values(by="genome_cov", ascending=False)

    def safe_reads(bc: str) -> int:
        return int(cells_reads.loc[bc, "reads_counts"]) if bc in cells_reads.index else 0

    max_genome_cov_cell = aligned_cells.iloc[0]
    max_genome_cov = max_genome_cov_cell["genome_cov"]
    max_cell_saturation = max_genome_cov_cell["cell_saturation"]
    max_total_cpg_number = max_genome_cov_cell["total_cpg_number"]
    max_reads_counts = safe_reads(str(max_genome_cov_cell.get("barcode", max_genome_cov_cell.name)))
    
    median_genome_cov_cell = aligned_cells.iloc[len(aligned_cells) // 2]
    median_genome_cov = median_genome_cov_cell["genome_cov"]
    median_cell_saturation = median_genome_cov_cell["cell_saturation"]
    median_total_cpg_number = median_genome_cov_cell["total_cpg_number"]
    median_reads_counts = safe_reads(str(median_genome_cov_cell.get("barcode", median_genome_cov_cell.name)))
    
    cell_info = {
        "Estimated Number of Cells": cells_reads.shape[0],
        "Genome Coverage rate of max cell": max_genome_cov,
        "Saturation of max cell": max_cell_saturation,
        "CPGs of max cell": int(max_total_cpg_number),
        "Reads of max cell": int(max_reads_counts),
        "Genome Coverage rate of median cell": median_genome_cov,
        "Saturation of median cell": median_cell_saturation,
        "CPGs of median cell": int(median_total_cpg_number),
        "Reads of median cell": int(median_reads_counts),
        "Reads in Cells": int(sum(cells_reads["reads_counts"]))  # Default value, should be calculated elsewhere
    }
    return cell_info

@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--outdir", help="outdir.")
@click.option("--samplename", help="samplename.")
@click.option("--summary_json", help="summary json.")
@click.option("--genome_info_json", help="genome info json.")
def outcsv(outdir, samplename, summary_json, genome_info_json):
    os.makedirs(outdir, exist_ok = True)
    with open(summary_json,'r') as fh:
        summary = json.load(fh)
    
    rawreads = summary["stat"]["total"]
    vaildreads = summary["stat"]["valid"]
    vaildratio = vaildreads / rawreads
    rate_7fratio = float(summary["stat"]["rate_7f"])
    rate_17lmeratio = float(summary["stat"]["rate_17lme"])
    rate_7f17lmeratio = float(summary["stat"]["rate_7f17lme"])
    conversion = f'{float(summary["stat"]["ct_mean"]):.2%}'
    cc_ratio = f'{float(summary["stat"]["cc_mean"]):.2%}'
    
    # Calculate new quality control metrics
    too_short = summary["stat"]["too_short"]
    forward = summary["stat"]["forward"]
    reverse = summary["stat"]["reverse"]
    
    # Get chimeric statistics, set to 0 if not present
    forward_chimeric = summary["stat"].get("forward_chimeric", 0)
    reverse_chimeric = summary["stat"].get("reverse_chimeric", 0)
    
    # Dropped_Too_Short: (valid-too_short)/total*100%
    dropped_too_short_ratio = too_short / rawreads
    
    # Dropped_Chimeric: (valid-too_short-forward_chimeric-reverse_chimeric)/total*100%
    dropped_chimeric_ratio = (forward_chimeric + reverse_chimeric) / (vaildreads - too_short)
    
    # Get methylation context metrics
    methylation_metrics = get_methylation_metrics(samplename, outdir)
    
    summary.update(methylation_metrics)
    raw = f'{summary["stat"]["total"]}'
    vaild = f'{vaildratio:.2%}'
    dropped_too_short = f'{dropped_too_short_ratio:.2%}'
    dropped_chimeric = f'{dropped_chimeric_ratio:.2%}'
    cpg_methylation = f'{methylation_metrics["cpg_methylation_rate"]:.1f}%'
    chg_methylation = f'{methylation_metrics["chg_methylation_rate"]:.1f}%'
    chh_methylation = f'{methylation_metrics["chh_methylation_rate"]:.1f}%'
    unknown_methylation = f'{methylation_metrics["unknown_methylation_rate"]:.1f}%'
    rate_7f = f'{rate_7fratio:.2%}'
    rate_17lme = f'{rate_17lmeratio:.2%}'
    rate_7f17lme = f'{rate_7f17lmeratio:.2%}'
    
    
    mapgenome = f'{summary["mapping"]["Reads Mapped to Genome"]:.2%}'
    confidently = f'{summary["mapping"]["Reads Mapped Confidently to Genome"]:.2%}'
    
    step3_dir = outdir
    if os.path.exists(f"{outdir}/{samplename}_methy/step3/"):
        step3_dir = f"{outdir}/{samplename}_methy/step3/"
    if summary.get("cells") is None:
        summary["cells"] = {}
    summary["cells"].update(get_total_cpg(f"{step3_dir}/{samplename}_.CGN-Merge.allc.tsv.gz", genome_info_json))
    total_cpgs = summary["cells"]["Total CPGs Detected"]
    cpg_coverage = f'{summary["cells"]["CpG Coverage rate"]:.2%}'
    
    cell_info = parse_cell_info(
        f"{step3_dir}/filtered_barcode_reads_counts.csv",
        f"{step3_dir}/{samplename}_cells.csv"
        )
    summary["cells"].update(cell_info)
    cov_max_cell = summary["cells"]["Genome Coverage rate of max cell"]
    cpgs_max_cell = summary["cells"]["CPGs of max cell"]
    reads_max_cell = summary["cells"]["Reads of max cell"]
    saturation_max_cell = summary["cells"]["Saturation of max cell"]
    cov_median_cell = summary["cells"]["Genome Coverage rate of median cell"]
    cpgs_median_cell = summary["cells"]["CPGs of median cell"]
    reads_median_cell = summary["cells"]["Reads of median cell"]
    saturation_median_cell = summary["cells"]["Saturation of median cell"]
    
    cov_max_cell = f'{summary["cells"]["Genome Coverage rate of max cell"]:.2%}'
    cpgs_max_cell = f'{summary["cells"]["CPGs of max cell"]}'
    
    reads_max_cell = f'{int(summary["cells"]["Reads of max cell"] / 2)}'
    saturation_max_cell = f'{summary["cells"]["Saturation of max cell"]:.2%}'
    cov_median_cell = f'{summary["cells"]["Genome Coverage rate of median cell"]:.2%}'
    cpgs_median_cell = f'{summary["cells"]["CPGs of median cell"]}'
    
    reads_median_cell = f'{int(summary["cells"]["Reads of median cell"] / 2)}'
    saturation_median_cell = f'{summary["cells"]["Saturation of median cell"]:.2%}'
    cellnum = f'{summary["cells"]["Estimated Number of Cells"]}'
    summary["cells"]["Fraction Reads in Cells"] = summary["cells"]["Reads in Cells"] / summary["mapping"]["uniquereads"]
    fraction = f'{summary["cells"]["Fraction Reads in Cells"]:.2%}' 
    with open(summary_json, 'w') as fh:
        json.dump(summary, fh, indent=4)

    header=('Samplename,Estimated_Number_of_Cells,Number_of_Reads,Valid_Barcode_Ratio,Dropped_Too_Short,Dropped_Chimeric,Valid_7F_Reads_Rate,Valid_17LME_Reads_Rate,Valid_7F17LME_Reads_Rate,C-T_Conversion,C-C_Ratio,'
            'Reads_Mapped_to_Genome,Reads_Mapped_Confidently_to_Genome,CpG_Methylation_Rate,CHG_Methylation_Rate,CHH_Methylation_Rate,Unknown_Methylation_Rate,CpG_Coverage_Rate,'
            'Total_CPGs_Detected,Genome_Coverage_Rate_of_Max_Cell,CPGs_of_Max_Cell,Reads_of_Max_Cell,Saturation_of_Max_Cell,'
            'Genome_Coverage_Rate_of_Median_Cell,CPGs_of_Median_Cell,Reads_of_Median_Cell,Saturation_of_Median_Cell,Fraction_Reads_in_Cells')

    summary_data = [
             samplename,
             cellnum,
             raw,
             vaild,
             dropped_too_short,
             dropped_chimeric,
             rate_7f,
             rate_17lme,
             rate_7f17lme,
             conversion,
             cc_ratio,
             mapgenome,
             confidently,
             cpg_methylation,
             chg_methylation,
             chh_methylation,
             unknown_methylation,
             cpg_coverage,
             total_cpgs,
             cov_max_cell,
             cpgs_max_cell,
             reads_max_cell,
             saturation_max_cell,
             cov_median_cell,
             cpgs_median_cell,
             reads_median_cell,
             saturation_median_cell,
             fraction
           ]

    with open(os.path.join(outdir, f'{samplename}_wgs_summary.csv'), 'w') as fh:
        fh.write(header + '\n')
        fh.write(','.join(str(_).replace(',', '') for _ in summary_data)+ '\n')
        
if __name__ == "__main__":
    outcsv()
