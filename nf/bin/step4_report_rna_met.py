#!/usr/bin/env python
import json
import os
from jinja2 import Environment, FileSystemLoader
import numpy as np
import pandas as pd
from loguru import logger
import re
import gzip
from scipy.io import mmread
import click
from typing import List, Dict, Any
from seeksoultools.utils.countUtil import calculate_metrics
from collections import defaultdict

CSV_KEY = {
    'Joint': [
        'Sample ID', 'Estimated number of cells', 'GEX Median genes per cell', 'MET CpG number of median cell',
        'Workflow version','Fastp version','SeekSoulTools version','Bismark version',
        'ALLCools version','Reference','Chemistry','Include introns'],
    'RNA': [
        'Number of Reads',
        'Valid Barcode',
        'Sequencing Saturation',
        'Too Short',
        'Q30 Bases in Barcode',
        'Q30 Bases in UMI',
        'Reads Mapped to Genome',
        'Reads Mapped Confidently to Genome',
        'Reads Mapped to Intergenic Regions',
        'Reads Mapped to Intronic Regions',
        'Reads Mapped to Exonic Regions',
        'Fraction Reads in Cells',
        'Mean Reads per Cell',
        'Median Genes per Cell',
        'Median UMI Counts per Cell',
        'Total Genes Detected'
    ],
    'MET':[
        'Number of Read Pairs',
        'Valid Barcodes',
        'Dropped Too Short',
        'Dropped Chimeric',
        'C-T Conversion',
        'C-C Ratio',
        'Read Pairs Mapped to Genome',
        'Read Pairs Mapped Confidently to Genome',
        'CpG Methylation Rate',
        'CHG Methylation Rate',
        'CHH Methylation Rate',
        'CpG Coverage Rate',
        'Total CpGs Detected',
        'Genome Coverage Rate of Median Cell',
        'Read Pairs of Median Cell',
        'Saturation of Median Cell',
        'Fraction Read Pairs in Cells'
    ]
}

PREFIX = {
    'Joint': '',
    'RNA':'GEX_',
    'MET':'MET_'
}

FIX_KEY = {
    'Joint': ['Estimated number of cells', 'GEX Median genes per cell', 'MET CpG number of median cell'],
    'RNA': ['Number of Reads', 'Too Short', 'Mean Reads per Cell', 'Median Genes per Cell',
            'Median UMI Counts per Cell', 'Total Genes Detected'],
    'MET': ['Number of Read Pairs','Total CpGs Detected', 'Read Pairs of Median Cell']
}

def find_key_in_nested_structure(data, target_key):
    """
    Find the value of a key in a nested dictionary or list.
    """
    if isinstance(data, dict):
        if target_key in data:
            return data[target_key]
        
        for value in data.values():
            result = find_key_in_nested_structure(value, target_key)
            if result is not None:
                return result
                
    elif isinstance(data, list):
        for item in data:
            result = find_key_in_nested_structure(item, target_key)
            if result is not None:
                return result
                
    return 

def fix_key_number(n_with_comma: str) -> int:
    """
    Convert a number string with a comma to an integer.
    """
    return int(n_with_comma.replace(',', ''))

def check_software_version() -> dict:
    sft_version = {}
    sft_version['Fastp'] = os.popen('fastp --version').read().strip().split(' ')[1]
    sft_version['SeekSoulTools'] = os.popen('seeksoultools --version').read().strip()
    sft_version['Bismark'] = re.search(
        r'Bismark Version:\s*(v[\d.]+)', 
        os.popen('bismark --version').read()).group(1)
    sft_version['Bowtie2'] = re.search(
        r'version\s*([\d.]+)', 
        os.popen('bowtie2 --version').read()).group(1)
    sft_version['ALLCools'] = os.popen('allcools --version').read().strip()
    return(sft_version)

def get_workflow_version_regex(config_path='nextflow.config'):
    try:
        with open(config_path, 'r') as f:
            content = f.read()
        
        # Match the version inside the manifest block
        pattern = r'manifest\s*\{[^}]*version\s*=\s*[\'"]([^\'"]+)[\'"]'
        match = re.search(pattern, content, re.DOTALL)
        if match:
            return match.group(1)
    except Exception as e:
        print(f"Failed to read config file: {e}")
    return ""

def try_get_value(dict_obj: dict, key: str, default: Any = None) -> Any:
    """
    Safely get a value from a dictionary. If the key is not found, return the default value.
    """
    return dict_obj.get(key, default)
   
def merge_rna_methylation_by_barcode(tsne_file: str, filtered_counts_file: str, cells_file: str, output_file: str, whitelist_file: str = None) -> pd.DataFrame:
    """
    Merge RNA and methylation tables into a single file based on barcode mapping.
    - tsne_file: GEX tSNE file (e.g., tsne_umi.xls; first column is the RNA barcode index)
    - filtered_counts_file: methylation barcode counts file (filtered_barcode_reads_counts.csv; columns: reads_counts, barcode)
    - cells_file: methylation cell metrics file (e.g., GYY2155663_cells.csv; includes cell_barcode where it looks like <m_cb>_allc.gz)
    - whitelist_file: whitelist mapping file (bUCB3_whitelist.csv; columns: gex_cb (RNA), m_cb (methylation))
    - output_file: path to the merged CSV output
    """
    # Read tSNE file (RNA barcode used as index)
    try:
        tsne_df = pd.read_table(tsne_file, index_col=0)
    except Exception as e:
        logger.error(f"Failed to read tSNE file: {tsne_file} -> {e}")
        raise
    # Use index as the RNA barcode column
    tsne_df = tsne_df.copy()
    if whitelist_file:
        tsne_df['gex_cb'] = tsne_df.index.astype(str)
    else:
        tsne_df['m_cb'] = tsne_df.index.astype(str)

    # Read filtered barcode counts (methylation barcodes)
    try:
        counts_df = pd.read_csv(filtered_counts_file)
    except Exception as e:
        logger.error(f"Failed to read methylation counts file: {filtered_counts_file} -> {e}")
        raise
    if 'barcode' not in counts_df.columns:
        raise ValueError(f"Missing 'barcode' column in {filtered_counts_file}")
    counts_df = counts_df.rename(columns={'barcode': 'm_cb'})

    # Read cells (methylation barcode is in cell_barcode; strip suffix)
    try:
        cells_df = pd.read_csv(cells_file)
    except Exception as e:
        logger.error(f"Failed to read cells file: {cells_file} -> {e}")
        raise
    if 'cell_barcode' not in cells_df.columns:
        raise ValueError(f"Missing 'cell_barcode' column in {cells_file}")
    # Remove suffix "_allc.gz" or possible ".gz"
    cells_df = cells_df.copy()
    cells_df['m_cb'] = (
        cells_df['cell_barcode']
        .astype(str)
        .str.replace('_allc.gz', '', regex=False)
        .str.replace('.gz', '', regex=False)
    )
    tsne_map_df = tsne_df.reset_index(drop=True)
    
    # Read whitelist (map gex_cb to m_cb)
    if whitelist_file:
        wl_df = pd.read_csv(whitelist_file)
        for col in ('gex_cb', 'm_cb'):
            if col not in wl_df.columns:
                raise ValueError(f"Whitelist file {whitelist_file} is missing column: {col}")
        # Deduplicate, keep the first mapping per gex_cb
        wl_df = wl_df[['gex_cb', 'm_cb']].drop_duplicates(subset=['gex_cb'])
        tsne_map_df = tsne_map_df.merge(wl_df, on='gex_cb', how='inner')
    # Merge counts by m_cb
    merged_df = tsne_map_df.merge(counts_df, on='m_cb', how='inner') 
    
    # Merge cells by m_cb; add suffix to avoid column conflicts
    merged_df = merged_df.merge(cells_df, on='m_cb', how='inner', suffixes=('', '_cells'))
    
    # comput cpg methylation levels and ch methylation levels
    cg_cov_col = [i for i in merged_df.columns[merged_df.columns.str.startswith('CG')].to_list() if i.endswith('_cov')]
    cg_mc_col = [i for i in merged_df.columns[merged_df.columns.str.startswith('CG')].to_list() if i.endswith('_mc')]
    ch_cov_col = [i for i in merged_df.columns[merged_df.columns.str.contains('^(CT|CC|CA)')].to_list() if i.endswith('_cov')]
    ch_mc_col = [i for i in merged_df.columns[merged_df.columns.str.contains('^(CT|CC|CA)')].to_list() if i.endswith('_mc')]

    merged_df['CpG_cov'] = merged_df[cg_cov_col].sum(axis=1)
    merged_df['CpG_mc'] = merged_df[cg_mc_col].sum(axis=1)
    merged_df['CpG%'] = round(merged_df['CpG_mc'] * 100 / merged_df['CpG_cov'], 2)
    
    merged_df['CH_cov'] = merged_df[ch_cov_col].sum(axis=1)
    merged_df['CH_mc'] = merged_df[ch_mc_col].sum(axis=1)
    merged_df['CH%'] = round(merged_df['CH_mc'] * 100 / merged_df['CH_cov'], 2)

    # Output
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        merged_df.to_csv(output_file, index=True, sep = '\t')
        logger.info(
            f"Merge completed, output file: {output_file} | tSNE rows: {len(tsne_df)} | counts matched: {merged_df['reads_counts'].notna().sum()} | cells matched: {merged_df['cell_barcode'].notna().sum()}"
        )
    except Exception as e:
        logger.error(f"Failed to write merged result: {output_file} -> {e}")
        raise
    return(merged_df)

def pre_barcode_rank_data(raw_dir, filtered_dir) -> list:
    with gzip.open(f"{raw_dir}/barcodes.tsv.gz", "rt") as fh:
        barcodes = [line.strip() for line in fh]

    with gzip.open(f"{filtered_dir}/barcodes.tsv.gz", "rt") as fh:
        cells = {line.strip() for line in fh if line.strip()}

    with gzip.open(f"{raw_dir}/matrix.mtx.gz") as fh:
        m = mmread(fh)

    umi_sum = m.sum(axis=0)

    df = pd.DataFrame({"barcode": barcodes, "UMI": umi_sum.A[0]})
    df["is_cell"] = False
    df.loc[df["barcode"].isin(cells), "is_cell"] = True
    df = (df.sort_values(["UMI", "is_cell"], ascending=[False, False])
        .reset_index(drop=True)
        .reset_index(names="idx"))
    return(df)

def barcode_rank_data(
    pre_barcode_rank_data_df: pd.DataFrame,
    rna_met_df: pd.DataFrame
) -> List[Dict[str, Any]]:
    """
    Segment UMI-ranked data by adjacent barcode categories (RNA+DNA, RNA-only, Background)
    to generate segment lists for plotting.

    Parameters:
    - pre_barcode_rank_data_df: DataFrame with columns [idx, barcode, UMI, is_cell], already sorted by UMI descending; idx is the rank index.
    - rna_met_df: DataFrame with column [gex_cb], defining the set of "RNA+DNA" barcodes.
    - tail_background_chunk_size: when the last segment is Background, split it by this size (default 1000).

    Returns:
    - List[Dict], each dict like {x: int, y: List[float], text: str, color: str}
      where:
        x: minimum idx of the segment (segment start)
        y: list of UMIs within the segment (original order)
        text: segment type ("RNA+DNA" | "RNA-only" | "Background")
        color: color mapping (RNA+DNA -> orange; RNA-only -> blue; Background -> blue)
    """
    # Validate required columns
    required_pre = {"idx", "barcode", "UMI", "is_cell"}
    required_rna = {"gex_cb"}
    missing_pre = required_pre - set(pre_barcode_rank_data_df.columns)
    missing_rna = required_rna - set(rna_met_df.columns)
    if missing_pre:
        raise ValueError(f"pre_barcode_rank_data_df missing columns: {sorted(missing_pre)}")
    if missing_rna:
        raise ValueError(f"rna_met_df missing columns: {sorted(missing_rna)}")

    # Preserve order: sort by idx ascending (ensure consistency with UMI ranking)
    df = pre_barcode_rank_data_df.copy()
    df = df.sort_values(by=["idx"], ascending=True)

    # Set of RNA+DNA barcodes
    rna_set = set(rna_met_df["gex_cb"].astype(str))

    # Color mapping
    color_map = {
        "RNA+DNA": "rgba(80, 80, 201, 1.0)",
        "RNA-only": "rgba(255, 211, 26, 0.5)",
        "Background": "rgba(221, 221, 221, 1.0)",
    }

    # Annotate cells_status for each row
    def _status(row) -> str:
        b = str(row["barcode"])  # Prevent type mismatch
        if b in rna_set:
            return "RNA+DNA"
        elif bool(row["is_cell"]):
            return "RNA-only"
        else:
            return "Background"

    df = df.assign(cells_status=df.apply(_status, axis=1))

    # Split adjacent rows with the same status into segments (keep idx and UMI lists)
    segments_data: List[Dict[str, Any]] = []  # Temporary: {"status": str, "idxs": List[int], "umis": List[float]}
    current = {"status": None, "idxs": [], "umis": []}

    for _, row in df.iterrows():
        s = row["cells_status"]
        if current["status"] is None:
            current["status"] = s
            current["idxs"] = [int(row["idx"])]
            current["umis"] = [row["UMI"]]
        elif s == current["status"]:
            current["idxs"].append(int(row["idx"]))
            current["umis"].append(row["UMI"])
        else:
            segments_data.append(current)
            current = {
                "status": s,
                "idxs": [int(row["idx"])],
                "umis": [row["UMI"]],
            }
    if current["status"] is not None:
        segments_data.append(current)
    
    # Convert to final output format
    result: List[Dict[str, Any]] = []
    for seg in segments_data:
        status = seg["status"]
        idxs = seg["idxs"]
        umis = seg["umis"]
        result.append({
            "x": int(min(idxs) if idxs else 0),
            "y": list(umis),
            "text": status,
            "color": color_map.get(status, "blue"),
        })

    return result


def pre_diff_data(diff_table, n=50):
    df = pd.read_table(diff_table)
    tmp = df.groupby('cluster').apply(lambda x: x.sort_values(by='avg_log2FC', ascending=False).head(n))
    diff_data = tmp.loc[:,['Ensembl', 'gene', 'avg_log2FC', 'p_val_adj','cluster']].to_dict('records')
    return diff_data

def get_gex_tsne(tsnefile):
    df = pd.read_table(tsnefile, index_col = 0, header = 0)
    # df = df.drop(['orig.ident', 'nFeature_RNA', 'percent.mito'], axis=1)
    tsne_d = df.to_dict(orient='list')
    tsne1 = tsne_d['tSNE_1']
    tsne2 = tsne_d['tSNE_2']
    nCount_RNA = tsne_d['nCount_RNA']
    nFeature_RNA = tsne_d['nFeature_RNA']
    mito = tsne_d['percent.mito']
    cluster = tsne_d['RNA_snn_res.0.8']
    Total_CpG_number = tsne_d['total_cpg_number']
    CpG_methylation_level = tsne_d['CpG%']
    CH_methylation_level = tsne_d['CH%']
    Genome_coverage = tsne_d['genome_cov']
    return tsne1, tsne2, nCount_RNA, nFeature_RNA,mito, cluster, Total_CpG_number, CpG_methylation_level, CH_methylation_level, Genome_coverage

@click.command()
@click.option('--gexjson', required=True, help='Path to the GEX JSON file.')
@click.option('--metjson', required=True, help='Path to the MET JSON file.')
@click.option('--tsne_file', required=True, help='Path to the tSNE file.')
@click.option('--filtered_counts_file', required=True, help='Path to the filtered counts file.')
@click.option('--cells_file', required=True, help='Path to the cells file.')
@click.option('--whitelist_file', required=False, help='Path to the whitelist file.')
@click.option('--outdir', required=True, help='Path to the output directory.')
@click.option('--samplename', required=True, help='Sample name.')
@click.option('--rawname', required=True, help='Raw name.')
@click.option('--nf_config', required=True, help='Path to the Nextflow config file.')
@click.option('--raw_dir', required=True, help='Path to the raw directory.')
@click.option('--filtered_dir', required=True, help='Path to the filtered directory.')
@click.option('--diff_data', required=True, help='Path to the diff data file.')
@click.option('--gtf', required=True, help='Path to the GTF file.')
@click.option('--counts_file', required=True, help='Path to the counts file.')
@click.option('--detail_file', required=True, help='Path to the detail file.')
def report(gexjson, metjson, 
           tsne_file, filtered_counts_file, cells_file, 
           outdir, samplename, rawname, nf_config, 
           raw_dir, filtered_dir, diff_data, gtf, counts_file, detail_file, 
           whitelist_file=None,
           **kwargs):
    os.makedirs(outdir, exist_ok=True)
    datajson = os.path.join(os.path.dirname(__file__), './utils/report_rna_met/sgrnamet.json')
    assert os.path.exists(datajson), f'{datajson} not found!'
    
    software_version = check_software_version()
    
    with open(gexjson) as gfh:
        gex_summary = json.load(gfh)
    with open(metjson) as mfh:
        met_summary = json.load(mfh)
    with open(datajson) as dfh:
        data_summary = json.load(dfh)
    if rawname == "rawname":
        rawname = samplename
    df = merge_rna_methylation_by_barcode(
        tsne_file = tsne_file, 
        filtered_counts_file = filtered_counts_file,
        cells_file = cells_file,
        whitelist_file = whitelist_file,
        output_file = os.path.join(outdir, f'{samplename}_rna_met_merge.xls'),
    )    
    # joint: title
    data_summary["Joint"][0]["left"][0]["data"]["Estimated number of cells"] = f'{df.shape[0]:,}'
    data_summary["Joint"][0]["right"][0]["data"]["GEX Median genes per cell"] = f'{int(df["nFeature_RNA"].median()):,}'
    # joint: left_Sample riget_Joint Metrics
    data_summary["Joint"][1]["left"][0]["data"]["Sample ID"] = samplename
    data_summary["Joint"][1]["left"][0]["data"]["Sample description"] = rawname
    data_summary["Joint"][1]["left"][0]["data"]["Workflow version"] = get_workflow_version_regex(nf_config)
    data_summary["Joint"][1]["left"][0]["data"]["Fastp version"] = software_version['Fastp']
    data_summary["Joint"][1]["left"][0]["data"]["SeekSoulTools version"] = software_version['SeekSoulTools']
    data_summary["Joint"][1]["left"][0]["data"]["Bismark version"] = software_version['Bismark']
    data_summary["Joint"][1]["left"][0]["data"]["ALLCools version"] = software_version['ALLCools']
    data_summary["Joint"][1]["left"][0]["data"]["Reference"] = gex_summary['reference']
    data_summary["Joint"][1]["left"][0]["data"]["Chemistry"] = met_summary['stat']['chemistry']
    data_summary["Joint"][1]["left"][0]["data"]["Include introns"] = gex_summary['include_introns']

    gex_tsnefile = os.path.join(outdir, f'{samplename}_rna_met_merge.xls')
    if not os.path.exists(gex_tsnefile):
        logger.info(f"Warning : The path of '{gex_tsnefile}' is not exists.")
        gex_tsne1 = []
        gex_tsne2 = []
        gex_nCount_RNA = []
        gex_nFeature_RNA = []
        gex_cluster = []
    else:
        gex_tsne1, gex_tsne2, gex_nCount_RNA, \
        gex_nFeature_RNA, mito, gex_cluster, MET_Total_CpG_number, \
        CpG_methylation_level, CH_methylation_level, Genome_coverage = get_gex_tsne(gex_tsnefile)
    
    data_summary["Joint"][1]["right"][0]["data"]["data"] = barcode_rank_data(
        pre_barcode_rank_data_df = pre_barcode_rank_data(raw_dir, filtered_dir),
        rna_met_df = df
    )
    
    data_summary["tsne"]["data"]["coordinate"]["tSNE1"] = gex_tsne1
    data_summary["tsne"]["data"]["coordinate"]["tSNE2"] = gex_tsne2
    data_summary["tsne"]["data"]["nCount_RNA"] = gex_nCount_RNA
    data_summary["tsne"]["data"]["nFeature_RNA"] = gex_nFeature_RNA
    data_summary["tsne"]["data"]["mito"] = [ round(i,2) for i in mito ]
    data_summary["tsne"]["data"]["cluster"] = gex_cluster
    data_summary["tsne"]["data"]["CpG_number"] = [ int(i) for i in MET_Total_CpG_number ]
    CpG_methylation_level = [ 0 if np.isnan(i) else i for i in CpG_methylation_level ]
    CH_methylation_level = [ 0 if np.isnan(i) else i for i in CH_methylation_level ]
    data_summary["tsne"]["data"]["CpG_methylation"] = [ round(i,2) for i in CpG_methylation_level ]
    data_summary["tsne"]["data"]["CH_methylation"] = [ round(i,2) for i in CH_methylation_level ]
    
    data_summary["tsne"]["data"]["range"]["nCount_RNA"] = [ int(min(gex_nCount_RNA)), int(max(gex_nCount_RNA)) ]
    data_summary["tsne"]["data"]["range"]["nFeature_RNA"] = [ int(min(gex_nFeature_RNA)), int(max(gex_nFeature_RNA)) ]
    data_summary["tsne"]["data"]["range"]["mito"] = [ round(min(mito),2), round(max(mito),2) ]
    data_summary["tsne"]["data"]["range"]["cluster"] = [ str(i) for i in list(set(sorted(gex_cluster))) ]
    
    data_summary["tsne"]["data"]["range"]["CpG_number"] = [ int(min(MET_Total_CpG_number)), int(max(MET_Total_CpG_number)) ]
    data_summary["tsne"]["data"]["range"]["CpG_methylation"] = [ round(min(CpG_methylation_level),2), round(max(CpG_methylation_level),2) ]
    data_summary["tsne"]["data"]["range"]["CH_methylation"] = [ round(min(CH_methylation_level),2), round(max(CH_methylation_level),2) ]

    # rna: title
    data_summary["RNA"][0]["left"][0]["data"]["Estimated number of cells"] = f'{df.shape[0]:,}'
    data_summary["RNA"][0]["right"][0]["data"]["GEX Median genes per cell"] = f'{int(df["nFeature_RNA"].median()):,}'
    
    
    # rna: left_Sequencing riget_Saturation
    data_summary["RNA"][1]["left"][0]["data"]["Number of Reads"] = f'{gex_summary["stat"]["total"]:,}'
    data_summary["RNA"][1]["left"][0]["data"]["Valid Barcode"] = f'{gex_summary["stat"]["valid"]/gex_summary["stat"]["total"]:.2%}'
    data_summary["RNA"][1]["left"][0]["data"]["Sequencing Saturation"] = f'{gex_summary["cells"]["Sequencing Saturation"]:.2%}'
    data_summary["RNA"][1]["left"][0]["data"]["Too Short"] = f'{gex_summary["stat"]["too_short"]:,}'
    b_total_base = sum([sum(v) for v in gex_summary["barcode_q"].values()])
    b30_base = sum([sum(v[30:]) for v in gex_summary["barcode_q"].values()])
    data_summary["RNA"][1]["left"][0]["data"]["Q30 Bases in Barcode"] = f'{b30_base/b_total_base:.2%}'
    u_total_base = sum([sum(v) for v in gex_summary["umi_q"].values()])
    u30_base = sum([sum(v[30:]) for v in gex_summary["umi_q"].values()])
    data_summary["RNA"][1]["left"][0]["data"]["Q30 Bases in UMI"] = f'{u30_base/u_total_base:.2%}'
    data_summary["RNA"][1]["right"][0]["data"]["Reads Mapped to Genome"] = f'{gex_summary["mapping"]["Reads Mapped to Genome"]:.2%}'
    data_summary["RNA"][1]["right"][0]["data"]["Reads Mapped Confidently to Genome"] = f'{gex_summary["mapping"]["Reads Mapped Confidently to Genome"]:.2%}'
    data_summary["RNA"][1]["right"][0]["data"]["Reads Mapped to Intergenic Regions"] = f'{gex_summary["mapping"]["Reads Mapped to Intergenic Regions"]:.2%}'
    data_summary["RNA"][1]["right"][0]["data"]["Reads Mapped to Intronic Regions"] = f'{gex_summary["mapping"]["Reads Mapped to Intronic Regions"]:.2%}'
    data_summary["RNA"][1]["right"][0]["data"]["Reads Mapped to Exonic Regions"] = f'{gex_summary["mapping"]["Reads Mapped to Exonic Regions"]:.2%}'
    
    df["gex_cb"].to_csv(os.path.join(outdir, f'tmp_filtered_barcode.xls'), header = False, index = False, sep = "\t")
    gex_summary_cells_update, downsample_dict = calculate_metrics(
        counts_file = counts_file,
        detail_file = detail_file,
        filterd_barcodes_file = os.path.join(outdir, f'tmp_filtered_barcode.xls'),
        gtf = gtf,
        basedir = os.path.join(outdir),
        filterd_features_file = None
    )
    os.remove(os.path.join(outdir, f'tmp_filtered_barcode.xls'))
    data_summary["RNA"][2]["left"][0]["data"]["Estimated Number of Cells"] = f'{df.shape[0]:,}'
    data_summary["RNA"][2]["left"][0]["data"]["Fraction Reads in Cells"] = f'{gex_summary_cells_update["Fraction Reads in Cells"]:.2%}'
    data_summary["RNA"][2]["left"][0]["data"]["Mean Reads per Cell"] = f'{int(gex_summary["stat"]["total"]/int(df.shape[0])):,}'
    data_summary["RNA"][2]["left"][0]["data"]["Median Genes per Cell"] = f'{int(gex_summary_cells_update["Median Genes per Cell"]):,}'
    data_summary["RNA"][2]["left"][0]["data"]["Median UMI Counts per Cell"] = f'{int(gex_summary_cells_update["Median UMI Counts per Cell"]):,}'
    data_summary["RNA"][2]["left"][0]["data"]["Total Genes Detected"] = f'{gex_summary_cells_update["Total Genes Detected"]:,}'
    
    data_summary["RNA"][3]["left"][0]["data"]["x"] = [0, ] + gex_summary["downsample"]["Reads"]
    data_summary["RNA"][3]["left"][0]["data"]["y"] = [0, ] + gex_summary["downsample"]["median"]
    data_summary["RNA"][3]["right"][0]["data"]["x"] = [0, ] + gex_summary["downsample"]["Reads"]
    data_summary["RNA"][3]["right"][0]["data"]["y"] = [0, ] + gex_summary["downsample"]["saturation"]
    
    # diff
    diff_table = diff_data
    if not os.path.exists(diff_table):
        logger.info(f"Warning : The path of '{diff_table}' is not exists.")
    else:
        diff_data = pre_diff_data(diff_table)
        data_summary["diff"] = diff_data

    
    data_summary["MET"][0]["left"][0]["data"]["Estimated number of cells"] = f'{df.shape[0]:,}'
    data_summary["MET"][0]["right"][0]["data"]["GEX Median genes per cell"] = f'{int(df["nFeature_RNA"].median()):,}'
    
    # met: left_Sequencing riget_median_fragments
    data_summary["MET"][1]["left"][0]["data"]["Number of Read Pairs"] = f'{met_summary["stat"]["total"]:,}'
    data_summary["MET"][1]["left"][0]["data"]["Valid Barcodes"] = f'{met_summary["stat"]["valid"]/met_summary["stat"]["total"]:.2%}'
    data_summary["MET"][1]["left"][0]["data"]["Dropped Too Short"] = f'{met_summary["stat"]["too_short"] / met_summary["stat"]["total"]:.2%}'
    
    forward_chimeric = try_get_value(met_summary["stat"], "forward_chimeric", 0)
    reverse_chimeric = try_get_value(met_summary["stat"], "reverse_chimeric", 0)
    too_short = try_get_value(met_summary["stat"], "too_short", 0)
    vaildreads = try_get_value(met_summary["stat"], "valid", 0)
    data_summary["MET"][1]["left"][0]["data"]["Dropped Chimeric"] = f'{(forward_chimeric + reverse_chimeric) / (vaildreads - too_short):.2%}'
    data_summary["MET"][1]["left"][0]["data"]["C-T Conversion"] = f'{met_summary["stat"]["ct_mean"]:.2%}'
    data_summary["MET"][1]["left"][0]["data"]["C-C Ratio"] = f'{met_summary["stat"]["cc_mean"]:.2%}'
    data_summary["MET"][1]["right"][0]["data"]["Read Pairs Mapped to Genome"] = f'{met_summary["mapping"]["Reads Mapped to Genome"]:.2%}'
    data_summary["MET"][1]["right"][0]["data"]["Read Pairs Mapped Confidently to Genome"] = f'{met_summary["mapping"]["Reads Mapped Confidently to Genome"]:.2%}'
    data_summary["MET"][1]["right"][0]["data"]["CpG Methylation Rate"] = f'{round(met_summary["cpg_methylation_rate"],2)}%'
    data_summary["MET"][1]["right"][0]["data"]["CHG Methylation Rate"] = f'{round(met_summary["chg_methylation_rate"],2)}%'
    data_summary["MET"][1]["right"][0]["data"]["CHH Methylation Rate"] = f'{round(met_summary["chh_methylation_rate"],2)}%'
    data_summary["MET"][1]["right"][0]["data"]["CpG Coverage Rate"] = f'{met_summary["cells"]["CpG Coverage rate"]:.2%}'
    data_summary["MET"][1]["right"][0]["data"]["Total CpGs Detected"] = f'{met_summary["cells"]["Total CPGs Detected"]:,}'
    
    data_summary["MET"][2]["left"][0]["data"]["Estimated Number of Cells"] = f'{met_summary["cells"]["Estimated Number of Cells"]:,}'
    data_summary["MET"][2]["left"][0]["data"]["Genome Coverage Rate of Median Cell"] = f'{met_summary["cells"]["Genome Coverage rate of median cell"]:.2%}'
    data_summary["MET"][2]["left"][0]["data"]["CpGs of Median Cell"] = f'{int(met_summary["cells"]["CPGs of median cell"]):,}'
    data_summary["MET"][2]["left"][0]["data"]["Read Pairs of Median Cell"] = f'{int(met_summary["cells"]["Reads of median cell"] / 2):,}'
    data_summary["MET"][2]["left"][0]["data"]["Saturation of Median Cell"] = f'{met_summary["cells"]["Saturation of median cell"]:.2%}'
    data_summary["MET"][2]["left"][0]["data"]["Fraction Read Pairs in Cells"] = f'{met_summary["cells"]["Fraction Reads in Cells"]:.2%}' 
    data_summary["MET"][2]["right"][0]["data"]["y"]['Cell Saturation'] = [ round(i*100, 2) for i in df["cell_saturation"].tolist()]
    data_summary["MET"][2]["right"][0]["data"]["y"]['Cell Genome Coverage'] = [ round(i*100,2) for i in Genome_coverage ]
    
    data_summary["MET"][3]["left"][0]["data"]["y"]['CG methylation'] = [ round(i,2) for i in CpG_methylation_level ]
    data_summary["MET"][3]["left"][0]["data"]["y"]['CH methylation'] = [ round(i,2) for i in CH_methylation_level ]
    data_summary["MET"][3]["right"][0]["data"]["x"] = np.log10((df["reads_counts"] / 2) + 1).round(4).tolist()
    data_summary["MET"][3]["right"][0]["data"]["y"] = np.log10(df["total_cpg_number"] + 1).round(4).tolist()
    
    data_summary["Joint"][0]["right"][0]["data"]["MET CpG number of median cell"] = f'{int(met_summary["cells"]["CPGs of median cell"]):,}'
    data_summary["RNA"][0]["right"][0]["data"]["MET CpG number of median cell"] = f'{int(met_summary["cells"]["CPGs of median cell"]):,}'
    data_summary["MET"][0]["right"][0]["data"]["MET CpG number of median cell"] = f'{int(met_summary["cells"]["CPGs of median cell"]):,}'
    
    
    
    with open(os.path.join(outdir, f'{samplename}_rna_met.json'), 'w') as fh:
        json.dump(data_summary, fh, indent = 4)
    # report
    template_dir_new = os.path.abspath(os.path.join(os.path.dirname(__file__), './utils/report_rna_met'))
    env = Environment(loader=FileSystemLoader(template_dir_new))
    template = env.get_template('base.html')
    with open(os.path.join(outdir, f'{samplename}_rna_methyl_report.html'), 'w') as fh:
        fh.write(template.render(websummary_json_data=json.dumps(data_summary).replace("5'", "5\\'").replace("3'", "3\\'")))
        
    rna_methy_summary_combine = defaultdict(dict)
    for k,v in CSV_KEY.items():
        for i in v:
            i_with_no_space = i.replace(" ","_")
            rna_methy_summary_combine[f'{PREFIX[k]}{i_with_no_space}'] = find_key_in_nested_structure(data_summary[k], i)
            if i in FIX_KEY[k]:
                rna_methy_summary_combine[f'{PREFIX[k]}{i_with_no_space}'] = fix_key_number(rna_methy_summary_combine[f'{PREFIX[k]}{i_with_no_space}'])
    pd.DataFrame([rna_methy_summary_combine]).to_csv(os.path.join(outdir, f'{samplename}_rna_methy_summary.csv'), index=False)



if __name__ == "__main__":
    report()