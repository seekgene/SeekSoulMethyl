import click
import json
import os
import gzip
import pandas as pd
import numpy as np
import subprocess
from collections import defaultdict
from scipy.sparse import coo_matrix
from scipy.io import mmwrite

def cblist(cbfile):
    cb_list = []
    if cbfile.endswith(".gz"):
        _open = gzip.open
    else:
        _open = open
    with _open(cbfile, "rt") as f:
        for line in f:
            cb = line.strip()
            cb_list.append(cb)
    return set(cb_list)

def ml_dict(infile):
    n = 0
    total = 0
    with open(infile, 'r') as fh:
        for line in fh:
            #lines = line.strip().split('\t')
            #if not lines[3].startswith('CG'): continue
            total += 1

    return total

def mk_dict(infile):
    d = {}
    n = 0
    
    df = pd.DataFrame(columns=['sizeID', 'type'], dtype=np.int16)

    with open(infile, 'r') as fh:
        genelist = [line.strip().split('\t')[-1] for line in fh]
        unique_genelist = list(set(genelist))
    
    for line in unique_genelist:
        sizename = line.strip()
        d[sizename] = n
        n += 1
        # sizelist.append(sizename)
        df = pd.concat([df, pd.DataFrame({'sizeID': [sizename], 'type': ['ML Expression']})], ignore_index=True)
    
    return n, d, df

def merged(cbmlfile):
    df = pd.read_csv(cbmlfile, delimiter='\t', names=['ml', 'gene'])
    df_merged = df.groupby('gene')['ml'].mean().reset_index()
    # print(df_merged)
    gene_ml_dict = df_merged.set_index('gene')['ml'].to_dict()
    
    return gene_ml_dict

def totalcpg(allcpgfile):
    df = pd.read_csv(allcpgfile, 
                     sep='\t',
                     usecols=[3],  # 只读取第4列
                     compression='gzip' if allcpgfile.endswith('.gz') else None)
    
    total = df[df.iloc[:, 0].str.startswith('CG')].shape[0]
    
    print(f'total CPG number: {total}')
    print(f'file: {allcpgfile}')
    return total

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--cbfile', required=True,help='filter cb list')
@click.option('--outdir', required=True,help='outdir')
@click.option('--samplename', required=True,help='samplename')
@click.option('--outcsv', required=True,help='output csv file')
@click.option('--allcpgfile', required=True,help='bulk cpg file')
@click.option('--indir', required=True,help='directory containing each barcode cpg info')
@click.option('--summary_json', required=True,help='summary json')
def main(cbfile, outdir, samplename, outcsv, allcpgfile, indir, summary_json):
    filter_cblist = cblist(cbfile)
    allcpgnum = totalcpg(allcpgfile=allcpgfile)

    df_cbMC = pd.DataFrame(columns=['cellID', 'CPG_num'], dtype=np.int16)
    for i in filter_cblist:
        cpgpath = f'{indir}/{i}_CpG.bed'
        cpg_num = ml_dict(infile=cpgpath)
        df_cbMC = pd.concat([df_cbMC, pd.DataFrame({'cellID': [i], 'CPG_num': [cpg_num]})], ignore_index=True)

    # cb mC num
    index_l = df_cbMC.index.values.tolist()
    index_l1 = np.array(index_l) + 1
    df_cbMC = df_cbMC.sort_values(by='CPG_num', ascending=False)
    df_cbMC['Orderid'] = index_l1
    df_cbMC.to_csv(outcsv, sep='\t', index=False)

    # 合并3个表格
    cb_genomecov_xls = os.path.join(outdir, samplename+'_cb_genomecov.xls')
    df0 = pd.read_csv(cb_genomecov_xls, sep='\t')

    filter_umi_xls = os.path.join(outdir, 'filter_umicounts.xls')
    df1 = pd.read_csv(filter_umi_xls, sep='\t')

    merged_tmp = pd.merge(df1, df0, on='Barcode')
    merged_df = pd.merge(df_cbMC, merged_tmp, left_on='cellID', right_on='Barcode')

    # 按cov1列由高到底排序
    sorted_df = merged_df.sort_values(by='genomecov_1x', ascending=False)
    sortxls = os.path.join(outdir, 'filter_gcov_umi_sort.xls')
    # sorted_df = pd.read_csv(sortxls, sep='\t')
    sorted_df.to_csv(sortxls, sep='\t', index=False)
    # 输出cov1最高的cb的umi列的值和reads列的值
    cell_max_cov = float(sorted_df.loc[sorted_df['genomecov_1x'].idxmax(), ['genomecov_1x']])
    cell_max_umi = int(sorted_df.loc[sorted_df['genomecov_1x'].idxmax(), ['Umi_num']])
    cell_max_reads = int(sorted_df.loc[sorted_df['genomecov_1x'].idxmax(), ['Reads_num']])
    cell_max_cpgs = int(sorted_df.loc[sorted_df['genomecov_1x'].idxmax(), ['CPG_num']])
    cell_max_saturation = 1 - (cell_max_umi / cell_max_reads)

    median_index = len(sorted_df) // 2
    # 输出对应索引位置的cb、umi和reads列的值
    cell_median_cov = float(sorted_df.iloc[median_index][['genomecov_1x']])
    cell_median_umi = int(sorted_df.iloc[median_index][['Umi_num']])
    cell_median_reads = int(sorted_df.iloc[median_index][['Reads_num']])
    cell_median_cpgs = int(sorted_df.iloc[median_index][['CPG_num']])
    cell_median_saturation = 1 - (cell_median_umi / cell_median_reads)

    with open(summary_json, "r") as fh:
        summary = json.load(fh)

    with open(summary_json, "w") as fhout:
        summary["cells"]["Genome Coverage rate of max cell"] = cell_max_cov
        summary["cells"]["CPGs of max cell"] = cell_max_cpgs
        summary["cells"]["UMIs of max cell"] = cell_max_umi
        summary["cells"]["Reads of max cell"] = cell_max_reads
        summary["cells"]["Saturation of max cell"] = cell_max_saturation
        summary["cells"]["Genome Coverage rate of median cell"] = cell_median_cov
        summary["cells"]["CPGs of median cell"] = cell_median_cpgs
        summary["cells"]["UMIs of median cell"] = cell_median_umi
        summary["cells"]["Reads of median cell"] = cell_median_reads
        summary["cells"]["Saturation of median cell"] = cell_median_saturation
        summary["cells"]["Total CPGs Detected"] = allcpgnum
        json.dump(
            summary,
            fhout,
            indent=4,
            default=lambda o: int(o) if isinstance(o, np.int64) else o
        )


if __name__ == '__main__':
    main()
