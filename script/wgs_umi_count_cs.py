import click
import os
import pandas as pd
import numpy as np
import gzip
import json
from collections import defaultdict

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
    return cb_list


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--infile', required=True,help='input, step3/counts.xls.')
@click.option('--rawcsv', required=True,help='output.')
@click.option('--filtercsv', required=True,help='output.')
@click.option('--gexcb', required=True, help='input.')
@click.option('--cbcsv', required=True, help='input.')
@click.option('--outdir', required=True, help='output.')
@click.option('--samplename', required=True, help='output.')

def count(infile, rawcsv, filtercsv, gexcb, cbcsv, outdir, samplename):
    d = {}
    print("开始读取文件...")
    with open(infile,"rt") as split_file:
        for line in split_file:
            if line.startswith("cellID"): continue
            ls = line.strip().split("\t")
            barcode = ls[0]
            gene = ls[1]
            umi = ls[2]
            reads = ls[3]
            d[barcode] = d.get(barcode, {'num_UMI':0, 'num_Reads':0})
            d[barcode]['num_UMI'] += int(umi)
            d[barcode]['num_Reads'] += int(reads)


    print("创建字典完成~")
    print("开始统计...")

    df = pd.DataFrame(columns=['Barcode', 'Umi_num', 'Reads_num'], dtype=np.int16)

    for item in d:
        num_umi = d[item]['num_UMI']
        num_reads = d[item]['num_Reads']
        df = pd.concat([df, pd.DataFrame({'Barcode': [item], 'Umi_num': [num_umi], 'Reads_num': [num_reads]})], ignore_index=True)

    df = df.sort_values(by='Umi_num', ascending=False)
    df.to_csv(rawcsv, sep='\t', index=False)
    print("output step3/raw_umicounts.xls")
    
    gexcblist = cblist(gexcb)
    rawdf = pd.read_csv(cbcsv)
    celldf = rawdf[rawdf['gex_cb'].isin(gexcblist)]
    mcb_list = celldf['m_cb'].tolist()

    # 筛选出UMI大于umicut的行
    df0 = df[df['Barcode'].isin(mcb_list)]
    # 保存为df0
    df0.to_csv(filtercsv, sep='\t', index=False)
    print("output step3/filter_umicounts.xls")

    filter_reads_sum = df0['Reads_num'].sum()
    raw_reads_sum = df['Reads_num'].sum()
    raw_umi_sum = df['Umi_num'].sum()
    umi_median = df0['Umi_num'].median()

    saturation = 1 - (raw_umi_sum / raw_reads_sum)
    fraction = filter_reads_sum / raw_reads_sum

    print(f'Estimated Number of Cells: {len(df0)}')
    print(f'Sequencing Saturation: {saturation:.2%}')
    print(f'Fraction Reads in Cells: {fraction:.2%}')
    print(f'Median UMI Counts per Cell: {int(umi_median)}')
    with open(os.path.join(outdir, samplename+"_summary.json"), "r") as fh:
        summary = json.load(fh)
        rawreads = summary["stat"]["total"]
    print(f'Mean Reads per Cell: {int(rawreads / len(df0))}')
    print("统计完成！")
    with open(os.path.join(outdir, samplename+"_summary.json"), "w") as fhout:
        summary_cells = defaultdict()
        summary_cells["Estimated Number of Cells"] = len(df0)
        summary_cells["Sequencing Saturation"] = saturation
        summary_cells["Fraction Reads in Cells"] = fraction
        summary_cells["Median UMI Counts per Cell"] = int(umi_median)
        summary_cells["Mean Reads per Cell"] = int(rawreads / len(df0))

        summary["cells"] = summary_cells
        json.dump(
            summary,
            fhout,
            indent=4,
            default=lambda o: int(o) if isinstance(o, np.int64) else o
        )
if __name__ == '__main__':
    count()
