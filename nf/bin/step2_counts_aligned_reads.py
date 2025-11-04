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
    统计双端比对文件中每个barcode比对上的reads数目
    
    Args:
        bam: 输入的BAM文件路径
        outdir: 输出目录
    """
    input_bam = pysam.AlignmentFile(bam, "rb")
    barcode_counts = defaultdict(int)
    
    # 遍历BAM文件中的每条reads
    for read in input_bam.fetch(until_eof=True):
        # 只处理主要比对（非辅助比对和非次级比对）
        if read.is_secondary or read.is_supplementary:
            continue
            
        # 从CB标签中获取barcode
        if read.has_tag("CB"):
            barcode = read.get_tag("CB")
            # 增加该barcode的计数
            barcode_counts[barcode] += 1
    
    # 关闭BAM文件
    input_bam.close()
    
    # 确保输出目录存在
    os.makedirs(outdir, exist_ok=True)
    
    # 将结果转换为DataFrame
    df = pd.DataFrame({
        'barcode': list(barcode_counts.keys()),
        'aligned_reads': list(barcode_counts.values())
    })
    
    # 按reads数量降序排序
    df = df.sort_values('aligned_reads', ascending=False)
    
    # 输出结果到文件
    output_file = os.path.join(outdir, os.path.basename(bam).replace('.bam', '_cb_aligned_reads_counts.csv'))
    df.to_csv(output_file, index=False)
    
    logger.info(f"统计完成，共有{len(barcode_counts)}个barcode，结果已保存到{output_file}")
    
    return output_file

@click.command()
@click.option('-b', '--bam', required=True, help='输入的BAM文件路径')
@click.option('-o', '--outdir', required=True, help='输出目录')
def main(bam, outdir):
    """统计双端比对文件中每个barcode比对上的reads数目"""
    logger.info(f"开始处理BAM文件: {bam}")
    output_file = counts_reads(bam, outdir)
    logger.info(f"处理完成，结果保存在: {output_file}")

if __name__ == '__main__':
    main()