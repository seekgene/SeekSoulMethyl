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
    估计细胞数量并过滤barcode
    如果raw_aligned_reads_counts是目录，则处理该目录下所有以_barcode_counts.csv结尾的文件
    并按barcode列进行汇总
    """
    os.makedirs(outdir, exist_ok=True)
    
    # 创建一个字典来存储每个barcode的总reads数
    barcode_total_reads = {}
    
    # 检查输入是文件还是目录
    if os.path.isdir(raw_aligned_reads_counts):
        print(f"处理目录: {raw_aligned_reads_counts}")
        # 获取目录下所有以_barcode_counts.csv结尾的文件
        csv_files = glob.glob(os.path.join(raw_aligned_reads_counts, "*_cb_aligned_reads_counts.csv"))
        
        if not csv_files:
            print(f"在{raw_aligned_reads_counts}中没有找到*_cb_aligned_reads_counts.csv文件")
            return
            
        print(f"找到{len(csv_files)}个文件需要处理")
        
        # 读取所有文件并汇总barcode的reads数
        for csv_file in csv_files:
            print(f"处理文件: {os.path.basename(csv_file)}")
            try:
                df = pd.read_csv(csv_file)
                
                # 确保文件包含必要的列
                if 'barcode' not in df.columns or 'aligned_reads' not in df.columns:
                    print(f"警告: 文件{csv_file}没有必要的列(barcode, aligned_reads)，跳过。")
                    continue
                
                # 遍历每一行，将reads数加到对应barcode的总数中
                for _, row in df.iterrows():
                    barcode = row['barcode']
                    reads = row['aligned_reads']
                    if barcode in barcode_total_reads:
                        barcode_total_reads[barcode] += reads
                    else:
                        barcode_total_reads[barcode] = reads
                        
            except Exception as e:
                print(f"处理文件{csv_file}时出错: {str(e)}")
                continue
    else:
        # 直接读取单个文件
        print(f"处理单个文件: {raw_aligned_reads_counts}")
        try:
            df = pd.read_csv(raw_aligned_reads_counts)
            # 确保文件包含必要的列
            if 'barcode' not in df.columns or 'aligned_reads' not in df.columns:
                print(f"警告: 文件{raw_aligned_reads_counts}没有必要的列(barcode, aligned_reads)。")
                return
                
            # 遍历每一行，将reads数加到对应barcode的总数中
            for _, row in df.iterrows():
                barcode = row['barcode']
                reads = row['aligned_reads']
                barcode_total_reads[barcode] = reads
                
        except Exception as e:
            print(f"读取文件{raw_aligned_reads_counts}时出错: {str(e)}")
            return
    
    if not barcode_total_reads:
        print("没有找到有效数据。退出。")
        return
        
    # 将汇总结果转换为DataFrame
    merged_df = pd.DataFrame({
        'barcode': list(barcode_total_reads.keys()),
        'aligned_reads': list(barcode_total_reads.values())
    })
    
    # 保存合并后的数据
    merged_file = os.path.join(outdir, 'merged_barcode_counts.csv')
    merged_df.to_csv(merged_file, index=False)
    print(f"合并后的数据已保存到{merged_file}")
    
    # 按reads数量降序排序
    merged_df = merged_df.sort_values(by='aligned_reads', ascending=False)
    
    # 计算阈值
    percentile = 99
    threshold_index = int(expected_cell_num * (1 - percentile / 100))
    if threshold_index >= len(merged_df):
        threshold_index = len(merged_df) - 1
        
    readscut = merged_df.iloc[threshold_index]['aligned_reads'] * 0.1
    
    # 筛选出reads大于readscut的行
    filtered_df = merged_df[merged_df['aligned_reads'] > readscut]
    
    # 保存过滤后的结果
    filtered_file = os.path.join(outdir, 'filtered_barcode_read_counts.csv')
    filtered_df.to_csv(filtered_file, sep=',', index=False, header=False)
    filtered_df["barcode"].to_csv('filtered_barcode', sep=',', index=False, header=False)
    
    print(f"过滤后的barcode已保存到{filtered_file}")
    print(f"原始barcode数量: {len(merged_df)}, 过滤后barcode数量: {len(filtered_df)}")

if __name__ == '__main__':
    estimated_cells()