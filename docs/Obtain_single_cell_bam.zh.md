# 如何获取单个细胞的 BAM 文件

为了加速甲基化数据分析，我们在 FASTQ 处理阶段会根据已纠错的细胞 barcode 的前若干个碱基（默认为 4 个，也可配置为 1 或 2 个）对文件进行拆分。实际拆分长度可直接从文件名中识别。

同时，由于 Bismark 比对时原始链（OT/OB）和互补链（CTOT/CTOB）需要不同的参数，我们会在拆分 barcode 的基础上，进一步将 reads 分流为 forward FASTQ 和 reverse FASTQ 文件。

最终交付的 BAM 文件即为上述拆分后的 FASTQ 经 Bismark 比对生成的。若您需要获取单细胞的 BAM 文件，请按照以下步骤操作：

>[!Note]
>如果您在自己的服务器上使用 SeekSoulMethyl 流程分析了数据，单细胞 BAM 文件已自动生成并存储于 `${sample}/${sample}_methy/step3/split_bams/` 目录下，无需手动执行以下操作。

**Step 1: 将每个 BAM 按照 read 名称进行排序（name-sort）**

```shell
# ....
samtools sort -n -o WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn.bam WTJW969_forward_TTTG_1_bismark_bt2_pe.bam
samtools sort -n -o WTJW969_reverse_TTTG_1_bismark_bt2_pe_sortn.bam WTJW969_reverse_TTTG_1_bismark_bt2_pe.bam
#....
```

**Step 2: 从每个排序后的 BAM 文件中拆分出单细胞 BAM 文件**

[step3_split_bams.py](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/nf/bin/step3_split_bams.py)

如果数据是 **DD-MET5**，对每个排序后的 BAM 文件进行如下处理：

```shell
# gex_barcodes为RNA数据的filtered_feature_bc_matrix/barcodes.tsv.gz文件
gex_barcodes="/path/to/filtered_feature_bc_matrix/barcodes.tsv.gz"

# 如果数据为DD-MET5
# ....
python step3_split_bams.py \
--bam WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn.bam \
--outdir ./ \
--samplename WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn \
--core 1 \
--gexcb ${gex_barcodes}

python step3_split_bams.py \
--bam WTJW969_reverse_TTTG_1_bismark_bt2_pe_sortn.bam \
--outdir ./ \
--samplename WTJW969_reverse_TTTG_1_bismark_bt2_pe_sortn \
--core 1 \
--gexcb ${gex_barcodes}
# ...
```

如果数据为 **DD-MET3**，需要使用 `--cbcsv` 参数，指定转录组与甲基化数据的 barcode 白名单对应表 [DD-M_bUCB3_whitelist.csv](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/nf/bin/barcodes/DD-M_bUCB3_whitelist.csv)，然后对每个排序后的 BAM 进行如下处理：

```shell
# gex_barcodes为RNA数据的filtered_feature_bc_matrix/barcodes.tsv.gz文件
gex_barcodes="/path/to/filtered_feature_bc_matrix/barcodes.tsv.gz"
# ...
python step3_split_bams.py \
--bam WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn.bam \
--outdir ./ \
--samplename WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn \
--core 1 \
--gexcb ${gex_barcodes} \
--cbcsv DD-M_bUCB3_whitelist.csv

python step3_split_bams.py \
--bam WTJW969_reverse_TTTG_1_bismark_bt2_pe_sortn.bam \
--outdir ./ \
--samplename WTJW969_reverse_TTTG_1_bismark_bt2_pe_sortn \
--core 1 \
--gexcb ${gex_barcodes} \
--cbcsv DD-M_bUCB3_whitelist.csv
#...
```

**Step 3: 合并单个细胞的 forward 和 reverse BAM 文件**

通常每个细胞同时包含 forward 和 reverse 的 reads，因此必须将两者合并，才能得到该细胞完整的 BAM 文件。

以上步骤完成后，以样本名 `WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn.bam` 为例，`step3_split_bams.py` 会在指定 `outdir` 下自动创建一个子目录：

- forward：`./WTJW969_forward_TTTG_1/`
- reverse：`./WTJW969_reverse_TTTG_1/`

单细胞拆分后的 BAM 和对应的 `*_filtered_barcode`、`*_filtered_barcode_reads_counts.csv` 文件均位于这些子目录中。

```shell
set -e
merged_fr_bam="./merged_fr_bam/"
mkdir -p $merged_fr_bam

# 根据实际样本名修改下面两个目录
split_forward_bams_dir="./WTJW969_forward_TTTG_1/"
split_reverse_bams_dir="./WTJW969_reverse_TTTG_1/"

forward_bams=($(find "$split_forward_bams_dir" -name "*.bam" | sort))
reverse_bams=($(find "$split_reverse_bams_dir" -name "*.bam" | sort))

declare -A fmap
declare -A rmap
declare -A seen

for fb in "${forward_bams[@]}"; do bn=$(basename "$fb" .bam); fmap["$bn"]="$fb"; done
for rb in "${reverse_bams[@]}"; do bn=$(basename "$rb" .bam); rmap["$bn"]="$rb"; done

for bc in "${!fmap[@]}"; do seen["$bc"]=1; done
for bc in "${!rmap[@]}"; do if [[ -z "${seen[$bc]}" ]]; then seen["$bc"]=1; fi; done
for bc in "${!seen[@]}"; do
    f="${fmap[$bc]:-}"
    r="${rmap[$bc]:-}"
    output_bam="${merged_fr_bam}/${bc}.bam"
    if [[ -n "$f" && -n "$r" ]]; then
        samtools merge -n -@ 2 -o "$output_bam" "$f" "$r"
    elif [[ -n "$f" ]]; then
        cp "$f" "$output_bam"
    elif [[ -n "$r" ]]; then
        cp "$r" "$output_bam"
    fi
done

echo "Merging and deduplicating barcode files..."
# 合并 forward 和 reverse 目录下的 filtered_barcode
cat \
  ${split_forward_bams_dir}/*_filtered_barcode \
  ${split_reverse_bams_dir}/*_filtered_barcode \
  | sort | uniq > ${merged_fr_bam}/WTJW969_forward_TTTG_1_merge_filtered_barcode

# Merge all filtered_barcode_reads_counts files and aggregate reads_counts by barcode
echo "Merging and aggregating reads_counts files..."

# 跳过各文件首行表头，按 barcode 汇总 reads_counts
awk -F ',' '{
    if (NR == 1) {
        print
        next
    }
    if (NF >= 2) {
        barcode = $2
        reads = $1
        total[barcode] += reads
    }
} END {
    for (barcode in total) {
        print total[barcode] "," barcode
    }
}' \
  ${split_forward_bams_dir}/*_filtered_barcode_reads_counts.csv \
  ${split_reverse_bams_dir}/*_filtered_barcode_reads_counts.csv \
  > ${merged_fr_bam}/WTJW969_forward_TTTG_1_merge_filtered_barcode_reads_counts.csv

echo "BAM file merging, reads_counts aggregation completed"
```
## 批量处理
如果您需要批量处理多个文件，推荐使用脚本 [batch_single_cell_bam.py](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/nf/bin/utils/batch_single_cell_bam.py)。该脚本已优化，支持并行合并（充分利用多核 CPU）及断点续跑功能。

```shell
python batch_single_cell_bam.py \
--assay_type DD-MET5 \
--bam_dir /path/to/bismark/ \
--outdir /path/to/output/ \
--gex_barcodes /path/to/filtered_feature_bc_matrix/barcodes.tsv.gz \
--threads 4 \
--parallel_jobs 8

# assay_type：数据类型， DD-MET5或者DD-MET3
# bam_dir：包含所有bam的目录
# outdir：结果输出目录
# gex_barcodes：转录组数据的filtered_feature_bc_matrix/barcodes.tsv.gz文件
# threads：每个任务的线程数，默认4
# parallel_jobs：并行任务数（影响排序、拆分及合并步骤的并发度），默认 8
```

`--bam_dir`为输入目录，结构如下：
```text
/path/to/bismark/
├── WTJW969_forward_AAAG_1_bismark_bt2_pe.bam
├── WTJW969_forward_AAAT_1_bismark_bt2_pe.bam
├── WTJW969_forward_AAGA_1_bismark_bt2_pe.bam
├── ...
├── WTJW969_reverse_TTGT_1_bismark_bt2_pe.bam
├── WTJW969_reverse_TTTA_1_bismark_bt2_pe.bam
└── WTJW969_reverse_TTTG_1_bismark_bt2_pe.bam
```
以上执行完成后，单个细胞的bam文件存储在`/path/to/output/split_bam`目录下

>[!Note]
>使用32CPU，32GB的服务器拆分2196个细胞的BAM文件，用时大概2h23m。
