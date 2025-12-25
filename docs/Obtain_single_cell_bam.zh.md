# 如何获取单个细胞的 BAM 文件

为了加快甲基化数据分析速度，我们在 FASTQ 阶段就按照细胞 barcode（barcode已完成纠错） 的前若干个碱基拆分 FASTQ 文件。默认使用前 4 个碱基，但具体长度可以根据流程参数设置为 1 个或 2 个碱基，实际使用的长度可以从文件名中判断。

在按 barcode 拆分 FASTQ 的同时，我们还会区分原始链和互补链的 reads。因为原始链和互补链在 bismark 比对过程中使用的参数不同，所以在上述拆分规则的基础上，会进一步将 FASTQ 拆分为 forward FASTQ 和 reverse FASTQ 文件。


最终我们释放的 BAM 文件，就是在上述拆分基础上，经 bismark 比对后得到的结果。如果用户需要获取单个细胞的 BAM 文件，可以按照以下步骤操作：

>[!Note]
>如果您自己使用我们的SeekSoulMethyl分析了数据，完整的单个细胞的bam文件存储在 `${sample}/${sample}_methy/step3/split_bams/` 目录下，无需进行以下操作。

**step1 将每个 BAM 按照 read 名称进行排序（name-sort）**

```shell
# ....
samtools sort -n -o WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn.bam WTJW969_forward_TTTG_1_bismark_bt2_pe.bam
samtools sort -n -o WTJW969_reverse_TTTG_1_bismark_bt2_pe_sortn.bam WTJW969_reverse_TTTG_1_bismark_bt2_pe.bam
#....
```

**step2 从每个排序后的 BAM 文件中拆分出单细胞 BAM 文件**

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

**step3 合并单个细胞的 forward 和 reverse BAM 文件**

一般每个细胞既有 forward 也有 reverse 的 reads，因此需要将 forward 和 reverse 的 BAM 文件合并成一个 BAM，才是完整的单细胞 BAM 文件。

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
