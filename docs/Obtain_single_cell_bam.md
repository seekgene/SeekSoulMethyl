# How to Obtain Single-Cell BAM Files

To accelerate methylation data analysis, we split FASTQ files at the FASTQ processing stage according to the first several bases of the error-corrected cell barcode. By default, we use the first 4 bases, but this length can be configured (e.g., 1 or 2 bases); the actual length used can be inferred from the file names.
 
In addition, since Bismark alignment requires different parameters for the original strand (OT/OB) and the complementary strand (CTOT/CTOB), we further split the reads into Forward FASTQ and Reverse FASTQ files based on the barcode splitting.
 
The final BAM files provided are the Bismark alignment results generated from these split FASTQ files. If you need to obtain single-cell BAM files, please follow the steps below:

>[!Note]
>If you have analyzed the data using the SeekSoulMethyl pipeline, the single-cell BAM files are automatically generated and stored in the directory `${sample}/${sample}_methy/step3/split_bams/`. You do not need to perform the following operations manually.

**Step 1: Sort each BAM file by read name (name-sort)**

```shell
# ....
samtools sort -n -o WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn.bam WTJW969_forward_TTTG_1_bismark_bt2_pe.bam
samtools sort -n -o WTJW969_reverse_TTTG_1_bismark_bt2_pe_sortn.bam WTJW969_reverse_TTTG_1_bismark_bt2_pe.bam
#....
```

**Step 2: Split single-cell BAM files from each name-sorted BAM**

[step3_split_bams.py](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/nf/bin/step3_split_bams.py)

If the data type is **DD-MET5**, run the following commands for each name-sorted BAM:

```shell
# gex_barcodes: RNA filtered_feature_bc_matrix/barcodes.tsv.gz
gex_barcodes="/path/to/filtered_feature_bc_matrix/barcodes.tsv.gz"

# If the data type is DD-MET5
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

If the data type is **DD-MET3**, you need to use the `--cbcsv` parameter to provide the barcode whitelist mapping between transcriptome and methylation data [DD-M_bUCB3_whitelist.csv](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/nf/bin/barcodes/DD-M_bUCB3_whitelist.csv), then run the following for each name-sorted BAM:

```shell
# gex_barcodes: RNA filtered_feature_bc_matrix/barcodes.tsv.gz
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

**Step 3: Merge Forward and Reverse BAM Files for Each Single Cell**
 
Typically, each cell contains both Forward and Reverse reads. Therefore, you must merge the two to obtain the complete BAM file for that cell.

After Step 2, using the sample `WTJW969_forward_TTTG_1_bismark_bt2_pe_sortn.bam` as an example, `step3_split_bams.py` will automatically create a subdirectory under the specified `outdir`:

- forward: `./WTJW969_forward_TTTG_1/`
- reverse: `./WTJW969_reverse_TTTG_1/`

The split single-cell BAM files and their corresponding `*_filtered_barcode` and `*_filtered_barcode_reads_counts.csv` files are all located in these directories.

```shell
set -e
merged_fr_bam="./merged_fr_bam/"
mkdir -p $merged_fr_bam

# Modify the following two directories according to your sample names
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
# Merge filtered_barcode files from forward and reverse directories
cat \
  ${split_forward_bams_dir}/*_filtered_barcode \
  ${split_reverse_bams_dir}/*_filtered_barcode \
  | sort | uniq > ${merged_fr_bam}/WTJW969_forward_TTTG_1_merge_filtered_barcode

# Merge all filtered_barcode_reads_counts files and aggregate reads_counts by barcode
echo "Merging and aggregating reads_counts files..."

# Skip header line, aggregate reads_counts by barcode
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

## Batch Processing

If you need to process multiple files in batches, we highly recommend using the script [batch_single_cell_bam.py](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/nf/bin/utils/batch_single_cell_bam.py).Before using it, please ensure that you have downloaded [step3_split_bams.py](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/nf/bin/step3_split_bams.py) and placed it in the same directory as `batch_single_cell_bam.py` or in its parent directory.

```shell
python batch_single_cell_bam.py \
--assay_type DD-MET5 \
--bam_dir /path/to/bismark/ \
--outdir /path/to/output/ \
--gex_barcodes /path/to/filtered_feature_bc_matrix/barcodes.tsv.gz \
--threads 4 \
--parallel_jobs 8

# assay_type: Data type, DD-MET5 or DD-MET3
# bam_dir: Directory containing all BAM files
# outdir: Output directory for results
# gex_barcodes: RNA filtered_feature_bc_matrix/barcodes.tsv.gz file
# threads: Number of threads per task, default is 4
# parallel_jobs: Number of parallel tasks (affects concurrency of sorting, splitting, and merging steps), default is 8
```

`--bam_dir` is the input directory with the following structure:
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

After execution, the single-cell BAM files will be stored in the `/path/to/output/split_bam` directory.

>[!Note]
>Using a server with 32 CPUs and 32GB RAM, splitting BAM files for 2196 cells takes approximately 2h 23m.
