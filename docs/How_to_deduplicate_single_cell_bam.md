# How to Deduplicate Single-Cell BAM Files

## Introduction

When performing Copy Number Variation (CNV) analysis or other analyses sensitive to read counts, ensuring input data accuracy is crucial. Random fragmentation and PCR amplification during single-cell methylation library construction introduce PCR bias. To eliminate this bias, we must deduplicate the BAM files.

This deduplication process consists of two key steps:
1. **Flag Correction**: Correcting Flag anomalies produced by Bismark alignment.
2. **UMI Deduplication**: Deduplicating based on UMI and alignment position.

## Preparation

### Input Data
Split single-cell BAM files.

### Dependencies
Please ensure the following tools are installed in your environment:
- **Python 3**
- **[umi_tools](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html)**
- **samtools** (depends on pysam)

Install `umi_tools`:
```shell
pip install umi_tools
```

### Required Scripts
You need to download the following two scripts:
- [`step1_swap_flags_on_reverse.py`](../nf/bin/utils/step1_swap_flags_on_reverse.py)
- [`step2_umi_tools_dedup.py`](../nf/bin/utils/step2_umi_tools_dedup.py)

## Detailed Steps

### Step 1: Flag Correction

**Background**:
When processing OT/OB (original top/bottom) and CTOT/CTOB (complementary to original top/bottom) reads, Bismark [automatically swaps the Flags of CTOT/CTOB reads](https://github.com/FelixKrueger/Bismark/issues/41). Specifically, R1 is flagged as "second in pair" and R2 as "first in pair". To ensure downstream deduplication tools (like umi_tools) correctly identify R1 and R2, we need to revert these Flags to normal.

**Command**:

```shell
# Configuration
outdir="./"
input_dir="/path/to/step3/split_bams/merged/"  # Directory containing single-cell BAMs
workers=16                                      # Number of parallel workers

# 1. Generate file list
find -L "$input_dir" -name "*.bam" | sort > "$outdir/bam_list.txt"

# 2. Create output directory
mkdir -p $outdir/flag_reverse/

# 3. Execute Flag Correction
python step1_swap_flags_on_reverse.py \
  --bam_list $outdir/bam_list.txt \
  --out_dir $outdir/flag_reverse \
  --max_workers $workers
```

**Output**:
The Flag-corrected BAM files will be generated in the `$outdir/flag_reverse/` directory, with the suffix `.swapflags.bam`.

### Step 2: UMI Deduplication

**Background**:
Since the ends of library fragments are randomly fragmented, we use **R1 alignment position + UMI sequence** as the Unique Identifier for deduplication.

**Command**:
Use `umi_tools` for deduplication with key parameters `--extract-umi-method tag --umi-tag UR --paired --ignore-tlen`.

```shell
# 1. Generate list of corrected BAMs
find "$outdir/flag_reverse" -name "*.bam" | sort > "$outdir/swapflags_bam_list.txt"

# 2. Create output directory
mkdir -p $outdir/umi_tools_swapflags_dedup_bam/

# 3. Execute Deduplication
python step2_umi_tools_dedup.py \
  --bam_list $outdir/swapflags_bam_list.txt \
  --out_dir $outdir/umi_tools_swapflags_dedup_bam/ \
  --max_workers $workers \
  --skip_existing
```

**Output**:
The deduplicated BAM files will be stored in the `$outdir/umi_tools_swapflags_dedup_bam/` directory, with the suffix `.swapflags.dedup.bam`. These files are the final deduplicated single-cell BAM files.
