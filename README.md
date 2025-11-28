# SeekSoulMethyl
SeekSoulMethyl is a single-cell transcriptome + methylation dual-omics analysis pipeline designed for analyzing data from Beijing SeekGene Biotechnology Co., Ltd. single-cell transcriptome + methylation dual-omics kit.

## Data Structure
RNA-MET data comes in two chemistry types. DD-MET3 (dual-label) means the RNA and DNA methylation data barcodes are different for the same cell, and the RNA library is 3′-end trancriptome. DD-MET5 (single-label) means the RNA and DNA methylation data barcodes are the same for the same cell, and the RNA library is 5′-end trancriptome. Below we describe the DNA methylation library structures for both chemistries.

**DD-MET3 Methylation Library Structure**
<figure style="text-align: center;">
<img src="./docs/DD-MET3_library_structure.png" alt="DD-MET3 Methylation Library" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 1. DD-MET3 methylation library structure</figcaption>
</figure>

Structure notes:
- SP1/SP2: Adapter sequences
- barcode: 17 bp cell barcode
- 7F: 7 bp linker sequence
- 17L and ME: 17 bp fixed sequence <span style="color:#43a047;">Y</span>gt<span style="color:#43a047;">Y</span><span style="color:#43a047;">Y</span>gt<span style="color:#43a047;">Y</span>gttg<span style="color:#43a047;">Y</span>t<span style="color:#43a047;">Y</span>gt
- ME: 19 bp fixed sequence AGATGTGTATAAGAGA<span style="color:#43a047;">Y</span>AG
- 9 bp: extension sequence from the Tn5 insertion fragment

**DD-MET5 Methylation Library Structure**
<figure style="text-align: center;">
<img src="./docs/DD-MET5_library_structure.png" alt="DD-MET5 Methylation Library" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 2. DD-MET5 methylation library structure</figcaption>
</figure>

Structure notes:
- SP1/SP2: Adapter sequences
- barcode: 17 bp cell barcode
- UMI: 12 bp UMI sequence
- TSO: 13 bp TSO sequence TTT<span style="color:#43a047;">Y</span>TTATATGGG
- 17L and ME: 17 bp fixed sequence <span style="color:#43a047;">Y</span>gt<span style="color:#43a047;">Y</span><span style="color:#43a047;">Y</span>gt<span style="color:#43a047;">Y</span>gttg<span style="color:#43a047;">Y</span>t<span style="color:#43a047;">Y</span>gt
- ME: 19 bp fixed sequence AGATGTGTATAAGAGA<span style="color:#43a047;">Y</span>AG
- 9 bp: extension sequence from the Tn5 insertion fragment


## Installation

1. Clone the repository:
```bash
git clone https://github.com/seekgene/SeekSoulMethyl.git
cd SeekSoulMethyl
```

2. Create and activate conda environment:

For users in China:
```bash
conda env create -n seeksoulmethyl -f conda_dependencies.zh.yml
conda activate seeksoulmethyl
```

For international users:
```bash
conda env create -n seeksoulmethyl -f conda_dependencies.yml
conda activate seeksoulmethyl
```

3. Install the package:
```bash
cd dependence
pip install . \
  simpleqc/target/wheels/simpleqc-0.1.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl \
  search-pattern/target/wheels/search_pattern-0.1.0-py3-none-manylinux_2_5_x86_64.manylinux1_x86_64.whl
cd ..

```
We will download our modified versions of Bismark and ALLCools for analysis.

- [Bismark](https://github.com/seekgene/Bismark.git) adds the CB (error-corrected barcode) tag and the UR (raw UMI) tag to BAM files.
- [ALLCools](https://github.com/seekgene/ALLCools.git) performs UMI deduplication and methylation level calculation based on the UR tag.

```shell
# Clone and install custom ALLCools
conda activate seeksoulmethyl
git clone https://github.com/seekgene/ALLCools.git && \
pip install ./ALLCools && \
rm -rf ./ALLCools

# Clone and install custom Bismark
git clone https://github.com/seekgene/Bismark.git && \
bin_path=$(dirname `which python`)
cp -r ./Bismark/* $bin_path/ && \
    chmod +x $bin_path/bismark* && \
    chmod +x $bin_path/deduplicate_bismark && \
    rm -rf ./Bismark

```

## Download Reference Database
```bash
# Download human reference genome (GRCh38)
wget -dc -O human-reference-GRCh38.tar.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/v1.1/human-reference-GRCh38.tar.gz"
wget -dc -O human-reference-GRCh38.tar.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/v1.1/human-reference-GRCh38.tar.gz.md5"

# Download mouse reference genome (GRCm39)
wget -dc -O mouse-reference-GRCm39.tar.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/v1.1/mouse-reference-GRCm39.tar.gz"
wget -dc -O mouse-reference-GRCm39.tar.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/v1.1/mouse-reference-GRCm39.tar.gz.md5"

# Extract reference genomes
tar -xzf human-reference-GRCh38.tar.gz
tar -xzf mouse-reference-GRCm39.tar.gz
```

## Download Test Data (Optional)
For testing the pipeline with a small dataset, you can download the tiny test data:

```bash
# Download transcriptome test data
wget -dc -O XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/fastq/XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz"
wget -dc -O XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/fastq/XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz.md5"
wget -dc -O XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/fastq/XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz"
wget -dc -O XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/fastq/XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz.md5"

# Download methylation test data
wget -dc -O XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz"
wget -dc -O XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz.md5"
wget -dc -O XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz"
wget -dc -O XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz.md5"


```

**Note**: This is a small test dataset for pipeline validation. For production analysis, use your own sequencing data.

## Repository Layout

After cloning, the key Nextflow entry points and modules are:

- `nf/main.nf`: Main workflow for RNA + methylation end-to-end processing (script/SeekSoulMethyl/nf/main.nf:1)
- `nf/methy_only.nf`: Workflow for methylation-only data (script/SeekSoulMethyl/nf/methy_only.nf:1)
- `nf/modules/`: Step-wise process modules:
  - `step1.nf` preprocessing, QC, barcode extraction, RNA analysis
  - `step2.nf` Bismark alignment and BAM sorting
  - `step3.nf` per-cell BAM splitting, ALLC generation/merge, multi-scale datasets
  - `step4.nf` summaries, dimensionality reduction, joint report
  - `utils.nf` helpers for methylation-only workflow (reads counting and cell estimation)
- `nf/bin/`: Helper scripts and resources (e.g., barcode whitelists)
- `nf/nextflow.config`: Executors and resource configuration (script/SeekSoulMethyl/nf/nextflow.config:1)

## Usage

### Activate Environment
```bash
conda activate seeksoulmethyl
```

### Run Dual-omics Analysis (Shell script)
```bash
# sc_methy_workflow.sh can be found in the SeekSoulMethyl directory you cloned
bash sc_methy_workflow.sh \
/path/to/expression_R1.fastq.gz \
/path/to/expression_R2.fastq.gz \
/path/to/methy_R1.fastq.gz \
/path/to/methy_R2.fastq.gz \
--sample WTJW880 \
--outdir /path/to/results \
--database_dir /path/to/human-reference-GRCh38 \
--chemistry DD-MET3 \
--core 64 \
--filter_ch 2
```

### Input Parameters:

- **$1**: Single-cell transcriptome Read1 fastq file
- **$2**: Single-cell transcriptome Read2 fastq file  
- **$3**: Single-cell methylation Read1 fastq file
- **$4**: Single-cell methylation Read2 fastq file
- **sample**: Sample name
- **outdir**: Output directory path
- **database_dir**: Reference genome database path
- **chemistry**: Methylation chemistry (DD-MET3 or DD-MET5). DD-MET3 is a dual-label dataset, meaning the RNA and DNA methylation data barcodes are different for the same cell, while DD-MET5 is single-label, meaning the RNA and DNA methylation data barcodes are the same for the same cell.
- **core**: Number of CPU cores
- **filter_ch**: Filter reads that contain > n CH methylation sites.

## Process Details
### RNA Processing Workflow
RNA data is analyzed using [SeekSoulTools](http://seeksoul.seekgene.com/en/v1.3.0/2.tutorial/1.rna/4.description.html). See the official Algorithms Overview for detailed steps. Cells used in the downstream methylation library are determined based on the RNA library cell barcodes.

### Methylation Processing Workflow
#### Step 1: Preprocessing and Barcode Parsing
**Barcode extraction and correction**
Based on the designed structure, we locate the barcode in the read and extract the corresponding sequence. If the extracted barcode is in the whitelist, it is counted as a valid barcode; otherwise, it is considered invalid.
Sequencing errors can occur. When a whitelist is provided, SeekSoulTools can attempt barcode correction. If correction is enabled and an invalid barcode has a one-base mismatch (Hamming distance = 1) from an entry in the whitelist:
* If exactly one whitelist candidate matches: correct the invalid barcode to that whitelist barcode.
* If multiple whitelist candidates match: correct to the candidate supported by the highest read count.

**UMI extraction**
UMI positions are read from the designed structure and extracted without correction.

**Forward and Reverse reads determination**
From the positions corresponding to 17L and ME, there are 7 bases that can be C or converted T (highlighted below). Use the final two C/T positions for forward and reverse reads determination: both C indicates a reverse read; otherwise it is forward read.
Forward: <span style="color:#43a047;">T</span>gt<span style="color:#43a047;">TT</span>gt<span style="color:#43a047;">T</span>gttg<span style="color:#43a047;">T</span>t<span style="color:#43a047;">T</span>gtAGATGTGTATAAGAGA<span style="color:#43a047;">T</span>
Reverse: <span style="color:#43a047;">C</span>gt<span style="color:#43a047;">CC</span>gt<span style="color:#43a047;">C</span>gttg<span style="color:#43a047;">C</span>t<span style="color:#43a047;">C</span>gtAGATGTGTATAAGAGA<span style="color:#43a047;">C</span>
Reverse reads correspond to CTOT/CTOB in methylation terminology; forward reads correspond to OT/OB.

**C–T conversion rate**
Aside from the two positions used for forward and reverse reads determination, the remaining five C/T bases are used to compute the C–T conversion rate:
<div style="font-size: 0.8em;">

$$
    \text{CT conversion} = \frac{\text{Total number of "T" across five positions in forward reads}}{\text{Number of forward reads} \times 5}
$$

</div>

**Filter reads (optional)**
Filter based on the number of non-CpG methylated C bases in a read pair. By default, pairs with > 2 non-CpG methylated Cs detected in read1/read2 are removed.

#### Step 2: Bismark alignment and sorting by name
**Alignment and tagging**
We use Bismark for methylation alignment. Our modified Bismark adds `--add_barcode` and `--add_umi` to tag BAM files by read name with CB (error-corrected barcode) and UR (raw UMI). For forward-strand reads, we use `-X 1000` to allow insert sizes up to 1000 bp; for reverse-strand reads, we use `--pbat` and `-X 1000`.
After alignment, sort BAMs by read name using `samtools sort -n`; the name-sorted BAMs serve as inputs for downstream analysis.

#### Step 3: ALLCools analysis
**Split by cell barcode**
Split name-sorted BAMs by RNA-derived cell barcodes into per-cell BAM files, each containing uniquely mapped reads for one cell.
**Generate ALLC files**
Sort each per-cell BAM by position and convert to ALLC using ALLCools `bam-to-allc`. Our modified ALLCools performs UR-tag-based UMI correction and deduplication per C site; see the [SeekGene ALLCools repository](https://github.com/seekgene/ALLCools) for details.
**Generate MCDS**
Run `allcools generate-dataset` to bin the genome (chrom10k/20k/50k/100k/500k/1M) and compute per-cell methylation matrices.

#### Step 4: Reduction and clustering
By default, perform dimensionality reduction with LSI on chrom20k bins using ALLCools, followed by UMAP visualization.

### System Requirements
If you use `sc_methy_workflow.sh`, the system requirements are as follows:
- **CPU**: 64 cores (recommended)
- **Memory**: 128GB RAM (recommended)
- **Storage**: At least 500GB available space
- **OS**: Linux (recommended Ubuntu 18.04+ or CentOS 7+)

## Quick Start with Nextflow (Recommended)

If you want to process samples in batch and get pipeline-level reports, use Nextflow:

```bash
conda activate seeksoulmethyl
conda install -n seeksoulmethyl -c bioconda nextflow

# Prepare the sample sheet (see example above)
nextflow run SeekSoulMethyl/nf/main.nf \
  --outdir /path/to/results \
  --samplesheet samplelist.csv \
  -w /path/to/work \
  -c SeekSoulMethyl/nf/nextflow.config \
  -profile aliyun_k8s \
  --database_dir /path/to/reference \
  --split_fastq 4 \
  --filter_ch 2 \
  --chemistry DD-MET3 > methy.log
```
<figure style="text-align: center;">
<img src="./docs/nf_SeekSoulMethyl_workflow.png" alt="SeekSoulMethyl Pipeline" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 3. SeekSoulMethyl nextflow pipeline workflow</figcaption>
</figure>

## Running with Nextflow

1. Install nextflow:
```bash
conda install -n seeksoulmethyl -c bioconda nextflow
```
2. Prepare your input samplesheet
```
cat > samplelist.csv << EOF
sample_id,expression_r1,expression_r2,methylation_r1,methylation_r2
XYRD-WTJW880,/path/to/XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz,/path/to/XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz,/path/to/XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz,/path/to/XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz
EOF

# expression_r1: Single-cell transcriptome Read1 fastq file
# expression_r2: Single-cell transcriptome Read2 fastq file
# methylation_r1: Single-cell methylation Read1 fastq file
# methylation_r2: Single-cell methylation Read2 fastq file
```

3. Run the pipeline:
```bash
nextflow run -bg SeekSoulMethyl/nf/main.nf \
--outdir /path/to/tiny_demo/results/ \
--samplesheet samplelist.csv \
-w /path/to/tiny_demo/results/work \
-c SeekSoulMethyl/nf/nextflow.config \
-profile aliyun_k8s \
--database_dir /path/to/human-reference-GRCh38/ \
--split_fastq 1 \
--filter_ch 2 \
--chemistry DD-MET3 > methy.log

# --outdir: final results directory
# --samplesheet: input samplesheet file
# -w: working directory for nextflow
# -c: nextflow configuration file. Must be modified according to your server configuration!!!!!
# -profile: nextflow profile for aliyun k8s
# --database_dir: reference genome database directory
# --split_fastq: To speed up the analysis process, split fastq according to first n bases of barcode. Default is 4.
# --filter_ch: filter reads that contain > 2 CH methylation sites.
# --chemistry: methylation chemistry (DD-MET3 or DD-MET5)
```

### Notes on nextflow.config

`nextflow.config` declares the execution environment, resource quotas, and run policies. You must customize it to your own server or cluster.

- Location: `SeekSoulMethyl/nf/nextflow.config` (specify with `-c SeekSoulMethyl/nf/nextflow.config` when running).
- Key items to tailor to your infrastructure:
  - Executor: `process.executor` (e.g., `local`, `slurm`, `pbs`, `k8s`, `awsbatch`).
  - Resources: `process.cpus`, `process.memory`, `process.time`, or fine-grained overrides via `withLabel`/`withName`.
  - Work directory: `workDir` (can also be set via `-w`); ensure it is writable and has sufficient space.
  - Environment: enable one of `conda.enabled`, `docker.enabled`, or `singularity.enabled` according to what your server supports.

Example configurations (replace paths and parameters with values valid on your servers):

```groovy
profiles {
  // Local machine
  local {
    process.executor = 'local'
    workDir          = '/path/to/work'
    process.cpus     = 8
    process.memory   = '32 GB'
    conda.enabled    = true
    // Or containers: singularity.enabled = true / docker.enabled = true
  }

  // Slurm cluster (HPC)
  slurm {
    process.executor   = 'slurm'
    workDir            = '/lustre/work/your_user'
    process.cpus       = 8
    process.memory     = '32 GB'
    process.queue      = 'normal'
    process.clusterOptions = '-A your_account --qos=normal'

    withLabel: 'high_mem' {
      cpus   = 16
      memory = '64 GB'
    }
  }

  // Kubernetes (e.g., Alibaba Cloud ACK)
  aliyun_k8s {
    process.executor   = 'k8s'
    workDir            = '/mnt/nf-work'    // persistent volume path
    k8s {
      namespace        = 'your-namespace'
      storageClaimName = 'your-pvc'
      cpu              = 4
      memory           = '16 GB'
    }
    // If using container images globally: docker.enabled = true
  }
}
```

Tips:
- Pick the `-profile` that matches your environment (e.g., `slurm`, `local`, `aliyun_k8s`), then adapt the parameters.
- Keep `workDir` on fast storage with ample capacity; the final results directory is controlled by `--outdir`.
- If you use the README’s conda environment, prefer `conda.enabled`; if your cluster favors containers, use Docker/Singularity.

References:
- Nextflow configuration & profiles: https://www.nextflow.io/docs/latest/config.html
- Executors (local, Slurm, K8s, etc.): https://www.nextflow.io/docs/latest/executor.html
- Kubernetes guide: https://www.nextflow.io/docs/latest/kubernetes.html

## Nextflow Step-by-Step Details

This section describes each Nextflow process with inputs, core logic, key parameters, and outputs for troubleshooting and interpretation. The workflow is defined in `nf/main.nf` (script/SeekSoulMethyl/nf/main.nf:1) and processes are implemented in `nf/modules/*.nf`.

### Step 1: Preprocessing and Barcode Parsing (step1.nf)
- Compute genome-wide CpG count: `COMPUTE_CPG_SITES` (script/SeekSoulMethyl/nf/modules/step1.nf:6)
  - Input: `params.genomefa`, `params.chrom_size_path`
  - Action: run `count_cg_sites.py` to count CpG sites
  - Output: `genome_cg_info.json`

- Expression FASTQ QC (multi-group): `FASTP_EXPRESSION_MULTI` (script/SeekSoulMethyl/nf/modules/step1.nf:22)
  - Input: paired FASTQs per sample (groups G1/G2/...)
  - Action: `fastp` trimming and QC
  - Output: cleaned `*_expression_clean_R1/2.fastq.gz`, `*.html`, `*.json`

- Methylation FASTQ QC (multi-group): `FASTP_METHYLATION_MULTI` (script/SeekSoulMethyl/nf/modules/step1.nf:88)
  - Input: paired FASTQs per sample (groups G1/G2/...)
  - Action: `fastp` QC (adapter detection disabled, trimming as per pipeline)
  - Output: cleaned `*_methylation_clean_R1/2.fastq.gz`, `*.html`, `*.json`

- RNA alignment and quantification: `SEEKSOULTOOLS_RNA` (script/SeekSoulMethyl/nf/modules/step1.nf:156)
  - Input: cleaned expression R1/R2 lists
  - Action: `seeksoultools rna run` (STAR mapping, counting, filtering, clustering, DE)
  - Output: `Analysis/step3/filtered_feature_bc_matrix/` etc.

- Methylation barcode parsing: `METHYLATION_BARCODE_EXTRACTION` (script/SeekSoulMethyl/nf/modules/step1.nf:206)
  - Input: cleaned methylation R1/R2 lists, whitelist `params.methy_barcode_wl`
  - Action: run `barcode_cs_multi.py` to parse/correct cell barcodes and UMIs by `params.chemistry`; optional `--split_fastq n` to shard reads by first n barcode bases
  - Output: `step1/${sample}_forward_*{1,2}.fq.gz`, `step1/${sample}_reverse_*{1,2}.fq.gz`, `${sample}_methy_summary.json`, optional `${sample}_barcode_stats.txt`

- Build forward/reverse pairing lists: `PARSE_FASTQ_FILES` (script/SeekSoulMethyl/nf/modules/step1.nf:255)
  - Input: forward/reverse FASTQ sets
  - Action: scan `step1/` and pair files by identifiers
  - Output: `forward_pairs.txt`, `reverse_pairs.txt`

- Post-barcode-extraction QC: `FASTP_METHYLATION_BARCODE_EXTRACT` (script/SeekSoulMethyl/nf/modules/step1.nf:368)
  - Input: paired sub-FASTQs
  - Action: `fastp` QC
  - Output: per-pair `*.html`, `*.json`

### Step 2: Bismark alignment and BAM sorting (step2.nf)
- Forward-strand alignment: `BISMARK_ALIGNMENT_FORWARD` (script/SeekSoulMethyl/nf/modules/step2.nf:3)
  - Key flags: `--add_barcode`, `--add_umi`; max insert size `-X 1000`
  - Output: `*_bismark_bt2_pe.bam`, `*_bismark_bt2_PE_report.txt`

- Reverse (PBAT) alignment: `BISMARK_ALIGNMENT_REVERSE` (script/SeekSoulMethyl/nf/modules/step2.nf:31)
  - Key flags: `--pbat`, `--add_barcode`, `--add_umi`
  - Output: same as above

- Sort by read name: `SORT_BAM_BY_NAME` (script/SeekSoulMethyl/nf/modules/step2.nf:61)
  - Action: `samtools sort -n`
  - Output: `*_sortbyname.bam`

### Step 3: Per-cell split, ALLC generation/merge, dataset building (step3.nf)
- Split BAMs by cell barcode: `SPLIT_BAM_FILES` (script/SeekSoulMethyl/nf/modules/step3.nf:1)
  - Input: name-sorted BAM and GEX barcodes `barcodes.tsv.gz`
  - Action: run `step3_split_bams.py` to split by cell barcode and keep shared cells only
  - Output: per-cell BAM dir, `*_filtered_barcode`, `*_filtered_barcode_reads_counts.csv`

- Merge forward/reverse per-cell BAMs and counts: `MERGE_BISMARK_BAM` (script/SeekSoulMethyl/nf/modules/step3.nf:27)
  - Action: `samtools merge -n` per matching cell; merge/deduplicate barcodes and read counts
  - Output: `*_merged_fr_bam/`, `*_merge_filtered_barcode`, `*_merge_filtered_barcode_reads_counts.csv`

- Generate per-cell ALLC: `ALLCOOLS_BAM_TO_ALLC` (script/SeekSoulMethyl/nf/modules/step3.nf:86)
  - Action: run `step3_bam_to_allc.py` (custom ALLCools), UR-based dedup and methylation quantification
  - Output: per-cell `*_allc.gz` and index

- Merge cell metrics: `MERGE_FILTERED_BARCODE_READS_COUNTS` (script/SeekSoulMethyl/nf/modules/step3.nf:114)
  - Action: merge barcodes and read counts, create `${sample}_cells.csv` and `.json`
  - Output: `filtered_barcode`, `filtered_barcode_reads_counts.csv`, `${sample}_cells.{csv,json}`

- Build multi-scale dataset: `ALLCOOLS_GENERATE_DATASETS` (script/SeekSoulMethyl/nf/modules/step3.nf:136)
  - Action: `allcools generate-dataset` for chrom10k/20k/50k/100k/500k/1M, metrics like `count` and `hypo-score`
  - Output: `${sample}.mcds`

- Merge ALLC (when sharded): `ALLCOOLS_SUBMERGE` (script/SeekSoulMethyl/nf/modules/step3.nf:188), `ALLCOOLS_MERGE` (script/SeekSoulMethyl/nf/modules/step3.nf:223)
  - Action: merge per-shard/per-sample ALLCs
  - Output: `${sample}_merge_allc.gz` and index

- Extract CG context ALLC: `ALLCOOLS_EXTRACT` (script/SeekSoulMethyl/nf/modules/step3.nf:252)
  - Action: `allcools extract-allc --mc_contexts CGN --strandness merge`
  - Output: `*.CGN-Merge*`

(Methylation-only workflow `methy_only.nf` additionally includes `COUNTS_MAPPED_READS` and `ESTIMATED_CELLS` for read-count-based cell estimation and barcode filtering, see script/SeekSoulMethyl/nf/modules/utils.nf:1 and :17)

### Step 4: Summary, DR & joint report (step4.nf)
- Methylation summary: `METHYLATION_SUMMARY` (script/SeekSoulMethyl/nf/modules/step4.nf:2)
  - Action: `step4_wgs_summary.py` aggregates Bismark reports, cell metrics and CpG stats to produce `${sample}_methy_summary.json` and `${sample}_wgs_summary.csv`

- LSI/PCA clustering and visualization: `METHYLATION_LSI_PCA_CLUSTERING` (script/SeekSoulMethyl/nf/modules/step4.nf:26)
  - Action: `step4_allcools_PCA_cluster.py --var_dim chrom20k --reduc lsi`
  - Output: `*.h5ad`, `*.pdf`, `*.png`

- RNA+MET joint report: `MULTI_REPORT` (script/SeekSoulMethyl/nf/modules/step4.nf:52)
  - Action: `step4_report_rna_met.py` integrates RNA and methylation outputs
  - Output: `${sample}_rna_methyl_report.html`, `${sample}_rna_met.json`

## Outputs (by outdir structure)
- `fastp/`: QC reports for raw and post-barcode FASTQs (html/json)
- `${sample}/${sample}_exp/`: RNA analysis directory with filtered matrix, clustering, DE results
- `${sample}/${sample}_methy/step1/`: barcode-parsed and sharded FASTQs
- `${sample}/${sample}_methy/step2/`: Bismark BAMs and reports
- `${sample}/${sample}_methy/step3/`:
  - `split_bams/` and `split_bams/merged/`: per-cell BAMs, merged BAMs, barcode counts
  - `allcools/` and `allcools_generate_datasets/`: per-cell ALLCs and `${sample}.mcds`
  - `${sample}_merge_allc.gz`, `*.CGN-Merge*`
- `${sample}/${sample}_methy/step4/`: clustering plots and `*.h5ad`
- `${sample}/`:
  - `${sample}_methy_summary.json`, `${sample}_wgs_summary.csv`
  - `${sample}_rna_methyl_report.html` (if running the main workflow)
- Nextflow run artifacts (as configured by `-c nf/nextflow.config`): `execution_report.html`, `execution_timeline.html`, `pipeline_dag.html`, `execution_trace.txt` (script/SeekSoulMethyl/nf/nextflow.config:17)

## Methylation-only workflow (test version, currently not recommended for use)
Use the simplified workflow when you only have methylation reads:
```bash
nextflow run SeekSoulMethyl/nf/methy_only.nf \
  --outdir /path/to/results \
  --samplesheet samplelist.csv \
  -w /path/to/work \
  -c SeekSoulMethyl/nf/nextflow.config \
  -profile aliyun_k8s \
  --database_dir /path/to/reference \
  --split_fastq 4 \
  --filter_ch 2 \
  --chemistry DD-MET3
```

## Key parameters and tips
- `--database_dir`: reference directory with `fasta/genome.fa`, `genes/genes.gtf`, `star/`, `bed/chr_nochrM.bed` (script/SeekSoulMethyl/nf/main.nf:19)
- `--chemistry`: `DD-MET3` or `DD-MET5`; also sets RNA chemistry and barcode whitelist (script/SeekSoulMethyl/nf/main.nf:27)
- `--split_fastq`: shard by the first n barcode bases (default 4; increases parallelism but adds scheduling/merge overhead) (script/SeekSoulMethyl/nf/main.nf:36)
- `--filter_ch`: filter reads with > n CH methylation sites (default 2) (script/SeekSoulMethyl/nf/modules/step1.nf:241)
- The samplesheet header must be `sample_id` (script/SeekSoulMethyl/nf/main.nf:116)

## Execution environment and resources
- Recommended container image: `registry-vpc.cn-beijing.aliyuncs.com/seekgene/seeksoulmethyl:1.1.0` (script/SeekSoulMethyl/nf/nextflow.config:69)
- Choose `-profile aliyun_k8s` or `-profile docker`, and edit `nf/nextflow.config` for your infra
- Heavy steps: Bismark/ALLCools need substantial CPU/RAM; defaults are set in `withName` blocks, scale up if needed (script/SeekSoulMethyl/nf/nextflow.config:313)


## FAQ
- Samplesheet parsing error: ensure first column is `sample_id`, use absolute paths (script/SeekSoulMethyl/nf/main.nf:111)
- Missing `${sample}.mcds`: check `ALLCOOLS_BAM_TO_ALLC` produced per-cell `*_allc.gz` and `chrom_size_path` is correct
- Stuck at Bismark: verify reference indices and that `params.bismark_ref` is visible in the container
- Resume runs: use `-resume` with the same `-w` work directory

## License

This project is licensed under the MIT License - see the LICENSE file for details.
