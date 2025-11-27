# SeekSoulMethyl
SeekSoulMethyl is a single-cell transcriptome + methylation dual-omics analysis pipeline designed for analyzing data from Beijing SeekGene Biotechnology Co., Ltd. single-cell transcriptome + methylation dual-omics kit.

## 数据结构
RNA-MET数据分为双标签数据和单标签数据。双标签数据对应的chemistry为DD-MET3，即同一个细胞的RNA和DNA甲基化数据的barcode标签不一致，且RNA转录组属于3‘转录组数据。单标签的数据对应的chemistry是DD-MET5，同一个细胞的RNA和DNA甲基化数据的barcode标签一致，且RNA转录组属于5‘转录组数据。下面详细描述一下双标签数据和单标签数据的DNA甲基化文库结构。

**DD-MET3 Methlytion Library Structure**
<figure style="text-align: center;">
<img src="./docs/DD-MET3_library_structure.png" alt="DD-MET3 Methylation Library" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">图 1. DD-MET3 甲基化文库结构示意图</figcaption>
</figure>

结构说明：
- SP1/SP2：接头序列；
- barcode：17bp的cell barcode；
- 7F：7bp的连接序列；
- 17L和ME：17bp的固定序列 <span style="color:#43a047;">Y</span>gt<span style="color:#43a047;">Y</span><span style="color:#43a047;">Y</span>gt<span style="color:#43a047;">Y</span>gttg<span style="color:#43a047;">Y</span>t<span style="color:#43a047;">Y</span>gt；
- ME: 19bp的固定序列AGATGTGTATAAGAGA<span style="color:#43a047;">Y</span>AG；
- 9bp：为Tn5酶打断后的插入片段的延伸序列。

**DD-MET5 Methlytion Library Structure**
<figure style="text-align: center;">
<img src="./docs/DD-MET5_library_structure.png" alt="DD-MET5 Methylation Library" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">图 2. DD-MET5 甲基化文库结构示意图</figcaption>
</figure>

结构说明：
- SP1/SP2：接头序列；
- barcode：17bp的cell barcode；
- UMI: 12bp的UMI序列；
- TSO：13bpTSO序列TTT<span style="color:#43a047;">Y</span>TTATATGGG；
- 17L和ME：17bp的固定序列 <span style="color:#43a047;">Y</span>gt<span style="color:#43a047;">Y</span><span style="color:#43a047;">Y</span>gt<span style="color:#43a047;">Y</span>gttg<span style="color:#43a047;">Y</span>t<span style="color:#43a047;">Y</span>gt；
- ME: 19bp的固定序列AGATGTGTATAAGAGA<span style="color:#43a047;">Y</span>AG；
- 9bp：为Tn5酶打断后的插入片段的延伸序列。


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
### RNA数据处理流程
RNA数据使用[SeekSoulTools](http://seeksoul.seekgene.com/en/v1.3.0/2.tutorial/1.rna/4.description.html)进行分析。具体的处理步骤请见官网的[Algorithms Overview](http://seeksoul.seekgene.com/en/v1.3.0/2.tutorial/1.rna/4.description.html#)。后续甲基化文库的细胞基于转录组文库判定的细胞barcode进行分析。

### 甲基化数据处理流程
#### Step 1: Preprocessing and Barcode Parsing
**Barcode Extraction and Correction**
根据结构设计，确定barcode所在序列位置，取出相应序列。当取出的barcode序列在白名单中时，我们认为它是有效barcode，计入有效barcode的reads数量；当barcode不在白名单中时，我们认为它是无效barcode。
测序过程中，有一定几率发生测序错误。在提供有白名单的情况下，SeekSoulTools可以尝试barcode矫正。在启用矫正时，当无效barcode一个碱基错配（一个hamming distance）的序列存在于白名单中：
* 只有唯一一个序列存在于白名单中：我们将这个无效barcode矫正为白名单中barcode；
* 有多个序列存在于白名单中：我们将这个无效barcode改为read支持数量最多的序列。

**UMI Extraction**
根据结构设计，确定UMI所在序列位置，取出相应序列，不进行任何矫正。

**Forward和Reverse链判定**
按照位置，提取reads位置上对应17L和ME序列，共有7个碱基包含C或者转化后的T（下面标绿色）。使用其中最后两个位置上的C/T来判定，如果这两个位置全部为C则为reverse链reads，否则为forward链reads。
Forward: <span style="color:#43a047;">T</span>gt<span style="color:#43a047;">TT</span>gt<span style="color:#43a047;">T</span>gttg<span style="color:#43a047;">T</span>t<span style="color:#43a047;">T</span>gtAGATGTGTATAAGAGA<span style="color:#43a047;">T</span>
Reverse: <span style="color:#43a047;">C</span>gt<span style="color:#43a047;">CC</span>gt<span style="color:#43a047;">C</span>gttg<span style="color:#43a047;">C</span>t<span style="color:#43a047;">C</span>gtAGATGTGTATAAGAGA<span style="color:#43a047;">C</span>
Reverse链就是我们在甲基化数据中常说的CTOT和CTOB链，forward链就是我们在甲基化数据中常说的OT和OB链。

**C-T转化率**
除了上述用于判定Forward和Reverse链判定的碱基，剩余的5个C/T碱基用于计算C-T转化率，公式如下：
<div style="font-size: 0.8em;">

$$
    CT conversion = （Forward reads中5个位置上“T”碱基的总数 / Forward reads数*5）
$$

</div>

**Filter Reads(Optional)**
是通过统计Reads中非CpG的甲基化C的数量来过滤的，默认为>2，即read1/read2中检测到的非CpG的甲基化C的的数量大于2，该read pair被去掉。

#### Step 2: Bismark Alignment and Sorting by name
**Bismark Alignment and Sorting by name**
使用Bismark进行甲基化数据的比对。我们修改了Bismark的代码，增加了--add_barcode和--add_umi参数，根据read name给bam加上CB 和 UR tag。CB为error-corrected barcode，UR为未纠错的 raw UMI。对于Forward链，Bismark中使用-X 1000参数进行比对，允许read1和read2 的insert size为1000bp。对于Reverse链，Bismark中使用--pbat 和-X 1000参数进行比对。
比对完成后，使用samtools sort -n对bam文件进行排序，排序后的bam文件为后续分析的输入。

#### Step 3: ALLCools analysis
**Split by Cell Barcode**
根据转录组识别到的cell barcode，将bam文件拆分成多个小文件，每个小文件包含一个cell的所有uniqly mapped reads。
**Generate allc file**
将每个单细胞的bam文件按照position排序，然后使用ALLCools bam-to-allc工具将每个cell的bam文件转换为allc文件。我们使用的是modified ALLCools，基于 UR-tag对每个C位点的reads进行UMI的矫正和去重，具体矫正和去重操作详见[seekgene仓库里的ALLCools](https://github.com/seekgene/ALLCools)。
**Generate MCDS**
使用allcools generate-dataset 将基因组按照chrom10k/20k/50k/100k/500k/1M bin划分，计算每个细胞在这些bins里的甲基化水平矩阵。

#### Step 4: Reduction and Clustering
使用allcools，默认基于chrom20k bins按照LSI算法对细胞进行降维分析。然后使用UMAP进行可视化。



## System Requirements

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
- Samplesheet parsing error: ensure first column is `sample_id`, use absolute paths or `oss://`, comma-separate multiple files (script/SeekSoulMethyl/nf/main.nf:111)
- Missing `${sample}.mcds`: check `ALLCOOLS_BAM_TO_ALLC` produced per-cell `*_allc.gz` and `chrom_size_path` is correct
- Stuck at Bismark: verify reference indices and that `params.bismark_ref` is visible in the container
- Resume runs: use `-resume` with the same `-w` work directory

## License

This project is licensed under the MIT License - see the LICENSE file for details.
