# SeekSoulMethyl

SeekSoulMethyl is a single-cell transcriptome + methylation analysis pipeline designed to analyze data generated using the Beijing SeekGene BioSciences Co., Ltd. SeekOne DD Single Cell Multiome Methylation + RNA kit.

## Data Structure

SeekOne DD Single Cell Multiome Methylation + RNA kit comes in two chemistries. DD-MET3 (dual-label) means the RNA and DNA methylation data barcodes are different for the same cell, and the RNA library is a 3′-end transcriptome library. DD-MET5 (single-label) means the RNA and DNA methylation data barcodes are the same for the same cell, and the RNA library is a 5′-end transcriptome library. Below we describe the DNA methylation library structures for both chemistries.

**DD-MET3 Methylation Library Structure**

<figure style="text-align: center;">
<img src="./docs/images/DD-MET3_library_structure.png" alt="DD-MET3 Methylation Library" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 1. DD-MET3 methylation library structure</figcaption>
</figure>

Structure notes:

- SP1/SP2: Adapter sequences
- barcode: 17 bp cell barcode
- 7F: 7 bp linker sequence
- 17L: 17 bp fixed sequence **C**gt**CC**gt**C**gttg**C**t**C**gt
- ME: 19 bp fixed sequence AGATGTGTATAAGAGA**C**AG
- 9 bp: extension sequence from the Tn5 insertion fragment

**DD-MET5 Methylation Library Structure**

<figure style="text-align: center;">
<img src="./docs/images/DD-MET5_library_structure.png" alt="DD-MET5 Methylation Library" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 2. DD-MET5 methylation library structure</figcaption>
</figure>

Structure notes:

- SP1/SP2: Adapter sequences
- barcode: 17 bp cell barcode
- UMI: 12 bp UMI sequence
- TSO: 13 bp TSO sequence TTT**C**TTATATGGG
- 17L: 17 bp fixed sequence **C**gt**CC**gt**C**gttg**C**t**C**gt
- ME: 19 bp fixed sequence AGATGTGTATAAGAGA**C**AG
- 9 bp: extension sequence from the Tn5 insertion fragment

Since the enzymatic treatment converts unmethylated cytosines (C) to thymines (T), the C bases in SP1 and SP2 are methylated to prevent errors in the sequencing adapters during this conversion. Furthermore, the barcodes used for methylation data do not contain any C bases. In contrast, the C bases in 7F, 17L, and ME are not methylated and will be converted to T during the enzymatic process; we use these fixed sequences to calculate the C-to-T conversion rate.

## Installation

1. Clone the repository:

```bash
git clone https://github.com/seekgene/SeekSoulMethyl.git
cd SeekSoulMethyl
```

1. Create and activate conda environment:

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

1. Install the package:

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

- `nf/main.nf`: Top-level entry. Select sub-workflow via `--workflow` (`rna_met`, `methy_only`, `force_cell`).
- `nf/subworkflows/`: Workflow definitions:
  - `rna_met/main.nf`: Transcriptome + methylation end-to-end processing.
  - `met_only/main.nf`: Methylation-only workflow.
  - `forcecell/main.nf`: Force-cell workflow (recomputes/updates results using previous outputs).
- `nf/modules/`: Step-wise process modules:
  - `step1.nf` preprocessing, QC, barcode extraction, transcriptome analysis.
  - `step2.nf` Bismark alignment and BAM sorting.
  - `step3.nf` per-cell BAM splitting, ALLC generation/merge, multi-scale datasets.
  - `step4.nf` summaries, dimensionality reduction, joint report.
  - `utils.nf` helpers for methylation-only workflow (reads counting and cell estimation).
- `nf/bin/`: Helper scripts and resources (e.g., barcode whitelists).
- `nf/nextflow.config`: Executors and resource configuration.
- `nf/nextflow_schema.json`: Pipeline parameter schema.
- `sc_methy_workflow.sh`: Shell script to run the dual-omics analysis pipeline.

We provide two methods for data analysis:

1. **Shell Script**: Run the analysis pipeline directly via the `sc_methy_workflow.sh` script.
2. **Nextflow Pipeline**: Run the Nextflow pipeline via `nf/main.nf`.

Details of both methods are provided below.

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

If a sample has multiple datasets, separate file paths with commas. Ensure that the FASTQ files are listed in the correct order for each dataset.

```shell
bash sc_methy_workflow.sh \
/path/to/WTJW969_E_L003_R1.fq.gz,/path/to/WTJW969_E_L004_R1.fq.gz \
/path/to/WTJW969_E_L003_R2.fq.gz,/path/to/WTJW969_E_L004_R2.fq.gz \
/path/to/WTJW969_Met_L000_R1.fq.gz,/path/to/WTJW969_Met_L001_R1.fq.gz,/path/to/WTJW969_Met_L002_R1.fq.gz,/path/to/WTJW969_Met_L003_R1.fq.gz,/path/to/WTJW969_Met_L004_R1.fq.gz \
/path/to/WTJW969_Met_L000_R2.fq.gz,/path/to/WTJW969_Met_L001_R2.fq.gz,/path/to/WTJW969_Met_L002_R2.fq.gz,/path/to/WTJW969_Met_L003_R2.fq.gz,/path/to/WTJW969_Met_L004_R2.fq.gz \
--sample WTJW969 \
--outdir /path/to/results \
--database_dir /path/to/human-reference-GRCh38 \
--chemistry DD-MET5 \
--core 64 \
--filter_ch 2
```

### Input Parameters:

- **$1**: Single-cell transcriptome Read1 fastq file path.
- **$2**: Single-cell transcriptome Read2 fastq file path.
- **$3**: Single-cell methylation Read1 fastq file path.
- **$4**: Single-cell methylation Read2 fastq file path.
- **sample**: Sample name.
- **outdir**: Output directory path.
- **database\_dir**: Reference genome database path.
- **chemistry**: Methylation chemistry (DD-MET3 or DD-MET5). DD-MET3 is a dual-label dataset, meaning the RNA and DNA methylation data barcodes are different for the same cell, while DD-MET5 is single-label, meaning the RNA and DNA methylation data barcodes are the same for the same cell.
- **core**: Number of CPU cores.
- **filter\_ch**: Filter reads that contain > n CH methylation sites. If you do not want to enable this filter, set filter\_ch to 0.

## Process Details

### Transcriptome Processing Workflow

Transcriptome data is analyzed using [SeekSoulTools](https://github.com/seekgene/SeekSoulMethyl/tree/nf_rna_methy/dependence/seeksoultools). See the official [Algorithms Overview](http://seeksoul.seekgene.com/en/v1.3.0/2.tutorial/1.rna/4.description.html)for detailed steps. Cells used in the downstream methylation library are determined based on the transcriptome library cell barcodes.

### Methylation Processing Workflow

#### Step 1: Preprocessing and Barcode Parsing

**Barcode extraction and correction**

Based on the designed structure, we locate the barcode in the read and extract the corresponding sequence. If the extracted barcode is in the whitelist, it is counted as a valid barcode; otherwise, SeekSoulMethyl attempts barcode correction, if the barcode has a one-base mismatch (Hamming distance = 1) from an entry in the whitelist:

- If exactly one whitelist candidate matches: correct the invalid barcode to that whitelist barcode.
- If multiple whitelist candidates match: correct to the candidate supported by the highest read count.

If correction fails, the read is discarded and considered a final invalid barcode read.

**UMI extraction**

UMI positions are read from the designed structure and extracted without correction.

**Forward and Reverse reads determination**

From the positions corresponding to 17L and ME, there are 7 bases that can be C or converted T (highlighted below). We use the first and the last two C/T positions to determine forward vs. reverse reads: if all three positions are C, it indicates a reverse read; otherwise, it is a forward read.

- Forward:  **T**gt**TT**gt**T**gttg**T**t**T**gtAGATGTGTATAAGAGA**T**
- Reverse:  **C**gt**CC**gt**C**gttg**C**t**C**gtAGATGTGTATAAGAGA**C**

Reverse reads correspond to CTOT/CTOB (reverse complement of the original strand) in methylation terminology; forward reads correspond to OT/OB (original strand).
The "forward" or "reverse" determination is annotated in the read name.

<figure style="text-align: center;">
<img src="./docs/images/fr_determinate.png" alt="Forward and reverse reads determination" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 3. Forward and reverse reads determination</figcaption>
</figure>

**C–T conversion rate**

We calculate the C-to-T conversion rate using the original C positions within 17L and ME sequences. Since these are fixed sequences prone to sequencing errors, we restrict the calculation to reads with verified structures:

- In DD-MET3, the 7F sequence must be `TTGCTGT` or `TTGTTGT`, the sequence spanning 17L and ME must be `GTAGATGTGTATAAGAGA`, and the bases at first and the last two original C positions must be T.
- In DD-MET5, the sequence spanning 17L and ME must be `GTAGATGTGTATAAGAGA`, and the bases at the first and last two original C positions must be T.

For the retained reads, we extract the bases at the corresponding positions to calculate the C-to-T conversion rate:

<figure style="text-align: center;">
<img src="./docs/images/CT_conversion.png" alt="CT conversion rate" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 4. CT conversion rate</figcaption>
</figure>

> \[!NOTE]
> The above filtering steps are used only for calculating the C-to-T conversion rate; reads that do not meet these criteria are not filtered out from the final output FASTQ files.

**Artifact removal**

Remove TSO/7F linker, 17L and ME sequences from Read1 according to their predefined positions in the library design.

**Adapter trimming**

Use cutadapt to remove ME adapter sequences introduced by R1/R2 read-through events (overlapping paired-end reads).

**Trim 9 bp gaps introduced by Tn5 transposase**

After removing adapters and other artificial sequences, we additionally trim the 9 bp gaps flanking the inserted fragment that are introduced by Tn5 transposition. These 9 bp regions can carry artificial methylation and spuriously elevate CH methylation adjacent to the insert boundaries, so they are removed prior to downstream analysis.

**Filter reads with too many non-CpG methylated C bases (optional)**

Filter based on the number of non-CpG methylated C bases in a read pair. By default, pairs with > 2 non-CpG methylated Cs detected in read1/read2 are removed. If you do not want to enable this filter, set filter\_ch to 0.

> \[!NOTE]
> This filtering strategy is based on findings by Lu et al. \[1], which suggest that nicks in synthesized adapters can trigger Bst polymerase nick translation. This activity incorporates 5-methyl-dCTPs, leading to completely unconverted reads that appear as artificial methylation signals.

**Filter too short reads**
After the preceding filtering and adapter trimming steps, if the length of R1 in a read pair is less than 20 bp or the length of R2 is less than 60 bp, the read pair is filtered out.

#### Step 2: Bismark alignment and sorting by name

**Alignment and tagging**

We use Bismark for methylation alignment. Our modified [Bismark](https://github.com/seekgene/Bismark.git) adds `--add_barcode` and `--add_umi` to tag BAM files by read name with CB (error-corrected barcode) and UR (raw UMI). For forward reads, we use `-X 1000` to allow insert sizes up to 1000 bp; for reverse reads, we use `--pbat` and `-X 1000`.
After alignment, sort BAMs by read name using `samtools sort -n`; the name-sorted BAMs serve as inputs for downstream analysis.

#### Step 3: ALLCools analysis

**Split by cell barcode**

Split name-sorted BAMs by RNA-derived cell barcodes into per-cell BAM files, each containing uniquely mapped reads for one cell.

**Generate ALLC files**

Sort each per-cell BAM by position and convert to ALLC using ALLCools `bam-to-allc`. Our modified ALLCools performs UR-tag-based UMI correction and deduplication per C site.

<figure style="text-align: center;">
<img src="./docs/images/umi_correction_detailed_diagram_en.png" alt="UMI correction detailed diagram" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 5. UMI correction detailed diagram</figcaption>
</figure>

See the [SeekGene ALLCools repository](https://github.com/seekgene/ALLCools) for details.

**Generate MCDS**

Run `allcools generate-dataset` to bin the genome (chrom10k/20k/50k/100k/500k/1M/geneslop2k) and compute per-cell methylation matrices. Geneslop2k bins are defined as 2k bp flanking each gene.

#### Step 4: Reduction and clustering

By default, perform dimensionality reduction with LSI on chrom20k bins using ALLCools, followed by UMAP visualization.

### System Requirements

If you use `sc_methy_workflow.sh`, the system requirements are as follows:

- **CPU**: 64 cores (recommended)
- **Memory**: 128GB RAM (recommended)
- **OS**: Linux (recommended Ubuntu 18.04+ or CentOS 7+)

## Running with Nextflow (Recommended)

If you want to process samples in batch and get pipeline-level reports, use Nextflow:

<figure style="text-align: center;">
<img src="./docs/images/nf_SeekSoulMethyl_workflow.png" alt="SeekSoulMethyl Pipeline" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">Figure 6. SeekSoulMethyl nextflow pipeline workflow</figcaption>
</figure>

1. Install nextflow:

```bash
conda install -n seeksoulmethyl -c bioconda nextflow
```

1. Prepare your input samplesheet

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

If a single sample has multiple datasets and the transcriptome and methylation FASTQ counts do not match (e.g., WTJW969), add multiple rows to the samplesheet, with each row representing one dataset.

```text
sample_id,expression_r1,expression_r2,methylation_r1,methylation_r2
WTJW969,/path/to/WTJW969_E_L003_R1.fq.gz,/path/to/WTJW969_E_L003_R2.fq.gz,/path/to/WTJW969_Met_L000_R1.fq.gz,/path/to/WTJW969_Met_L000_R2.fq.gz
WTJW969,/path/to/WTJW969_E_L004_R1.fq.gz,/path/to/WTJW969_E_L004_R2.fq.gz,/path/to/WTJW969_Met_L001_R1.fq.gz,/path/to/WTJW969_Met_L001_R2.fq.gz
WTJW969,,,/path/to/WTJW969_Met_L002_R1.fq.gz,/path/to/WTJW969_Met_L002_R2.fq.gz
WTJW969,,,/path/to/WTJW969_Met_L003_R1.fq.gz,/path/to/WTJW969_Met_L003_R2.fq.gz
WTJW969,,,/path/to/WTJW969_Met_L004_R1.fq.gz,/path/to/WTJW969_Met_L004_R2.fq.gz
```

1. Run the pipeline:

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
# -c: nextflow configuration file. Must be modified according to your server configuration.
# -profile: nextflow profile for aliyun k8s
# --database_dir: reference genome database directory
# --split_fastq: To speed up the analysis process, split fastq according to first n bases of barcode. Default is 4.
# --filter_ch: filter reads that contain > 2 CH methylation sites.
# --chemistry: methylation chemistry (DD-MET3 or DD-MET5)
```

### Notes on nextflow\.config

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

- Nextflow configuration & profiles: <https://www.nextflow.io/docs/latest/config.html>
- Executors (local, Slurm, K8s, etc.): <https://www.nextflow.io/docs/latest/executor.html>
- Kubernetes guide: <https://www.nextflow.io/docs/latest/kubernetes.html>

## Nextflow Step-by-Step Details

This section describes each Nextflow process with inputs, core logic, key parameters, and outputs for troubleshooting and interpretation. The workflow is defined in `nf/main.nf` and processes are implemented in `nf/modules/*.nf`.

### Step 1: Preprocessing and Barcode Parsing (step1.nf)

- Compute genome-wide CpG count: `COMPUTE_CPG_SITES`
  - Input: `params.genomefa`, `params.chrom_size_path`
  - Action: run `count_cg_sites.py` to count CpG sites
  - Output: `genome_cg_info.json`
- Expression FASTQ QC (multi-group): `FASTP_EXPRESSION_MULTI`
  - Input: paired FASTQs per sample (groups G1/G2/...)
  - Action: `fastp` trimming and QC
  - Output: cleaned `*_expression_clean_R1/2.fastq.gz`, `*.html`, `*.json`
- Methylation FASTQ QC (multi-group): `FASTP_METHYLATION_MULTI`
  - Input: paired FASTQs per sample (groups G1/G2/...)
  - Action: `fastp` QC (adapter detection disabled, trimming as per pipeline)
  - Output: cleaned `*_methylation_clean_R1/2.fastq.gz`, `*.html`, `*.json`
- RNA alignment and quantification: `SEEKSOULTOOLS_RNA`
  - Input: cleaned expression R1/R2 lists
  - Action: `seeksoultools rna run` (STAR mapping, counting, filtering, clustering, DE)
  - Output: `Analysis/step3/filtered_feature_bc_matrix/` etc.
- Methylation barcode parsing: `METHYLATION_BARCODE_EXTRACTION`
  - Input: cleaned methylation R1/R2 lists, whitelist `params.methy_barcode_wl`
  - Action: run `barcode_cs_multi.py` to parse/correct cell barcodes and UMIs by `params.chemistry`; optional `--split_fastq n` to shard reads by first n barcode bases
  - Output: `step1/${sample}_forward_*{1,2}.fq.gz`, `step1/${sample}_reverse_*{1,2}.fq.gz`, `${sample}_methy_summary.json`, optional `${sample}_barcode_stats.txt`
- Build forward/reverse pairing lists: `PARSE_FASTQ_FILES`
  - Input: forward/reverse FASTQ sets
  - Action: scan `step1/` and pair files by identifiers
  - Output: `forward_pairs.txt`, `reverse_pairs.txt`
- Post-barcode-extraction QC: `FASTP_METHYLATION_BARCODE_EXTRACT`
  - Input: paired sub-FASTQs
  - Action: `fastp` QC
  - Output: per-pair `*.html`, `*.json`

### Step 2: Bismark alignment and BAM sorting (step2.nf)

- Forward-strand alignment: `BISMARK_ALIGNMENT_FORWARD`
  - Key flags: `--add_barcode`, `--add_umi`; max insert size `-X 1000`
  - Output: `*_bismark_bt2_pe.bam`, `*_bismark_bt2_PE_report.txt`
- Reverse (PBAT) alignment: `BISMARK_ALIGNMENT_REVERSE`
  - Key flags: `--pbat`, `--add_barcode`, `--add_umi`
  - Output: same as above
- Sort by read name: `SORT_BAM_BY_NAME`
  - Action: `samtools sort -n`
  - Output: `*_sortbyname.bam`

### Step 3: Per-cell split, ALLC generation/merge, dataset building (step3.nf)

- Split BAMs by cell barcode: `SPLIT_BAM_FILES`
  - Input: name-sorted BAM and GEX barcodes `barcodes.tsv.gz`
  - Action: run `step3_split_bams.py` to split by cell barcode and keep shared cells only
  - Output: per-cell BAM dir, `*_filtered_barcode`, `*_filtered_barcode_reads_counts.csv`
- Merge forward/reverse per-cell BAMs and counts: `MERGE_BISMARK_BAM`
  - Action: `samtools merge -n` per matching cell; merge/deduplicate barcodes and read counts
  - Output: `*_merged_fr_bam/`, `*_merge_filtered_barcode`, `*_merge_filtered_barcode_reads_counts.csv`
- Generate per-cell ALLC: `ALLCOOLS_BAM_TO_ALLC`
  - Action: run `step3_bam_to_allc.py` (custom ALLCools), UR-based dedup and methylation quantification
  - Output: per-cell `*_allc.gz` and index
- Merge cell metrics: `MERGE_FILTERED_BARCODE_READS_COUNTS`
  - Action: merge barcodes and read counts, create `${sample}_cells.csv` and `.json`
  - Output: `filtered_barcode`, `filtered_barcode_reads_counts.csv`, `${sample}_cells.{csv,json}`
- Build multi-scale dataset: `ALLCOOLS_GENERATE_DATASETS`
  - Action: `allcools generate-dataset` for chrom10k/20k/50k/100k/500k/1M, metrics like `count` and `hypo-score`
  - Output: `${sample}.mcds`
- Merge ALLC (when sharded): `ALLCOOLS_SUBMERGE`, `ALLCOOLS_MERGE`
  - Action: merge per-shard/per-sample ALLCs
  - Output: `${sample}_merge_allc.gz` and index
- Extract CG context ALLC: `ALLCOOLS_EXTRACT`
  - Action: `allcools extract-allc --mc_contexts CGN --strandness merge`
  - Output: `*.CGN-Merge*`

(Methylation-only workflow `methy_only.nf` additionally includes `COUNTS_MAPPED_READS` and `ESTIMATED_CELLS` for read-count-based cell estimation and barcode filtering)

### Step 4: Summary, DR & joint report (step4.nf)

- Methylation summary: `METHYLATION_SUMMARY`
  - Action: `step4_wgs_summary.py` aggregates Bismark reports, cell metrics and CpG stats to produce `${sample}_methy_summary.json` and `${sample}_wgs_summary.csv`
- LSI/PCA clustering and visualization: `METHYLATION_LSI_PCA_CLUSTERING`
  - Action: `step4_allcools_PCA_cluster.py --var_dim chrom20k --reduc lsi`
  - Output: `*.h5ad`, `*.pdf`, `*.png`
- Transcriptome+Methylation joint report: `MULTI_REPORT`
  - Action: `step4_report_rna_met.py` integrates transcriptome and methylation outputs
  - Output: `${sample}_rna_methyl_report.html`, `${sample}_rna_met.json`

## Outputs (by outdir structresourceure)

- `fastp/`: QC reports for raw and post-barcode FASTQs (html/json)
- `${sample}/${sample}_exp/`: Transcriptome analysis directory with filtered matrix, clustering, DE results
- `${sample}/${sample}_methy/step1/`: Barcode-parsed and sharded FASTQs
- `${sample}/${sample}_methy/step2/`: Bismark BAMs and reports
- `${sample}/${sample}_methy/step3/`:
  - `split_bams/` and `split_bams/merged/`: Per-cell BAMs, merged BAMs, barcode counts
  - `allcools/` and `allcools_generate_datasets/`: Per-cell ALLCs and `${sample}.mcds`
  - `${sample}_merge_allc.gz`, `*.CGN-Merge*`
- `${sample}/${sample}_methy/step4/`: Clustering plots and `*.h5ad`
- `${sample}/`:
  - `${sample}_methy_summary.json`, `${sample}_wgs_summary.csv`
  - `${sample}_rna_methyl_report.html` (if running the main workflow)
- Nextflow run artifacts (as configured by `-c nf/nextflow.config`): `execution_report.html`, `execution_timeline.html`, `pipeline_dag.html`, `execution_trace.txt`

## Methylation-only workflow (test version, currently not recommended for use)

Use the simplified workflow when you only have methylation reads:

```bash
nextflow run SeekSoulMethyl/nf/main.nf \
  --workflow methy_only \
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

- `--database_dir`: reference directory with `fasta/genome.fa`, `genes/genes.gtf`, `star/`, `bed/chr_nochrM.bed`
- `--chemistry`: `DD-MET3` or `DD-MET5`; also sets transcriptome chemistry and barcode whitelist
- `--split_fastq`: shard by the first n barcode bases (default 4; increases parallelism but adds scheduling/merge overhead)
- `--filter_ch`: filter reads with > n CH methylation sites (default 2).  If you do not want to enable this filter, set `filter_ch` to 0.
- The samplesheet header must be `sample_id`

## Execution environment and resources

- Recommended container image: `registry.cn-beijing.aliyuncs.com/seekgene/seeksoulmethyl:1.1.2`
- Choose `-profile aliyun_k8s` or `-profile docker`, and edit `nf/nextflow.config` for your infra
- Heavy steps: Bismark/ALLCools need substantial CPU/RAM; defaults are set in `withName` blocks, scale up if needed

## FAQ

- Samplesheet parsing error: ensure first column is `sample_id`, use absolute paths
- Missing `${sample}.mcds`: check `ALLCOOLS_BAM_TO_ALLC` produced per-cell `*_allc.gz` and `chrom_size_path` is correct
- Stuck at Bismark: verify reference indices and that `params.bismark_ref` is visible in the container
- Resume runs: use `-resume` with the same `-w` work directory

## Documentation

Tutorials and reference materials are available under the [docs/](docs) directory:

- [How to build the reference genome database (`--database_dir`)](docs/How_to_build_reference_genome.md)
- [How to obtain single-cell BAM files](docs/Obtain_single_cell_bam.md) ([中文](docs/Obtain_single_cell_bam.zh.md))
- [How to deduplicate single-cell BAM files](docs/How_to_deduplicate_single_cell_bam.md) ([中文](docs/How_to_deduplicate_single_cell_bam.zh.md))

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## References

\[1] Lu X, Yuan Y, et al. Improved tagmentation-based whole-genome bisulfite sequencing for input DNA from less than 100 mammalian cells. Epigenomics. 2015;7(1):47-56. doi:10.2217/epi.14.76.

> "Furthermore, by manually checking the reads, we found a part of the reads were completely unconverted. We suspected that a nick in the synthesized adapters will cause the whole fragment displaced with incorporation of 5-methyl-dCTPs due to nick translation activity of Bst polymerase."

