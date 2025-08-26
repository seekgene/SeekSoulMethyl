# SeekSoulMethyl
SeekSoulMethyl is a single-cell transcriptome + methylation dual-omics analysis pipeline designed for analyzing data from Beijing SeekGene Biotechnology Co., Ltd. single-cell transcriptome + methylation dual-omics kit.

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
# Download methylation test data
wget -dc -O XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz"
wget -dc -O XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz.md5"
wget -dc -O XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz"
wget -dc -O XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz.md5"


```

**Note**: This is a small test dataset for pipeline validation. For production analysis, use your own sequencing data.

## Usage

### Activate Environment
```bash
conda activate seeksoulmethyl
```

### Run Dual-omics Analysis
```bash
# sc_methy_workflow.sh can be found in the SeekSoulMethyl directory you cloned
bash sc_methy_workflow.sh \
    sample \
    /path/to/expression_R1.fastq.gz \
    /path/to/expression_R2.fastq.gz \
    /path/to/methy_R1.fastq.gz \
    /path/to/methy_R2.fastq.gz \
    /path/to/output_dir \
    /path/to/database/refdata-cellranger-arc-GRCh38-2024-A \
    64
```

### Input Parameters:
- **$1**: Sample name
- **$2**: Single-cell transcriptome Read1 fastq file
- **$3**: Single-cell transcriptome Read2 fastq file  
- **$4**: Single-cell methylation Read1 fastq file
- **$5**: Single-cell methylation Read2 fastq file
- **$6**: Output directory path
- **$7**: Reference genome database path
- **$8**: Number of CPU cores



## System Requirements

- **CPU**: 64 cores (recommended)
- **Memory**: 128GB RAM (recommended)
- **Storage**: At least 500GB available space
- **OS**: Linux (recommended Ubuntu 18.04+ or CentOS 7+)

## System Configuration

The pipeline automatically optimizes system settings for optimal performance:

- **File Descriptor Limit**: Automatically sets `ulimit -n 99999` to handle high concurrency file operations
- **Samtools Memory**: Unified memory allocation of 1.5G for samtools operations across all system configurations
- **Bismark Parallelization**: Automatically calculates optimal thread allocation (1/8 of total CPU cores)

These optimizations ensure consistent performance across different system configurations and prevent resource-related issues during analysis.

## License

This project is licensed under the MIT License - see the LICENSE file for details.


