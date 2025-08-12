# SeekSoulMethyl
SeekSoulMethyl is a single-cell transcriptome + methylation dual-omics analysis pipeline designed for analyzing data from Beijing SeekGene Biotechnology Co., Ltd. single-cell transcriptome + methylation dual-omics kit.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/seekgenebio/SeekSoulMethyl.git
cd SeekSoulMethyl/seeksoultools-1.3.0
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

Note: If you encounter slow download speeds or installation issues with pip packages, you can try using alternative PyPI mirrors in your region:
```bash
# Use PyPI mirrors
pip config set global.index-url https://pypi.org/simple
```

3. Install the package:
```bash
pip install . \
  src/simpleqc/target/wheels/simpleqc-0.1.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl \
  src/search-pattern/target/wheels/search_pattern-0.1.0-py3-none-manylinux_2_5_x86_64.manylinux1_x86_64.whl
```

## Download Reference Database
```bash
wget -dc -O database.tar.gz "http://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo%2Fmethy_exp%2Fv1.1%2Fdatabase.tar.gz?Expires=37754989669&OSSAccessKeyId=LTAIxk81dBP5kdTu&Signature=VO1VlXN97Be4bWlTGbZS28uQ9ZQ%3D"
wget -dc -O database.tar.gz.md5 "http://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo%2Fmethy_exp%2Fv1.1%2Fdatabase.tar.gz.md5?Expires=37754989711&OSSAccessKeyId=LTAIxk81dBP5kdTu&Signature=9EmaZZEwSagnKewUunRJtPQ%2F9fE%3D"
tar -xzf database.tar.gz
```

## Usage

```bash
# Set the path to your SeekSoulMethyl conda installation
export PATH=/path/to/conda/envs/seeksoulmethyl/bin:$PATH
export LD_LIBRARY_PATH=/path/to/conda/envs/seeksoulmethyl/lib:$LD_LIBRARY_PATH
# sc_methy_workflow.sh can be found in the SeekSoulMethyl directory you cloned
sh /path/to/SeekSoulMethyl/script/sc_methy_workflow.sh \
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

- **CPU**: 64 cores
- **Memory**: 128GB RAM
- **Storage**: At least 500GB available space

## Support

For technical support or questions, please contact SeekGene support team or open an issue on GitHub.


