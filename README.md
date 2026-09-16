# SeekSoulMethyl
SeekSoul™ Methyl Tools  (SeekSoulMethyl) is a single-cell transcriptome + methylation analysis pipeline designed to analyze data generated using the Beijing [SeekGene](https://www.seekgene.com/) BioSciences Co., Ltd. SeekOne DD Single Cell Multiome Methylation + RNA kit.

## Release Notes


### v2.2.1 (2026-09) — MCDS merge data-loss fix

v2.2.1 fixes a silent data-loss bug in the methylation `MERGE_MCDS` step. In v2.0.0–v2.2.0, a random subset of cells was silently zeroed out in the large-scale bins — **~9% of cells in `chrom1M` and ~8% in `chrom500k`** — with no error or warning. The default clustering bin `chrom20k` was unaffected.

**Who should upgrade**: if your analysis uses the large-scale bins `chrom1M` / `chrom500k`, upgrade to v2.2.1 and re-run. If you only use the default `chrom20k` clustering, your existing results are intact and no re-run is needed.

Full report: [merge_mcds_data_loss_fix.md](docs/changelogs_detail/2.2.1/merge_mcds_data_loss_fix.md) · [中文](docs/changelogs_detail/2.2.1/merge_mcds_data_loss_fix.zh.md)

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

2. Create and activate conda environment:

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
wget -c -O human-reference-GRCh38.tar.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/v1.1/human-reference-GRCh38.tar.gz"
wget -c -O human-reference-GRCh38.tar.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/v1.1/human-reference-GRCh38.tar.gz.md5"

# Download mouse reference genome (GRCm39)
wget -c -O mouse-reference-GRCm39.tar.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/v1.1/mouse-reference-GRCm39.tar.gz"
wget -c -O mouse-reference-GRCm39.tar.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/v1.1/mouse-reference-GRCm39.tar.gz.md5"

# Extract reference genomes
tar -xzf human-reference-GRCh38.tar.gz
tar -xzf mouse-reference-GRCm39.tar.gz
```

## Download Test Data (Optional)

For testing the pipeline with a small dataset, you can download the tiny test data:

```bash
# Download transcriptome test data
wget -c -O XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/fastq/XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz"
wget -c -O XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/fastq/XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz.md5"
wget -c -O XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/fastq/XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz"
wget -c -O XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/fastq/XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz.md5"

# Download methylation test data
wget -c -O XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz"
wget -c -O XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz.md5"
wget -c -O XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz"
wget -c -O XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz.md5 "https://seekgene-public.oss-cn-beijing.aliyuncs.com/methy_demo/methy_exp/tiny_fastq/XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz.md5"


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

## Documentation


Tutorials and reference materials are available under the [docs/](docs) directory:

- [How to build the reference genome database (`--database_dir`)](docs/tutorials/How_to_build_reference_genome.md)
- [How to obtain single-cell BAM files](docs/tutorials/Obtain_single_cell_bam.md) ([中文](docs/tutorials/Obtain_single_cell_bam.zh.md))
- [How to deduplicate single-cell BAM files](docs/tutorials/How_to_deduplicate_single_cell_bam.md) ([中文](docs/tutorials/How_to_deduplicate_single_cell_bam.zh.md))
- [Usage — dual-omics (shell) & Nextflow](docs/guide/usage.md) ([中文](docs/guide/usage.zh.md))
- [Algorithm and processing details](docs/guide/algorithm.md) ([中文](docs/guide/algorithm.zh.md))
- [Output files](docs/guide/outputs.md) ([中文](docs/guide/outputs.zh.md))
- [FAQ](docs/guide/faq.md) ([中文](docs/guide/faq.zh.md))
- [Changelog details — v2.2.1 MCDS merge data-loss fix](docs/changelogs_detail/2.2.1/merge_mcds_data_loss_fix.md) ([中文](docs/changelogs_detail/2.2.1/merge_mcds_data_loss_fix.zh.md))

## License


This project is licensed under the MIT License - see the LICENSE file for details.

## References


[1] Lu X, Yuan Y, et al. Improved tagmentation-based whole-genome bisulfite sequencing for input DNA from less than 100 mammalian cells. Epigenomics. 2015;7(1):47-56. doi:10.2217/epi.14.76.
> "Furthermore, by manually checking the reads, we found a part of the reads were completely unconverted. We suspected that a nick in the synthesized adapters will cause the whole fragment displaced with incorporation of 5-methyl-dCTPs due to nick translation activity of Bst polymerase."
