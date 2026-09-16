# 输出文件

## 输出文件（按 outdir 结构）

- `fastp/`：原始和 barcode 提取后 FASTQ 的质控报告（html/json）
- `${sample}/${sample}_exp/`：转录组分析目录，包含过滤后的矩阵、聚类、差异表达结果
- `${sample}/${sample}_methy/step1/`：Barcode 解析和分片后的 FASTQ
- `${sample}/${sample}_methy/step2/`：Bismark BAM 和报告
- `${sample}/${sample}_methy/step3/`：
  - `split_bams/` 和 `split_bams/merged/`：单细胞 BAM、合并后的 BAM、barcode 计数
  - `allcools/` 和 `allcools_generate_datasets/`：单细胞 ALLC 和 `${sample}.mcds`
  - `${sample}_merge_allc.gz`、`*.CGN-Merge*`
- `${sample}/${sample}_methy/step4/`：聚类图和 `*.h5ad`
- `${sample}/`：
  - `${sample}_methy_summary.json`、`${sample}_wgs_summary.csv`
  - `${sample}_rna_methyl_report.html`（运行主工作流时生成）
- Nextflow 运行产物（由 `-c nf/nextflow.config` 配置）：`execution_report.html`、`execution_timeline.html`、`pipeline_dag.html`、`execution_trace.txt`
