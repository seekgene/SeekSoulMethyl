# Output Files

## Outputs (by outdir structure)

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
