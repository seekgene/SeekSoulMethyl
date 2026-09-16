# FAQ

## FAQ

- Samplesheet parsing error: ensure first column is `sample_id`, use absolute paths
- Missing `${sample}.mcds`: check `ALLCOOLS_BAM_TO_ALLC` produced per-cell `*_allc.gz` and `chrom_size_path` is correct
- Stuck at Bismark: verify reference indices and that `params.bismark_ref` is visible in the container
- Resume runs: use `-resume` with the same `-w` work directory
