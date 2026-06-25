# Changelog

All notable changes to this project will be documented in this file.

## [2.2.0] 2026-06-18
### Changed
- Default `chemistry` changed from `DD-MET3` to `DD-MET5` across all workflows (Nextflow and Shell), aligning with SeekOne DD kit defaults.
- Added fail-fast validation for expression/methylation FASTQ R1/R2 pair counts to avoid silently reusing the last mate file when sample inputs are mismatched.
- `sc_methy_workflow.sh`: Forward and reverse barcode fastp QC now runs in parallel (lightweight task, safe for large datasets).
- `sc_methy_workflow.sh`: Multi-batch FASTQ merge now concatenates gzip members directly instead of decompressing and recompressing.
- `step3_split_bams.py`: Reworked BAM splitting from multi-process repeated full-file traversal to single-pass barcode-group traversal, reducing redundant I/O.
- `barcode_cs_multi.py`: Added trie-backed Hamming search for `distance > 1` barcode correction and fixed `QcStat.update()` nested counter merging.
- `step3_bam_to_allc.py`: Replaced shell-string subprocess calls with argument-list calls and guaranteed cleanup of intermediate sorted BAM files.
- `step4_wgs_summary.py`: Replaced repeated full-report regex scans with line-by-line Bismark report parsing.
- Extracted common parameter initialization code (~70 lines) from `rna_met`, `met_only`, and `forcecell` into shared module `nf/modules/params_init.nf`, reducing code duplication.
- Added `withName` resource configurations for all force_cell processes and `MERGE_BISMARK_BAM`.
- Reduced default `METHYLATION_BARCODE_EXTRACTION` memory from 120 GB to 48 GB with retry-based scaling.
- Aligned Docker profile container image with README documentation.
- `nf/modules/step2.nf`: Bismark `--parallel` now scales with `task.cpus` (clamped to 1–16) instead of a hard-coded 8, so the configured `cpus = 32` for `BISMARK_ALIGNMENT_FORWARD`/`REVERSE` is actually used.
- `nf/nextflow.config`: Dynamic-memory processes now double requested memory on retry only when the previous attempt died from an OOM-class exit code (134, 137–140); other failures (e.g. 143 SIGTERM from node preemption) retry at the same memory level. Consolidated into a shared `oomAwareMemory` closure used by all 16 dynamic-memory processes. Requires Nextflow >= 24.10.0 (`task.previousTrace`).

### Fixed
- **`barcode_cs_multi.py`: Multi-read chain direction was re-determined from the post-adapter-trim R1 sequence in the multi rescue second pass, which could flip `forward`/`reverse` relative to the pre-trim direction used in the first pass.** This sent affected multi reads to the wrong Bismark branch and produced direction-inconsistent BAM/ALLC results across runs. Fixed by carrying the pre-trim `chain_direction` through the multi temp FASTQ read name (`barcode_old_candidate_list_umi_chain_direction_original_read_name`) and reusing it in the second pass.
- `barcode_cs_multi.py`: `get_new_bc()` one-mismatch correction for barcodes containing `N` overwrote prior candidates with `mm_dict = {...}`; switched to `.update(...)` so the `N`-position substitutions are merged with the rest of the one-base mismatch set rather than replacing them.
- `barcode_cs_multi.py`: `ct_mean` / `cc_mean` divided by zero when `line_A` / `line_B` was zero (no qualifying reads). Added `>0` guards that emit `0.0` plus a `*_missing_evidence` flag in the summary JSON so the downstream report can mark the metric as missing instead of crashing.
- `barcode_cs_multi.py`: Reads with average quality below 30 in the CT/CC window were silently dropped from conversion-rate accounting; added a `low_quality_skipped` counter so the loss is now visible in the summary.
- `step3_split_bams.py` / `step3.nf`: `methy_only` workflow passed methylation barcodes to `--gexcb` but the splitter still applied the `gex_cb → m_cb` whitelist mapping, producing an empty cell list. Added `--gexcb-is-methylation` and switched `methy_only` to use it.
- `step3_estimated_cells.py`: `filtered_barcode` was written with a bare relative path and landed in the process work directory instead of `outdir`; pinned to `os.path.join(outdir, 'filtered_barcode')`.
- `step3_merge_sc_metrics.py`: Barcode extraction matched `_all.gz.count.csv` but actual files are `_allc.gz.count.csv`; the mismatch fell through to a fallback that left an `_allc.gz` suffix on the barcode, breaking downstream joins. Added the correct suffix as the primary match.
- `step3_mcds_manager.py`: MCDS shard merge silently accepted duplicate `cell` values within or across shards, producing a corrupted merged MCDS. Added pre-merge `_validate_unique_obs_across_shards()` that fails fast with an explicit error if duplicates are detected.
- `step4_report_rna_met.py`: Shard-level methylation cell metrics with duplicate `m_cb` rows were merged as-is, producing duplicated/cartesian rows in the RNA+methylation join. Added `deduplicate_cells_by_m_cb()` that aggregates additive counters with `sum`, recomputes `weighted_mc_rate` / `cell_saturation` / `*_mc_rate` from the summed counters, and validated `reads_counts` exists and is numeric.
- `step4_wgs_summary.py`: `parse_cell_info` set the dataframe index to `barcode` without deduping, producing a non-unique index when reads-count CSVs had duplicates; group by `barcode` and sum `reads_counts` before indexing.
- `step4_wgs_summary.py`: `Fraction Reads in Cells` divided by `summary["mapping"]["uniquereads"]` without a key/zero guard; falls back to `0.0` when missing or zero.
- `step4_wgs_summary.py`: Per-cell summary CSV header renamed `Reads_of_Max/Median_Cell` → `Read_Pairs_of_Max/Median_Cell` to match the actual unit (paired-end reads).
- `nf/modules/step1.nf`: `CREATE_FORWARD_PAIRS` / `CREATE_REVERSE_PAIRS` now emit FASTQ pairs as Nextflow `path` outputs (linked into a per-process output directory) rather than relying on `publishDir` side effects; the previous design risked downstream tasks racing against the publish step.
- `nf/modules/step1.nf`: `METHYLATION_BARCODE_EXTRACTION` now writes empty gzipped placeholders when one strand produces no reads, so downstream outer-join logic can explicitly handle the missing direction instead of hitting a glob with zero matches.
- `nf/modules/step3.nf`: `MERGE_BISMARK_BAM` guarded against empty forward/reverse arrays (`for f in "${arr[@]}"` under `set -u` errored when the array was empty), and emits a warning if both sides are empty.
- `nf/modules/step3.nf`: `ALLCOOLS_MERGE` cleared `merge_list_real.txt` (`: > merge_list_real.txt`) before appending, so resumed runs do not see stale entries from a previous attempt.
- `nf/modules/force_cell_custom.nf`: Pre-analysis `pre_outdir` and `addtag_dir` are now declared as `path` inputs (instead of `val`), so Nextflow actually stages the directories into each consumer's work dir. The previous `val` form caused downstream processes to read from non-existent paths.
- `nf/modules/force_cell_custom.nf`: `FORCE_CELL_APPLY_ALLOCOOLS_CHANGES` switched from `find -maxdepth 1` to a recursive find (per-cell ALLC files live in nested directories) and added a post-delete verification that fails the process if a dropped barcode's files remain.
- `nf/modules/force_cell_custom.nf`: Replaced hard-coded `oss-cn-beijing-internal.aliyuncs.com` with `${params.oss_endpoint}` for `ossutil` calls.
- `nf/modules/force_cell_custom.nf`: Replaced `oss_cp ... || true` silent fallbacks for `allcools.tar.gz` / `allcools_generate_datasets.tar.gz` with an explicit fallback to download the underlying directory (and a warning), so a tar miss no longer silently leaves the consumer with no data.
- `nf/modules/step4.nf`: `METHYLATION_SUMMARY` consolidated the separate forward/reverse Bismark report inputs into a single `bismark_reports` collection and symlinks `*_bismark_bt2_PE_report.txt` into `step2/bismark/`, so the summary step no longer breaks when one strand is missing.
- `nf/modules/utils.nf`: `ESTIMATED_CELLS` input tuple dropped the obsolete `pair_id` element so forward and reverse mapped-read count CSVs are both received for cell estimation.
- `docs/How_to_build_reference_genome.md`: Fixed example config path `nextflow.new.config` → `nextflow.config`.
- Fixed `methy_only` barcode-extracted FASTQ path reconstruction to include the sample-level publish directory.
- Fixed `methy_only` cell estimation input shape so forward and reverse mapped-read count CSVs are both staged for `ESTIMATED_CELLS`.
- Fixed `rna_met` downstream grouping to use the actual merged forward/reverse pair count from `MERGE_BISMARK_BAM` instead of `Math.max(forward_count, reverse_count)`.
- Fixed barcode-extracted FASTQ dataflow so `rna_met` and `methy_only` consume Nextflow `path` outputs instead of reconstructing files from asynchronous `publishDir` locations.
- Added empty gzipped forward/reverse placeholders when barcode extraction has no reads for one strand, allowing downstream outer-join logic to handle missing directions explicitly.
- Fixed OSS `force_cell` staging to write directly into the stable pre-analysis staging directory referenced by downstream processes, avoiding publishDir symlink races.

### Added
- `nf/modules/params_init.nf`: Shared parameter initialization module with `initCommonParams()` and `sanitizeSpaces()` functions.
- `nf/bin/stage_file.sh`: Shared file staging helper for local and OSS inputs.
- `nf/bin/lint_config_schema.py`: Lightweight config/schema/version consistency lint script.
- `local_test` profile for lightweight local syntax/config checks.
- `params.oss_endpoint` (default `oss-cn-beijing-internal.aliyuncs.com`) registered in `nextflow_schema.json` so OSS endpoints are configurable.

## [2.1.2] 2026-04-30
### Fixed
- Fixed version number inconsistency in subworkflow help messages (`rna_met` and `met_only` showed `2.1.1`/`v1.0.0` instead of `2.1.2`).
- Fixed missing space before "for" in `README.md` line 206.
- Fixed several configuration issues.
### Removed
- Removed Fastp parameters `--max-len1 60 --max-len2 0`.
### Added
- Added `.gitignore` entry for `CLAUDE.md`.
- Added trailing newline to `.gitignore`.

## [2.1.1] 2026-03-24
### Fixed
- Resolved the input error "* Chromosome chrM doesn't present in the .genome file. *" in sc_methy_workflow.sh.

## [2.1.0] 2026-03-13
### Fixed
- Fixed an issue where batch FASTQ results were lost when forward and reverse fastq counts were unequal.

## [2.0.0] 2026-02-26
### Added
- Added a guide for building the reference genome database directory (`--database_dir`).
- Added a `force_cell` workflow for updating results based on previous outputs.
- Added validation for `sample_id` and `outdir`; whitespace is now disallowed and normalized to underscores.
### Changed
- Optimized MCDS generation by splitting into smaller tasks and merging results.
- Updated repository layout documentation to match the current DSL2 structure.

## [1.0.7] 2026-01-23
### Fixed
- Fixed a bug in the `step3` module where the `allc_file_path.txt` file was not being generated correctly.
- Updated `sc_methy_workflow.sh` to support multiple datasets per sample.

## [1.0.6] 2025-12-30
### Changed
- Updated `nextflow.config`.
- Updated `methy_only.nf`.
- Updated `main.nf`.

## [1.0.5] 2025-12-22
### Changed
- Updated barcode whitelist.

## [1.0.4] 2025-12-18
### Changed
- Updated the C-T conversion rate algorithm.

## [1.0.3] 2025-12-09
### Fixed
- Fixed a metric reporting issue in `rna_methyl_report.html` and `rna_methy_summary.csv` where values were incorrectly displayed as read numbers instead of read pairs (persisted in v1.0.1 and v1.0.2).

## [1.0.2] 2025-12-02

### Added
- Added a `project` parameter to tag each process with project-specific information.

### Fixed
- Fixed a bug where cells containing only reverse reads or only forward reads were incorrectly filtered out.

## [1.0.1] 2025-11-21

### Changed
- Updated workflow version from 1.0.0 to 1.0.1 in `nextflow.config`.
- Updated the MET section in HTML reports to label and report values as "read pairs" instead of "reads".

### Added
- Added a summary CSV file (`rna_methy_summary.csv`) containing the same metrics as the HTML report, saved in the same output directory.

## [1.0.0] 2025-10-31

### Added
- Initial release of the SeekSoulMethyl pipeline for single-cell transcriptome and methylation dual-omics analysis.
