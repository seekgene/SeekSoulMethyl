# Changelog

All notable changes to this project will be documented in this file.

## [2.1.2] 2026-04-30
### Fixed
- Removed Fastp parameters `--max-len1 60 --max-len2 0`.
- Fixed several configuration issues.

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
