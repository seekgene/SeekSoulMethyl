# Changelog

All notable changes to this project will be documented in this file.

## [1.0.4]
### Changed
- Updated the C-T conversion rate algorithm.

## [1.0.3]

### Fixed
- Fixed a metric reporting issue in `rna_methyl_report.html` and `rna_methy_summary.csv` where values were incorrectly displayed as read numbers instead of read pairs (persisted in v1.0.1 and v1.0.2).

## [1.0.2]

### Added
- Added a `project` parameter to tag each process with project-specific information.

### Fixed
- Fixed a bug where cells containing only reverse reads or only forward reads were incorrectly filtered out.

## [1.0.1]

### Changed
- Updated workflow version from 1.0.0 to 1.0.1 in `nextflow.config`.
- Updated the MET section in HTML reports to label and report values as "read pairs" instead of "reads".

### Added
- Added a summary CSV file (`rna_methy_summary.csv`) containing the same metrics as the HTML report, saved in the same output directory.

## [1.0.0]

### Added
- Initial release of the SeekSoulMethyl pipeline for single-cell transcriptome and methylation dual-omics analysis.

