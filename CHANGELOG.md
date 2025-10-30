# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased] - 2025-10-28
### Changed
- SeekSoulTools supports analysis of RNA–methylation multi-omics data with ME5 chemistry
- Added a Nextflow-based workflow to analyze RNA–methylation multi-omics data
- Bismark (https://github.com/seekgene/Bismark) BAM output now includes CB (error-corrected barcode) and UR (raw UMI) tags
- Removed Bismark deduplication and featureCounts quantification steps
- Adopted ALLCools (https://github.com/seekgene/ALLCools) for UMI-based deduplication and methylation level estimation
- Methylation analysis now defaults to 20K bins with LSI for dimensionality reduction and clustering
- SeekSoulTools add TSO trimming for RNA Read2 data.
- Change chemistry DD-M to DD-MET3, ME5 to DD-MET5.

## [Unreleased] - 2024-08-25

### Changed
- Updated bowtie2 version to 2.5.4 in conda_dependencies.yml
- Updated fastp parameters, added `--unqualified_percent_limit 80 --n_base_limit 10  --length_required 60 --max_len1 60 --max_len2 0 ` for compatibility with GeneMind sequencer
- Fixed git clone path
- Reorganized directory structure, moved conda_dependencies.yml and script to the root directory
- Changed database file download permissions to public, removed ossutil signature
- Modified environment activation method to use `conda activate`
- Converted Chinese comments to English
- Removed ossutil download functionality from Alibaba Cloud
- Added mouse reference genome support
- Optimized samtools memory allocation to 1.5G for consistent performance across all system configurations
- Added automatic file descriptor limit setting (`ulimit -n 99999`) to handle high concurrency file operations
- Fixed conda dependencies installation by separating wheel file installation from environment creation (following SeekSoulTools recommendations)
- Renamed src directory to dependence for better organization
- Renamed script directory to src and moved sc_methy_workflow.sh to root directory for better accessibility
- Moved setup files (setup.py, setup.cfg, pyproject.toml) to dependence directory for better organization


