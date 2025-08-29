#!/bin/bash
# SeekSoul Single-cell Transcriptome + Methylation Dual-omics Analysis Pipeline
# Applicable for analyzing single-cell transcriptome + methylation dual-omics reagent kit products from Beijing SeekSoul Biotechnology Co., Ltd.
# Input parameter description:
# $1: Sample name
# $2: Single-cell transcriptome Read1 fastq file
# $3: Single-cell transcriptome Read2 fastq file  
# $4: Single-cell methylation Read1 fastq file
# $5: Single-cell methylation Read2 fastq file
# $6: Result output path
# $7: Reference genome database path
# $8: Number of CPUs

# Set file descriptor limit

sample=${1}
exp_fq1=${2}
exp_fq2=${3}
methy_fq1=${4}
methy_fq2=${5}
outdir=${6}
database_dir=${7}
core=${8}



# Set the number of CPUs for bismark, recommended to be 1/8 of total CPUs
# Note: If bismark --parallel parameter is set to 8, ensure the machine has at least 64 CPUs
# For 500G data volume, --parallel set to 8, 64 CPU machine runs approximately 30 hours
bismark_core=$((core/8))

# Get script path
shell_path=`dirname "$0"`;
script_path=${shell_path}/src

# Set reference genome related file paths
genomeDir=${database_dir}/star          # STAR index directory
genomefa=$database_dir/fasta/genome.fa  # Reference genome fasta file
gtf=$database_dir/genes/genes.gtf       # Gene annotation file
genomebed=$database_dir/bed/chr_len.bed # Chromosome length file
chrom_size_path=$database_dir/bed/chr_nochrM.bed  # Chromosome length file without mitochondria
annofile=$database_dir/bed/chr_100kbins_anno.bed  # 100kb bins annotation file
bismark_genome=$database_dir/fasta      # bismark genome index directory

# Set barcode whitelist files
U3CB_methylation=${script_path}/barcodes/U3CB_methylation.txt  # Methylation library barcode whitelist
cbcsv=${script_path}/barcodes/bUCB3_whitelist.csv             # Correspondence between methylation and transcriptome library barcodes

# Set output paths
exp_outdir=${outdir}/${sample}_exp    # Transcriptome analysis result directory
methy_dir=${outdir}/${sample}_methy   # Methylation analysis result directory
# Cell barcode file obtained from transcriptome analysis, used for cell determination in methylation data
gexcb=$exp_outdir/${sample}/Analysis/step3/filtered_feature_bc_matrix/barcodes.tsv.gz

########################## Data Quality Control (fastp) ##########################
# Use fastp software to perform quality control on raw data, filter out low-quality reads to obtain clean data
mkdir -p ${outdir}/fastp
exp_fq1_clean=${outdir}/fastp/$(basename ${exp_fq1})
exp_fq2_clean=${outdir}/fastp/$(basename ${exp_fq2})
fastp \
    -i ${exp_fq1} \
    -I ${exp_fq2} \
    -o ${exp_fq1_clean} \
    -O ${exp_fq2_clean} \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3  --unqualified_percent_limit 80 --n_base_limit 10  --length_required 60 --max_len1 60 --max_len2 0 \
    -j ${outdir}/fastp/$(basename ${exp_fq1})_fastp.json \
    -h ${outdir}/fastp/$(basename ${exp_fq1})_fastp.html \
    --thread ${core}

# Perform quality control on methylation data
methy_fq1_clean=${outdir}/fastp/$(basename ${methy_fq1})
methy_fq2_clean=${outdir}/fastp/$(basename ${methy_fq2})
fastp \
    -i ${methy_fq1} \
    -I ${methy_fq2} \
    -o ${methy_fq1_clean} \
    -O ${methy_fq2_clean} \
    -j ${outdir}/fastp/$(basename ${methy_fq1})_fastp.json \
    -h ${outdir}/fastp/$(basename ${methy_fq1})_fastp.html \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3  --unqualified_percent_limit 80 --n_base_limit 10  --length_required 60 --disable_adapter_trimming \
    --thread ${core}

########################## Transcriptome Library Analysis (SeekSoulTools) ##########################
# Use SeekSoulTools to analyze transcriptome library
# Subsequent methylation library cells are based on cell barcodes determined by transcriptome library
mkdir -p ${exp_outdir}
seeksoultools rna run \
	--fq1 ${exp_fq1_clean} \
	--fq2 ${exp_fq2_clean} \
	--samplename ${sample} \
	--genomeDir ${genomeDir} \
	--gtf $gtf \
	--chemistry DDV2 \
	--core ${core} \
	--include-introns \
	--outdir ${exp_outdir}

########################## Methylation Library Analysis ##########################
########################## step1: Barcode Recognition and Adapter Removal ##########################
# Extract and process barcode/UMI according to Read1 structure design and parameters
# Main functions include three parts:
# 1. Identify and output reads containing barcodes in the whitelist
# 2. Identify and remove adapters
# 3. Remove 9bp bases related to TN5 enzyme, including 5' and 3' ends
# Considering the impact of reads length on alignment, only output reads with read1 length >= 20bp and read2 length >= 60bp
mkdir -p ${methy_dir}
mkdir -p ${methy_dir}/step1
python ${script_path}/barcode_cs_multi.py \
	--fq1 ${methy_fq1_clean} \
	--fq2 ${methy_fq2_clean} \
	--samplename ${sample} \
	--outdir ${methy_dir} \
	--barcode ${U3CB_methylation} \
	--chemistry DD-M \
	--core ${core}

# Perform quality control on fastq files after barcode recognition
fastp \
    -i ${methy_dir}/step1/${sample}_1.fq.gz \
    -I ${methy_dir}/step1/${sample}_2.fq.gz \
    -j ${methy_dir}/step1/${sample}_barcode_fastp.json \
    -h ${methy_dir}/step1/${sample}_barcode_fastp.html \
	--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3  --unqualified_percent_limit 80 --n_base_limit 10  --length_required 60 --disable_adapter_trimming \
    --thread ${core} 

########################## step2: Bismark Alignment and Deduplication ##########################
mkdir -p ${methy_dir}/step2/bismark

# Add barcode/umi information to read names for subsequent deduplication using deduplicate_bismark with barcode/umi information
python ${script_path}/add_umi_to_fastq_end_of_read_name.py \
    --fq1 ${methy_dir}/step1/${sample}_1.fq.gz \
    --fq2 ${methy_dir}/step1/${sample}_2.fq.gz \
    --samplename ${sample} \
    --outdir ${methy_dir}/step2/bismark

# bismark alignment, actual CPU usage should be 4 times the set core
# Parameter description:
# --parallel ${bismark_core}: Number of parallel threads
# --genome ${bismark_genome}: Reference genome path
# -X 1000: Maximum insert fragment length
# --temp_dir: Temporary file directory
# --non_directional: Non-strand-specific library
bismark \
	--parallel ${bismark_core} \
	--genome ${bismark_genome} \
	-1 ${methy_dir}/step2/bismark/${sample}_1_rename.fq.gz \
	-2 ${methy_dir}/step2/bismark/${sample}_2_rename.fq.gz \
	-o ${methy_dir}/step2/bismark/ \
	-X 1000 \
	--temp_dir ${methy_dir}/step2/bismark/ \
	--non_directional > ${methy_dir}/step2/bismark/bismark.log 2>&1

# deduplicate_bismark requires bam files to be sorted by read names, use samtools for sorting
samtools sort \
	-@ ${core} \
	-n \
	-o ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam \
	${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe.bam;

# deduplicate_bismark deduplication, use barcode and UMI information for PCR duplicate removal
deduplicate_bismark \
	-p \
	--output_dir ${methy_dir}/step2/bismark/ \
	--barcode \
	--bam ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam

# allcools requires bam files to be sorted by alignment position, use samtools for sorting
samtools sort \
    -@ ${core} \
    -o ${methy_dir}/step2/bismark/${sample}_1_deduplicated_sort.bam \
    ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.deduplicated.bam

# allcools counts methylation reads and coverage reads at CNN sites
# Output file format:
# Column 1: chromosome
# Column 2: position(1-based)
# Column 3: strand
# Column 4: sequence context
# Column 5: count of reads supporting methylation
# Column 6: read coverage
# Column 7: indicator of significant methylation (1 if no test is performed)
allcools bam-to-allc \
	--bam_path ${methy_dir}/step2/bismark/${sample}_1_deduplicated_sort.bam \
    --convert_bam_strandness \
    --reference_fasta ${genomefa} \
    --output_path ${methy_dir}/step2/bismark/${sample}_allc

########################## step2: Gene Annotation and Coverage Statistics ##########################
mkdir -p ${methy_dir}/step2/wgs

# Use featureCounts for gene annotation of bam files
# Parameter description:
# -T ${core}: Number of threads
# -t gene: Feature type is gene
# -s 0: Strand specificity is 0 (non-strand-specific)
# -M -O: Multi-mapping and overlapping reads processing
# -g gene_id: Use gene_id as gene identifier
# --fracOverlap 0.5: Overlap ratio threshold
# -p: Paired-end reads
# --countReadPairs: Count read pairs
featureCounts \
	-T ${core} -t gene -s 0 -M -O -g gene_id --fracOverlap 0.5 \
	-a ${gtf} \
	-p \
	--countReadPairs \
	-o ${methy_dir}/step2/wgs/wgs_gene_counts.txt \
	-R BAM \
	${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam

# Sort bam files by position
samtools sort \
	-@ ${core} \
	-o ${methy_dir}/step2/wgs/${sample}_sort.bam \
	${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam

# Use bedtools to calculate genome coverage
bedtools genomecov \
	-ibam ${methy_dir}/step2/wgs/${sample}_sort.bam > ${methy_dir}/step2/wgs/bedtools_genomecoverage.txt

# Calculate sample coverage in genome and write alignment and coverage information to json file
python ${script_path}/wgs_coverage_depth.py \
	--outdir ${methy_dir}/step2/wgs/ \
	--samplename ${sample} \
	--align_summary ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_PE_report.txt \
	--summary_json ${methy_dir}/${sample}_summary.json \
	--coveragefile ${methy_dir}/step2/wgs/bedtools_genomecoverage.txt

########################## step3: Cell Determination and Single-cell Analysis ##########################
# Sort featureCounts annotated bam files by read names
samtools sort \
        -n -O BAM -@ ${core} \
        -o ${methy_dir}/step2/wgs/${sample}_featureCounts_SortByName.bam \
        ${methy_dir}/step2/wgs/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam.featureCounts.bam

# Count the number of genes, umi, and reads for each barcode
# Output files:
# counts.xls: Column 1 is barcode, Column 2 is feature, Column 3 is umi count, Column 4 is reads count
python ${script_path}/step3dnam3.py \
	--bam ${methy_dir}/step2/wgs/${sample}_featureCounts_SortByName.bam \
	--outdir ${methy_dir}/step3

# Filter methylation data based on cell barcodes determined by transcriptome to obtain cell umicounts
# Output files:
# raw_umicounts.xls: All barcodes sorted by umi count
# filter_umicounts.xls: Barcodes determined as cells sorted by umi count
python ${script_path}/wgs_umi_count_cs.py  \
	--infile ${methy_dir}/step3/counts.xls \
	--rawcsv ${methy_dir}/step3/raw_umicounts.xls \
	--filtercsv ${methy_dir}/step3/filter_umicounts.xls \
	--gexcb ${gexcb} \
	--cbcsv ${cbcsv} \
	--outdir ${methy_dir} \
	--samplename ${sample}

# Output barcodes determined as cells to filtered_barcode file
cat ${methy_dir}/step3/filter_umicounts.xls |cut -f 1 |grep -v 'Barcode' > ${methy_dir}/step3/filtered_barcode

# Split deduplicated bam files into individual cell bam files according to filtered_barcode
python ${script_path}/split_bams.py \
    --bam ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.deduplicated.bam \
    --outdir ${methy_dir}/step3 \
	--samplename ${sample} \
	--filtered_barcode ${methy_dir}/step3/filtered_barcode

# For each cell's bam file, use samtools for sorting, use ALLCools to count methylation levels at each CNN site for each cell
# Output directories:
# allcools: Methylation level files for each cell
# allcools_generate_datasets: mcds files for downstream analysis
python ${script_path}/step3_run_allcools_and_generate_datasets.py \
	--indir ${methy_dir}/step3/split_bams \
	--samplename ${sample} \
	--outdir ${methy_dir}/step3/ \
	--genomefa ${genomefa} \
	--chrom_size_path ${chrom_size_path} \
	--filtered_barcode ${methy_dir}/step3/filtered_barcode

# Calculate coverage for each cell
# Output files:
# bed directory: bed files for each cell
# cb_cov directory: coverage files for each cell
# ${sample}_cb_genomecov.xls: Column 1 is barcode, Column 2 is coverage at 1X depth
python ${script_path}/wgs_sort_ml_anno_dev.py \
	--filtercb ${methy_dir}/step3/filtered_barcode \
	--core ${bismark_core} \
	--outdir ${methy_dir}/step3/ \
	--samplename ${sample} \
	--genomefa ${genomefa} \
	--allcpath ${methy_dir}/step3/allcools

# Count CpG number for each cell
# Output file: filter_gcov_umi_sort.xls
python ${script_path}/wgs_mk_martix_bins.py \
	--cbfile ${methy_dir}/step3/filtered_barcode \
	--outdir ${methy_dir}/step3 \
	--samplename ${sample} \
	--outcsv ${methy_dir}/step3/${sample}_CPGnum_cb.xls \
	--allcpgfile ${methy_dir}/step2/bismark/${sample}_allc.gz \
	--indir ${methy_dir}/step3/bed \
	--summary_json ${methy_dir}/${sample}_summary.json

# Generate final quality control table
python ${script_path}/wgs_summary_yf.py \
	--outdir ${methy_dir} \
	--samplename ${sample} \
	--summary_json ${methy_dir}/${sample}_summary.json

########################## step4: Single-cell Methylation Data Dimensionality Reduction and Clustering ##########################
# Use the results generated by ALLCools for single-cell methylation data dimensionality reduction and clustering
# Output directory step4 contains:
# tsne, umap dimensionality reduction clustering plots
# h5ad files
mkdir -p ${methy_dir}/step4
python ${script_path}/step4_allcools_PCA_cluster.py \
	--mcds_path ${methy_dir}/step3/allcools_generate_datasets/${sample}.mcds \
	--samplename ${sample} \
	--var_dim chrom1M \
	--filtered_barcode_file ${methy_dir}/step3/filtered_barcode \
	--outdir ${methy_dir}/step4
