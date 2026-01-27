#!/bin/bash
# SeekOne DD Single Cell Multiome Methylation + RNA Analysis Pipeline
# Applicable for analyzing SeekOne DD Single Cell Multiome Methylation + RNA kit products from Beijing SeekGene BioSciences Co., Ltd.
# 
# Usage:
#   sc_methy_workflow.sh <exp_fq1> <exp_fq2> <methy_fq1> <methy_fq2> \
#       --sample <sample_name> \
#       --outdir <output_directory> \
#       --database_dir <database_directory> \
#       [--chemistry DD-MET5|DD-MET3] \
#       [--core <cpu_count>] \
#       [--filter_ch <filter_channel>] \
#
# Positional parameters (required):
#   exp_fq1         Single-cell transcriptome Read1 fastq file
#   exp_fq2         Single-cell transcriptome Read2 fastq file
#   methy_fq1       Single-cell methylation Read1 fastq file
#   methy_fq2       Single-cell methylation Read2 fastq file
#
# Named parameters:
#   --sample        Sample name (required)
#   --outdir        Result output path (required)
#   --database_dir  Reference genome database path (required)
#   --chemistry     Chemistry type: DD-MET5 or DD-MET3 (default: DD-MET5)
#   --core          Number of CPUs (default: 64)
#   --filter_ch     Filter channel (default: 2)
#

# Set error handling
set -e

# Logging functions (pure text output for clean log files)
log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warning() {
    echo "[WARNING] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Help function
show_help() {
    cat << EOF
SeekOne DD Single Cell Multiome Methylation + RNA Analysis Pipeline

Usage:
  $0 <exp_fq1> <exp_fq2> <methy_fq1> <methy_fq2> \\
      --sample <sample_name> \\
      --outdir <output_directory> \\
      --database_dir <database_directory> \\
      [--chemistry DD-MET5|DD-MET3] \\
      [--core <cpu_count>] \\
      [--filter_ch <filter_channel>] \\

Positional parameters (required):
  exp_fq1         Single-cell transcriptome Read1 fastq file (comma-separated for multiple files)
  exp_fq2         Single-cell transcriptome Read2 fastq file (comma-separated for multiple files)
  methy_fq1       Single-cell methylation Read1 fastq file (comma-separated for multiple files)
  methy_fq2       Single-cell methylation Read2 fastq file (comma-separated for multiple files)

Named parameters:
  --sample        Sample name (required)
  --outdir        Result output path (required)
  --database_dir  Reference genome database path (required)
  --chemistry     Chemistry type: DD-MET5 or DD-MET3 (default: DD-MET5)
  --core          Number of CPUs (default: 64)
  --filter_ch     Filter CH sites (default: 2)
  --help, -h      Show this help message

Example:
  $0 exp_R1.fq.gz exp_R2.fq.gz methy_R1.fq.gz methy_R2.fq.gz \\
      --sample sample1 \\
      --outdir /path/to/output \\
      --database_dir /path/to/database \\
      --chemistry DD-MET5 \\
      --core 64
EOF
    exit 0
}

# Parse positional parameters
if [ $# -lt 4 ]; then
    log_error "Insufficient parameters provided"
    show_help
fi

exp_fq1=$1
exp_fq2=$2
methy_fq1=$3
methy_fq2=$4
shift 4

# Set default values
sample=""
outdir=""
database_dir=""
chemistry="DD-MET5"
core=64
filter_ch=2

# Parse named parameters
while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)
            sample="$2"
            shift 2
            ;;
        --outdir)
            outdir="$2"
            shift 2
            ;;
        --database_dir)
            database_dir="$2"
            shift 2
            ;;
        --chemistry)
            chemistry="$2"
            shift 2
            ;;
        --core)
            core="$2"
            shift 2
            ;;
        --filter_ch)
            filter_ch="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            ;;
        *)
            log_error "Unknown parameter: $1"
            show_help
            ;;
    esac
done

# Validate required parameters
log_info "Validating input parameters..."

if [ -z "$sample" ]; then
    log_error "Parameter --sample is required"
    exit 1
fi

if [ -z "$outdir" ]; then
    log_error "Parameter --outdir is required"
    exit 1
fi

if [ -z "$database_dir" ]; then
    log_error "Parameter --database_dir is required"
    exit 1
fi

# Validate chemistry parameter
if [ "$chemistry" != "DD-MET5" ] && [ "$chemistry" != "DD-MET3" ]; then
    log_error "Invalid chemistry parameter: $chemistry. Must be DD-MET5 or DD-MET3"
    exit 1
fi

# Validate input files
for file_list in "$exp_fq1" "$exp_fq2" "$methy_fq1" "$methy_fq2"; do
    IFS=',' read -ra files <<< "$file_list"
    for file in "${files[@]}"; do
        if [ ! -f "$file" ]; then
            log_error "Input file does not exist: $file"
            exit 1
        fi
    done
done

# Validate database directory
if [ ! -d "$database_dir" ]; then
    log_error "Database directory does not exist: $database_dir"
    exit 1
fi

# Display parameters
log_info "==================== Parameters ===================="
log_info "Sample name:             $sample"
log_info "Transcriptome R1:        $exp_fq1"
log_info "Transcriptome R2:        $exp_fq2"
log_info "Methylation R1:          $methy_fq1"
log_info "Methylation R2:          $methy_fq2"
log_info "Output directory:        $outdir"
log_info "Database directory:      $database_dir"
log_info "Chemistry type:          $chemistry"
log_info "CPU cores:               $core"
log_info "Filter channel:          $filter_ch"
log_info "===================================================="



# Set the number of CPUs for bismark, recommended to be 1/4 of total CPUs
# Note: If bismark --parallel parameter is set to 8, ensure the machine has at least 64 CPUs
# For 500G data volume, --parallel set to 8, 64 CPU machine runs approximately 30 hours
bismark_core=$((core/4))

# Get script path
shell_path=`dirname "$0"`;
script_path=${shell_path}/nf/bin

# Set reference genome related file paths
genomeDir=${database_dir}/star          # STAR index directory
genomefa=$database_dir/fasta/genome.fa  # Reference genome fasta file
gtf=$database_dir/genes/genes.gtf       # Gene annotation file
genomebed=$database_dir/bed/chr_len.bed # Chromosome length file
chrom_size_path=$database_dir/bed/chr_nochrM.bed  # Chromosome length file without mitochondria
bismark_genome=$database_dir/fasta      # bismark genome index directory

# Set barcode whitelist files
U3CB_methylation=${script_path}/barcodes/U3CB_methylation.txt  # Methylation library barcode whitelist
# Set correspondence file between methylation and transcriptome library barcodes based on chemistry
if [ "${chemistry}" = "DD-MET3" ]; then
    cbcsv=${script_path}/barcodes/DD-M_bUCB3_whitelist.csv
    log_info "Using DD-MET3 barcode whitelist: ${cbcsv}"
elif [ "${chemistry}" = "DD-MET5" ]; then
    cbcsv=${script_path}/barcodes/ME5_bUCB3_whitelist.csv
    log_info "Using DD-MET5 barcode whitelist: ${cbcsv}"
else
    log_warning "Unknown chemistry '${chemistry}', using default DD-M whitelist"
    cbcsv=${script_path}/barcodes/DD-M_bUCB3_whitelist.csv
    log_info "Using default barcode whitelist: ${cbcsv}"
fi

# Set output paths
exp_outdir=${outdir}/${sample}_exp    # Transcriptome analysis result directory
methy_dir=${outdir}/${sample}_methy   # Methylation analysis result directory
# Cell barcode file obtained from transcriptome analysis, used for cell determination in methylation data
gexcb=$exp_outdir/${sample}/Analysis/step3/filtered_feature_bc_matrix/barcodes.tsv.gz
########################## Step 1: Data Quality Control (fastp) ##########################
# Use fastp software to perform quality control on raw data, filter out low-quality reads to obtain clean data
log_info "Step 1: Starting data quality control with fastp..."
mkdir -p ${outdir}/fastp
exp_fq1_clean=${outdir}/fastp/${sample}_exp_R1_clean.fq.gz
exp_fq2_clean=${outdir}/fastp/${sample}_exp_R2_clean.fq.gz
# Prepare input files (handle multiple comma-separated files)
IFS=',' read -ra exp_fq1_files <<< "$exp_fq1"
IFS=',' read -ra exp_fq2_files <<< "$exp_fq2"
if [ ${#exp_fq1_files[@]} -ne ${#exp_fq2_files[@]} ]; then
    log_error "Transcriptome R1/R2 file count mismatch: ${#exp_fq1_files[@]} vs ${#exp_fq2_files[@]}"
    exit 1
fi
if [ ${#exp_fq1_files[@]} -eq 1 ]; then
    exp_fq1_merge=${exp_fq1_files[0]}
    exp_fq2_merge=${exp_fq2_files[0]}
else
    exp_fq1_merge=${outdir}/fastp/${sample}_exp_R1_merged.fq.gz
    exp_fq2_merge=${outdir}/fastp/${sample}_exp_R2_merged.fq.gz
    if ! zcat "${exp_fq1_files[@]}" | gzip -c > "${exp_fq1_merge}"; then
        log_error "Transcriptome R1 merge failed: ${exp_fq1_merge}"
        exit 1
    fi
    if ! zcat "${exp_fq2_files[@]}" | gzip -c > "${exp_fq2_merge}"; then
        log_error "Transcriptome R2 merge failed: ${exp_fq2_merge}"
        exit 1
    fi
fi

fastp \
    -i ${exp_fq1_merge} \
    -I ${exp_fq2_merge} \
    -o ${exp_fq1_clean} \
    -O ${exp_fq2_clean} \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3  --unqualified_percent_limit 80 --n_base_limit 10  --length_required 60 --max_len1 60 --max_len2 0 \
    -j ${outdir}/fastp/${sample}_exp_fastp.json \
    -h ${outdir}/fastp/${sample}_exp_fastp.html \
    --thread ${core}
if [ ${#exp_fq1_files[@]} -gt 1 ]; then
    rm -f "${exp_fq1_merge}" "${exp_fq2_merge}"
fi
log_info "Transcriptome QC completed"

# Perform quality control on methylation data
methy_fq1_clean=${outdir}/fastp/${sample}_methy_R1_clean.fq.gz
methy_fq2_clean=${outdir}/fastp/${sample}_methy_R2_clean.fq.gz

# Prepare input files (handle multiple comma-separated files)
IFS=',' read -ra methy_fq1_files <<< "$methy_fq1"
IFS=',' read -ra methy_fq2_files <<< "$methy_fq2"
if [ ${#methy_fq1_files[@]} -ne ${#methy_fq2_files[@]} ]; then
    log_error "Methylation R1/R2 file count mismatch: ${#methy_fq1_files[@]} vs ${#methy_fq2_files[@]}"
    exit 1
fi
if [ ${#methy_fq1_files[@]} -eq 1 ]; then
    methy_fq1_merge=${methy_fq1_files[0]}
    methy_fq2_merge=${methy_fq2_files[0]}
else
    methy_fq1_merge=${outdir}/fastp/${sample}_methy_R1_merged.fq.gz
    methy_fq2_merge=${outdir}/fastp/${sample}_methy_R2_merged.fq.gz
    if ! zcat "${methy_fq1_files[@]}" | gzip -c > "${methy_fq1_merge}"; then
        log_error "Methylation R1 merge failed: ${methy_fq1_merge}"
        exit 1
    fi
    if ! zcat "${methy_fq2_files[@]}" | gzip -c > "${methy_fq2_merge}"; then
        log_error "Methylation R2 merge failed: ${methy_fq2_merge}"
        exit 1
    fi
fi

fastp \
    -i ${methy_fq1_merge} \
    -I ${methy_fq2_merge} \
    -o ${methy_fq1_clean} \
    -O ${methy_fq2_clean} \
    -j ${outdir}/fastp/${sample}_methy_fastp.json \
    -h ${outdir}/fastp/${sample}_methy_fastp.html \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3  --unqualified_percent_limit 80 --n_base_limit 10  --length_required 60 --disable_adapter_trimming \
    --thread ${core}
if [ ${#methy_fq1_files[@]} -gt 1 ]; then
    rm -f "${methy_fq1_merge}" "${methy_fq2_merge}"
fi
log_info "Methylation QC completed"

########################## Step 2: Transcriptome and Methylation Barcode Processing ##########################
########################## Step 2.1: Transcriptome Library Analysis (SeekSoulTools) ##########################
# Use SeekSoulTools to analyze transcriptome library
# Subsequent methylation library cells are based on cell barcodes determined by transcriptome library
log_info "Step 2.1: Starting transcriptome analysis with SeekSoulTools..."
mkdir -p ${exp_outdir}

# Set chemistry parameter for transcriptome analysis based on input chemistry
if [ "${chemistry}" = "DD-MET5" ]; then
    exp_chemistry="DD-MET5"
elif [ "${chemistry}" = "DD-MET3" ]; then
    exp_chemistry="DDV2"
else
    exp_chemistry="${chemistry}"
fi

log_info "Input chemistry parameter: ${chemistry}"
log_info "Chemistry used for transcriptome analysis: ${exp_chemistry}"

seeksoultools rna run \
	--fq1 ${exp_fq1_clean} \
	--fq2 ${exp_fq2_clean} \
	--samplename ${sample} \
	--genomeDir ${genomeDir} \
	--gtf $gtf \
	--chemistry ${exp_chemistry} \
	--core ${core} \
	--include-introns \
	--outdir ${exp_outdir}
log_info "Transcriptome analysis completed"

########################## Step 2.2: Methylation Barcode Recognition and Adapter Removal ##########################
# Extract and process barcode/UMI according to Read1 structure design and parameters
# Main functions include three parts:
# 1. Identify and output reads containing barcodes in the whitelist
# 2. Identify and remove adapters
# 3. Remove 9bp bases related to TN5 enzyme, including 5' and 3' ends
# Considering the impact of reads length on alignment, only output reads with read1 length >= 20bp and read2 length >= 60bp
log_info "Step 2.2: Starting methylation barcode recognition and adapter removal..."
mkdir -p ${methy_dir}
mkdir -p ${methy_dir}/step1
python ${script_path}/barcode_cs_multi.py \
	--fq1 ${methy_fq1_clean} \
	--fq2 ${methy_fq2_clean} \
	--samplename ${sample} \
	--outdir ${methy_dir} \
	--barcode ${U3CB_methylation} \
	--chemistry ${chemistry} \
	--split_fastq 0 \
	--filter_ch ${filter_ch} \
	--core ${core}
log_info "Barcode recognition completed"

# Perform quality control on fastq files after barcode recognition
log_info "Performing QC on barcode-recognized reads..."
fastp \
    -i ${methy_dir}/step1/${sample}_forward_1.fq.gz \
    -I ${methy_dir}/step1/${sample}_forward_2.fq.gz \
    -j ${methy_dir}/step1/${sample}_forward_fastp.json \
    -h ${methy_dir}/step1/${sample}_forward_fastp.html \
	--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3  --unqualified_percent_limit 80 --n_base_limit 10  --length_required 60 --disable_adapter_trimming \
    --thread ${core} 
fastp \
    -i ${methy_dir}/step1/${sample}_reverse_1.fq.gz \
    -I ${methy_dir}/step1/${sample}_reverse_2.fq.gz \
    -j ${methy_dir}/step1/${sample}_reverse_fastp.json \
    -h ${methy_dir}/step1/${sample}_reverse_fastp.html \
	--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3  --unqualified_percent_limit 80 --n_base_limit 10  --length_required 60 --disable_adapter_trimming \
    --thread ${core} 
log_info "Post-barcode QC completed"

########################## Step 3: Methylation Data Processing and Single-Cell Analysis ##########################
########################## Step 3.1: Bismark Alignment ##########################
log_info "Step 3.1: Starting Bismark alignment..."
mkdir -p ${methy_dir}/step2/bismark

# bismark alignment, actual CPU usage should be 4 times the set core
# Parameter description:
# --parallel ${bismark_core}: Number of parallel threads
# --genome ${bismark_genome}: Reference genome path
# -X 1000: Maximum insert fragment length
# --temp_dir: Temporary file directory
# forward: Directional library alignment
# reverse: PBAT library alignment with --pbat parameter
# forward
log_info "Step 3.1.1: Running Bismark alignment for forward reads..."
bismark \
	--parallel ${bismark_core} \
	--genome ${bismark_genome} \
	-1 ${methy_dir}/step1/${sample}_forward_1.fq.gz \
	-2 ${methy_dir}/step1/${sample}_forward_2.fq.gz \
	-o ${methy_dir}/step2/bismark/ \
	-X 1000 \
    --add_barcode \
    --add_umi \
	--temp_dir ${methy_dir}/step2/bismark/ > ${methy_dir}/step2/bismark/bismark_forward.log 2>&1
log_info "Forward alignment completed"

# reverse
log_info "Step 3.1.2: Running Bismark alignment for reverse reads..."
bismark \
	--parallel ${bismark_core} \
	--genome ${bismark_genome} \
	-1 ${methy_dir}/step1/${sample}_reverse_1.fq.gz \
	-2 ${methy_dir}/step1/${sample}_reverse_2.fq.gz \
	-o ${methy_dir}/step2/bismark/ \
	-X 1000 \
	--pbat \
    --add_barcode \
    --add_umi \
	--temp_dir ${methy_dir}/step2/bismark/ > ${methy_dir}/step2/bismark/bismark_reverse.log 2>&1
log_info "Bismark alignment completed"

########################## Step 3.2: Sort BAM files by read name ##########################
# Sort forward and reverse BAM files by read name for splitting
log_info "Step 3.2: Sorting BAM files by read name..."
log_info "Step 3.2.1: Sorting forward BAM file..."
samtools sort \
	-@ ${core} \
	-n \
	-o ${methy_dir}/step2/bismark/${sample}_forward_1_bismark_bt2_pe_sortbyname.bam \
	${methy_dir}/step2/bismark/${sample}_forward_1_bismark_bt2_pe.bam
log_info "Forward BAM sorted"

log_info "Step 3.2.2: Sorting reverse BAM file..."
samtools sort \
	-@ ${core} \
	-n \
	-o ${methy_dir}/step2/bismark/${sample}_reverse_1_bismark_bt2_pe_sortbyname.bam \
	${methy_dir}/step2/bismark/${sample}_reverse_1_bismark_bt2_pe.bam
log_info "BAM sorting completed"

########################## Step 3.3: Split BAM files by cell barcode ##########################
log_info "Step 3.3: Splitting BAM files by cell barcodes..."
mkdir -p ${methy_dir}/step3/split_bams
# Split forward BAM file by cell barcodes
log_info "Splitting forward BAM file..."
python ${script_path}/step3_split_bams.py \
	--bam ${methy_dir}/step2/bismark/${sample}_forward_1_bismark_bt2_pe_sortbyname.bam \
	--samplename ${sample}_forward \
	--outdir ${methy_dir}/step3/split_bams \
	--gexcb ${gexcb} \
	--cbcsv ${cbcsv} \
	--core ${core}
log_info "Forward BAM splitting completed"
# Split reverse BAM file by cell barcodes
log_info "Splitting reverse BAM file..."
python ${script_path}/step3_split_bams.py \
	--bam ${methy_dir}/step2/bismark/${sample}_reverse_1_bismark_bt2_pe_sortbyname.bam \
	--samplename ${sample}_reverse \
	--outdir ${methy_dir}/step3/split_bams \
	--gexcb ${gexcb} \
	--cbcsv ${cbcsv} \
	--core ${core}
log_info "BAM splitting completed"

########################## Step 3.4: Merge forward and reverse single-cell BAM files ##########################
log_info "Step 3.4: Merging forward and reverse single-cell BAM files..."
mkdir -p ${methy_dir}/step3/split_bams/merged
mkdir -p ${methy_dir}/step3/split_bams/merged/merged_fr_bam

# Find all forward BAM files
# step3_split_bams.py creates directories by removing '_bismark_.*' from BAM filename
# e.g., XYRD-WTJW880_forward_1_bismark_bt2_pe_sortbyname.bam -> XYRD-WTJW880_forward_1
forward_bams_dir="${methy_dir}/step3/split_bams/${sample}_forward_1"
reverse_bams_dir="${methy_dir}/step3/split_bams/${sample}_reverse_1"

# Check if directories exist
if [ ! -d "${forward_bams_dir}" ]; then
    log_error "Forward BAMs directory not found: ${forward_bams_dir}"
    log_info "Looking for available directories..."
    ls -la ${methy_dir}/step3/split_bams/ || true
    exit 1
fi

if [ ! -d "${reverse_bams_dir}" ]; then
    log_error "Reverse BAMs directory not found: ${reverse_bams_dir}"
    log_info "Looking for available directories..."
    ls -la ${methy_dir}/step3/split_bams/ || true
    exit 1
fi

# Merge forward and reverse BAM files for each barcode
log_info "Merging BAM files for each barcode..."
merged_count=0

declare -A fmap
declare -A rmap
declare -A all_barcodes

# Collect forward BAMs
for f in "${forward_bams_dir}"/*.bam; do
    if [[ -f "$f" ]]; then
        bc=$(basename "$f" .bam)
        fmap["$bc"]="$f"
        all_barcodes["$bc"]=1
    fi
done

# Collect reverse BAMs
for r in "${reverse_bams_dir}"/*.bam; do
    if [[ -f "$r" ]]; then
        bc=$(basename "$r" .bam)
        rmap["$bc"]="$r"
        all_barcodes["$bc"]=1
    fi
done

for barcode in "${!all_barcodes[@]}"; do
    f_bam="${fmap[$barcode]}"
    r_bam="${rmap[$barcode]}"
    output_bam="${methy_dir}/step3/split_bams/merged/merged_fr_bam/${barcode}.bam"
    
    if [[ -n "$f_bam" && -n "$r_bam" ]]; then
        samtools merge -n -@ 4 -o "$output_bam" "$f_bam" "$r_bam"
    elif [[ -n "$f_bam" ]]; then
        cp "$f_bam" "$output_bam"
    elif [[ -n "$r_bam" ]]; then
        cp "$r_bam" "$output_bam"
    fi
    
    merged_count=$((merged_count + 1))
    
    # Log progress every 100 barcodes
    if [ $((merged_count % 100)) -eq 0 ]; then
        log_info "Processed $merged_count barcodes..."
    fi
done
log_info "BAM merging/copying completed. Total processed: $merged_count barcodes"

# Merge and deduplicate filtered_barcode files
log_info "Merging filtered barcode lists..."
cat ${forward_bams_dir}/${sample}_forward_1_filtered_barcode \
    ${reverse_bams_dir}/${sample}_reverse_1_filtered_barcode | \
    sort | uniq > ${methy_dir}/step3/split_bams/merged/merged_filtered_barcode
log_info "Filtered barcode list created"

# Merge and aggregate reads_counts files
log_info "Merging and aggregating reads counts..."
awk -F ',' 'BEGIN {OFS=","} 
    FNR==1 && NR==1 {print; next}
    FNR==1 {next}
    NF >= 2 {
        barcode = $2
        reads = $1
        total[barcode] += reads
    } 
    END {
        for (barcode in total) {
            print total[barcode], barcode
        }
    }' ${forward_bams_dir}/${sample}_forward_1_filtered_barcode_reads_counts.csv \
       ${reverse_bams_dir}/${sample}_reverse_1_filtered_barcode_reads_counts.csv \
    > ${methy_dir}/step3/split_bams/merged/merged_filtered_barcode_reads_counts.csv
log_info "Reads counts aggregation completed"

########################## Step 3.5: Run ALLCools bam-to-allc for single cells ##########################
log_info "Step 3.5: Running ALLCools bam-to-allc for single-cell methylation extraction..."

python ${script_path}/step3_bam_to_allc.py \
	--indir ${methy_dir}/step3/split_bams/merged/merged_fr_bam \
	--samplename ${sample} \
	--outdir ${methy_dir}/step3 \
	--genomefa ${genomefa} \
	--chrom_size_path ${chrom_size_path} \
	--filtered_barcode ${methy_dir}/step3/split_bams/merged/merged_filtered_barcode \
	--core ${core} \
	--tag UR
log_info "ALLCools bam-to-allc completed"

# Move allcools output to allcools/ directory (to match Nextflow structure)
log_info "Organizing output directory structure..."
mkdir -p ${methy_dir}/step3/allcools
mv ${methy_dir}/step3/merged_fr_bam_allcools ${methy_dir}/step3/allcools/
log_info "Output organized to: ${methy_dir}/step3/allcools/merged_fr_bam_allcools/"

########################## Step 3.6: Merge single-cell metrics ##########################
log_info "Step 3.6: Merging single-cell methylation metrics..."

python ${script_path}/step3_merge_sc_metrics.py \
	${methy_dir}/step3/allcools \
	-o ${methy_dir}/step3/${sample}_cells \
	--cbcsv ${cbcsv}
log_info "Single-cell metrics merging completed"

# Copy filtered_barcode_reads_counts.csv to step3 root for step4_wgs_summary.py
log_info "Organizing files for summary report..."
cp ${methy_dir}/step3/split_bams/merged/merged_filtered_barcode_reads_counts.csv \
   ${methy_dir}/step3/filtered_barcode_reads_counts.csv
log_info "File organization completed"

########################## Step 3.7: Generate ALLCools datasets ##########################
log_info "Step 3.7: Generating ALLCools datasets..."
mkdir -p ${methy_dir}/step3/allcools_generate_datasets

# Create allc_file_path.txt
log_info "Creating allc file path table..."
ls ${methy_dir}/step3/allcools/merged_fr_bam_allcools/*_allc.gz 2>/dev/null | \
    awk '{
        n=split($0, parts, "/");
        filename=parts[n];
        gsub(/_allc.gz/, "", filename);
        print filename "\t" $0
    }' > ${methy_dir}/step3/allcools_generate_datasets/allc_file_path.txt

# Count allc files
allc_count=$(wc -l < ${methy_dir}/step3/allcools_generate_datasets/allc_file_path.txt)
log_info "Found $allc_count allc files"

# Validate that we have allc files
if [ "$allc_count" -eq 0 ]; then
    log_error "No allc files found! Check if step3_bam_to_allc.py ran successfully."
    log_error "Expected location: ${methy_dir}/step3/allcools/merged_fr_bam_allcools/"
    exit 1
fi

# Generate ALLCools datasets with multiple bin sizes
log_info "Running allcools generate-dataset..."
allcools generate-dataset \
	--allc_table ${methy_dir}/step3/allcools_generate_datasets/allc_file_path.txt \
	--output_path ${methy_dir}/step3/allcools_generate_datasets/${sample}.mcds \
	--chrom_size_path ${chrom_size_path} \
	--obs_dim cell \
	--cpu ${core} \
	--regions chrom1M 1000000 \
	--regions chrom500k 500000 \
	--regions chrom100k 100000 \
	--regions chrom50k 50000 \
	--regions chrom20k 20000 \
	--regions chrom10k 10000 \
	--quantifiers chrom1M count CGN \
	--quantifiers chrom1M hypo-score CGN cutoff=0.9 \
	--quantifiers chrom500k count CGN \
	--quantifiers chrom500k hypo-score CGN cutoff=0.9 \
	--quantifiers chrom100k count CGN \
	--quantifiers chrom100k hypo-score CGN cutoff=0.9 \
	--quantifiers chrom50k count CGN \
	--quantifiers chrom50k hypo-score CGN cutoff=0.9 \
	--quantifiers chrom20k count CGN \
	--quantifiers chrom20k hypo-score CGN cutoff=0.9 \
	--quantifiers chrom10k count CGN \
	--quantifiers chrom10k hypo-score CGN cutoff=0.9
log_info "ALLCools datasets generation completed"

########################## Step 3.8: Merge single-cell allc to bulk allc ##########################
log_info "Step 3.8: Merging single-cell allc files to bulk allc..."

# Create merge list
ls ${methy_dir}/step3/allcools/merged_fr_bam_allcools/*_allc.gz 2>/dev/null > ${methy_dir}/step3/merge_list.txt

# Validate merge list
merge_count=$(wc -l < ${methy_dir}/step3/merge_list.txt)
log_info "Found $merge_count allc files to merge"

if [ "$merge_count" -eq 0 ]; then
    log_error "No allc files found for merging!"
    log_error "Check directory: ${methy_dir}/step3/allcools/merged_fr_bam_allcools/"
    exit 1
fi

# Run allcools merge
allcools merge \
	--cpu ${core} \
	--allc_paths ${methy_dir}/step3/merge_list.txt \
	--output_path ${methy_dir}/step3/${sample}_merge_allc.gz \
	--chrom_size_path ${chrom_size_path}
log_info "ALLCools merge completed"

########################## Step 3.9: Extract CGN context from merged allc ##########################
log_info "Step 3.9: Extracting CGN context from merged allc..."

allcools extract-allc \
	--cpu 1 \
	--allc_path ${methy_dir}/step3/${sample}_merge_allc.gz \
	--output_prefix ${methy_dir}/step3/${sample}_ \
	--chrom_size_path ${chrom_size_path} \
	--mc_contexts CGN \
	--strandness merge
log_info "ALLCools extraction completed"

########################## Step 4: Quality Control Report and Visualization ##########################
########################## Step 4.1: Generate methylation summary report ##########################
log_info "Step 4.1: Generating methylation summary report..."
mkdir -p ${methy_dir}/step4

# Compute total CpG sites in genome (if not already computed)
if [ ! -f "${methy_dir}/genome_cpg_sites.json" ]; then
    log_info "Computing total CpG sites in genome..."
    python ${script_path}/count_cg_sites.py \
        ${genomefa} \
        ${chrom_size_path} \
        > ${methy_dir}/genome_cpg_sites.json
    log_info "CpG sites computation completed"
fi
# Prepare files for step4_wgs_summary.py
# The script expects files in step3_dir root, but we have them in subdirectories
# Create symbolic links to make them accessible
log_info "Preparing files for summary report..."

# Link CGN-Merge file from step3/ to methy_dir root
if [ ! -L "${methy_dir}/${sample}_.CGN-Merge.allc.tsv.gz" ]; then
    ln -sf ${methy_dir}/step3/${sample}_.CGN-Merge.allc.tsv.gz \
           ${methy_dir}/${sample}_.CGN-Merge.allc.tsv.gz
fi

# Link cells.csv from step3/ to methy_dir root
if [ ! -L "${methy_dir}/${sample}_cells.csv" ]; then
    ln -sf ${methy_dir}/step3/${sample}_cells.csv \
           ${methy_dir}/${sample}_cells.csv
fi

# Link filtered_barcode_reads_counts.csv from step3/ to methy_dir root
if [ ! -L "${methy_dir}/filtered_barcode_reads_counts.csv" ]; then
    ln -sf ${methy_dir}/step3/filtered_barcode_reads_counts.csv \
           ${methy_dir}/filtered_barcode_reads_counts.csv
fi

# Generate summary report
# Note: step4_wgs_summary.py only accepts 4 parameters
python ${script_path}/step4_wgs_summary.py \
	--outdir ${methy_dir} \
	--samplename ${sample} \
	--summary_json ${methy_dir}/${sample}_summary.json \
	--genome_info_json ${methy_dir}/genome_cpg_sites.json
log_info "Summary report generated"

########################## Step 4.2: Generate Multi Report ##########################
log_info "Step 4.2: Generating multi-report..."

python ${script_path}/step4_report_rna_met.py \
    --gexjson ${exp_outdir}/${sample}/Analysis/${sample}_summary.json \
    --metjson ${methy_dir}/${sample}_summary.json \
    --tsne_file ${exp_outdir}/${sample}/Analysis/step4/tsne_umi.xls \
    --filtered_counts_file ${methy_dir}/filtered_barcode_reads_counts.csv \
    --cells_file ${methy_dir}/${sample}_cells.csv \
    --raw_dir ${exp_outdir}/${sample}/Analysis/step3/raw_feature_bc_matrix \
    --counts_file ${exp_outdir}/${sample}/Analysis/step3/counts.xls \
    --detail_file ${exp_outdir}/${sample}/Analysis/step3/detail.xls \
    --gtf ${gtf} \
    --filtered_dir ${exp_outdir}/${sample}/Analysis/step3/filtered_feature_bc_matrix \
    --diff_data ${exp_outdir}/${sample}/Analysis/step4/FindAllMarkers.xls \
    --whitelist_file ${cbcsv} \
    --outdir ${outdir} \
    --samplename ${sample} \
    --rawname ${sample} \
    --nf_config ${script_path}/../nextflow.config 
log_info "Multi-report generated"

########################## Step 4.3: LSI/PCA dimensionality reduction and clustering ##########################
log_info "Step 4.3: Running dimensionality reduction and clustering analysis..."

python ${script_path}/step4_allcools_PCA_cluster.py \
	--mcds_path ${methy_dir}/step3/allcools_generate_datasets/${sample}.mcds \
	--samplename ${sample} \
	--var_dim chrom20k \
	--filtered_barcode_file ${methy_dir}/step3/split_bams/merged/merged_filtered_barcode \
	--outdir ${methy_dir}/step4 \
	--reduc lsi
log_info "Dimensionality reduction and clustering completed"

log_info "=========================================="
log_info "Pipeline completed successfully!"
log_info "Results are saved in: ${outdir}"
log_info "=========================================="

log_info "\nOutput files:"
log_info "  - Transcriptome results: ${exp_outdir}"
log_info "  - Methylation results: ${methy_dir}"
log_info "  - Single-cell BAMs: ${methy_dir}/step3/split_bams/merged/merged_fr_bam/"
log_info "  - ALLCools datasets: ${methy_dir}/step3/allcools_generate_datasets/${sample}.mcds"
log_info "  - Summary report: ${methy_dir}/${sample}_summary.json"
log_info "  - Clustering results: ${methy_dir}/step4/"
