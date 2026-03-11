// split bams
process SPLIT_BAM_FILES {
    tag "$sample-$pair_id-SPLIT_BAM_FILES"
    resourceLabels label: "SPLIT_BAM_FILES_${params.project}_${sample}"
    
    input:
    tuple val(sample), val(pair_id), path(bismark_sortn_bam), path(gex_barcodes)
    
    output:
    tuple val(sample), val(pair_id), path("${bismark_sortn_bam.baseName.replaceAll(/_bismark_.*/, '')}/"), emit: split_bams_dir
    tuple val(sample), val(pair_id), path("${bismark_sortn_bam.baseName.replaceAll(/_bismark_.*/, '')}/${bismark_sortn_bam.baseName.replaceAll(/_bismark_.*/, '')}_filtered_barcode"), emit: filtered_barcode
    tuple val(sample), val(pair_id), path("${bismark_sortn_bam.baseName.replaceAll(/_bismark_.*/, '')}/${bismark_sortn_bam.baseName.replaceAll(/_bismark_.*/, '')}_filtered_barcode_reads_counts.csv"), emit: filtered_barcode_reads_counts
    
    script:
    """
    set -e
    # Split BAM files
    step3_split_bams.py \
        --bam ${bismark_sortn_bam} \
        --outdir . \
        --samplename ${bismark_sortn_bam.baseName} \
        --core ${task.cpus} \
        --gexcb ${gex_barcodes} \
        --cbcsv ${params.cbcsv}
    """
}
// merge single cell forward and reverse bam
process MERGE_BISMARK_BAM {
    tag "$sample-$pair_id-BISMARK_ALIGNMENT_MERGE"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/split_bams/merged/"
    resourceLabels label: "MERGE_BISMARK_BAM_${params.project}_${sample}"

    input:
    tuple val(sample), val(pair_id),
        path(forward_split_bams_dir),
        path(forward_filtered_barcodes),
        path(forward_filtered_barcode_reads_counts),
        path(reverse_split_bams_dir),
        path(reverse_filtered_barcodes),
        path(reverse_filtered_barcode_reads_counts)

    output:
    tuple val(sample),val(pair_id), path("${(forward_split_bams_dir.baseName == 'empty_split_bams_dir' ? reverse_split_bams_dir.baseName : forward_split_bams_dir.baseName)}_merge_filtered_barcode"), emit: merged_filtered_barcode
    tuple val(sample), val(pair_id), path("${(forward_split_bams_dir.baseName == 'empty_split_bams_dir' ? reverse_split_bams_dir.baseName : forward_split_bams_dir.baseName)}_merge_filtered_barcode_reads_counts.csv"), emit: merged_filtered_barcode_reads_counts
    tuple val(sample), val(pair_id), path("${(forward_split_bams_dir.baseName == 'empty_split_bams_dir' ? reverse_split_bams_dir.baseName : forward_split_bams_dir.baseName)}_merged_fr_bam"), emit: sc_merged_bam_dir
    
    script:
    def prefix = (forward_split_bams_dir.baseName == 'empty_split_bams_dir' ? reverse_split_bams_dir.baseName : forward_split_bams_dir.baseName)
    """
    set -e
    shopt -s nullglob

    mkdir -p ${prefix}_merged_fr_bam

    forward_bams=()
    if [[ -n "${forward_split_bams_dir}" && -d "${forward_split_bams_dir}" ]]; then
        mapfile -t forward_bams < <(find "${forward_split_bams_dir}/" -name "*.bam" | sort)
    fi

    reverse_bams=()
    if [[ -n "${reverse_split_bams_dir}" && -d "${reverse_split_bams_dir}" ]]; then
        mapfile -t reverse_bams < <(find "${reverse_split_bams_dir}/" -name "*.bam" | sort)
    fi

    declare -A fmap
    declare -A rmap
    for fb in "\${forward_bams[@]}"; do bn=\$(basename "\$fb" .bam); fmap["\$bn"]="\$fb"; done
    for rb in "\${reverse_bams[@]}"; do bn=\$(basename "\$rb" .bam); rmap["\$bn"]="\$rb"; done
    declare -A seen
    for bc in "\${!fmap[@]}"; do seen["\$bc"]=1; done
    for bc in "\${!rmap[@]}"; do seen["\$bc"]=1; done
    for bc in "\${!seen[@]}"; do
        f="\${fmap[\$bc]:-}"
        r="\${rmap[\$bc]:-}"
        output_bam="${prefix}_merged_fr_bam/\${bc}.bam"
        if [[ -n "\$f" && -n "\$r" ]]; then
            samtools merge -n -@ ${task.cpus} -o "\$output_bam" "\$f" "\$r"
        elif [[ -n "\$f" ]]; then
            cp "\$f" "\$output_bam"
        elif [[ -n "\$r" ]]; then
            cp "\$r" "\$output_bam"
        fi
    done
    
    echo "Merging and deduplicating barcode files..."
    barcode_files=( *_filtered_barcode )
    if (( \${#barcode_files[@]} > 0 )); then
        cat "\${barcode_files[@]}" | sort | uniq > ${prefix}_merge_filtered_barcode
    else
        : > ${prefix}_merge_filtered_barcode
    fi
    
    # Merge all filtered_barcode_reads_counts files and aggregate reads_counts by barcode
    echo "Merging and aggregating reads_counts files..."
    
    # Skip header line, merge all files, group by barcode and sum
    reads_count_files=( *_filtered_barcode_reads_counts.csv )
    if (( \${#reads_count_files[@]} > 0 )); then
        awk -F ',' '{
            if (NR==FNR && FNR==1) { print }
            if (NF >= 2 && FNR != 1 ) {
                barcode = \$2
                reads = \$1
                total[barcode] += reads
            }
        } END {
            for (barcode in total) {
                print total[barcode] "," barcode
            }
        }' "\${reads_count_files[@]}" > ${prefix}_merge_filtered_barcode_reads_counts.csv
    else
        echo "reads_counts,barcode" > ${prefix}_merge_filtered_barcode_reads_counts.csv
    fi
    
    echo "BAM file merging, reads_counts aggregation completed"
    """
}

// run allcools bam-to-allc
process ALLCOOLS_BAM_TO_ALLC {
    tag "$sample-$pair_id-ALLCOOLS_BAM_TO_ALLC"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/allcools/", enabled: !params.force_cell_number
    resourceLabels label: "ALLCOOLS_BAM_TO_ALLC_${params.project}_${sample}"
    
    input:
    tuple val(sample), val(pair_id), path(sc_merged_bam_dir), path(filtered_barcode)
    
    output:
    tuple val(sample), val(pair_id), path("${sc_merged_bam_dir.baseName}_allcools"), emit: allcools_allc_output
    
    script:
    def cores = Math.max(1, task.cpus - 2)
    """
    set -e     
    # Run allcools to generate datasets
    step3_bam_to_allc.py \
        --indir ${sc_merged_bam_dir} \
        --samplename ${sample} \
        --outdir . \
        --genomefa ${params.genomefa} \
        --chrom_size_path ${params.chrom_size_path} \
        --filtered_barcode ${filtered_barcode} \
        --core ${cores} \
        --tag UR
    
    """
}
// merge all single cell metrics
process MERGE_FILTERED_BARCODE_READS_COUNTS {
    tag "$sample-MERGE_FILTERED_BARCODE_READS_COUNTS"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/split_bams/merged/"
    resourceLabels label: "MERGE_FILTERED_BARCODE_READS_COUNTS_${params.project}_${sample}"

    input:
    tuple val(sample), path(merged_filtered_barcode), path(merged_filtered_barcode_reads_counts), path(allcools_allc_output)

    output:
    tuple val(sample), path("filtered_barcode"), path("filtered_barcode_reads_counts.csv"), emit: merged_filtered_barcode_reads_counts
    tuple val(sample), path("${sample}_cells.csv"), emit: allcools_cells_csv_output
    tuple val(sample), path("${sample}_cells.json"), emit: allcools_cells_json_output

    script:
    """
    set -e  
    cat *_merge_filtered_barcode > filtered_barcode
    awk '(NR==FNR){print}(NR!=FNR&&FNR!=1){print}' *_merge_filtered_barcode_reads_counts.csv > filtered_barcode_reads_counts.csv
    step3_merge_sc_metrics.py ./ -o ./${sample}_cells --cbcsv ${params.cbcsv}
    """
}

// run allcools generate-datasets
process ALLCOOLS_GENERATE_DATASETS {
    tag "$sample-ALLCOOLS_GENERATE_DATASETS"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/allcools_generate_datasets/"
    resourceLabels label: "ALLCOOLS_GENERATE_DATASETS_${params.project}_${sample}"
    
    input:
    tuple val(sample), path(allcools), path(filtered_barcode), path(gene_bed)
    
    output:
    tuple val(sample), path("${sample}.mcds"), emit: allcools_generate_datasets
    
    script:
    def cores = Math.max(1, task.cpus - 2)
    """
    set -e  
    # Optimized: Use find to avoid "Argument list too long" and pure bash for speed
    find -L . -name "*_allc.gz" | sort | while read id; do
        tmp=\${id%%_allc.gz}
        barcode=\${tmp##*/}
        echo "\${barcode}\t\${id}"
    done > allc_file_path.txt

    if [ ! -s allc_file_path.txt ]; then
        exit 1
    fi

    n_lines=\$(wc -l < allc_file_path.txt | awk '{print \$1}')
    
    # Generate regions and quantifiers parameters dynamically
    REGIONS_PARAMS=""
    QUANTIFIERS_PARAMS=""
    
    # Define bin sizes
    declare -A BINS_SIZES=(
        ["chrom1M"]=1000000
        ["chrom500k"]=500000
        ["chrom100k"]=100000
        ["chrom50k"]=50000
        ["chrom20k"]=20000
        ["chrom10k"]=10000
        ["geneslop2k"]=${gene_bed}
    )
    
    # Generate parameters for each bin size
    for region in "\${!BINS_SIZES[@]}"; do
        size=\${BINS_SIZES[\$region]}
        REGIONS_PARAMS="\$REGIONS_PARAMS --regions \$region \$size"
        QUANTIFIERS_PARAMS="\$QUANTIFIERS_PARAMS --quantifiers \$region count CGN --quantifiers \$region hypo-score CGN cutoff=0.9"
    done
    
    if [[ "\$n_lines" -eq 1 ]]; then
        real_bc=\$(cut -f1 allc_file_path.txt | head -n 1)
        echo "\$real_bc" > keep_barcodes.txt
        awk -F \$'\\t' 'NR==1{print \$1\"__DUMMY\\t\"\$2}' allc_file_path.txt > allc_file_path.dummy.txt
        cat allc_file_path.dummy.txt >> allc_file_path.txt

        tmp_out="./${sample}.mcds.__tmp"
        allcools generate-dataset \
            --allc_table allc_file_path.txt \
            --output_path "\$tmp_out" \
            --chrom_size_path ${params.chrom_size_path} \
            --obs_dim cell \
            --cpu ${cores} \
            \$REGIONS_PARAMS \
            \$QUANTIFIERS_PARAMS
        step3_mcds_manager.py subset --input_path "\$tmp_out" --barcode_path keep_barcodes.txt --output_path ./${sample}.mcds --mode stream
        rm -rf "\$tmp_out"
    else
        allcools generate-dataset \
            --allc_table allc_file_path.txt \
            --output_path ./${sample}.mcds \
            --chrom_size_path ${params.chrom_size_path} \
            --obs_dim cell \
            --cpu ${cores} \
            \$REGIONS_PARAMS \
            \$QUANTIFIERS_PARAMS
    fi
    """
}

process ALLCOOLS_GENERATE_DATASETS_PART {
    tag "$sample-$pair_id-ALLCOOLS_GENERATE_DATASETS_PART"
    resourceLabels label: "ALLCOOLS_GENERATE_DATASETS_PART_${params.project}_${sample}"
    
    input:
    tuple val(sample), val(pair_id), path(allcools_allc_output), path(gene_bed)
    
    output:
    tuple val(sample), val(pair_id), path("${sample}_${pair_id}.mcds"), optional: true, emit: allcools_generate_datasets_part
    
    script:
    def cores = Math.max(1, task.cpus - 2)
    """
    set -e  
    find -L . -name "*_allc.gz" | sort | while read id; do
        tmp=\${id%%_allc.gz}
        barcode=\${tmp##*/}
        echo "\${barcode}\t\${id}"
    done > allc_file_path.txt
    
    if [ ! -s allc_file_path.txt ]; then
        exit 0
    fi

    n_lines=\$(wc -l < allc_file_path.txt | awk '{print \$1}')
    
    REGIONS_PARAMS=""
    QUANTIFIERS_PARAMS=""
    
    declare -A BINS_SIZES=(
        ["chrom1M"]=1000000
        ["chrom500k"]=500000
        ["chrom100k"]=100000
        ["chrom50k"]=50000
        ["chrom20k"]=20000
        ["chrom10k"]=10000
        ["geneslop2k"]=${gene_bed}
    )
    
    for region in "\${!BINS_SIZES[@]}"; do
        size=\${BINS_SIZES[\$region]}
        REGIONS_PARAMS="\$REGIONS_PARAMS --regions \$region \$size"
        QUANTIFIERS_PARAMS="\$QUANTIFIERS_PARAMS --quantifiers \$region count CGN --quantifiers \$region hypo-score CGN cutoff=0.9"
    done
    
    if [[ "\$n_lines" -eq 1 ]]; then
        real_bc=\$(cut -f1 allc_file_path.txt | head -n 1)
        echo "\$real_bc" > keep_barcodes.txt
        awk -F \$'\\t' 'NR==1{print \$1\"__DUMMY\\t\"\$2}' allc_file_path.txt > allc_file_path.dummy.txt
        cat allc_file_path.dummy.txt >> allc_file_path.txt

        tmp_out="./${sample}_${pair_id}.mcds.__tmp"
        allcools generate-dataset \
            --allc_table allc_file_path.txt \
            --output_path "\$tmp_out" \
            --chrom_size_path ${params.chrom_size_path} \
            --obs_dim cell \
            --cpu ${cores} \
            \$REGIONS_PARAMS \
            \$QUANTIFIERS_PARAMS
        step3_mcds_manager.py subset --input_path "\$tmp_out" --barcode_path keep_barcodes.txt --output_path ./${sample}_${pair_id}.mcds --mode stream
        rm -rf "\$tmp_out"
    else
        allcools generate-dataset \
            --allc_table allc_file_path.txt \
            --output_path ./${sample}_${pair_id}.mcds \
            --chrom_size_path ${params.chrom_size_path} \
            --obs_dim cell \
            --cpu ${cores} \
            \$REGIONS_PARAMS \
            \$QUANTIFIERS_PARAMS
    fi
    """
}

process MERGE_MCDS {
    tag "$sample-MERGE_MCDS"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/allcools_generate_datasets/"
    resourceLabels label: "MERGE_MCDS_${params.project}_${sample}"
    
    input:
    tuple val(sample), path(mcds_parts)
    
    output:
    tuple val(sample), path("${sample}.mcds"), emit: allcools_generate_datasets
    
    script:
    """
    set -e
    step3_mcds_manager.py merge ${mcds_parts} --output_path ${sample}.mcds --merge_mode stream
    """
}

process ALLCOOLS_SUBMERGE {
    tag "$sample-$pair_id-ALLCOOLS_SUBMERGE"
    //publishDir "${params.outdir}/${sample}_methy/step3/"
    resourceLabels label: "ALLCOOLS_SUBMERGE_${params.project}_${sample}"

    input:
    tuple val(sample), val(pair_id), path(allcools_allc_output)

    output:
    tuple val(sample), val(pair_id), path("${sample}_${pair_id}/"), emit: allcools_submerge_allc

    script:
    def cores = Math.max(1, task.cpus - 2)
    """
    set -e     
    shopt -s nullglob
    find -L ./*/ -name "*_allc.gz" > merge_list.txt || true
    mkdir -p ${sample}_${pair_id}
    n_lines=`wc -l merge_list.txt | awk '{print \$1}'`
    if [[ \$n_lines -gt 1 ]]; then
        allcools merge \
        --cpu ${cores} \
        --allc_paths merge_list.txt \
        --output_path ${sample}_${pair_id}/${sample}_${pair_id}_merge_allc.gz \
        --chrom_size_path ${params.chrom_size_path}
    elif [[ \$n_lines -eq 1 ]]; then
        source_file=`head -n 1 merge_list.txt`
        cp \${source_file} ${sample}_${pair_id}/${sample}_${pair_id}_merge_allc.gz
        if [ -f "\${source_file}.tbi" ]; then
          cp \${source_file}.tbi ${sample}_${pair_id}/${sample}_${pair_id}_merge_allc.gz.tbi
        else
          touch ${sample}_${pair_id}/${sample}_${pair_id}_merge_allc.gz.tbi
        fi
    else
        touch ${sample}_${pair_id}/${sample}_${pair_id}_merge_allc.gz
        touch ${sample}_${pair_id}/${sample}_${pair_id}_merge_allc.gz.tbi
    fi
    """

}
// merge single cell allc to bulk allc
process ALLCOOLS_MERGE {
    tag "$sample-ALLCOOLS_MERGE"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/"
    resourceLabels label: "ALLCOOLS_MERGE_${params.project}_${sample}"

    input:
    // Accept a collection of per-pair merged allc files for a sample
    tuple val(sample), path(allcools_submerge_allc)

    output:
    tuple val(sample), path("${sample}_merge_allc.gz"), path("${sample}_merge_allc.gz.tbi"), emit: allcools_merge_allc

    script:
    def cores = Math.max(1, task.cpus - 2)
    """
    set -e     
    # Run allcools to merge datasets, about 12h
    set -e
    ls */*_allc.gz > merge_list.txt
    cat merge_list.txt | while read id; do if [[ -s \${id} ]]; then echo "\${id}" >> merge_list_real.txt;fi;done
    allcools merge \
    --cpu ${cores} \
    --allc_paths merge_list_real.txt \
    --output_path ${sample}_merge_allc.gz \
    --chrom_size_path ${params.chrom_size_path}
    """


}
// extract cg context allc
process ALLCOOLS_EXTRACT {
    tag "$sample-ALLCOOLS_EXTRACT"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/"
    resourceLabels label: "ALLCOOLS_EXTRACT_${params.project}_${sample}"

    input:
    tuple val(sample), path(allcools_merge_allc), path(allcools_merge_allc_tbi)

    output:
    tuple val(sample), path("*.CGN-Merge*"), emit: allcools_extract_allc_output

    script:
    """
    set -e      
    # Run allcools to extract datasets, about 16 min
    allcools extract-allc \
    --cpu 1 \
    --allc_path ${sample}_merge_allc.gz \
    --output_prefix ${sample}_ \
    --chrom_size_path ${params.chrom_size_path} \
    --mc_contexts CGN \
    --strandness merge
    """
}
