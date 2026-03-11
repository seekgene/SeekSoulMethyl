/*
 * Process definitions
 */

// compute cpg sites
process COMPUTE_CPG_SITES {
    resourceLabels label: "COMPUTE_CPG_SITES_${params.project}"
    tag "COMPUTE_CPG_SITES"
    publishDir "${params.outdir}/"

    output:
    path("genome_cg_info.json"), emit: cpg_sites
    
    script:
    """
    set -e
    # Calculate CpG sites
    count_cg_sites.py ${params.genomefa} ${params.chrom_size_path} > genome_cg_info.json
    """
}

// Expression data quality control (multi-group)
process FASTP_EXPRESSION_MULTI {
    tag "$sample-FASTP_EXPRESSION_MULTI-$group"
    publishDir "${params.outdir}/fastp"
    resourceLabels label: "FASTP_EXPRESSION_MULTI_${params.project}_$sample"
    
    
    input:
    tuple val(sample), val(group), val(exp_r1), val(exp_r2)
    
    output:
    tuple val(sample), val(group), path("${sample}_${group}_expression_clean_R1.fastq.gz"), path("${sample}_${group}_expression_clean_R2.fastq.gz"), emit: rna_fastp_multi_data
    path "*.{html,json}"
    
    script:
    def cores = Math.max(1, task.cpus - 2)
    """
    set -e
    stage_file() {
        local src="\$1"
        local name=\$(basename "\$src")
        if [[ "\$src" == oss://* ]]; then
            echo "Downloading \$src" >&2
            cfg="${projectDir}/bin/.ossutilconfig"
            if [ -f "\$cfg" ]; then
                ossutil cp --sign-version v4 --region cn-beijing -c "\$cfg" "\$src" "\$name" 1>&2
            else
                ossutil cp --sign-version v4 --region cn-beijing -e oss-cn-beijing-internal.aliyuncs.com -i "\$AccessKeyId" -k "\$AccessKeySecret" "\$src" "\$name" 1>&2
            fi
            if [ \$? -ne 0 ]; then
                echo "ERROR: ossutil failed for \$src" >&2
                exit 1
            fi
        else
            if [ -f "\$src" ]; then
                ln -s "\$src" "\$name"
            else
                cp "\$src" "\$name"
            fi
        fi
        echo "\$name"
    }
    R1_LOCAL=\$(stage_file "${exp_r1}")
    R2_LOCAL=\$(stage_file "${exp_r2}")

    if [ ! -f "\$R1_LOCAL" ]; then
        echo "ERROR: R1 file not found: \$R1_LOCAL" >&2
        exit 1
    fi
    if [ ! -f "\$R2_LOCAL" ]; then
        echo "ERROR: R2 file not found: \$R2_LOCAL" >&2
        exit 1
    fi

    fastp \
        -i \$R1_LOCAL \
        -I \$R2_LOCAL \
        -o ${sample}_${group}_expression_clean_R1.fastq.gz \
        -O ${sample}_${group}_expression_clean_R2.fastq.gz \
        -w ${cores} \
        -h ${sample}_${group}_expression_fastp.html \
        -j ${sample}_${group}_expression_fastp.json \
        --cut_tail --cut_tail_window_size 1 \
        --cut_tail_mean_quality 3  --unqualified_percent_limit 80 \
        --n_base_limit 10  --length_required 60 --max_len1 60 --max_len2 0
    """
}

// Methylation data quality control (multi-group)
process FASTP_METHYLATION_MULTI {
    tag "$sample-FASTP_METHYLATION_MULTI-$group"
    publishDir "${params.outdir}/fastp"
    resourceLabels label: "FASTP_METHYLATION_MULTI_${params.project}_$sample"
    
    input:
    tuple val(sample), val(group), val(methy_r1), val(methy_r2)
    
    output:
    tuple val(sample), val(group), path("${sample}_${group}_methylation_clean_R1.fastq.gz"), path("${sample}_${group}_methylation_clean_R2.fastq.gz"), emit: methy_fastp_multi_data
    path "*.{html,json}"
    
    script:
    def cores = Math.max(1, task.cpus - 2)
    """
    set -e
    stage_file() {
        local src="\$1"
        local name=\$(basename "\$src")
        if [[ "\$src" == oss://* ]]; then
            echo "Downloading \$src" >&2
            cfg="${projectDir}/bin/.ossutilconfig"
            if [ -f "\$cfg" ]; then
                ossutil cp --sign-version v4 --region cn-beijing -c "\$cfg" "\$src" "\$name" 1>&2
            else
                ossutil cp --sign-version v4 --region cn-beijing -e oss-cn-beijing-internal.aliyuncs.com -i "\$AccessKeyId" -k "\$AccessKeySecret" "\$src" "\$name" 1>&2
            fi
            if [ \$? -ne 0 ]; then
                echo "ERROR: ossutil failed for \$src" >&2
                exit 1
            fi
        else
            if [ -f "\$src" ]; then
                ln -s "\$src" "\$name"
            else
                cp "\$src" "\$name"
            fi
        fi
        echo "\$name"
    }
    R1_LOCAL=\$(stage_file "${methy_r1}")
    R2_LOCAL=\$(stage_file "${methy_r2}")

    if [ ! -f "\$R1_LOCAL" ]; then
        echo "ERROR: R1 file not found: \$R1_LOCAL" >&2
        exit 1
    fi
    if [ ! -f "\$R2_LOCAL" ]; then
        echo "ERROR: R2 file not found: \$R2_LOCAL" >&2
        exit 1
    fi

    fastp \
        -i \$R1_LOCAL \
        -I \$R2_LOCAL \
        -o ${sample}_${group}_methylation_clean_R1.fastq.gz \
        -O ${sample}_${group}_methylation_clean_R2.fastq.gz \
        -w ${cores} \
        -h ${sample}_${group}_methylation_fastp.html \
        -j ${sample}_${group}_methylation_fastp.json \
        --disable_adapter_trimming \
        --cut_tail --cut_tail_window_size 1 \
        --cut_tail_mean_quality 3  --unqualified_percent_limit 80 \
        --n_base_limit 10  --length_required 60
    """
}


// RNA expression analysis
process SEEKSOULTOOLS_RNA {
    tag "$sample-SEEKSOULTOOLS_RNA"
    publishDir "${params.outdir}/${sample}/${sample}_exp/"
    resourceLabels label: "SEEKSOULTOOLS_RNA_${params.project}_$sample"
    
    input:
    tuple val(sample), path(exp_r1_list), path(exp_r2_list)
    
    output:
    tuple val(sample), path("${sample}/Analysis/step3/filtered_feature_bc_matrix/barcodes.tsv.gz"), emit: gex_barcodes
    tuple val(sample), path("${sample}/Analysis/${sample}_gex_summary.json"), emit: gex_summary_json
    tuple val(sample), path("${sample}/Analysis/step3/filtered_feature_bc_matrix"), emit: filtered_dir
    tuple val(sample), path("${sample}/Analysis/step3/raw_feature_bc_matrix"), emit: raw_dir
    tuple val(sample), path("${sample}/Analysis/step3/counts.xls"), emit: counts_xls
    tuple val(sample), path("${sample}/Analysis/step3/detail.xls"), emit: detail_xls
    tuple val(sample), path("${sample}/Analysis/step4/tsne_umi.xls"), emit: tsne_umi
    tuple val(sample), path("${sample}/Analysis/step4/FindAllMarkers.xls"), emit:diff_data
    path("${sample}/")
    
    script:
    def cores = Math.max(1, task.cpus - 2)
    // Normalize inputs to lists to handle single-file vs multi-file cases
    def r1s = (exp_r1_list instanceof java.util.List) ? exp_r1_list : [exp_r1_list]
    def r2s = (exp_r2_list instanceof java.util.List) ? exp_r2_list : [exp_r2_list]
    if (r1s.size() != r2s.size()) {
        throw new IllegalArgumentException("Mismatched expression R1/R2 inputs for ${sample}: R1=${r1s.size()} R2=${r2s.size()}")
    }
    def fq_args = (0..<r1s.size()).collect { i ->
        def r1n = r1s[i].toString().tokenize('/')[-1]
        def r2n = r2s[i].toString().tokenize('/')[-1]
        "--fq1 \"${r1n}\" --fq2 \"${r2n}\""
    }.join(' ')
    """
    set -e
    seeksoultools rna run \
        --samplename ${sample} \
        ${fq_args} \
        --outdir . \
        --core ${cores} \
        --chemistry ${params.exp_chemistry} \
        --include-introns \
        --gtf ${params.gtf} \
        --genomeDir ${params.genomeDir}
    rm -rf "${sample}/Analysis/.test"
    rm -rf "${sample}/Analysis/step2/STAR/"*__STARtmp
    mv "${sample}/Analysis/${sample}_summary.json" "${sample}/Analysis/${sample}_gex_summary.json"

    """
}

// Methylation barcode extraction
process METHYLATION_BARCODE_EXTRACTION {
    tag "$sample-METHYLATION_BARCODE_EXTRACTION"
    publishDir "${params.outdir}/${sample}/${sample}_methy/"
    resourceLabels label: "METHYLATION_BARCODE_EXTRACTION_${params.project}_$sample"
    
    input:
    tuple val(sample), path(methy_r1_list), path(methy_r2_list)
    
    output:
    tuple val(sample), path("step1/${sample}_forward_*1.fq.gz"), path("step1/${sample}_forward_*2.fq.gz"), path("step1/${sample}_reverse_*1.fq.gz"), path("step1/${sample}_reverse_*2.fq.gz"), path("${sample}_methy_summary.json"), emit: methy_barcode_output
    tuple val(sample), path("step1/${sample}_barcode_stats.txt"), emit: all_cb_step1_reads_counts, optional: true
    
    script:
    def cores = Math.max(1, task.cpus - 2)
    // Normalize inputs to lists to handle single-file vs multi-file cases
    def mr1s = (methy_r1_list instanceof java.util.List) ? methy_r1_list : [methy_r1_list]
    def mr2s = (methy_r2_list instanceof java.util.List) ? methy_r2_list : [methy_r2_list]
    if (mr1s == null || mr1s.size()==0) throw new IllegalArgumentException("No methylation R1 inputs for ${sample}")
    if (mr2s == null || mr2s.size()==0) throw new IllegalArgumentException("No methylation R2 inputs for ${sample}")
    if (mr1s.size() != mr2s.size()) {
        throw new IllegalArgumentException("Mismatched methylation R1/R2 inputs for ${sample}: R1=${mr1s.size()} R2=${mr2s.size()}")
    }
    def fq_args_m = (0..<mr1s.size()).collect { i ->
        def r1n = mr1s[i].toString().tokenize('/')[-1]
        def r2n = mr2s[i].toString().tokenize('/')[-1]
        "--fq1 \"${r1n}\" --fq2 \"${r2n}\""
    }.join(' ')
    """
    set -e
    barcode_cs_multi.py \
        ${fq_args_m} \
        --barcode ${params.methy_barcode_wl} \
        --outdir . \
        --samplename ${sample} \
        --core ${cores} \
        --chemistry ${params.chemistry} \
        --filter_ch ${params.filter_ch} \
        --split_fastq ${params.split_fastq}
    mv "${sample}_summary.json" "${sample}_methy_summary.json"
    if [[ ${params.split_fastq} -gt 0 ]]; then
        rm step1/${sample}_forward_1.fq.gz 
        rm step1/${sample}_forward_2.fq.gz 
        rm step1/${sample}_reverse_1.fq.gz 
        rm step1/${sample}_reverse_2.fq.gz
    fi
    """
}    

// Parse and group fastq files - pair based on identifiers
process PARSE_FASTQ_FILES {
    tag "$sample-PARSE_FASTQ_FILES"
    publishDir "${params.outdir}/${sample}/${sample}_methy/"
    resourceLabels label: "PARSE_FASTQ_FILES_${params.project}_$sample"
    
    input:
    tuple val(sample), path(forward_r1_files), path(forward_r2_files), path(reverse_r1_files), path(reverse_r2_files), path(summary_json)
    
    output:
    tuple val(sample), path("forward_pairs.txt"), emit: forward_pairs
    tuple val(sample), path("reverse_pairs.txt"), emit: reverse_pairs
    
    script:
    """
    set -e
    # Create forward file pair list - pair based on identifiers
    echo "Creating forward pairs based on identifiers..."
    
    # Extract identifiers from all forward R1 files
    for r1_file in ${sample}_forward_*1.fq.gz; do
        if [ -f "\$r1_file" ]; then
            # Extract identifier (e.g., extract AA from sample_forward_AA_1.fq.gz)
            identifier=\$(echo "\$r1_file" | sed 's/.*_forward_\\(.*\\)1\\.fq\\.gz/\\1/')
            r2_file="${sample}_forward_\${identifier}2.fq.gz"
            
            # Check if corresponding R2 file exists
            if [ -f "\$r2_file" ]; then
                echo "\$r1_file,\$r2_file" >> forward_pairs.txt
                echo "Paired: \$r1_file <-> \$r2_file"
            else
                echo "Warning: No matching R2 file found for \$r1_file"
            fi
        fi
    done
    
    # Create reverse file pair list - pair based on identifiers
    echo "Creating reverse pairs based on identifiers..."
    
    # Extract identifiers from all reverse R1 files
    for r1_file in ${sample}_reverse_*1.fq.gz; do
        if [ -f "\$r1_file" ]; then
            # Extract identifier (e.g., extract AA from sample_reverse_AA_1.fq.gz)
            identifier=\$(echo "\$r1_file" | sed 's/.*_reverse_\\(.*\\)1\\.fq\\.gz/\\1/')
            r2_file="${sample}_reverse_\${identifier}2.fq.gz"
            
            # Check if corresponding R2 file exists
            if [ -f "\$r2_file" ]; then
                echo "\$r1_file,\$r2_file" >> reverse_pairs.txt
                echo "Paired: \$r1_file <-> \$r2_file"
            else
                echo "Warning: No matching R2 file found for \$r1_file"
            fi
        fi
    done
    
    # Ensure output files exist (even if empty)
    touch forward_pairs.txt reverse_pairs.txt
    
    echo "Forward pairs:"
    cat forward_pairs.txt || echo "No forward pairs found"
    echo "Reverse pairs:"
    cat reverse_pairs.txt || echo "No reverse pairs found"
    """
}

// Create forward file pair distribution process
process CREATE_FORWARD_PAIRS {
    tag "$sample-CREATE_FORWARD_PAIRS"
    resourceLabels label: "CREATE_FORWARD_PAIRS_${params.project}_$sample"
    
    input:
    tuple val(sample), path(pairs_file), path(r1_fastq), path(r2_fastq)
    
    output:
    tuple val(sample), stdout
    
    script:
    """
    set -e
    # Read pairs file and output each pair as tab-separated values
    while IFS=',' read -r r1_file r2_file; do
        if [ -n "\$r1_file" ] && [ -n "\$r2_file" ]; then
            # Extract pair_id from filename (e.g., AA from HC2_4_forward_AA_1.fq.gz)
            pair_id=\$(echo "\$r1_file" | sed 's/.*_forward\\(.*\\)1\\.fq\\.gz/\\1/')
            echo "${sample},\${pair_id},\${r1_file},\${r2_file}"
        fi
    done < ${pairs_file}
    """
}

// Create reverse file pair distribution process
process CREATE_REVERSE_PAIRS {
    tag "$sample-CREATE_REVERSE_PAIRS"
    resourceLabels label: "CREATE_REVERSE_PAIRS_${params.project}_$sample"
    
    input:
    tuple val(sample), path(pairs_file), path(r1_fastq), path(r2_fastq)
    
    output:
    tuple val(sample), stdout
    
    script:
    """
    set -e
    # Read pairs file and output each pair as tab-separated values
    while IFS=',' read -r r1_file r2_file; do
        if [ -n "\$r1_file" ] && [ -n "\$r2_file" ]; then
            # Extract pair_id from filename (e.g., AA from HC2_4_reverse_AA_1.fq.gz)
            pair_id=\$(echo "\$r1_file" | sed 's/.*_reverse\\(.*\\)1\\.fq\\.gz/\\1/')
            echo -e "${sample},\${pair_id},\${r1_file},\${r2_file}"
        fi
    done < ${pairs_file}
    """
}

// Fastp quality control after methylation barcode extraction
process FASTP_METHYLATION_BARCODE_EXTRACT {
    tag "$sample-$pair_id-FASTP_METHYLATION_BARCODE_EXTRACT"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step1/"
    resourceLabels label: "FASTP_METHYLATION_BARCODE_EXTRACT_${params.project}_$sample"
   
    input:
    tuple val(sample), val(pair_id),path(methy_barcode_r1), path(methy_barcode_r2)
    
    output:
    path "*.{html,json}"
    
    script:
    """
    set -e
    fastp \
        -i ${methy_barcode_r1} \
        -I ${methy_barcode_r2} \
        -w ${task.cpus} \
        -h ${methy_barcode_r1.baseName}_fastp.html \
        -j ${methy_barcode_r1.baseName}_fastp.json \
        --disable_adapter_trimming \
        --cut_tail --cut_tail_window_size 1 \
        --cut_tail_mean_quality 3  --unqualified_percent_limit 80 \
        --n_base_limit 10  --length_required 60
    """

}
