
// Bismark forward strand alignment
process BISMARK_ALIGNMENT_FORWARD {
    tag "$sample-BISMARK_ALIGNMENT_FORWARD"
    publishDir "${params.outdir}/${sample}_methy/step2/bismark"
    
    input:
    tuple val(sample), val(pair_id), path(methy_forward_r1), path(methy_forward_r2)
    
    output:
    tuple val(sample), val(pair_id), path("${methy_forward_r1.baseName.replaceAll(/\.fq$/, '')}_bismark_bt2_pe.bam"), emit: bismark_forward_bam
    tuple val(sample), val(pair_id), path("${methy_forward_r1.baseName.replaceAll(/\.fq$/, '')}_bismark_bt2_PE_report.txt"), emit: bismark_forward_report
    
    script:
    """
    set -e   
    # Bismark alignment
    bismark \
        --genome ${params.bismark_ref} \
        --parallel 8 \
        -1 ${methy_forward_r1} \
        -2 ${methy_forward_r2} \
        -o . \
        -X 1000 \
        --add_barcode \
        --add_umi \
        --temp_dir . > bismark.log 2>&1
    """
}
// Bismark reverse strand alignment
process BISMARK_ALIGNMENT_REVERSE {
    tag "$sample-BISMARK_ALIGNMENT_REVERSE"
    publishDir "${params.outdir}/${sample}_methy/step2/bismark"
    
    input:
    tuple val(sample), val(pair_id), path(methy_reverse_r1), path(methy_reverse_r2)
    
    output:
    tuple val(sample), val(pair_id), path("${methy_reverse_r1.baseName.replaceAll(/\.fq$/, '')}_bismark_bt2_pe.bam"), emit: bismark_reverse_bam
    tuple val(sample), val(pair_id), path("${methy_reverse_r1.baseName.replaceAll(/\.fq$/, '')}_bismark_bt2_PE_report.txt"), emit: bismark_reverse_report
    
    script:
    """
    set -e    
    # Bismark alignment
    bismark \
        --genome ${params.bismark_ref} \
        --parallel 8 \
        -1 ${methy_reverse_r1} \
        -2 ${methy_reverse_r2} \
        -o . \
        -X 1000 \
        --pbat \
        --add_barcode \
        --add_umi \
        --temp_dir . > bismark.log 2>&1
    """
}

// BAM file sorting (by read name)
process SORT_BAM_BY_NAME {
    tag "$sample-SORT_BAM_BY_NAME"
    //publishDir "${params.outdir}/${sample}_methy/step2/bismark/"
    
    input:
    tuple val(sample), val(pair_id), path(bismark_bam)
    
    output:
    tuple val(sample), val(pair_id), path("${bismark_bam.baseName}_sortbyname.bam"), emit: bismark_sortn_bam
    
    script:
    def cores = Math.max(1, task.cpus - 2)
    """
    set -e
    # Sort by read name
    samtools sort \
        -n -@ ${cores} \
        -o ${bismark_bam.baseName}_sortbyname.bam \
        ${bismark_bam}
    """
}
