process COUNTS_MAPPED_READS {
    tag "$sample-COUNTS_MAPPED_READS"
    //publishDir "${params.outdir}/${sample}_methy/step2/"
    
    input:
    tuple val(sample), val(pair_id), path(bismark_bam)

    output:
    tuple val(sample), val(pair_id), path("*_cb_aligned_reads_counts.csv"), emit: cb_aligned_reads_counts

    script:
    """
    step2_counts_aligned_reads.py -b ${bismark_bam} -o .
    """
}

process ESTIMATED_CELLS {
    tag "$sample-ESTIMATED_CELLS"
    publishDir "${params.outdir}/${sample}_methy/step3/"
    
    input:
    tuple val(sample),val(pair_id), path(cb_aligned_reads_counts)
    
    output:
    tuple val(sample), path("filtered_barcode"), emit: filtered_barcode
    
    script:
    """
    echo "estimated cells"
    step3_estimated_cells.py -r ./ --outdir .
    """
}