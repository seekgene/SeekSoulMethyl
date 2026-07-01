process COUNTS_MAPPED_READS {
    tag "$sample-COUNTS_MAPPED_READS"
    //publishDir "${params.outdir}/${sample}_methy/step2/"
    resourceLabels label: "COUNTS_MAPPED_READS_${params.project}_${sample}"
    
    input:
    tuple val(sample), val(pair_id), path(bismark_bam)

    output:
    tuple val(sample), val(pair_id), path("*_cb_aligned_reads_counts.csv"), emit: cb_aligned_reads_counts

    script:
    """
    set -e
    step2_counts_aligned_reads.py -b ${bismark_bam} -o .
    """
}

process ESTIMATED_CELLS {
    tag "$sample-ESTIMATED_CELLS"
    publishDir "${params.outdir}/${sample}_methy/step3/"
    resourceLabels label: "ESTIMATED_CELLS_${params.project}_${sample}"

    input:
    tuple val(sample), path(cb_aligned_reads_counts)

    output:
    tuple val(sample), path("filtered_barcode"), emit: filtered_barcode

    script:
    """
    set -e
    echo "estimated cells"
    step3_estimated_cells.py -r ./ --outdir . -e ${params.expected_cell_num}
    """
}

process GTF_TO_GENE_BED {
    tag "GTF_TO_GENE_BED"
    publishDir "${params.outdir}/"
    time "30m"
    resourceLabels label: "GTF_TO_GENE_BED_${params.project}"
    
    output:
    path("geneslop2k.bed"), emit: gene_bed
    
    script:
    """
    set -e
    echo "gtf to gene bed"
    gtf_to_gene_bed.py --gtf ${params.gtf} --out genes.bed;
    bedtools slop -b 2000 -i genes.bed -g ${params.chrom_size_path_full} > geneslop2k.bed
    """
}
