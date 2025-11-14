// Generate summary
process METHYLATION_SUMMARY {
    tag "$sample-METHYLATION_SUMMARY"
    publishDir "${params.outdir}/${sample}/${sample}_methy/"
    
    input:
    tuple val(sample),path(bismark_forward_report), path(bismark_reverse_report), path(allcools_cells_csv_output), 
    path(filtered_barcode),path(filtered_barcode_reads_counts), path(summary_json), path(allcools_extract_allc), path(cpg_sites)
    
    output:
    tuple val(sample), path("${sample}_methy_summary.json"), path("${sample}_wgs_summary.csv"), emit: methy_summary
    
    script:
    """
    set -e
    # Generate summary report
    step4_wgs_summary.py \
        --outdir . \
        --samplename ${sample} \
        --summary_json ${sample}_methy_summary.json \
        --genome_info_json ${cpg_sites}
    """
}

// LSI or PCA reduction and clustering analysis
process METHYLATION_LSI_PCA_CLUSTERING {
    tag "METHYLATION_LSI_PCA_CLUSTERING-$sample"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step4/", mode: "copy"
    
    input:
    tuple val(sample), path(mcds_file), path(filtered_barcode)
    
    output:
    path "*.h5ad"
    path "*.pdf"
    path "*.png"
    
    script:
    """
    set -e 
    step4_allcools_PCA_cluster.py \
        --mcds_path ${mcds_file} \
        --samplename ${sample} \
        --var_dim chrom20k \
        --filtered_barcode_file ${filtered_barcode} \
        --outdir . \
        --reduc lsi
    """
}

// generate RNA-MET multi report
process MULTI_REPORT {
    tag "MULTI_REPORT-$sample"
    publishDir "${params.outdir}/${sample}/", mode: "copy"
    
    input:
    tuple val(sample), path(gexjson), path(filtered_dir), 
    path(raw_dir), path(counts_xls), path(detail_xls), path(tsne_file), path(rna_diff_data),
    path(methyjson), path(methy_filtered_counts_file), path(methy_cells_info_file)
    

    output:
    tuple val(sample), path("${sample}_rna_methyl_report.html"), path("${sample}_rna_met.json"), emit: rna_methy_report
    
    script:
    """
    set -e
    step4_report_rna_met.py \
    --gexjson ${gexjson} \
    --metjson ${methyjson} \
    --tsne_file ${tsne_file} \
    --filtered_counts_file ${methy_filtered_counts_file} \
    --cells_file ${methy_cells_info_file} \
    --raw_dir ${raw_dir} \
    --counts_file ${counts_xls} \
    --detail_file ${detail_xls} \
    --gtf ${params.gtf} \
    --filtered_dir ${filtered_dir} \
    --diff_data ${rna_diff_data} \
    --whitelist_file ${params.cbcsv} \
    --outdir . \
    --samplename ${sample} \
    --rawname ${sample} \
    --nf_config ${projectDir}/nextflow.config 
    """
}