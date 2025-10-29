#!/usr/bin/env nextflow

/*
 * Single-cell RNA-seq and methylation analysis pipeline
 * Version: 4.0.0
 */

nextflow.enable.dsl=2

// Parameter definitions

// Input parameters
params.samplesheet = null

// Output directory
params.outdir = "./results"
    
// Database file paths
params.database_dir = ""
params.genomeDir = "${params.database_dir}/star"
params.genomefa = "${params.database_dir}/fasta/genome.fa"
params.gtf = "${params.database_dir}/genes/genes.gtf"
params.bismark_ref = "${params.database_dir}/fasta/"
params.chrom_size_path = "${params.database_dir}/bed/chr_nochrM.bed"
params.methy_barcode_wl = "${projectDir}/bin/barcodes/U3CB_methylation.txt"
params.chemistry = "DD-M"
if (params.chemistry == "DD-M") {
    params.exp_chemistry = "DDV2"
    params.cbcsv = "${projectDir}/bin/barcodes/DD-M_bUCB3_whitelist.csv"
    params.methy_barcode_wl = "${projectDir}/bin/barcodes/U3CB_methylation.txt"
}else if (params.chemistry == "ME5") {
    params.exp_chemistry = "ME5"
    params.cbcsv = "${projectDir}/bin/barcodes/ME5_bUCB3_whitelist.csv"
    params.methy_barcode_wl = "${projectDir}/bin/barcodes/U3CB_methylation.txt"
}
params.split_fastq = 4
params.filter_ch = 2
// Help information
params.help = false


// Help message
def helpMessage() {
    log.info"""
    Single-cell RNA-seq and methylation analysis pipeline - v4.0.0
    
    Usage:
    Batch sample analysis:
        nextflow run sc_methy_workflow_v4.nf --samplesheet samples.csv --outdir results
    """.stripIndent()
}

// Display help information
if (params.help) {
    helpMessage()
    exit 0
}

// Parameter validation
if (!params.samplesheet) {
    error "Error: --samplesheet parameter must be provided"
}

include {
   COMPUTE_CPG_SITES;
   FASTP_EXPRESSION_MULTI;
   FASTP_METHYLATION_MULTI;
   SEEKSOULTOOLS_RNA;
   METHYLATION_BARCODE_EXTRACTION;
   
   PARSE_FASTQ_FILES;
   CREATE_FORWARD_PAIRS;
   CREATE_REVERSE_PAIRS;
   FASTP_METHYLATION_BARCODE_EXTRACT as FASTP_METHYLATION_BARCODE_EXTRACT_F;
   FASTP_METHYLATION_BARCODE_EXTRACT as FASTP_METHYLATION_BARCODE_EXTRACT_R
  
} from './modules/step1'

include {
   BISMARK_ALIGNMENT_FORWARD;
   BISMARK_ALIGNMENT_REVERSE;
   SORT_BAM_BY_NAME as SORT_BAM_BY_NAME_F;
   SORT_BAM_BY_NAME as SORT_BAM_BY_NAME_R;
   
} from './modules/step2'

include {
    SPLIT_BAM_FILES as SPLIT_BAM_FILES_F;
    SPLIT_BAM_FILES as SPLIT_BAM_FILES_R;
    MERGE_BISMARK_BAM;
    MERGE_FILTERED_BARCODE_READS_COUNTS;
    ALLCOOLS_BAM_TO_ALLC;
    ALLCOOLS_GENERATE_DATASETS;
    ALLCOOLS_MERGE;
    ALLCOOLS_EXTRACT
} from './modules/step3'

include {
    METHYLATION_SUMMARY;
    METHYLATION_LSI_PCA_CLUSTERING;
    MULTI_REPORT
} from './modules/step4'

// Create input channel
def create_input_channel() {
    if (params.samplesheet) {
        // Batch sample mode - intelligent grouping and merging based on sample_id
        return Channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row ->
                // Parse each row data, collect all file paths
                def sample_id = row.sample_id
                def row_files = [:]
                
                // Collect all non-empty file paths
                row.each { key, value ->
                    if (key != 'sample_id' && value && value.trim() != '') {
                        if (!row_files[key]) {
                            row_files[key] = []
                        }
                        
                        // Handle comma-separated file paths
                        def file_paths = value.split(',')
                        file_paths.each { path ->
                            def trimmed_path = path.trim()
                            if (trimmed_path != '') {
                                // For OSS paths, handle directly as strings
                                if (trimmed_path.startsWith('oss://')) {
                                    row_files[key].add(trimmed_path)
                                } else {
                                    // Support wildcard paths
                                    if (trimmed_path.contains('*')) {
                                        row_files[key].addAll(file(trimmed_path).collect())
                                    } else {
                                        row_files[key].add(file(trimmed_path))
                                    }
                                }
                            }
                        }
                    }
                }
                
                return tuple(sample_id, row_files)
            }
            .groupTuple(by: 0) // Group by sample_id
            .map { sample_id, file_groups ->
                // Merge all files with the same sample_id
                def merged_files = [:]
                file_groups.each { file_group ->
                    file_group.each { key, file_list ->
                        if (!merged_files[key]) {
                            merged_files[key] = []
                        }
                        merged_files[key].addAll(file_list)
                    }
                }
                
                return tuple(sample_id, merged_files)
            }
    } else {
        // Single sample mode
        def exp_r1 = (params.expression_r1 && params.expression_r1.startsWith('oss://')) ? params.expression_r1 : file(params.expression_r1)
        def exp_r2 = (params.expression_r2 && params.expression_r2.startsWith('oss://')) ? params.expression_r2 : file(params.expression_r2)
        def methy_r1 = (params.methylation_r1 && params.methylation_r1.startsWith('oss://')) ? params.methylation_r1 : file(params.methylation_r1)
        def methy_r2 = (params.methylation_r2 && params.methylation_r2.startsWith('oss://')) ? params.methylation_r2 : file(params.methylation_r2)
        
        return Channel.of(tuple(
            params.sample_name,
            [
                expression_r1: [exp_r1],
                expression_r2: [exp_r2],
                methylation_r1: [methy_r1],
                methylation_r2: [methy_r2],
                download_url: [params.download_url ?: 'local']
            ]
        ))
    }
}

/*
 * Workflow definition
 */
workflow {
    // Create input channel
    input_ch = create_input_channel()
    
    // Total CG sites in genome
    cpg_sites = COMPUTE_CPG_SITES()
    
    // Build per-group tuples (no merging; each group goes through fastp)
    exp_groups = input_ch.map { sample_id, files ->
        def r1s = files['expression_r1'] ?: []
        def r2s = files['expression_r2'] ?: []
        def pairs = []
        def n = Math.max(r1s.size(), r2s.size())
        (0..<n).each { idx ->
            def r1 = idx < r1s.size() ? r1s[idx] : r1s.last()
            def r2 = idx < r2s.size() ? r2s[idx] : r2s.last()
            pairs.add(tuple(sample_id, "G${idx+1}", r1, r2))
        }
        return pairs
    }.flatMap { it }

    methy_groups = input_ch.map { sample_id, files ->
        def r1s = files['methylation_r1'] ?: []
        def r2s = files['methylation_r2'] ?: []
        def pairs = []
        def n = Math.max(r1s.size(), r2s.size())
        (0..<n).each { idx ->
            def r1 = idx < r1s.size() ? r1s[idx] : r1s.last()
            def r2 = idx < r2s.size() ? r2s[idx] : r2s.last()
            pairs.add(tuple(sample_id, "G${idx+1}", r1, r2))
        }
        return pairs
    }.flatMap { it }

    // Run fastp per group
    exp_clean_multi = FASTP_EXPRESSION_MULTI(exp_groups)
    methy_clean_multi = FASTP_METHYLATION_MULTI(methy_groups)

    // Assemble cleaned file lists for downstream multi-input tools (stageable paths)
    exp_clean_pairs = exp_clean_multi.rna_fastp_multi_data
        .groupTuple(by: 0)
        .map { t ->
            // t is [sample, groups, r1s, r2s]
            tuple(t[0], t[2], t[3])
        }

    methy_clean_pairs = methy_clean_multi.methy_fastp_multi_data
        .groupTuple(by: 0)
        .map { t ->
            // t is [sample, groups, r1s, r2s]
            tuple(t[0], t[2], t[3])
        }

    // RNA expression analysis with multi-group inputs
    rna_results = SEEKSOULTOOLS_RNA(exp_clean_pairs)

    // Methylation barcode extraction with multi-group inputs
    methy_barcode = METHYLATION_BARCODE_EXTRACTION(methy_clean_pairs)

    // Parse and group fastq files
    parsed_files = PARSE_FASTQ_FILES(methy_barcode.methy_barcode_output)
    
    // Create file pairs and distribute to Bismark alignment
    forward_pairs_raw = CREATE_FORWARD_PAIRS(parsed_files.forward_pairs
    .combine(methy_barcode.methy_barcode_output
    .map {it -> tuple(it[0], it[1], it[2])}, by:0))
    reverse_pairs_raw = CREATE_REVERSE_PAIRS(parsed_files.reverse_pairs.
    combine(methy_barcode.methy_barcode_output
    .map {it -> tuple(it[0], it[1], it[2])}, by:0))
    //forward_pairs_raw.view()
    //reverse_pairs_raw.view()
    // Process stdout output, convert multi-line strings to individual tuples
    // Calculate pair counts per sample for dynamic size in groupTuple operations
    pair_counts_per_sample = forward_pairs_raw
        .map { sample, stdout_content ->
            def count = 0
            stdout_content.split('\n').each { line ->
                if (line.trim()) {
                    def parts = line.split(',')
                    if (parts.size() == 4) {
                        count++
                    }
                }
            }
            return tuple(sample, count)
        }
    
    forward_pairs = forward_pairs_raw
        .map { sample, stdout_content ->
            def pairs = []
            stdout_content.split('\n').each { line ->
                if (line.trim()) {
                    def parts = line.split(',')
                    if (parts.size() == 4) {
                        def sample_id = parts[0]
                        def pair_id = parts[1]
                        def r1_file = file("${params.outdir}/${sample}_methy/step1/${parts[2]}")
                        def r2_file = file("${params.outdir}/${sample}_methy/step1/${parts[3]}")
                        pairs.add(tuple(sample_id, pair_id, r1_file, r2_file))
                    }
                }
            }
            return pairs
        }
        .flatMap { it }
    
    reverse_pairs = reverse_pairs_raw
        .map { sample, stdout_content ->
            def pairs = []
            stdout_content.split('\n').each { line ->
                if (line.trim()) {
                    def parts = line.split(',')
                    if (parts.size() == 4) {
                        def sample_id = parts[0]
                        def pair_id = parts[1]
                        def r1_file = file("${params.outdir}/${sample}_methy/step1/${parts[2]}")
                        def r2_file = file("${params.outdir}/${sample}_methy/step1/${parts[3]}")
                        pairs.add(tuple(sample_id, pair_id, r1_file, r2_file))
                    }
                }
            }
            return pairs
        }
        .flatMap { it }
    
    // Quality control again for fastq after barcode and adapter removal - process forward and reverse pairs
    //forward_pairs.view { "Forward pairs: $it" }
    
    methy_barcode_fastp_forward = FASTP_METHYLATION_BARCODE_EXTRACT_F(forward_pairs)
    methy_barcode_fastp_reverse = FASTP_METHYLATION_BARCODE_EXTRACT_R(reverse_pairs)
    
    
    // Bismark alignment - parallel processing of multiple file pairs
    bismark_aligned_forward = BISMARK_ALIGNMENT_FORWARD(forward_pairs)
    bismark_aligned_reverse = BISMARK_ALIGNMENT_REVERSE(reverse_pairs)
    
    // BAM file sorting (by read name)
    sorted_by_name_f = SORT_BAM_BY_NAME_F(bismark_aligned_forward.bismark_forward_bam)
    sorted_by_name_r = SORT_BAM_BY_NAME_R(bismark_aligned_reverse.bismark_reverse_bam)

    split_bams_f = SPLIT_BAM_FILES_F(sorted_by_name_f.bismark_sortn_bam.combine(rna_results.gex_barcodes,by:0))
    split_bams_r = SPLIT_BAM_FILES_R(sorted_by_name_r.bismark_sortn_bam.combine(rna_results.gex_barcodes,by:0))
    
        
    // merge single cell forward and reverse bismark bam
    sc_bismark_merge_bam = MERGE_BISMARK_BAM(split_bams_f.split_bams_dir
    .combine(split_bams_f.filtered_barcode, by:[0,1])
    .combine(split_bams_f.filtered_barcode_reads_counts, by:[0,1])
    .combine(
        split_bams_r.split_bams_dir
        .combine(split_bams_r.filtered_barcode, by:[0,1])
        .combine(split_bams_r.filtered_barcode_reads_counts, by:[0,1]), 
        by: [0,1]))
    
    // Run allcools bam-to-allc - requires multiple inputs
    allc_generated = ALLCOOLS_BAM_TO_ALLC(
        sc_bismark_merge_bam.sc_merged_bam_dir
        .combine(sc_bismark_merge_bam.merged_filtered_barcode, by: [0,1]))
    
    
    // Merge filtered barcode reads counts - 按样本聚合，样本间并行
    merged_counts = MERGE_FILTERED_BARCODE_READS_COUNTS(
    sc_bismark_merge_bam.merged_filtered_barcode
    .combine(pair_counts_per_sample, by: 0)
    .map { sample_id, pair_id, barcode_data, pair_count -> 
        tuple(groupKey(sample_id, pair_count), barcode_data)
    }
    .groupTuple()
    .map { group_key, barcode_list -> 
        tuple(group_key.toString(), barcode_list)
    }
    .combine(
        sc_bismark_merge_bam.merged_filtered_barcode_reads_counts
        .combine(pair_counts_per_sample, by: 0)
        .map { sample_id, pair_id, reads_data, pair_count -> 
            tuple(groupKey(sample_id, pair_count), reads_data)
        }
        .groupTuple()
        .map { group_key, reads_list -> 
            tuple(group_key.toString(), reads_list)
        },
    by: 0)
    .combine(
        allc_generated.allcools_allc_output
        .combine(pair_counts_per_sample, by: 0)
        .map { sample_id, pair_id, allc_data, pair_count -> 
            tuple(groupKey(sample_id, pair_count), allc_data)
        }
        .groupTuple()
        .map { group_key, allc_list -> 
            tuple(group_key.toString(), allc_list)
        }, by: 0)
    )
    
    
    // Run allcools generate datasets
    allcools_datasets = ALLCOOLS_GENERATE_DATASETS(
        allc_generated.allcools_allc_output
        .combine(pair_counts_per_sample, by: 0)
        .map { sample_id, pair_id, allc_data, pair_count -> 
            tuple(groupKey(sample_id, pair_count), allc_data)
        }
        .groupTuple()
        .map { group_key, allc_list -> 
            tuple(group_key.toString(), allc_list)
        }
        .combine(
        merged_counts.merged_filtered_barcode_reads_counts
        .map{it -> tuple(it[0], it[1])}, by: 0))
    
    // Run allcools merge datasets
    allcools_merged = ALLCOOLS_MERGE(allc_generated.allcools_allc_output
        .combine(pair_counts_per_sample, by: 0)
        .map { sample_id, pair_id, allc_data, pair_count -> 
            tuple(groupKey(sample_id, pair_count), allc_data)
        }
        .groupTuple()
        .map { group_key, allc_list -> 
            tuple(group_key.toString(), allc_list)
        })
    
    // Run allcools extract datasets
    allcools_extracted = ALLCOOLS_EXTRACT(allcools_merged.allcools_merge_allc)
    // Generate methylation summary report
    methylation_summary = METHYLATION_SUMMARY(
        bismark_aligned_forward.bismark_forward_report
        .combine(pair_counts_per_sample, by: 0)
        .map { sample_id, pair_id, report_data, pair_count -> 
            tuple(groupKey(sample_id, pair_count), report_data)
        }
        .groupTuple()
        .map { group_key, report_list -> 
            tuple(group_key.toString(), report_list)
        }
        .combine(
        bismark_aligned_reverse.bismark_reverse_report
        .combine(pair_counts_per_sample, by: 0)
        .map { sample_id, pair_id, report_data, pair_count -> 
            tuple(groupKey(sample_id, pair_count), report_data)
        }
        .groupTuple()
        .map { group_key, report_list -> 
            tuple(group_key.toString(), report_list)
        }, by: 0)
        .combine(merged_counts.allcools_cells_csv_output, by: 0)
        .combine(merged_counts.merged_filtered_barcode_reads_counts, by: 0)
        .combine(methy_barcode.methy_barcode_output.map {it -> tuple(it[0], it.last())}, by: 0)
        .combine(allcools_extracted.allcools_extract_allc_output, by: 0)
        .combine(cpg_sites.cpg_sites))
    
    // PCA clustering analysis
    pca_clustering = METHYLATION_LSI_PCA_CLUSTERING(
        allcools_datasets.allcools_generate_datasets
        .combine(
            merged_counts.merged_filtered_barcode_reads_counts
            .map{it -> tuple(it[0], it[1])}, by: 0))
    
    // Multi-report
    /*
    multi_report = MULTI_REPORT(
       rna_results.gex_summary_json
       .combine(rna_results.filtered_dir, by: 0)
       .combine(rna_results.raw_dir, by: 0)
       .combine(rna_results.counts_xls, by: 0)
       .combine(rna_results.detail_xls, by: 0)
       .combine(rna_results.tsne_umi, by: 0)
       .combine(rna_results.diff_data, by: 0)
       .combine(methylation_summary.methy_summary.map {it -> tuple(it[0], it[1])}, by: 0)
       .combine(merged_counts.merged_filtered_barcode_reads_counts.map {it -> tuple(it[0], it[2])}, by: 0)
       .combine(merged_counts.allcools_cells_csv_output, by: 0)
    )
    */
}

/*
* Workflow completion information
*/
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
    log.info "Total samples processed: ${workflow.success ? 'All samples completed successfully' : 'Some samples failed'}"
}

workflow.onError {
    log.error "Pipeline execution stopped with the following message: $workflow.errorMessage"
}
