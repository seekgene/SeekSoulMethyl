#!/usr/bin/env nextflow

/*
 * SeekOne DD Single Cell Methylation analysis pipeline
 * Version: 2.1.2
 */

nextflow.enable.dsl=2

// Parameter definitions

params.samplesheet = params.samplesheet ?: null
params.outdir = params.outdir ?: "./results"

def sanitizeSpaces(v) {
    if (v == null) return null
    return v.toString().replaceAll(/\s+/, '_')
}

def __outdir_raw = params.outdir?.toString()
if (__outdir_raw != null && (__outdir_raw =~ /\s/)) {
    def __outdir_new = sanitizeSpaces(__outdir_raw)
    log.warn "outdir contains whitespace; replaced with '_' : '${__outdir_raw}' -> '${__outdir_new}'"
    params.outdir = __outdir_new
}
if (params.outdir?.toString() =~ /\s/) {
    error "Error: outdir still contains whitespace after normalization: '${params.outdir}'"
}
    
// Database file paths
params.database_dir = params.database_dir ?: ""

def _isNullLike = { v ->
    if (v == null) return true
    def s = v.toString().trim()
    return s == '' || s.equalsIgnoreCase('null')
}

if (_isNullLike(params.genomeDir)) params['genomeDir'] = "${params.database_dir}/star"
if (_isNullLike(params.genomefa)) params['genomefa'] = "${params.database_dir}/fasta/genome.fa"
if (_isNullLike(params.gtf)) params['gtf'] = "${params.database_dir}/genes/genes.gtf"
if (_isNullLike(params.bismark_ref)) params['bismark_ref'] = "${params.database_dir}/fasta/"
if (_isNullLike(params.chrom_size_path_full)) params['chrom_size_path_full'] = "${params.database_dir}/bed/chr_len.bed"
if (_isNullLike(params.chrom_size_path)) {
    def __nochrM = "${params.database_dir}/bed/chr_nochrM.bed"
    params['chrom_size_path'] = file(__nochrM).exists() ? __nochrM : params.chrom_size_path_full
}
if (_isNullLike(params.methy_barcode_wl)) params['methy_barcode_wl'] = "${projectDir}/bin/barcodes/U3CB_methylation.txt"

params.chemistry = params.chemistry ?: "DD-MET3"
if (params.chemistry == "DD-MET3") {
    if (_isNullLike(params.exp_chemistry)) params['exp_chemistry'] = "DDV2"
    if (_isNullLike(params.cbcsv)) params['cbcsv'] = "${projectDir}/bin/barcodes/DD-M_bUCB3_whitelist.csv"
}else if (params.chemistry == "DD-MET5") {
    if (_isNullLike(params.exp_chemistry)) params['exp_chemistry'] = "DD-MET5"
    if (_isNullLike(params.cbcsv)) params['cbcsv'] = "${projectDir}/bin/barcodes/ME5_bUCB3_whitelist.csv"
}
params.split_fastq = params.split_fastq ?: 4
params.filter_ch = params.filter_ch ?: 2

if (!params.project || params.project.toString().trim() == '') {
    def __outdir_clean = params.outdir?.toString()?.replaceAll('/+$','')
    def __tokens = __outdir_clean.split('/').findAll { it != '' }
    def __proj = 'project'
    if (__tokens) {
        def __idxProj = __tokens.indexOf('proj')
        if (__idxProj != -1 && __tokens.size() > __idxProj + 1) {
            __proj = __tokens[__idxProj + 1]
        } else if (__tokens[-1] == 'results') {
            def __prev = __tokens.size() >= 2 ? __tokens[-2] : null
            if (__prev && (__prev ==~ /\d{4}-\d{2}-\d{2}(?:-\d{2}-\d{2}(?:-\d{2})?)?/)) {
                __proj = __tokens.size() >= 3 ? __tokens[-3] : 'project'
            } else {
                __proj = __prev ?: 'project'
            }
        } else {
            __proj = __tokens[-1]
        }
    }
    params.project = __proj
}
params.help = params.help ?: false
// Help message
def helpMessage() {
    log.info"""
    Single-cell RNA-seq and methylation analysis pipeline - v2.1.2
    
    Usage:
    Batch sample analysis:
        nextflow run main.nf --samplesheet samples.csv --outdir results --database_dir refdata-cellranger-arc-GRCh38-2024-A
    """.stripIndent()
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
  
} from '../../modules/step1'

include {
   BISMARK_ALIGNMENT_FORWARD;
   BISMARK_ALIGNMENT_REVERSE;
   SORT_BAM_BY_NAME as SORT_BAM_BY_NAME_F;
   SORT_BAM_BY_NAME as SORT_BAM_BY_NAME_R;
   
} from '../../modules/step2'

include {
    SPLIT_BAM_FILES as SPLIT_BAM_FILES_F;
    SPLIT_BAM_FILES as SPLIT_BAM_FILES_R;
    MERGE_BISMARK_BAM;
    MERGE_FILTERED_BARCODE_READS_COUNTS;
    ALLCOOLS_BAM_TO_ALLC;
    ALLCOOLS_GENERATE_DATASETS;
    ALLCOOLS_SUBMERGE;
    ALLCOOLS_MERGE;
    ALLCOOLS_EXTRACT
} from '../../modules/step3'

include {
    METHYLATION_SUMMARY;
    METHYLATION_LSI_PCA_CLUSTERING;
    MULTI_REPORT
} from '../../modules/step4'

include {
    COUNTS_MAPPED_READS as COUNTS_MAPPED_READS_F;
    COUNTS_MAPPED_READS as COUNTS_MAPPED_READS_R;
    ESTIMATED_CELLS;
    GTF_TO_GENE_BED
} from '../../modules/utils'

// Helper: build a stable per-sample grouping key (avoid name clash with Nextflow's groupKey aggregator)
def sampleGroupKey(sample_id, pair_count) {
    return sample_id?.toString() ?: 'unknown'
}


// Create input channel
def create_input_channel() {
    if (params.samplesheet) {
        // Batch sample mode - intelligent grouping and merging based on sample_id
        return Channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row ->
                // Parse each row data, collect all file paths
                def __sample_raw = row.sample_id
                def sample_id = __sample_raw == null ? '' : __sample_raw.toString().trim()
                if (sample_id == '') {
                    error "Error: samplesheet contains empty sample_id"
                }
                def __sample_new = sanitizeSpaces(sample_id)
                if (__sample_new != sample_id) {
                    log.warn "samplesheet sample_id contains whitespace; replaced with '_' : '${sample_id}' -> '${__sample_new}'"
                    sample_id = __sample_new
                }
                if (sample_id =~ /\s/) {
                    error "Error: sample_id still contains whitespace after normalization: '${sample_id}'"
                }
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
workflow METHY_ONLY {
    if (params.help) {
        helpMessage()
        exit 0
    }
    def __workdir_raw = workflow.workDir?.toString()
    if (__workdir_raw != null && (__workdir_raw =~ /\s/)) {
        error "Error: work directory contains whitespace: '${__workdir_raw}'. Use --workdir (auto '_' normalization) or -w with no spaces."
    }
    if (!params.samplesheet) {
        error "Error: --samplesheet parameter must be provided"
    }
    // Create input channel
    input_ch = create_input_channel()
    
    // Total CG sites in genome
    cpg_sites = COMPUTE_CPG_SITES()

    // Generate gene bed file
    gene_bed = GTF_TO_GENE_BED()
    
    // Build per-group tuples (no merging; each group goes through fastp)
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
     .view { "METHY_GROUPS: ${it}" }

    // Run fastp per group
    methy_clean_multi = FASTP_METHYLATION_MULTI(methy_groups)

    // Compute expected fastp group counts per sample
    methy_group_counts = methy_groups
        .groupTuple(by: 0)
        .map { t -> tuple(t[0], t[1].size()) }

    // Assemble cleaned file lists per sample using groupTuple with groupKey
    methy_clean_pairs = methy_clean_multi.methy_fastp_multi_data
        .combine(methy_group_counts, by: 0)
        .map { sample_id, group_id, r1, r2, count ->
            tuple(groupKey(sample_id, count), r1, r2)
        }
        .groupTuple()
        .map { sample_key, r1_list, r2_list ->
            tuple(sample_key.toString(), r1_list, r2_list)
        }

    

    // Methylation barcode extraction with multi-group inputs
    methy_barcode = METHYLATION_BARCODE_EXTRACTION(methy_clean_pairs)

    // Parse and group fastq files
    parsed_files = PARSE_FASTQ_FILES(methy_barcode.methy_barcode_output)
    
    // Create file pairs and distribute to Bismark alignment
    forward_pairs_raw = CREATE_FORWARD_PAIRS(
        parsed_files.forward_pairs
        .combine(
            methy_barcode.methy_barcode_output
            .map {it -> tuple(it[0], it[1], it[2])}, by:0))
    reverse_pairs_raw = CREATE_REVERSE_PAIRS(
        parsed_files.reverse_pairs
        .combine(
            methy_barcode.methy_barcode_output
            .map {it -> tuple(it[0], it[1], it[2])}, by:0))
    //forward_pairs_raw.view()
    //reverse_pairs_raw.view()
    // Process stdout output, convert multi-line strings to individual tuples
    // Calculate pair counts per sample for dynamic size in groupTuple operations
    pair_counts_per_sample_f = forward_pairs_raw
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

    pair_counts_per_sample_r = reverse_pairs_raw
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

    pair_counts_per_sample = pair_counts_per_sample_f
        .join(pair_counts_per_sample_r, by: 0, remainder: true)
        .map { it ->
            def sample = it[0]
            def f_count = 0
            def r_count = 0
            if (it.size() == 3) {
                f_count = (it[1] ?: 0) as Integer
                r_count = (it[2] ?: 0) as Integer
            } else if (it.size() == 2) {
                if (it[1] instanceof List && it[1].size() >= 2) {
                    f_count = ((it[1][0] ?: 0) as Integer)
                    r_count = ((it[1][1] ?: 0) as Integer)
                } else {
                    f_count = (it[1] ?: 0) as Integer
                    r_count = 0
                }
            } else {
                throw new IllegalStateException("Unexpected join tuple size=${it.size()} value=${it}")
            }
            tuple(sample, Math.max(f_count, r_count))
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

    aligned_reads_counts_f = COUNTS_MAPPED_READS_F(bismark_aligned_forward.bismark_forward_bam)
    aligned_reads_counts_r = COUNTS_MAPPED_READS_R(bismark_aligned_reverse.bismark_reverse_bam)

    estiamted_cell = ESTIMATED_CELLS(
        aligned_reads_counts_f.cb_aligned_reads_counts
        .combine(pair_counts_per_sample, by: 0)
        .map { sample_id, pair_id, reads_counts, pair_count -> 
            tuple(groupKey(sample_id, pair_count), reads_counts)
        }
        .groupTuple()
        .map { group_key, reads_counts_list -> 
            tuple(group_key.toString(), reads_counts_list)
        }
        .combine(
            aligned_reads_counts_r.cb_aligned_reads_counts
            .combine(pair_counts_per_sample, by: 0)
            .map { sample_id, pair_id, reads_counts, pair_count -> 
                tuple(groupKey(sample_id, pair_count), reads_counts)
            }
            .groupTuple()
            .map { group_key, reads_counts_list -> 
                tuple(group_key.toString(), reads_counts_list)
            }, by : 0)
    )
    split_bams_f = SPLIT_BAM_FILES_F(sorted_by_name_f.bismark_sortn_bam.combine(estiamted_cell.filtered_barcode,by:0))
    split_bams_r = SPLIT_BAM_FILES_R(sorted_by_name_r.bismark_sortn_bam.combine(estiamted_cell.filtered_barcode,by:0))
    
        
    // merge single cell forward and reverse bismark bam
    forward_data = split_bams_f.split_bams_dir
        .combine(split_bams_f.filtered_barcode, by: [0, 1])
        .combine(split_bams_f.filtered_barcode_reads_counts, by: [0, 1])
        .map { sample, pair_id, bam_dir, barcode, reads_counts ->
            tuple(sample, pair_id, bam_dir, barcode, reads_counts)
        }

    reverse_data = split_bams_r.split_bams_dir
        .combine(split_bams_r.filtered_barcode, by: [0, 1])
        .combine(split_bams_r.filtered_barcode_reads_counts, by: [0, 1])
        .map { sample, pair_id, bam_dir, barcode, reads_counts ->
            tuple(sample, pair_id, bam_dir, barcode, reads_counts)
        }

    empty_bam_dir = file("${projectDir}/assets/empty_split_bams_dir")
    empty_barcode = file("${projectDir}/assets/empty_filtered_barcode")
    empty_reads_counts = file("${projectDir}/assets/empty_filtered_barcode_reads_counts.csv")

    combined_data = forward_data
        .join(reverse_data, by: [0, 1], remainder: true)
        .map { it ->
            def sample = it[0]
            def pair_id = it[1]

            def f_bam_dir = null
            def f_barcode = null
            def f_reads_counts = null
            def r_bam_dir = null
            def r_barcode = null
            def r_reads_counts = null

            if (it.size() == 8) {
                f_bam_dir = it[2]; f_barcode = it[3]; f_reads_counts = it[4]
                r_bam_dir = it[5]; r_barcode = it[6]; r_reads_counts = it[7]
            } else if (it.size() == 6) {
                if (it[2] == null) {
                    r_bam_dir = it[3]; r_barcode = it[4]; r_reads_counts = it[5]
                } else if (it[5] == null) {
                    f_bam_dir = it[2]; f_barcode = it[3]; f_reads_counts = it[4]
                } else if (it[5] instanceof List && it[5].size() >= 3) {
                    f_bam_dir = it[2]; f_barcode = it[3]; f_reads_counts = it[4]
                    r_bam_dir = it[5][0]; r_barcode = it[5][1]; r_reads_counts = it[5][2]
                } else {
                    throw new IllegalStateException("Unexpected join tuple shape: ${it}")
                }
            } else {
                throw new IllegalStateException("Unexpected join tuple size=${it.size()} value=${it}")
            }

            tuple(sample, pair_id,
                f_bam_dir ?: empty_bam_dir, f_barcode ?: empty_barcode, f_reads_counts ?: empty_reads_counts,
                r_bam_dir ?: empty_bam_dir, r_barcode ?: empty_barcode, r_reads_counts ?: empty_reads_counts)
        }

    sc_bismark_merge_bam = MERGE_BISMARK_BAM(combined_data)
    
    // Run allcools bam-to-allc - requires multiple inputs
    allc_generated = ALLCOOLS_BAM_TO_ALLC(
        sc_bismark_merge_bam.sc_merged_bam_dir
        .combine(sc_bismark_merge_bam.merged_filtered_barcode, by: [0,1]))

    pair_counts_effective = sc_bismark_merge_bam.sc_merged_bam_dir
        .map { sample_id, pair_id, bam_dir -> tuple(sample_id, pair_id) }
        .groupTuple(by: 0)
        .map { sample_id, pair_ids -> tuple(sample_id, pair_ids.size()) }

    // Merge filtered barcode reads counts 
    merged_counts = MERGE_FILTERED_BARCODE_READS_COUNTS(
    sc_bismark_merge_bam.merged_filtered_barcode
    .combine(pair_counts_effective, by: 0)
    .map { sample_id, pair_id, barcode_data, pair_count -> 
        tuple(groupKey(sample_id, pair_count), barcode_data)
    }
    .groupTuple()
    .map { group_key, barcode_list -> 
        tuple(group_key.toString(), barcode_list)
    }
    .combine(
        sc_bismark_merge_bam.merged_filtered_barcode_reads_counts
        .combine(pair_counts_effective, by: 0)
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
        .combine(pair_counts_effective, by: 0)
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
        .combine(pair_counts_effective, by: 0)
        .map { sample_id, pair_id, allc_data, pair_count -> 
            tuple(groupKey(sample_id, pair_count), allc_data)
        }
        .groupTuple()
        .map { group_key, allc_list -> 
            tuple(group_key.toString(), allc_list)
        }
        .combine(
        merged_counts.merged_filtered_barcode_reads_counts
        .map{it -> tuple(it[0], it[1])}, by: 0)
        .combine(gene_bed.gene_bed))
    if (params.split_fastq > 0) {
        // Run allcools merge for split dataset
        allc_submerge = ALLCOOLS_SUBMERGE(allc_generated.allcools_allc_output)
        // Run allcools merge datasets
        allcools_merged = ALLCOOLS_MERGE(
            allc_submerge.allcools_submerge_allc
            .combine(pair_counts_effective, by: 0)
            .map { sample_id, pair_id, allc_data , pair_count -> 
                // pass only the gz file paths as a collected list per sample
                tuple(groupKey(sample_id, pair_count), allc_data)
            }
            .groupTuple()
            .map { group_key, allc_list -> 
                tuple(group_key.toString(), allc_list)
            }
        )
    }else{
        allcools_merged = ALLCOOLS_MERGE(
            allc_generated.allcools_allc_output
            .combine(pair_counts_effective, by: 0)
            .map { sample_id, pair_id, allc_data, pair_count -> 
                // pass only the gz file paths as a collected list per sample
                tuple(groupKey(sample_id, pair_count), allc_data)
            }
            .groupTuple()
            .map { group_key, allc_list -> 
                tuple(group_key.toString(), allc_list)
            }
        )
    }
    
    // Run allcools extract datasets
    allcools_extracted = ALLCOOLS_EXTRACT(allcools_merged.allcools_merge_allc)
    // Generate methylation summary report
    methylation_summary = METHYLATION_SUMMARY(
        bismark_aligned_forward.bismark_forward_report
            .join(bismark_aligned_reverse.bismark_reverse_report, by: [0, 1], remainder: true)
            .map { it ->
                def sample_id = it[0]
                def pair_id = it[1]
                def forward_report = null
                def reverse_report = null
                if (it.size() == 4) {
                    forward_report = it[2]
                    reverse_report = it[3]
                } else if (it.size() == 3 && it[2] instanceof List && it[2].size() >= 2) {
                    forward_report = it[2][0]
                    reverse_report = it[2][1]
                } else {
                    throw new IllegalStateException("Unexpected join tuple shape: ${it}")
                }
                def empty_bismark_report = file("${projectDir}/assets/empty_bismark_report.txt")
                tuple(sample_id, pair_id, forward_report ?: empty_bismark_report, reverse_report ?: empty_bismark_report)
            }
            .combine(pair_counts_effective, by: 0)
            .map { sample_id, pair_id, forward_report, reverse_report, pair_count ->
                tuple(groupKey(sample_id, pair_count), forward_report, reverse_report)
            }
            .groupTuple()
            .map { group_key, forward_reports, reverse_reports ->
                tuple(group_key.toString(), forward_reports, reverse_reports)
            }
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
