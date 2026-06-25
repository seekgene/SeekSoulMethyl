#!/usr/bin/env nextflow

/*
 * SeekOne DD Force Cell analysis pipeline
 * Version: 2.2.0
 */

nextflow.enable.dsl=2

import groovy.json.JsonSlurper

params.pre_analysis_path = params.pre_analysis_path ? params.pre_analysis_path.toString().replaceAll('/+$','') : null
params.outdir = params.outdir ? params.outdir.toString().replaceAll('/+$','') : "./results"
params.pre_outdir = params.pre_outdir ?: (params.pre_analysis_path ? "${params.pre_analysis_path}/nextflowYaml/results" : null)

params.sample = params.sample ?: null
params.force_cell_number = params.force_cell_number ?: null

params.database_dir = params.database_dir ?: ""

// Import shared parameter initialization
include { initCommonParams } from '../../modules/params_init'

include {
    PRECHECK_SAMPLE;
    RESOLVE_PRE_ANALYSIS_ROOT;
    RUN_RNA_FORCE;
    STAGE_METHY_ASSETS;
    STAGE_BISMARK_ASSETS;
    FORCE_CELL_BARCODE_DIFF;
    FORCE_CELL_APPLY_ALLOCOOLS_CHANGES;
    FORCE_CELL_DISCOVER_BUCKET_DIRS;
    FORCE_CELL_ADAPT_BUCKET_FOR_SUBMERGE;
    FORCE_CELL_RECOMPUTE_CELLS_METRICS;
    FORCE_CELL_UPDATE_FILTERED_READS_COUNTS;
    FORCE_CELL_MERGE_AND_SUBSET_MCDS;
} from '../../modules/force_cell_custom'

include { COMPUTE_CPG_SITES } from '../../modules/step1'
include { SORT_BAM_BY_NAME as SORT_BAM_BY_NAME_F; SORT_BAM_BY_NAME as SORT_BAM_BY_NAME_R } from '../../modules/step2'
include {
    SPLIT_BAM_FILES as SPLIT_BAM_FILES_F;
    SPLIT_BAM_FILES as SPLIT_BAM_FILES_R;
    MERGE_BISMARK_BAM;
    ALLCOOLS_BAM_TO_ALLC;
    MERGE_FILTERED_BARCODE_READS_COUNTS;
    ALLCOOLS_GENERATE_DATASETS_PART;
    MERGE_MCDS;
    ALLCOOLS_SUBMERGE;
    ALLCOOLS_MERGE;
    ALLCOOLS_EXTRACT
} from '../../modules/step3'
include { METHYLATION_SUMMARY; METHYLATION_LSI_PCA_CLUSTERING; MULTI_REPORT } from '../../modules/step4'
include { GTF_TO_GENE_BED } from '../../modules/utils'  


def _parseBismarkStrand = { String baseName ->
    if (baseName.contains('_forward_')) return 'forward'
    if (baseName.contains('_reverse_')) return 'reverse'
    return 'unknown'
}

def _parseBismarkPairId = { String baseName ->
    def m = (baseName =~ /_(forward|reverse)_([A-Za-z0-9]+)_(\\d+)_bismark_bt2_pe$/)
    if (m) return "${m[0][2]}_${m[0][3]}"
    def m2 = (baseName =~ /_(forward|reverse)_(\\d+)_bismark_bt2_pe$/)
    if (m2) return "${m2[0][2]}"
    return baseName
        .replaceFirst(/_(forward|reverse)_/, '_')
        .replaceFirst(/_bismark_bt2_pe$/, '')
}

def _fileNotEmpty = { filePath ->
    def f = new File(filePath.toString())
    if (!f.exists()) return false
    return f.length() > 0
}

// Help message
def helpMessage() {
    log.info"""
    SeekOne DD Force Cell analysis pipeline - v2.2.0

    Usage:
        nextflow run main.nf --workflow force_cell \
            --pre_analysis_path /path/to/previous_run \
            --outdir /path/to/force_results \
            --database_dir /path/to/reference \
            --sample sampleA,sampleB \
            --force_cell_number '{"sampleA":5000,"sampleB":6000}'
    """.stripIndent()
}

workflow FORCE_CELL {
    initCommonParams(params, projectDir)

    if (params.help) {
        helpMessage()
        exit 0
    }

    if (!params.pre_analysis_path || !params.pre_outdir) {
        error "Error: --pre_analysis_path (and/or --pre_outdir) must be provided"
    }
    if (!params.sample) {
        error "Error: --sample must be provided (comma-separated sample names)"
    }
    if (!params.force_cell_number) {
        error "Error: --force_cell_number must be provided (json map)"
    }

    def __outdir_clean = (params.outdir ?: '').toString().replaceAll('/+$','')
    def __pre_outdir_clean = (params.pre_outdir ?: '').toString().replaceAll('/+$','')
    if (__outdir_clean && __pre_outdir_clean && __outdir_clean == __pre_outdir_clean) {
        error "Error: force_cell --outdir must not be the same as the previous results directory (--pre_outdir)"
    }

    def samples = params.sample.toString().split(',').collect { it.trim() }.findAll { it }
    def forceMap = new JsonSlurper().parseText(params.force_cell_number.toString())
    def sampleTuples = samples.collect { s ->
        def n = forceMap[s]
        if (n == null) {
            throw new IllegalArgumentException("Missing force_cell_number for sample: ${s}")
        }
        tuple(s, n as Integer)
    }

    sample_force = Channel.fromList(sampleTuples)
    sample_only = Channel.fromList(samples)

    sample_ctx = RESOLVE_PRE_ANALYSIS_ROOT(sample_force).ctx
        .map { sample, force_cell, pre_outdir, addtag_dir ->
            tuple(sample, force_cell, pre_outdir, pre_outdir, addtag_dir)
        }

    precheck = PRECHECK_SAMPLE(sample_ctx)

    cpg_sites = COMPUTE_CPG_SITES()
    gene_bed = GTF_TO_GENE_BED()

    rna_forced = RUN_RNA_FORCE(precheck.ok)
    methy_assets = STAGE_METHY_ASSETS(precheck.ok)
    bismark_assets = STAGE_BISMARK_ASSETS(precheck.ok)

    bismark_bams_flat = bismark_assets.bismark_bams
        .flatMap { sample, bams ->
            def xs = (bams instanceof List) ? bams : [bams]
            xs.collect { tuple(sample, it) }
        }

    bismark_reports_flat = bismark_assets.bismark_reports
        .flatMap { sample, reports ->
            def xs = (reports instanceof List) ? reports : [reports]
            xs.collect { tuple(sample, it) }
        }

    barcode_diff = FORCE_CELL_BARCODE_DIFF(
        rna_forced.gex_barcodes
            .combine(methy_assets.filtered_counts_csv, by: 0)
            .map { sample, gex_barcodes, filtered_counts_csv ->
                tuple(sample, gex_barcodes, filtered_counts_csv)
            }
    )

    has_add_by_sample = barcode_diff.diff
        .map { sample, drop, add, add_gex, target -> tuple(sample, _fileNotEmpty(add)) }

    drop_barcodes = barcode_diff.diff.map { sample, drop, add, add_gex, target -> tuple(sample, drop) }
    add_barcodes = barcode_diff.diff.map { sample, drop, add, add_gex, target -> tuple(sample, add) }
    add_gex_barcodes = barcode_diff.diff.map { sample, drop, add, add_gex, target -> tuple(sample, add_gex) }
    target_methy_barcodes = barcode_diff.diff.map { sample, drop, add, add_gex, target -> tuple(sample, target) }

    bismark_bams_for_add = bismark_bams_flat
        .combine(has_add_by_sample, by: 0)
        .filter { sample, bam, has_add -> has_add }
        .map { sample, bam, has_add -> tuple(sample, bam) }

    bismark_tagged = bismark_bams_for_add
        .map { sample, bam ->
            def bn = bam.baseName.toString()
            def strand = _parseBismarkStrand(bn)
            def pair_id = _parseBismarkPairId(bn)
            tuple(sample, strand, pair_id, bam)
        }

    bismark_forward = bismark_tagged
        .filter { sample, strand, pair_id, bam -> strand == 'forward' }
        .map { sample, strand, pair_id, bam -> tuple(sample, pair_id, bam) }

    bismark_reverse = bismark_tagged
        .filter { sample, strand, pair_id, bam -> strand == 'reverse' }
        .map { sample, strand, pair_id, bam -> tuple(sample, pair_id, bam) }

    sorted_forward = SORT_BAM_BY_NAME_F(bismark_forward)
    sorted_reverse = SORT_BAM_BY_NAME_R(bismark_reverse)

    split_forward = SPLIT_BAM_FILES_F(
        sorted_forward.bismark_sortn_bam
            .combine(add_gex_barcodes, by: 0)
            .map { sample, pair_id, sortn_bam, add_gex ->
                tuple(sample, pair_id, sortn_bam, add_gex)
            }
    )

    split_reverse = SPLIT_BAM_FILES_R(
        sorted_reverse.bismark_sortn_bam
            .combine(add_gex_barcodes, by: 0)
            .map { sample, pair_id, sortn_bam, add_gex ->
                tuple(sample, pair_id, sortn_bam, add_gex)
            }
    )

    forward_data = split_forward.split_bams_dir
        .combine(split_forward.filtered_barcode, by: [0, 1])
        .combine(split_forward.filtered_barcode_reads_counts, by: [0, 1])
        .map { sample, pair_id, bam_dir, barcode, reads_counts ->
            tuple(sample, pair_id, bam_dir, barcode, reads_counts)
        }

    reverse_data = split_reverse.split_bams_dir
        .combine(split_reverse.filtered_barcode, by: [0, 1])
        .combine(split_reverse.filtered_barcode_reads_counts, by: [0, 1])
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

    add_bismark_merged = MERGE_BISMARK_BAM(combined_data)

    add_allcools = ALLCOOLS_BAM_TO_ALLC(
        add_bismark_merged.sc_merged_bam_dir
            .combine(add_bismark_merged.merged_filtered_barcode, by: [0, 1])
            .map { sample, pair_id, sc_merged_bam_dir, merged_filtered_barcode ->
                tuple(sample, pair_id, sc_merged_bam_dir, merged_filtered_barcode)
            }
    )

    add_allc_dirs_by_sample = sample_only
        .map { s -> tuple(s, []) }
        .mix(add_allcools.allcools_allc_output.map { sample, pair_id, dir -> tuple(sample, [dir]) })
        .groupTuple(by: 0)
        .map { sample, lists ->
            def out = []
            lists.each { v ->
                if (v instanceof List) out.addAll(v)
            }
            tuple(sample, out)
        }

    add_reads_counts_by_sample = sample_only
        .map { s -> tuple(s, []) }
        .mix(add_bismark_merged.merged_filtered_barcode_reads_counts.map { sample, pair_id, csv -> tuple(sample, [csv]) })
        .groupTuple(by: 0)
        .map { sample, lists ->
            def out = []
            lists.each { v ->
                if (v instanceof List) out.addAll(v)
            }
            tuple(sample, out)
        }

    allcools_forced = FORCE_CELL_APPLY_ALLOCOOLS_CHANGES(
        methy_assets.allcools_dir
            .combine(drop_barcodes, by: 0)
            .combine(add_allc_dirs_by_sample, by: 0)
            .map { sample, allcools_dir, drop_file, add_dirs ->
                tuple(sample, allcools_dir, drop_file, add_dirs)
            }
    )

    buckets = FORCE_CELL_DISCOVER_BUCKET_DIRS(allcools_forced.allcools_force)

    bucket_tuples = buckets.buckets_dir
        .flatMap { sample, buckets_dir ->
            def fs = buckets_dir.toFile().listFiles()
            if (fs == null) return []
            fs.toList().findAll { it.exists() }.collect { f ->
                tuple(sample, f.getName(), f.toPath())
            }
        }

    adapted_buckets = FORCE_CELL_ADAPT_BUCKET_FOR_SUBMERGE(bucket_tuples)

    allc_submerge = ALLCOOLS_SUBMERGE(adapted_buckets.adapted)

    allcools_merged = ALLCOOLS_MERGE(
        allc_submerge.allcools_submerge_allc
            .map { sample, bucket_id, allc_data -> tuple(sample, allc_data) }
            .groupTuple(by: 0)
            .map { sample, allc_list -> tuple(sample, allc_list) }
    )

    allcools_extracted = ALLCOOLS_EXTRACT(allcools_merged.allcools_merge_allc)

    cells_metrics = FORCE_CELL_RECOMPUTE_CELLS_METRICS(allcools_forced.allcools_force)

    updated_reads_counts = FORCE_CELL_UPDATE_FILTERED_READS_COUNTS(
        methy_assets.filtered_counts_csv
            .combine(drop_barcodes, by: 0)
            .combine(add_barcodes, by: 0)
            .combine(add_reads_counts_by_sample, by: 0)
            .map { sample, old_counts_csv, drop_file, add_file, add_counts_csvs ->
                tuple(sample, old_counts_csv, drop_file, add_file, add_counts_csvs)
            }
    )

    add_mcds_parts = ALLCOOLS_GENERATE_DATASETS_PART(
        add_allcools.allcools_allc_output
            .combine(gene_bed.gene_bed)
            .map { sample, pair_id, allc_dir, gene_bed_path ->
                tuple(sample, pair_id, allc_dir, gene_bed_path)
            }
    )

    add_mcds_parts_by_sample = sample_only
        .map { s -> tuple(s, []) }
        .mix(add_mcds_parts.allcools_generate_datasets_part.map { sample, pair_id, mcds_part -> tuple(sample, [mcds_part]) })
        .groupTuple(by: 0)
        .map { sample, lists ->
            def out = []
            lists.each { v ->
                if (v instanceof List) out.addAll(v)
            }
            tuple(sample, out)
        }

    mcds_final = FORCE_CELL_MERGE_AND_SUBSET_MCDS(
        methy_assets.mcds_dir
            .combine(add_mcds_parts_by_sample, by: 0)
            .combine(target_methy_barcodes, by: 0)
            .map { sample, mcds_orig, mcds_parts, target_file ->
                tuple(sample, mcds_orig, mcds_parts, target_file)
            }
    )

    bismark_reports_tagged = bismark_reports_flat
        .map { sample, report ->
            def bn = report.baseName.toString()
            def strand = _parseBismarkStrand(bn)
            tuple(sample, strand, report)
        }

    forward_reports_by_sample = sample_only
        .map { s -> tuple(s, []) }
        .mix(bismark_reports_tagged.filter { sample, strand, report -> strand == 'forward' }.map { sample, strand, report -> tuple(sample, [report]) })
        .groupTuple(by: 0)
        .map { sample, lists ->
            def out = []
            lists.each { v ->
                if (v instanceof List) out.addAll(v)
            }
            tuple(sample, out)
        }

    reverse_reports_by_sample = sample_only
        .map { s -> tuple(s, []) }
        .mix(bismark_reports_tagged.filter { sample, strand, report -> strand == 'reverse' }.map { sample, strand, report -> tuple(sample, [report]) })
        .groupTuple(by: 0)
        .map { sample, lists ->
            def out = []
            lists.each { v ->
                if (v instanceof List) out.addAll(v)
            }
            tuple(sample, out)
        }

    methy_summary = METHYLATION_SUMMARY(
        forward_reports_by_sample
            .combine(reverse_reports_by_sample, by: 0)
            .combine(cells_metrics.cells.map { sample, cells_csv, cells_json -> tuple(sample, cells_csv) }, by: 0)
            .combine(updated_reads_counts.updated.map { sample, filtered_barcode, filtered_counts_csv -> tuple(sample, filtered_barcode, filtered_counts_csv) }, by: 0)
            .combine(methy_assets.methy_summary_json, by: 0)
            .combine(allcools_extracted.allcools_extract_allc_output, by: 0)
            .combine(cpg_sites.cpg_sites)
            .map { sample, forward_reports, reverse_reports, cells_csv, filtered_barcode, filtered_counts_csv, summary_json, allc_extract, cpg_sites_json ->
                tuple(sample, forward_reports, reverse_reports, cells_csv, filtered_barcode, filtered_counts_csv, summary_json, allc_extract, cpg_sites_json)
            }
    )

    pca_cluster = METHYLATION_LSI_PCA_CLUSTERING(
        mcds_final.mcds_final
            .combine(updated_reads_counts.updated.map { sample, filtered_barcode, filtered_counts_csv -> tuple(sample, filtered_barcode) }, by: 0)
            .map { sample, mcds_file, filtered_barcode ->
                tuple(sample, mcds_file, filtered_barcode)
            }
    )

    multi_report = MULTI_REPORT(
        rna_forced.gex_summary_json
            .combine(rna_forced.filtered_dir, by: 0)
            .combine(rna_forced.raw_dir, by: 0)
            .combine(rna_forced.counts_xls, by: 0)
            .combine(rna_forced.detail_xls, by: 0)
            .combine(rna_forced.tsne_umi, by: 0)
            .combine(rna_forced.diff_data, by: 0)
            .combine(methy_summary.methy_summary.map { sample, json, wgs -> tuple(sample, json) }, by: 0)
            .combine(updated_reads_counts.updated.map { sample, filtered_barcode, filtered_counts_csv -> tuple(sample, filtered_counts_csv) }, by: 0)
            .combine(cells_metrics.cells.map { sample, cells_csv, cells_json -> tuple(sample, cells_csv) }, by: 0)
            .map { sample, gexjson, filtered_dir, raw_dir, counts_xls, detail_xls, tsne_file, rna_diff, methyjson, methy_filtered_counts, methy_cells ->
                tuple(sample, gexjson, filtered_dir, raw_dir, counts_xls, detail_xls, tsne_file, rna_diff, methyjson, methy_filtered_counts, methy_cells)
            }
    )
}
