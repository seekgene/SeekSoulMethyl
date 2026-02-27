process RESOLVE_PRE_ANALYSIS_ROOT {
    tag "${sample}-RESOLVE_PRE_ANALYSIS_ROOT"
    resourceLabels label: "FORCE_CELL_RESOLVE_PRE_${params.project}_${sample}"
    publishDir "${params.outdir}/pre_analysis", mode: 'symlink', saveAs: { filename ->
        if (filename.startsWith("addtagYaml/") || filename.startsWith("nextflowYaml/")) {
            return filename
        }
        return null
    }

    input:
    tuple val(sample), val(force_cell)

    output:
    tuple val(sample), val(force_cell), stdout, emit: ctx
    path("addtagYaml/${sample}_exp/${sample}_exp_SortedByCoordinate_withTag.bam"),optional:true
    path("addtagYaml/${sample}_exp/${sample}_exp_SortedByCoordinate_withTag.bam.bai"),optional:true
    path("nextflowYaml/results/${sample}/"),optional:true

    script:
    """
    set -e
    oss_cp() {
        local src="\$1"
        local dst="\$2"
        mkdir -p "\$(dirname "\$dst")"
        ossutil cp --sign-version v4 --region cn-beijing -e oss-cn-beijing-internal.aliyuncs.com -i "\$AccessKeyId" -k "\$AccessKeySecret" "\$src" "\$dst" 1>&2
        
    }
    oss_cp_r() {
        local src="\$1"
        local dst="\$2"
        mkdir -p "\$dst"
        ossutil cp --sign-version v4 --region cn-beijing -e oss-cn-beijing-internal.aliyuncs.com -i "\$AccessKeyId" -k "\$AccessKeySecret" -r "\$src" "\$dst" 1>&2
    }
    if [[ "${params.pre_analysis_path}" == oss://* ]]; then
        pre_path="${params.pre_analysis_path}"
        pre_local="./"
        pre_outdir_oss="${params.pre_outdir}"
        sample="${sample}"

        oss_cp "\${pre_outdir_oss}/\${sample}/\${sample}_rna_methy_summary.csv" "\${pre_local}/nextflowYaml/results/\${sample}/\${sample}_rna_methy_summary.csv"

        oss_cp "\${pre_path}/addtagYaml/\${sample}_exp/\${sample}_exp_SortedByCoordinate_withTag.bam" "\${pre_local}/addtagYaml/\${sample}_exp/\${sample}_exp_SortedByCoordinate_withTag.bam"
        oss_cp "\${pre_path}/addtagYaml/\${sample}_exp/\${sample}_exp_SortedByCoordinate_withTag.bam.bai" "\${pre_local}/addtagYaml/\${sample}_exp/\${sample}_exp_SortedByCoordinate_withTag.bam.bai"

        oss_cp_r "\${pre_outdir_oss}/\${sample}/\${sample}_exp/\${sample}/Analysis/" "\${pre_local}/nextflowYaml/results/\${sample}/\${sample}_exp/\${sample}/Analysis/" 
        oss_cp "\${pre_outdir_oss}/\${sample}/\${sample}_methy/step3/split_bams/merged/\${sample}_cells.csv" "\${pre_local}/nextflowYaml/results/\${sample}/\${sample}_methy/step3/split_bams/merged/\${sample}_cells.csv"
        oss_cp "\${pre_outdir_oss}/\${sample}/\${sample}_methy/step3/split_bams/merged/filtered_barcode_reads_counts.csv" "\${pre_local}/nextflowYaml/results/\${sample}/\${sample}_methy/step3/split_bams/merged/filtered_barcode_reads_counts.csv"
        oss_cp "\${pre_outdir_oss}/\${sample}/\${sample}_methy/${sample}_methy_summary.json" "\${pre_local}/nextflowYaml/results/\${sample}/\${sample}_methy/${sample}_methy_summary.json"

        oss_cp "\${pre_outdir_oss}/\${sample}/\${sample}_methy/step3/allcools.tar.gz" "\${pre_local}/nextflowYaml/results/\${sample}/\${sample}_methy/step3/allcools.tar.gz" || true
        oss_cp "\${pre_outdir_oss}/\${sample}/\${sample}_methy/step3/allcools_generate_datasets.tar.gz" "\${pre_local}/nextflowYaml/results/\${sample}/\${sample}_methy/step3/allcools_generate_datasets.tar.gz" || true

        oss_cp_r "\${pre_outdir_oss}/\${sample}/\${sample}_methy/step2/bismark/" "\${pre_local}/nextflowYaml/results/\${sample}/\${sample}_methy/step2/bismark"

        echo -e "${params.outdir}/pre_analysis\t${params.outdir}/pre_analysis/nextflowYaml/results\t${params.outdir}/pre_analysis/addtagYaml" 
    else
        pre_path="${params.pre_analysis_path}"
        pre_local="./"
        sample="${sample}"
        if [ -d "\${pre_path}/nextflowYaml/results/\${sample}" ]; then
            pre_path="\${pre_path}/nextflowYaml/results/"
        fi
        echo -e "\${pre_path}\t\${pre_path}\t\${pre_path}/../../addtagYaml" 
    fi
    """
}

process PRECHECK_SAMPLE {
    tag "${sample}-PRECHECK_SAMPLE"
    resourceLabels label: "FORCE_CELL_PRECHECK_${params.project}_${sample}"

    input:
    tuple val(sample), val(force_cell), val(pre_root), val(pre_outdir), val(addtag_dir)

    output:
    tuple val(sample), val(force_cell), val(pre_root), val(pre_outdir), val(addtag_dir), emit: ok

    script:
    """
    set -e

    test -n "${pre_root}"
    test -n "${pre_outdir}"

    test -f "${pre_outdir}/${sample}/${sample}_rna_methy_summary.csv"

    fc_bam="${pre_root}/${sample}/${sample}_exp/${sample}/Analysis/step2/featureCounts/${sample}_SortedByName.bam"
    if [[ "${params.pre_analysis_path}" != oss://* ]] && [ -f "\$fc_bam" ]; then
      :
    else
      test -f "${addtag_dir}/${sample}_exp/${sample}_exp_SortedByCoordinate_withTag.bam"
      test -f "${addtag_dir}/${sample}_exp/${sample}_exp_SortedByCoordinate_withTag.bam.bai"
    fi

    test -f "${pre_outdir}/${sample}/${sample}_methy/step3/split_bams/merged/${sample}_cells.csv"
    test -f "${pre_outdir}/${sample}/${sample}_methy/step3/split_bams/merged/filtered_barcode_reads_counts.csv"
    test -f "${pre_outdir}/${sample}/${sample}_methy/${sample}_methy_summary.json"

    if [ ! -d "${pre_outdir}/${sample}/${sample}_methy/step3/allcools" ] && [ ! -f "${pre_outdir}/${sample}/${sample}_methy/step3/allcools.tar.gz" ]; then
      exit 1
    fi
    if [ ! -d "${pre_outdir}/${sample}/${sample}_methy/step3/allcools_generate_datasets" ] && [ ! -f "${pre_outdir}/${sample}/${sample}_methy/step3/allcools_generate_datasets.tar.gz" ]; then
      exit 1
    fi

    ls "${pre_outdir}/${sample}/${sample}_methy/step2/bismark/"*_bismark_bt2_PE_report.txt >/dev/null 2>&1
    ls "${pre_outdir}/${sample}/${sample}_methy/step2/bismark/"*_bismark_bt2_pe.bam >/dev/null 2>&1
    """
}

process RUN_RNA_FORCE {
    tag "${sample}-RUN_RNA_FORCE"
    publishDir "${params.outdir}/${sample}/${sample}_exp/"
    resourceLabels label: "FORCE_CELL_RNA_${params.project}_${sample}"

    input:
    tuple val(sample), val(force_cell), val(pre_root), val(pre_outdir), val(addtag_dir)

    output:
    tuple val(sample), path("${sample}/Analysis/step3/filtered_feature_bc_matrix/barcodes.tsv.gz"), emit: gex_barcodes
    tuple val(sample), path("${sample}/Analysis/${sample}_gex_summary.json"), emit: gex_summary_json
    tuple val(sample), path("${sample}/Analysis/step3/filtered_feature_bc_matrix"), emit: filtered_dir
    tuple val(sample), path("${sample}/Analysis/step3/raw_feature_bc_matrix"), emit: raw_dir
    tuple val(sample), path("${sample}/Analysis/step3/counts.xls"), emit: counts_xls
    tuple val(sample), path("${sample}/Analysis/step3/detail.xls"), emit: detail_xls
    tuple val(sample), path("${sample}/Analysis/step4/tsne_umi.xls"), emit: tsne_umi
    tuple val(sample), path("${sample}/Analysis/step4/FindAllMarkers.xls"), emit: diff_data
    path("${sample}/"), emit: rna_root

    script:
    def bam_addtag = "${addtag_dir}/${sample}_exp/${sample}_exp_SortedByCoordinate_withTag.bam"
    def counts_xls = "${pre_outdir}/${sample}/${sample}_exp/${sample}/Analysis/step3/counts.xls"
    def detail_xls = "${pre_outdir}/${sample}/${sample}_exp/${sample}/Analysis/step3/detail.xls"
    def raw_dir = "${pre_outdir}/${sample}/${sample}_exp/${sample}/Analysis/step3/raw_feature_bc_matrix"
    def fc_bam = "${pre_outdir}/${sample}/${sample}_exp/${sample}/Analysis/step2/featureCounts/${sample}_SortedByName.bam"
    def cores = Math.max(1, task.cpus - 2)
    def samtools_cores = Math.max(1, Math.min(4, cores))
    """
    set -e
    mkdir -p ${sample}/Analysis
    cp ${pre_outdir}/${sample}/${sample}_exp/${sample}/Analysis/${sample}_gex_summary.json ${sample}/Analysis/${sample}_summary.json
    if [[ "${params.pre_analysis_path}" != oss://* ]] && [ -f "${fc_bam}" ]; then 
      seeksoultools rna step3 \
      --bam ${fc_bam} \
      --outdir . \
      --samplename ${sample} \
      --gtf ${params.gtf} \
      --forceCell ${force_cell}
    else
      samtools view -h "${bam_addtag}" | awk 'BEGIN {OFS="\\t"} {
        if (\$0 ~ /^@/) {
            print
        } else {
            cb=""; ur="";
            for(i=12;i<=NF;i++) {
                if(\$i ~ /^CB:Z:/) cb=substr(\$i, 6);
                if(\$i ~ /^UR:Z:/) ur=substr(\$i, 6);
            }
            if(cb!="" && ur!="") {
                \$1 = cb "_" ur "_" \$1
            }
            print
        }
      }' | samtools sort -n -@ ${samtools_cores} -o input.bam -
      seeksoultools rna step3 \
      --bam input.bam \
      --outdir . \
      --samplename ${sample} \
      --gtf ${params.gtf} \
      --forceCell ${force_cell}
    fi

    seeksoultools rna step4 \
      --matrix ${sample}/Analysis/step3/filtered_feature_bc_matrix \
      --outdir . \
      --samplename ${sample}

    seeksoultools rna report --samplename ${sample} --rna_wd ${sample}/ --outdir ${sample}/ --organism human

    if [ -f "${sample}/Analysis/${sample}_summary.json" ] && [ ! -f "${sample}/Analysis/${sample}_gex_summary.json" ]; then
      mv "${sample}/Analysis/${sample}_summary.json" "${sample}/Analysis/${sample}_gex_summary.json"
    fi
    """
}

process STAGE_METHY_ASSETS {
    tag "${sample}-STAGE_METHY_ASSETS"
    resourceLabels label: "FORCE_CELL_STAGE_METHY_${params.project}_${sample}"

    input:
    tuple val(sample), val(force_cell), val(pre_root), val(pre_outdir), val(addtag_dir)

    output:
    tuple val(sample), path("allcools"), emit: allcools_dir
    tuple val(sample), path("allcools_generate_datasets"), emit: datasets_dir
    tuple val(sample), path("allcools_generate_datasets/${sample}.mcds"), emit: mcds_dir
    tuple val(sample), path("merged/${sample}_cells.csv"), emit: cells_csv
    tuple val(sample), path("merged/filtered_barcode_reads_counts.csv"), emit: filtered_counts_csv
    tuple val(sample), path("${sample}_methy_summary.json"), emit: methy_summary_json

    script:
    """
    set -e
    mkdir -p merged

    cp "${pre_outdir}/${sample}/${sample}_methy/step3/split_bams/merged/${sample}_cells.csv" merged/${sample}_cells.csv
    cp "${pre_outdir}/${sample}/${sample}_methy/step3/split_bams/merged/filtered_barcode_reads_counts.csv" merged/filtered_barcode_reads_counts.csv
    cp "${pre_outdir}/${sample}/${sample}_methy/${sample}_methy_summary.json" ${sample}_methy_summary.json

    if [ -d "${pre_outdir}/${sample}/${sample}_methy/step3/allcools" ]; then
      ln -s "${pre_outdir}/${sample}/${sample}_methy/step3/allcools" allcools
    else
      ln -s "${pre_outdir}/${sample}/${sample}_methy/step3/allcools.tar.gz" allcools.tar.gz
      tar -xvzf allcools.tar.gz
    fi

    if [ -d "${pre_outdir}/${sample}/${sample}_methy/step3/allcools_generate_datasets" ]; then
      ln -s "${pre_outdir}/${sample}/${sample}_methy/step3/allcools_generate_datasets" allcools_generate_datasets
    else
      ln -s "${pre_outdir}/${sample}/${sample}_methy/step3/allcools_generate_datasets.tar.gz" allcools_generate_datasets.tar.gz
      tar -xvzf allcools_generate_datasets.tar.gz
    fi
    """
}

process STAGE_BISMARK_ASSETS {
    tag "${sample}-STAGE_BISMARK_ASSETS"
    resourceLabels label: "FORCE_CELL_STAGE_BISMARK_${params.project}_${sample}"

    input:
    tuple val(sample), val(force_cell), val(pre_root), val(pre_outdir), val(addtag_dir)

    output:
    tuple val(sample), path("bismark/*_bismark_bt2_pe.bam"), emit: bismark_bams
    tuple val(sample), path("bismark/*_bismark_bt2_PE_report.txt"), emit: bismark_reports

    script:
    """
    set -e
    mkdir -p bismark
    ln -s "${pre_outdir}/${sample}/${sample}_methy/step2/bismark/"* bismark/
    """
}

process FORCE_CELL_BARCODE_DIFF {
    tag "${sample}-FORCE_CELL_BARCODE_DIFF"
    resourceLabels label: "FORCE_CELL_BARCODE_DIFF_${params.project}_${sample}"
    publishDir "${params.outdir}/force_cell_barcode_diff/${sample}"

    input:
    tuple val(sample), path(gex_barcodes), path(methy_counts_csv)

    output:
    tuple val(sample), path("drop_barcodes.tsv"), path("add_barcodes.tsv"), path("add_gex_barcodes.tsv"), path("target_methy_barcodes.tsv"), emit: diff
    tuple val(sample), path("barcode_diff.meta.tsv"), emit: meta

    script:
    """
    set -e
    force_cell_barcode_diff.py \
      --gex_barcodes ${gex_barcodes} \
      --methy_counts_csv ${methy_counts_csv} \
      --chemistry ${params.chemistry} \
      --cbcsv ${params.cbcsv} \
      --outdir .
    """
}

process FORCE_CELL_APPLY_ALLOCOOLS_CHANGES {
    tag "${sample}-FORCE_CELL_APPLY_ALLOCOOLS_CHANGES"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/"
    resourceLabels label: "FORCE_CELL_APPLY_ALLOCOOLS_CHANGES_${params.project}_${sample}"

    input:
    tuple val(sample), path(allcools_work), path(drop_barcodes), path(add_allc_dirs)

    output:
    tuple val(sample), path("allcools"), emit: allcools_force

    script:
    def split_n = (params.split_fastq ?: 4) as Integer
    """
    set -e
    in_allcools="${allcools_work}"
    if [ -e allcools ] && [ "\$(readlink allcools 2>/dev/null || true)" != "" -o "\$(test -d allcools && echo dir || true)" != "" ]; then
      if [ "\$(basename "${allcools_work}")" = "allcools" ]; then
        mv allcools allcools_in
        in_allcools="allcools_in"
      fi
    fi
    mkdir -p allcools
    cp -Lr "\${in_allcools}"/* allcools/

    removed_total=0
    while read bc; do
      bc="\${bc//\$'\\r'/}"
      bc="\${bc//[[:space:]]/}"
      [ -z "\$bc" ] && continue
      n=\$(find -L -maxdepth 1 allcools/ -type f \\( \
        -name "\${bc}_allc.gz" -o \
        -name "\${bc}_allc.gz.tbi" -o \
        -name "\${bc}_allc.gz.count.csv" \\
      \\) -print -delete | wc -l)
      removed_total=\$((removed_total+n))
    done < ${drop_barcodes}
    echo "Removed files: \$removed_total" > drop_removed_files.txt

    for d in ${add_allc_dirs}; do
      [ -d "\$d" ] || continue
      find -L "\$d" -maxdepth 1 -type f -name "*_allc.gz" | while read f; do
        bc=\$(basename "\$f" "_allc.gz")
        bucket=\${bc:0:${split_n}}
        dest=\$(find -L allcools/ -maxdepth 1 -type d -name "*_\${bucket}_*_allcools" -print | head -n 1 || true)
        if [ -z "\$dest" ]; then
          dest="allcools/${sample}_forward_\${bucket}_1_merged_fr_bam_allcools"
          mkdir -p "\$dest"
        fi
        cp -f "\$f" "\$dest/\${bc}_allc.gz"
        if [ -f "\$f.tbi" ]; then cp -f "\$f.tbi" "\$dest/\${bc}_allc.gz.tbi"; fi
        if [ -f "\$f.count.csv" ]; then cp -f "\$f.count.csv" "\$dest/\${bc}_allc.gz.count.csv"; fi
      done
    done
    """
}

process FORCE_CELL_DISCOVER_BUCKET_DIRS {
    tag "${sample}-FORCE_CELL_DISCOVER_BUCKET_DIRS"
    resourceLabels label: "FORCE_CELL_DISCOVER_BUCKET_DIRS_${params.project}_${sample}"

    input:
    tuple val(sample), path(allcools_force)

    output:
    tuple val(sample), path("buckets"), emit: buckets_dir

    script:
    """
    set -e
    mkdir -p buckets
    shopt -s nullglob
    n=0
    for d in ${allcools_force}/*_allcools; do
      if [ -d "\$d" ]; then
        target=\$(readlink -f "\$d")
        ln -s "\${target}" "buckets/\$(basename "\$d")"
        n=\$((n+1))
      fi
    done
    if [ "\$n" -eq 0 ]; then
      echo "ERROR: no bucket directories found under allcools_force (expected *_allcools). Check staging/allcools copy." 1>&2
      find -L ${allcools_force} -maxdepth 1 -type d -print 1>&2 || true
      exit 1
    fi
    """
}

process FORCE_CELL_ADAPT_BUCKET_FOR_SUBMERGE {
    tag "${sample}-${bucket_id}-FORCE_CELL_ADAPT_BUCKET_FOR_SUBMERGE"
    resourceLabels label: "FORCE_CELL_ADAPT_BUCKET_FOR_SUBMERGE_${params.project}_${sample}"

    input:
    tuple val(sample), val(bucket_id), path(bucket_dir)

    output:
    tuple val(sample), val(bucket_id), path("adapted_bucket"), emit: adapted

    script:
    """
    set -e
    mkdir -p adapted_bucket
    shopt -s nullglob
    for f in ${bucket_dir}/*_allc.gz; do
      bc=\$(basename "\$f" "_allc.gz")
      mkdir -p "adapted_bucket/\$bc"
      t=\$(readlink -f "\$f")
      ln -s "\$t" "adapted_bucket/\$bc/\${bc}_allc.gz"
      if [ -f "\$f.tbi" ]; then ln -s "\$(readlink -f "\$f.tbi")" "adapted_bucket/\$bc/\${bc}_allc.gz.tbi"; fi
      if [ -f "\$f.count.csv" ]; then ln -s "\$(readlink -f "\$f.count.csv")" "adapted_bucket/\$bc/\${bc}_allc.gz.count.csv"; fi
    done
    """
}

process FORCE_CELL_RECOMPUTE_CELLS_METRICS {
    tag "${sample}-FORCE_CELL_RECOMPUTE_CELLS_METRICS"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/split_bams/merged/"
    resourceLabels label: "FORCE_CELL_RECOMPUTE_CELLS_METRICS_${params.project}_${sample}"

    input:
    tuple val(sample), path(allcools_force)

    output:
    tuple val(sample), path("${sample}_cells.csv"), path("${sample}_cells.json"), emit: cells

    script:
    """
    set -e
    step3_merge_sc_metrics.py ${allcools_force} -o ./${sample}_cells --cbcsv ${params.cbcsv}
    """
}

process FORCE_CELL_UPDATE_FILTERED_READS_COUNTS {
    tag "${sample}-FORCE_CELL_UPDATE_FILTERED_READS_COUNTS"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/split_bams/merged/"
    resourceLabels label: "FORCE_CELL_UPDATE_FILTERED_READS_COUNTS_${params.project}_${sample}"

    input:
    tuple val(sample), path(old_counts_csv), path(drop_barcodes), path(add_barcodes), path(add_reads_counts_csvs)

    output:
    tuple val(sample), path("filtered_barcode"), path("filtered_barcode_reads_counts.csv"), emit: updated

    script:
    """
    set -e
    args=""
    for f in ${add_reads_counts_csvs}; do
      args="\${args} --add_reads_counts_csv \${f}"
    done
    mv ${old_counts_csv} old_filtered_barcode_reads_counts.csv
    old_counts_csv="old_filtered_barcode_reads_counts.csv"
    force_cell_update_reads_counts.py \
      --old_counts_csv \${old_counts_csv} \
      --drop_barcodes ${drop_barcodes} \
      --add_barcodes ${add_barcodes} \
      \${args} \
      --outdir .
    """
}

process FORCE_CELL_MERGE_AND_SUBSET_MCDS {
    tag "${sample}-FORCE_CELL_MERGE_AND_SUBSET_MCDS"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/allcools_generate_datasets/"
    resourceLabels label: "FORCE_CELL_MERGE_AND_SUBSET_MCDS_${params.project}_${sample}"

    input:
    tuple val(sample), path(mcds_orig), path(add_mcds_parts), path(target_methy_barcodes)

    output:
    tuple val(sample), path("${sample}.mcds"), emit: mcds_final

    script:
    """
    set -e
    mv ${mcds_orig} input_orig.mcds
    if [ -z "${add_mcds_parts}" ]; then
      python ${projectDir}/bin/step3_mcds_manager.py subset \
        --input_path input_orig.mcds \
        --barcode_path ${target_methy_barcodes} \
        --output_path ${sample}.mcds \
        --mode stream
    else
      python ${projectDir}/bin/step3_mcds_manager.py merge \
        input_orig.mcds ${add_mcds_parts} \
        --output_path merged_all.mcds \
        --merge_mode stream
      python ${projectDir}/bin/step3_mcds_manager.py subset \
        --input_path merged_all.mcds \
        --barcode_path ${target_methy_barcodes} \
        --output_path ${sample}.mcds \
        --mode stream
    fi
    """
}
