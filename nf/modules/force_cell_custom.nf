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
    path("addtagYaml/${sample}_exp/${sample}_exp_SortedByCoordinate_withTag.bam")
    path("addtagYaml/${sample}_exp/${sample}_exp_SortedByCoordinate_withTag.bam.bai")
    path("nextflowYaml/results/${sample}/")

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
        echo -e "\${pre_path}\t\${pre_path}\t\$\${pre_path}/../../addtagYaml" 
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
    test -n "${addtag_dir}"

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
    def cores = Math.max(1, task.cpus - 2)
    def samtools_cores = Math.max(1, Math.min(4, cores))
    """
    set -e
    mkdir -p ${sample}/Analysis
    cp ${pre_outdir}/${sample}/${sample}_exp/${sample}/Analysis/${sample}_gex_summary.json ${sample}/Analysis/${sample}_summary.json
    if [[ "${params.pre_analysis_path}" != oss://* ]] && [ -f "${counts_xls}" ] && [ -f "${detail_xls}" ] && [ -d "${raw_dir}" ]; then 
      seeksoultools rna callcell \
      --raw_matrix ${raw_dir} \
      --outdir . \
      --samplename ${sample} \
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

    ln -s "${pre_outdir}/${sample}/${sample}_methy/step3/split_bams/merged/${sample}_cells.csv" merged/${sample}_cells.csv
    ln -s "${pre_outdir}/${sample}/${sample}_methy/step3/split_bams/merged/filtered_barcode_reads_counts.csv" merged/filtered_barcode_reads_counts.csv
    ln -s "${pre_outdir}/${sample}/${sample}_methy/${sample}_methy_summary.json" ${sample}_methy_summary.json

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

