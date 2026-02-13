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

process FORCE_CELL_PREPARE_WORKDIR {
    tag "${sample}-FORCE_CELL_PREPARE_WORKDIR"
    resourceLabels label: "FORCE_CELL_PREPARE_WORKDIR_${params.project}_${sample}"

    input:
    tuple val(sample), path(allcools_dir), path(datasets_dir), path(mcds_dir)

    output:
    tuple val(sample), path("allcools_work"), path("datasets_work"), path("datasets_work/${sample}.mcds"), emit: work

    script:
    """
    set -e
    mkdir -p allcools_work datasets_work
    (cp -al ${allcools_dir}/. allcools_work/ 2>/dev/null) || cp -a ${allcools_dir}/. allcools_work/
    (cp -al ${datasets_dir}/. datasets_work/ 2>/dev/null) || cp -a ${datasets_dir}/. datasets_work/
    test -d datasets_work/${sample}.mcds
    """
}

process FORCE_CELL_APPLY_ALLOCOOLS_CHANGES {
    tag "${sample}-FORCE_CELL_APPLY_ALLOCOOLS_CHANGES"
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/", mode: "symlink", overwrite: true, saveAs: { "allcools" }
    resourceLabels label: "FORCE_CELL_APPLY_ALLOCOOLS_CHANGES_${params.project}_${sample}"

    input:
    tuple val(sample), path(allcools_work), path(drop_barcodes), path(add_allc_dirs)

    output:
    tuple val(sample), path("allcools_force"), emit: allcools_force

    script:
    def split_n = (params.split_fastq ?: 4) as Integer
    """
    set -e
    mkdir -p allcools_force
    (cp -al ${allcools_work}/. allcools_force/ 2>/dev/null) || cp -a ${allcools_work}/. allcools_force/

    removed_total=0
    while read bc; do
      bc="\${bc//\$'\\r'/}"
      bc="\${bc//[[:space:]]/}"
      [ -z "\$bc" ] && continue
      n=\$(find allcools_force -type f \\( \
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
        dest=\$(find allcools_force -maxdepth 1 -type d -name "*_\${bucket}_*_allcools" -print | head -n 1 || true)
        if [ -z "\$dest" ]; then
          dest="allcools_force/force_added_\${bucket}_allcools"
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
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/"
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
    publishDir "${params.outdir}/${sample}/${sample}_methy/step3/"
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
