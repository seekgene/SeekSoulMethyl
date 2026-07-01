/*
 * Shared parameter initialization module
 * Used by rna_met and met_only subworkflows to avoid duplicate code
 */

// Sanitize whitespace in string values (replace spaces with underscores)
def sanitizeSpaces(v) {
    if (v == null) return null
    return v.toString().replaceAll(/\s+/, '_')
}

// Check if a value is null-like (null, empty string, or string "null")
def isNullLike(v) {
    if (v == null) return true
    def s = v.toString().trim()
    return s == '' || s.equalsIgnoreCase('null')
}

// Initialize common parameters from database_dir and other defaults
def initCommonParams(params, projectDir) {
    // Sanitize outdir whitespace
    def __outdir_raw = params.outdir?.toString()
    if (__outdir_raw != null && (__outdir_raw =~ /\s/)) {
        def __outdir_new = sanitizeSpaces(__outdir_raw)
        log.warn "outdir contains whitespace; replaced with '_' : '${__outdir_raw}' -> '${__outdir_new}'"
        params.outdir = __outdir_new
    }
    if (params.outdir?.toString() =~ /\s/) {
        error "Error: outdir still contains whitespace after normalization: '${params.outdir}'"
    }

    // Database file paths derived from database_dir
    if (isNullLike(params.genomeDir)) params['genomeDir'] = "${params.database_dir}/star"
    if (isNullLike(params.genomefa)) params['genomefa'] = "${params.database_dir}/fasta/genome.fa"
    if (isNullLike(params.gtf)) params['gtf'] = "${params.database_dir}/genes/genes.gtf"
    if (isNullLike(params.bismark_ref)) params['bismark_ref'] = "${params.database_dir}/fasta/"
    if (isNullLike(params.chrom_size_path_full)) params['chrom_size_path_full'] = "${params.database_dir}/bed/chr_len.bed"
    if (isNullLike(params.chrom_size_path)) {
        def __nochrM = "${params.database_dir}/bed/chr_nochrM.bed"
        params['chrom_size_path'] = file(__nochrM).exists() ? __nochrM : params.chrom_size_path_full
    }
    if (isNullLike(params.methy_barcode_wl)) params['methy_barcode_wl'] = "${projectDir}/bin/barcodes/U3CB_methylation.txt"

    // Chemistry-dependent defaults
    params.chemistry = params.chemistry ?: "DD-MET5"
    if (params.chemistry == "DD-MET3") {
        if (isNullLike(params.exp_chemistry)) params['exp_chemistry'] = "DDV2"
        if (isNullLike(params.cbcsv)) params['cbcsv'] = "${projectDir}/bin/barcodes/DD-M_bUCB3_whitelist.csv"
    } else if (params.chemistry == "DD-MET5") {
        if (isNullLike(params.exp_chemistry)) params['exp_chemistry'] = "DD-MET5"
        if (isNullLike(params.cbcsv)) params['cbcsv'] = "${projectDir}/bin/barcodes/ME5_bUCB3_whitelist.csv"
    }
    if (isNullLike(params.split_fastq)) params['split_fastq'] = 4
    if (isNullLike(params.filter_ch)) params['filter_ch'] = 2

    // Infer project name from outdir if not provided
    if (!params.containsKey('project') || !params.project || params.project.toString().trim() == '') {
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

    return params
}
