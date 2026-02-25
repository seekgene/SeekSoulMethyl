nextflow.enable.dsl=2

params.workflow = params.workflow ?: 'RNA_MET'

include { METHY_ONLY } from './subworkflows/met_only/main'
include { MAIN as RNA_MET } from './subworkflows/rna_met/main'
include { FORCE_CELL } from './subworkflows/forcecell/main'

workflow {
    def w = params.workflow.toString().trim().toLowerCase()
    if (w == 'rna_met' || w == 'rna-met' || w == 'rnamet') {
        RNA_MET()
    } else if (w == 'methy_only' || w == 'methy-only' || w == 'methy') {
        METHY_ONLY()
    } else if (w == 'force_cell' || w == 'force-cell' || w == 'forcecell' || w == 'force') {
        FORCE_CELL()
    } else {
        error "Error: unknown workflow '${params.workflow}'. Use one of: rna_met | methy_only | force_cell"
    }
}
