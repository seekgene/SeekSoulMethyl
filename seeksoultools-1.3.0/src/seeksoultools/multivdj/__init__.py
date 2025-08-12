from pathlib import Path
import click
from ..utils.helper import logger, include_introns_callback, str2path, validate_cores, check_dependencies, check_external_tools, check_all
from typing import List, Optional
import os
import shutil 

@click.group(help="5' RNA and vdj ")
@click.pass_context
def multivdj(obj):
    pass
@multivdj.command()
@click.option("--samplename", required=True, type=str, help="Sample name.")
@click.option("--rnafq1", required=True, multiple=True, type=click.Path(exists=True), help="Path(s) to read1 file(s) of rna.")
@click.option("--rnafq2", required=True, multiple=True, type=click.Path(exists=True), help="Path(s) to read2 file(s) of rna.")
@click.option("--tcrfq1", required=False, multiple=True, type=click.Path(exists=True), help="Path(s) to read1 file(s) of tcr.")
@click.option("--tcrfq2", required=False, multiple=True, type=click.Path(exists=True), help="Path(s) to read2 file(s) of tcr.")
@click.option("--bcrfq1", required=False, multiple=True, type=click.Path(exists=True), help="Path(s) to read1 file(s) of bcr.")
@click.option("--bcrfq2", required=False, multiple=True, type=click.Path(exists=True), help="Path(s) to read2 file(s) of bcr.")
@click.option("--star_path", "star_path", default="STAR",
              help="Path to executable STAR aligner.")
@click.option("--organism", type=click.Choice(["human", "mouse",]), required=False, help="Organism type.")
@click.option("--cfg", type=click.Path(exists=True), callback=str2path, required=False, help="Path to ref config.")
@click.option("--outdir", required=True, type=click.Path(), callback=str2path, help="Working directory.")
@click.option("--read_pair", default=False, is_flag=True, help="Use r1 and r2.")
@click.option("--core", default=8, type=int, help="Number of threads.")
@click.option("--reads_limit", hidden=True, default=80000, type=int, show_default=True, help="Read number limit for each barcode, negative means no limit.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True,
              help="include introns or not.")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to GTF file.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True,
              help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", required=False,
              help="Force pipeline to use this number of cells, skipping cell calling algorithm.")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, hidden=True,
              help="cluster, adjacency, directional")
@validate_cores
def run(
    samplename: str,
    rnafq1: List[Path],
    rnafq2: List[Path],
    tcrfq1: List[Path],
    tcrfq2: List[Path],
    bcrfq1: List[Path],
    bcrfq2: List[Path],
    genomeDir: str,
    star_path: str,
    gtf: str,
    organism: str,
    outdir: Path,
    read_pair: bool = False,
    cfg: Optional[Path] = None,
    core: int = 8,
    reads_limit: int = 80000,
    **kwargs
    ):
    ### check 
    has_tcr = bool(tcrfq1 and tcrfq2)
    has_bcr = bool(bcrfq1 and bcrfq2)
    if not (has_tcr or has_bcr):
        raise click.BadParameter("At least one set of TCR (tcrfq1 and tcrfq2) or BCR (bcrfq1 and bcrfq2) data")
    if (tcrfq1 and not tcrfq2 and  bcrfq2) or (bcrfq1  and not bcrfq2 and tcrfq2):
        raise click.BadParameter("Must be a complete TCR data (tcrfq1 and tcrfq2) or a BCR datat (bcrfq1 and bcrfq2)")
    ### check lib
    check_all()
    ### rna
    wd = outdir / samplename
    rna = wd / "rna"
    rnadir = wd / "rna" / "Analysis"
    kwargs.update({'chemistry': 'DD5V1'})
    from ..utils.barcode import check_rna_options
    chemistry_kwargs = check_rna_options(fq1=rnafq1, fq2=rnafq2, outdir=rna, **kwargs)
    kwargs.update(chemistry_kwargs)
    from ..utils.barcode import barcode_main
    barcode_main(
        fq1=rnafq1,
        fq2=rnafq2,
        samplename=samplename,
        outdir=rnadir,
        **kwargs)
    fq = os.path.join(rnadir, "step1", f"{samplename}_2.fq.gz")
    from ..rna.step2 import align
    align(
        fq=[fq], 
        genomeDir=genomeDir, gtf=gtf,
        samplename=samplename,
        outdir=rnadir, 
        core=core, 
        star_path=star_path, 
        **kwargs)
    bam = os.path.join(rnadir, "step2", "featureCounts",  f"{samplename}_SortedByName.bam")
    from ..rna.step3 import count
    count(bam, rnadir, gtf, **kwargs)
    raw_matrix = os.path.join(rnadir, "step3", "raw_feature_bc_matrix")
    from ..rna.step3 import cell_calling
    cell_calling(raw_matrix, rnadir, samplename, gtf, **kwargs)
    matrix = rnadir / "step3" / "filtered_feature_bc_matrix"
    matrix = os.path.join(rnadir, "step3", "filtered_feature_bc_matrix")
    from ..rna.step4 import do_seurat
    do_seurat(matrix=matrix,
              samplename=samplename,
              outdir=rnadir,
              dims=15, minpct=0.1,
              logfc=0.25, rscript_path="Rscript",
               **kwargs)
    from ..rna.report import report
    report(samplename=samplename, rna_wd=wd / "rna" )
    
    ### tcr:
    if tcrfq1 and tcrfq2:
        chain = "TR"
        tcrdir = wd / "tcr"
        from ..vdj.utils import get_config, get_white_list
        white_list = get_white_list()
        ref, fa, leader = get_config(organism, chain, cfg)
        from ..vdj.run import run_trust4
        run_trust4(
            "tcr", tcrfq1, tcrfq2, fa, ref, leader, organism, white_list,
            chain, tcrdir, core, read_pair, reads_limit
            )
        from ..vdj.report import report
        _ref, _fa, leader = get_config(organism, chain, cfg)
        report(tcrdir, "tcr", chain, leader )
    if bcrfq1 and bcrfq2:
        chain = "IG"
        bcrdir = wd / "bcr"
        from ..vdj.utils import get_config, get_white_list
        white_list = get_white_list()
        ref, fa, leader = get_config(organism, chain, cfg)
        from ..vdj.run import run_trust4
        run_trust4(
            "bcr", bcrfq1, bcrfq2, fa, ref, leader, organism, white_list, 
            chain, bcrdir, core, read_pair, reads_limit
            )
        from ..vdj.report import report
        _ref, _fa, leader = get_config(organism, chain, cfg)
        report(bcrdir, "bcr", chain, leader )
    ## all report
    from ..rna.report import report
    if bcrfq1 and bcrfq2:
       bcr_wd = bcrdir
    else:
       bcr_wd = None
    if tcrfq1 and tcrfq2:
       tcr_wd = tcrdir
    else:
       tcr_wd = None
    
    report(samplename=samplename, 
           rna_wd=rna,
           tcr_wd=tcr_wd,
           bcr_wd=bcr_wd,
           outdir=wd,
           organism=organism, 
           **kwargs)
    
    ### merge csv
    ana_dir = wd / 'Outs' / 'Analysis'
    shutil.rmtree(ana_dir, ignore_errors=True)
