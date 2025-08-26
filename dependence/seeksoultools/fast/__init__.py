import os
import json
import click
from ..utils.helper import logger, include_introns_callback, validate_cores
from ..utils.helper import check_all
_steps = {
    "step1": [],
    "step2": ["STAR", "SortByPos", "SortByName", "qualimap", "FeatureCounts", ],
    "step3": [],
    "step4": [],
}

@click.group(help="quantifies singlecell gene expression")
@click.option("--steps", default=None, type=click.Path(), help="json format.")
@click.pass_obj
def fast(obj, steps):
    if  steps:
        with open(steps) as fh:
            obj["steps"] = json.load(fh)
    else:
        obj["steps"] = _steps

@fast.command(help="extract cell barcode and umi.")
@click.pass_obj
@click.option("--fq1", "fq1", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True,
              help="Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True,
              help="Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(resolve_path=True),
              help="Output dir.")
@click.option("--chemistry", type=click.Choice(["DD-Q",]),
              help="DD-Q.")
@click.option("--shift", is_flag=True, default=False, hidden=True,
              help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", hidden=True,
              help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, hidden=True,
              help="Barcode white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--structure", hidden=True,
              help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, hidden=True,
              help="Linker white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--skip_len", "use_short_read", is_flag=True, default=False, show_default=True, help="Skip filtering short reads after adapter filter, short reads will be used.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@validate_cores
def step1(obj, **kwargs):
    if "step1" in obj["steps"]:
        kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
        os.makedirs(kwargs["outdir"], exist_ok=True)
        from ..utils.barcode import check_rna_options
        chemistry_kwargs = check_rna_options(**kwargs)
        kwargs.update(chemistry_kwargs)
        with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
            json.dump(kwargs, fh, indent=4)

        from ..utils.barcode import barcode_main
        barcode_main(**kwargs)

@fast.command(help="align reads to genome.")
@click.pass_obj
@click.option("--fq", "ofq", required=True, multiple=True, help="Read2 fq file")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(), help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--star_path", "star_path", default="STAR", help="Path to executable STAR aligner.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True, help="include introns or not.")
@click.option("--scoremin",  help="outFilterScoreMinOverLread", required=False)
@click.option("--matchnmin", help="outFilterMatchNminOverLread", required=False)
@click.option("--sc5p",is_flag=True,default=False,show_default=True,help="If set, the single cell data is considered as 5' data.")
@click.option('--fra', default=0.1, show_default=True,help="downsample fraction")
@click.option("--rRNAgenomeDir", "rRNAgenomeDir",required=False,
              help="rRNA genomeDir ")
@click.option("--rRNAgtf", "rRNAgtf", required=False,
              help="")
@validate_cores
def step2(obj, **kwargs):
    from .step2 import align
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    align(**kwargs)

@fast.command(help="quantifies.")
@click.option("--bam", required=True, help="Bam file which contain alignment info.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, help="cluster, adjacency, directional")
@click.option("--expectNum", "expectNum", default=3000, show_default=True, help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--scoremin", type=float, help="outFilterScoreMinOverLread", required=False)
@click.option("--matchnmin", type=float, help="outFilterMatchNminOverLread", required=False)
@click.option("--forceCell", "forceCell", help="Force pipeline to use this number of cells, skipping cell calling algorithm.",required=False)
@click.pass_obj
def step3(obj, **kwargs):
    from .step3 import count, cell_calling
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    count(**kwargs)
    raw_matrix = os.path.join(kwargs["outdir"], "step3", "raw_feature_bc_matrix")
    kwargs["raw_matrix"] = raw_matrix
    cell_calling(**kwargs)

@fast.command(help="cell calling")
@click.option("--raw_matrix", required=True, help="Raw Feature-barcode matrix.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True, help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", help="Force pipeline to use this number of cells, skipping cell calling algorithm.",required=False)
@click.pass_obj
def callcell(obj, **kwargs):
    from .step3 import cell_calling
    cell_calling(**kwargs)

@fast.command(help="seurat.")
@click.option("--matrix", type=click.Path(), help="Feature-barcode matrix.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--dims", default=15, show_default=True, help="Number of dimension used as input in dimensional reduction.")
@click.option("--minpct", default=0.1, show_default=True, help="Minimum percentage that a feature to be detected in either of the two groups of cells.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--logfc", default=0.25, show_default=True, help="Limit testing to genes to this number of fold difference(log_scale) between the two groups of cells")
@click.pass_obj
def step4(obj, **kwargs):
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    from .step4 import do_seurat
    do_seurat(**kwargs)

@fast.command(help="report.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--scoremin", help="outFilterScoreMinOverLread", required=False)
@click.option("--matchnmin", help="outFilterMatchNminOverLread", required=False)
@click.pass_obj
def report(obj, **kwargs):
    from .report import report
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"])
    os.makedirs(kwargs["outdir"], exist_ok=True)
    report(**kwargs)

@fast.command(help="run all steps.")
@click.pass_obj
@click.option("--fq1", "fq1", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True,
              help="Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True,
              help="Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(resolve_path=True),
              help="Output dir.")
@click.option("--shift", is_flag=True, default=False, hidden=True,
              help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", hidden=True,
              help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, hidden=True,
              help="Barcode white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--structure", hidden=True,
              help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, hidden=True,
              help="Linker white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--skip_len", "use_short_read", is_flag=True, default=False, show_default=True, 
              help="Skip filtering short reads after adapter filter, short reads will be used.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True, help="include introns or not.")
@click.option("--scoremin", type=float, help="outFilterScoreMinOverLread", required=False)
@click.option("--matchnmin", type=float, help="outFilterMatchNminOverLread", required=False)
@click.option('--fra', default=0.1, show_default=True,help="downsample fraction")
@click.option("--sc5p", is_flag=True, default=False, show_default=True, hidden=True,
              help="If set, the single cell data is considered as 5' data.")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to GTF file.")
@click.option("--star_path", "star_path", default="STAR", 
              help="Path to executable STAR aligner.")
@click.option("--chemistry", type=click.Choice(["DD-Q",]),
              help="DD-Q.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True,
              help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", required=False,
              help="Force pipeline to use this number of cells, skipping cell calling algorithm.")
@click.option("--rRNAgenomeDir", "rRNAgenomeDir",required=False, type=click.Path(exists=True, resolve_path=True),
              help="")
@click.option("--rRNAgtf", "rRNAgtf", required=False, type=click.Path(exists=True, resolve_path=True),
              help="rRNA annotation gtf ,required=False")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, hidden=True,
              help="cluster, adjacency, directional")
@validate_cores
def run(obj, **kwargs):

    report_out = os.path.join(kwargs["outdir"], kwargs["samplename"])
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    check_all()
    if "step1" in obj["steps"]:
        os.makedirs(kwargs["outdir"], exist_ok=True)
        from ..utils.barcode import check_rna_options
        chemistry_kwargs = check_rna_options(**kwargs)
        kwargs.update(chemistry_kwargs)
        with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
            json.dump(kwargs, fh, indent=4)

        from ..utils.barcode import barcode_main
        barcode_main(**kwargs)

    fq1 = os.path.join(kwargs['outdir'], 'step1', kwargs['samplename'] + '_1.fq.gz')
    fq2 = os.path.join(kwargs['outdir'], 'step1', kwargs['samplename'] + '_2.fq.gz')
    kwargs['ofq'] = [fq1, fq2]

    if "step2" in obj["steps"]:
        from .step2 import align
        align(**kwargs)

    bam = os.path.join(kwargs['outdir'],"step2", 'featureCounts', f'{kwargs["samplename"]}_SortedByName.bam')
    kwargs['bam'] = bam

    if "step3" in obj["steps"]:
        from .step3 import count
        count(**kwargs)

    raw_matrix = os.path.join(kwargs["outdir"], "step3", "raw_feature_bc_matrix")
    kwargs["raw_matrix"] = raw_matrix
    if "step3" in obj["steps"]:
        from .step3 import cell_calling
        cell_calling(**kwargs)

    matrix = os.path.join(kwargs["outdir"], "step3", "filtered_feature_bc_matrix")
    kwargs["matrix"] = matrix
    if "step4" in obj["steps"]:
        from .step4 import do_seurat
        do_seurat(**kwargs)

    from .report import report
    kwargs['outdir'] = report_out
    report(**kwargs)
