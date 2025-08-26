import os
import json
import click
from ..utils.chemistry import CHEMISTRY
from ..utils.helper import logger, include_introns_callback, str2path, validate_cores
from ..utils.helper import check_all
_steps = {
    "step1": [],
    "step2": ["STAR", "SortByPos", "FeatureCounts", "SortByName"],
    "step3": [],
    "step4": [],
}

@click.group(help="quantifies singlecell gene expression")
@click.option("--steps", default=None, type=click.Path(), help="json format.")
@click.pass_obj
def rna(obj, steps):
    if  steps:
        with open(steps) as fh:
            obj["steps"] = json.load(fh)
    else:
        obj["steps"] = _steps



@rna.command(help="extract cell barcode and umi.")
@click.option("--fq1", "fq1", required=True, type=click.Path(), multiple=True, help="Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(), multiple=True, help="Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--shift", is_flag=True, default=False, show_default=True, help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, help="Barcode white list file, can specify multiple times.")
@click.option("--structure", help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, help="Linker white list file, can specify multiple times.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--skip_len", "use_short_read", is_flag=True, default=False, show_default=True, help="Skip filtering short reads after adapter filter, short reads will be used.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--chemistry", type=click.Choice(["DDV1", "DDV2", "DDVS", "DD5V1", "MM", "MM-D", "DD-Q", "custom"]), help="DDV1, DDV2, DDVS, DD5V1, MM, MM-D, DD-Q.")
@validate_cores
@click.pass_obj
def step1(obj, **kwargs):
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    from ..utils.barcode import check_rna_options
    chemistry_kwargs = check_rna_options(**kwargs)
    kwargs.update(chemistry_kwargs)
    with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
        json.dump(kwargs, fh, indent=4)
    from ..utils.barcode import barcode_main
    barcode_main(**kwargs)

@rna.command(help="align reads to genome.")
@click.option("--fq", required=True, multiple=True, help="Read2 fq file")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(), help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--star_path", "star_path", default="STAR", help="Path to executable STAR aligner.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True, help="include introns or not.")
@click.option("--sc5p",is_flag=True,default=False,show_default=True,help="If set, the single cell data is considered as 5' data.")
@validate_cores
@click.pass_obj

def step2(obj, **kwargs):
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    from .step2 import align
    align(**kwargs)

@rna.command(help="quantifies.")
@click.option("--bam", required=True, help="Bam file which contain alignment info.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, help="cluster, adjacency, directional")
@click.option("--expectNum", "expectNum", default=3000, show_default=True, help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", help="Force pipeline to use this number of cells, skipping cell calling algorithm.",required=False)
@click.pass_obj
def step3(obj, **kwargs):
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    from .step3 import count, cell_calling
    count(**kwargs)
    raw_matrix = os.path.join(kwargs["outdir"], "step3", "raw_feature_bc_matrix")
    kwargs["raw_matrix"] = raw_matrix
    cell_calling(**kwargs)

@rna.command(help="cell calling")
@click.option("--raw_matrix", required=True, help="Raw Feature-barcode matrix.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True, help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", help="Force pipeline to use this number of cells, skipping cell calling algorithm.",required=False)
@click.pass_obj
def callcell(obj, **kwargs):
    from .step3 import cell_calling
    cell_calling(**kwargs)


@rna.command(help="seurat.")
@click.option("--matrix", type=click.Path(), help="Feature-barcode matrix.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--dims", default=15, show_default=True, help="Number of dimension used as input in dimensional reduction.")
@click.option("--minpct", default=0.1, show_default=True, help="Minimum percentage that a feature to be detected in either of the two groups of cells.")
@click.option("--logfc", default=0.25, show_default=True, help="Limit testing to genes to this number of fold difference(log_scale) between the two groups of cells")
@click.pass_obj
def step4(obj, **kwargs):
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    from .step4 import do_seurat
    do_seurat(**kwargs)

@rna.command(help="report.")
@click.option("--samplename", required=True, type=str, help="Sample name.")
@click.option("--rna_wd", required=True, type=click.Path(), callback=str2path, help="Working directory.")
@click.option("--tcr_wd", type=click.Path(), callback=str2path, help="TCR working directory.")
@click.option("--bcr_wd", type=click.Path(), callback=str2path, help="BCR working directory.")
@click.option("--cfg", type=click.Path(exists=True), callback=str2path, required=False, help="Path to ref config for VDJ analysis.")
@click.option("--outdir", type=click.Path(), callback=str2path, help="Output dir.")
@click.option("--organism", type=click.Choice(["human", "monkey", "mouse", "rabbit", "rat"]), required=False, help="Organism type.")
@click.pass_obj
def report(obj, **kwargs):
    from .report import report
    report(**kwargs)


@rna.command(help="run all steps.")
@click.pass_obj
@click.option("--fq1", "fq1", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True,
              help="Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(exists=True, resolve_path=True),
              multiple=True, help="Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True,
              help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(resolve_path=True),
              help="Output dir.")
@click.option("--shift", is_flag=True, default=False, hidden=True,
              help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", hidden=True,
              help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, hidden=True,
              help="Barcode white list file, can specify multiple times, used to describe read1 structure.")
# @click.option("--match_type", multiple=True, hidden=True,
#               help="Match_types of barcode during barcode correction. Default is 1")
@click.option("--structure", hidden=True,
              help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, hidden=True,
              help="Linker white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True,
              help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True,
              help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True,
              help="Do not rescue barcode match multi when do correction.")
@click.option("--skip_len", "use_short_read", is_flag=True, default=False, show_default=True, 
              help="Skip filtering short reads after adapter filter, short reads will be used.")
@click.option("--core", default=4, show_default=True,
              help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True,
              help="include introns or not.")
@click.option("--sc5p", is_flag=True, default=False, show_default=True, hidden=True,
              help="If set, the single cell data is considered as 5' data.")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to GTF file.")
@click.option("--star_path", "star_path", default="STAR",
              help="Path to executable STAR aligner.")
@click.option("--chemistry", type=click.Choice(["DDV1", "DDV2", "DDVS", "DD5V1", "MM", "MM-D", "DD-Q", "Auto", "custom"]),
              help="DDV1, DDV2, DDVS, DD5V1, MM, MM-D, DD-Q, Auto.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True,
              help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", required=False,
              help="Force pipeline to use this number of cells, skipping cell calling algorithm.")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, hidden=True,
              help="cluster, adjacency, directional")
@validate_cores
def run(obj, **kwargs):
    report_out = os.path.join(kwargs["outdir"], kwargs["samplename"])
    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"], "Analysis")
    os.makedirs(kwargs["outdir"], exist_ok=True)
    check_all()
    if "step1" in obj["steps"]:
        from ..utils.barcode import check_rna_options
        chemistry_kwargs = check_rna_options(**kwargs)
        kwargs.update(chemistry_kwargs)
        with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
            json.dump(kwargs, fh, indent=4)
        from ..utils.barcode import barcode_main
        barcode_main(**kwargs)
    fq = os.path.join(kwargs["outdir"], "step1", f"{kwargs['samplename']}_2.fq.gz")

    # paired: [r1.fq.gz, r2.fq.gz]
    kwargs["fq"] = [fq, ]

    if "step2" in obj["steps"]:
        from .step2 import align
        align(stpes=_steps["step2"], **kwargs)
    bam = os.path.join(kwargs["outdir"], "step2", "featureCounts",  f"{kwargs['samplename']}_SortedByName.bam")
    kwargs["bam"] = bam
   
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
    report(samplename=kwargs["samplename"], rna_wd=report_out)
