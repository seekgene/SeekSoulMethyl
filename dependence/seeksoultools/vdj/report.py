import os
import json
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from ..utils.report.websummaryVDJ import websummaryBCR, websummaryTCR
from ..utils.helper import logger

def report(vdj_wd:Path, samplename:str, chain:str, leader_ref:Path|None,
            matrix:Path|None=None, rna_wd:Path|None=None, outdir:Path|None=None, pa:str="Outs"):
    if not outdir:
        outdir = vdj_wd
    # find clone
    turst_wd = vdj_wd / "Analysis" / "trust4"
    raw_airr_file = turst_wd / f"{samplename}_barcode_airr.tsv"
    final_out = turst_wd / f"{samplename}_final.out"
    annofile = turst_wd / f"{samplename}_annot.fa"
    if (rna_wd is not None) and (matrix is None):
        matrix = rna_wd / "Analysis" / "step3" / "filtered_feature_bc_matrix"
    from .find_clones import complete_airr
    ddl_wd = outdir / "Analysis" / "dandelion"
    ddl_wd.mkdir(exist_ok=True, parents=True)
    
    logger.info("run dandelion")
    try:
        complete_airr(samplename, raw_airr_file, final_out, annofile, chain, ddl_wd, matrix, 3)
    except IndexError as e:
        logger.error(f"{str(e)}")

    # summary
    from .summary import get_summary
    get_summary(samplename, vdj_wd, chain, leader_ref, outdir, pa)

    # report
    if chain == "TR":
        object = websummaryTCR(
            summary_csv=f"{outdir}/{pa}/metrics_summary.csv",
            samplename = samplename,
        )
    elif chain == "IG":
        object = websummaryBCR(
            summary_csv=f"{outdir}/{pa}/metrics_summary.csv",
            samplename = samplename,
        )
    
    bccount =  turst_wd/f"{samplename}_bc_count.tsv"
    _bccount = turst_wd/f"_{samplename}_bc_count.tsv"
    data_json = object.to_json(
        clonotypes_tsv=f"{outdir}/{pa}/{samplename}_clonotypes.tsv",
        airr_tsv=f"{outdir}/{pa}/{samplename}_airr_rearrangement.tsv",
        bc_count_tsv=_bccount if _bccount.exists() else bccount
    )
    if rna_wd is not None:
        from ..utils.report.websummaryRNA import websummaryRNA
        object = websummaryRNA(
            summary_json=f"{rna_wd}/Analysis/{samplename}_summary.json",
        )

        rna_json = object.to_json(
            diff_table=f"{rna_wd}/Analysis/step4/FindAllMarkers.xls",
            dim_table=f"{rna_wd}/Analysis/step4/tsne_umi.xls",
            filtered_dir=f"{rna_wd}/Analysis/step3/filtered_feature_bc_matrix",
            raw_dir=f"{rna_wd}/Analysis/step3/raw_feature_bc_matrix",
            T_airr=f"{vdj_wd}/{pa}/{samplename}_airr_rearrangement.tsv" if chain=="TR" else None,
            B_airr=f"{vdj_wd}/{pa}/{samplename}_airr_rearrangement.tsv" if chain=="IG" else None,
        )
        data_json.update(rna_json)

    template_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/report'))
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('base.html')
    with open(os.path.join(outdir, pa,  f'{samplename}_report.html'), 'w') as fh:
        fh.write(template.render(websummary_json_data=json.dumps(data_json).replace("5'", "5\\'").replace("3'", "3\\'")))
    return data_json
