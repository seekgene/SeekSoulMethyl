import os
import json
from jinja2 import Environment, FileSystemLoader
from functools import reduce
from ..utils.report.websummaryRNA import websummaryRNA
from ..utils.report.websummaryVDJ import websummaryBCR, websummaryTCR

def report(outdir, samplename, rnadir=None, tcrdir=None, bcrdir=None):
    os.makedirs(outdir, exist_ok=True)

    data_json_string = []
    if rnadir:
        object = websummaryRNA(
            summary_json=f"{rnadir}/{samplename}_summary.json",
        )
        data_json_string.append(object.to_json(
                diff_table=f"{rnadir}/step4/FindAllMarkers.xls",
                dim_table=f"{rnadir}/step4/tsne_umi.xls",
                filtered_dir=f"{rnadir}/step3/filtered_feature_bc_matrix",
                raw_dir=f"{rnadir}/step3/raw_feature_bc_matrix",
                T_csv=f"{tcrdir}/filtered_contig_annotations.csv" if tcrdir else None,
                B_csv=f"{bcrdir}/filtered_contig_annotations.csv" if bcrdir else None,
            )
        )
    if tcrdir:
        object = websummaryTCR(
            summary_csv=f"{tcrdir}/metrics_summary.csv",
            samplename=samplename,
        )
        data_json_string.append(object.to_json(
                clonotypes_csv=f"{tcrdir}/clonotypes.csv",
                all_contig_annotations_csv=f"{tcrdir}/all_contig_annotations.csv",
                filtered_dir=f"{rnadir}/step3/filtered_feature_bc_matrix" if rnadir else None, 
                outputdir=f"{outdir}/filtered",
            )
        )
    if bcrdir:
        object = websummaryBCR(
            summary_csv=f"{bcrdir}/metrics_summary.csv",
            samplename=samplename,
        )
        data_json_string.append(object.to_json(
                clonotypes_csv=f"{bcrdir}/clonotypes.csv",
                all_contig_annotations_csv=f"{bcrdir}/all_contig_annotations.csv",
                filtered_dir=f"{rnadir}/step3/filtered_feature_bc_matrix" if rnadir else None, 
                outputdir=f"{outdir}/filtered",
            )
        )

    data = reduce(lambda x, y: {**x, **y}, [json.loads(j.replace("\\'", "'")) for j in data_json_string], {})
    if data:
        template_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'report'))
        env = Environment(loader=FileSystemLoader(template_dir))
        template = env.get_template('base.html')
        with open(os.path.join(outdir, f'{samplename}_report.html'), 'w') as fh:
            fh.write(template.render(websummary_json_data=json.dumps(data).replace("5'", "5\\'").replace("3'", "3\\'")))
