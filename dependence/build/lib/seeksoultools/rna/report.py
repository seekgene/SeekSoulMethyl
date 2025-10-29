import os
import pandas as pd
from pathlib import Path
import lzstring
from ..utils.plotUtil import plotUtils
from ..utils.report.websummaryRNA import websummaryRNA
import shutil
data = {}

def line_color(p):
    _color = ['#dbdbf4',
             '#c5c5f5',
             '#afaff6',
             '#9999f7',
             '#8383f8',
             '#6d6df9',
             '#5858fb',
             '#4141fc',
             '#2b2bfd',
             '#1616fe']
    if p==1:
        return '#0000ff'
    elif p==0:
        return '#dddddd'
    else:
        return _color[int(p/0.1)]

def diff_table(f, N=20):
    df = (pd.read_table(f)
         .assign(cluster=lambda df: df['cluster'].map(lambda x: f'cluster{x}'))
    )
    # table_id, classes
    fold_change_title = 'avg_logFC'
    if 'avg_log2FC' in df.columns:
        fold_change_title = 'avg_log2FC'
    df2 = (df.drop(columns=['p_val',  'pct.1',  'pct.2'])
        .groupby('cluster')
        .apply(lambda x: x.sort_values(['p_val_adj', fold_change_title], ascending=[True, False]).head(N))
        .reset_index(drop=True)
        .pivot_table(values=[fold_change_title, 'p_val_adj'], index=['Ensembl', 'gene'], columns='cluster')
        .swaplevel(axis=1)
        .sort_index(axis=1, level=0, key=lambda x: x.str.replace('cluster', '').astype(int))
    )
    # index_names
    return (df2.to_html(table_id='marker_table', classes='display', na_rep='-', index_names=False)
        .replace('border="1"', ''))

def reduction_data(reduction_umi):
    df = pd.read_table(reduction_umi, index_col=0)
    df = df.drop(['orig.ident', 'nFeature_RNA', 'percent.mito'], axis=1)
    data['reduction'] = df.to_dict(orient='list')
    data['params'] = [ _ for _ in df.columns if _.startswith(('RNA_snn_res','seurat_clusters')) ]
    data['labs'] = df.columns[0:2].to_list()


def report(samplename, rna_wd, tcr_wd=None, bcr_wd=None, outdir=None, organism=None, cfg=None, **kwargs):
    import json
    from jinja2 import Environment, FileSystemLoader
    if not outdir:
        outdir = rna_wd
    outdir = os.path.join(outdir, 'Outs')
    os.makedirs(outdir, exist_ok=True)
    rna_wd = os.path.join(rna_wd, 'Analysis')
    summary_file = os.path.join(rna_wd, f'{samplename}_summary.json')
    assert os.path.exists(summary_file), f'{summary_file} not found!'
    with open(summary_file) as fh:
        summary = json.load(fh)

    sequencing_table = {}
    sequencing_table['Number of Reads'] = f'{summary["stat"]["total"]:,}'
    sequencing_table['Valid Barcodes'] = f'{summary["stat"]["valid"]/summary["stat"]["total"]:.2%}'
    sequencing_table['Sequencing Saturation'] = f'{summary["cells"]["Sequencing Saturation"]:.2%}'
    del summary["cells"]["Sequencing Saturation"]

    if 'no_anchor' in summary["stat"]:
        sequencing_table['Without Anchor'] = f'{summary["stat"]["no_anchor"]:,}'

    if 'trimmed' in summary["stat"]:
        sequencing_table['Trimmed'] =  f'{summary["stat"]["trimmed"]:,}'
    if 'too_short' in summary["stat"]:
        sequencing_table['Too Short'] =  f'{summary["stat"]["too_short"]:,}'
    b_total_base = sum([sum(v) for v in summary["barcode_q"].values()])
    b30_base = sum([sum(v[30:]) for v in summary["barcode_q"].values()])
    sequencing_table['Q30 Bases in Barcode'] = f'{b30_base/b_total_base:.2%}'
    u_total_base = sum([sum(v) for v in summary["umi_q"].values()])
    u30_base = sum([sum(v[30:]) for v in summary["umi_q"].values()])
    sequencing_table['Q30 Bases in UMI'] = f'{u30_base/u_total_base:.2%}'

    mapping_table = {k: f'{v:.2%}' for k, v in summary["mapping"].items()}

    cells_table = dict([(k, f'{v:,}') if isinstance(v, int) else (k,f'{v:.2%}') for k,v in summary["cells"].items()])
    
    sample_table = {
        'Name': samplename, 
        'Description': '',
        'Transcriptome': summary["reference"],
        'Chemistry': summary["stat"].get("chemistry", "custom"),
        'Include introns': summary["include_introns"],
        'Seeksoul tools Version': summary["__version__"]
    }
    
    reduction_xls = os.path.join(rna_wd, 'step4', 'tsne_umi.xls')
    assert os.path.exists(reduction_xls), f'{reduction_xls} not found!'
    reduction_data(reduction_xls)

    cells_gz = os.path.join(rna_wd, 'step3/filtered_feature_bc_matrix/barcodes.tsv.gz')
    barcodes_gz = os.path.join(rna_wd, 'step3/raw_feature_bc_matrix/barcodes.tsv.gz')
    mtx_gz = os.path.join(rna_wd, 'step3/raw_feature_bc_matrix/matrix.mtx.gz')
    barcod_rank_data = plotUtils.barcode_rank_plot_rna(cells_gz, barcodes_gz, mtx_gz)
    data["barcode_rank_data"]=barcod_rank_data

    f =  os.path.join(rna_wd, 'step4', 'FindAllMarkers.xls')
    marker_table = diff_table(f)

    header=('Samplename,Estimated_Number_of_Cells,Mean_Reads_per_Cell,Median_Genes_per_Cell,Number_of_Reads,'
            'Valid_Barcodes,Sequencing_Saturation,Reads_Mapped_Confidently_to_Genome,Fraction_Reads_in_Cells,'
            'Total_Genes_Detected,Median_UMI_Counts_per_Cell')

    summary_data = [
             samplename,
             cells_table['Estimated Number of Cells'],
             cells_table['Mean Reads per Cell'],
             cells_table['Median Genes per Cell'],
             sequencing_table['Number of Reads'],
             sequencing_table['Valid Barcodes'],
             sequencing_table['Sequencing Saturation'],
             mapping_table['Reads Mapped Confidently to Genome'],
             cells_table['Fraction Reads in Cells'],
             cells_table['Total Genes Detected'],
             cells_table['Median UMI Counts per Cell']
           ]

    with open(os.path.join(outdir, f'{samplename}_summary.csv'), 'w') as fh:
        fh.write(header + '\n')
        fh.write(','.join(str(_).replace(',', '') for _ in summary_data)+ '\n')

    matrix = Path(os.path.join(rna_wd, 'step3/filtered_feature_bc_matrix'))
    T_airr, tcr_json = None, None
    B_airr, bcr_json = None, None

    if tcr_wd:
        from ..vdj.report import report
        from ..vdj.utils import get_config
        outdir = Path(outdir)
        _ref, _fa, leader = get_config(organism, "TR", cfg)
        tcr_wd = Path(tcr_wd).expanduser().resolve()
        samplename_t = tcr_wd.name
        tcr_json = report(tcr_wd, samplename_t, chain="TR", matrix=matrix, outdir=outdir, leader_ref=leader, pa="vdj_t")
        T_metrics_summary = outdir / "vdj_t" / "metrics_summary.csv"
        T_metrics_summary = T_metrics_summary.rename(T_metrics_summary.parent/f"{samplename}_t_metrics_summary.csv")
        T_report_html= outdir / "vdj_t" / f"{samplename_t}_report.html"
        T_report_html = T_report_html.rename(T_report_html.parent/f"{samplename}_t_report.html")
        T_airr = outdir / "vdj_t" / f"{samplename_t}_airr_rearrangement.tsv"
        T_airr = T_airr.rename(T_airr.parent/f"{samplename}_t_airr_rearrangement.tsv")
        T_clonotypes = outdir / "vdj_t" / f"{samplename_t}_clonotypes.tsv"
        T_clonotypes = T_clonotypes.rename(T_clonotypes.parent/f"{samplename}_t_clonotypes.tsv")
        T_annotations = outdir / "vdj_t" / f"{samplename_t}_filtered_contig_annotations.tsv"
        T_annotations = T_annotations.rename(T_annotations.parent/f"{samplename}_t_filtered_contig_annotations.tsv")
        T_json = outdir / "vdj_t" / f"{samplename_t}_cell_barcodes.json"
        T_json = T_json.rename(T_json.parent/f"{samplename}_t_cell_barcodes.json")
        T_contig = outdir / "vdj_t" / f"{samplename_t}_filtered_contig.fasta"
        T_contig = T_contig.rename(T_contig.parent/f"{samplename}_t_filtered_contig.fasta")
        T_consensus = outdir / "vdj_t" / f"{samplename_t}_consensus.fasta"
        T_consensus = T_consensus.rename(T_consensus.parent/f"{samplename}_t_consensus.fasta")

    if bcr_wd:
        from ..vdj.report import report
        from ..vdj.utils import get_config
        outdir = Path(outdir)
        _ref, _fa, leader = get_config(organism, "IG", cfg)
        bcr_wd = Path(bcr_wd).expanduser().resolve()
        samplename_b = bcr_wd.name
        if not samplename_b:
            samplename_b = samplename
        bcr_json = report(bcr_wd, samplename_b, chain="IG", matrix=matrix, outdir=outdir, leader_ref=leader, pa="vdj_b")
        B_metrics_summary = outdir / "vdj_b" / "metrics_summary.csv"
        B_metrics_summary = B_metrics_summary.rename(B_metrics_summary.parent/f"{samplename}_b_metrics_summary.csv")
        B_report_html= outdir / "vdj_b" / f"{samplename_b}_report.html"
        B_report_html = B_report_html.rename(B_report_html.parent/f"{samplename}_b_report.html")
        B_airr = outdir / "vdj_b" / f"{samplename_b}_airr_rearrangement.tsv"
        B_airr = B_airr.rename(B_airr.parent/f"{samplename}_b_airr_rearrangement.tsv")
        B_clonotypes = outdir / "vdj_b" / f"{samplename_b}_clonotypes.tsv"
        B_clonotypes = B_clonotypes.rename(B_clonotypes.parent/f"{samplename}_b_clonotypes.tsv")
        B_annotations = outdir / "vdj_b" / f"{samplename_b}_filtered_contig_annotations.tsv"
        B_annotations = B_annotations.rename(B_annotations.parent/f"{samplename}_b_filtered_contig_annotations.tsv")
        B_json = outdir / "vdj_b" / f"{samplename_b}_cell_barcodes.json"
        B_json = B_json.rename(B_json.parent/f"{samplename}_b_cell_barcodes.json")
        B_contig = outdir / "vdj_b" / f"{samplename_b}_filtered_contig.fasta"
        B_contig = B_contig.rename(B_contig.parent/f"{samplename}_b_filtered_contig.fasta")
        B_consensus = outdir / "vdj_b" / f"{samplename_b}_consensus.fasta"
        B_consensus = B_consensus.rename(B_consensus.parent/f"{samplename}_b_consensus.fasta")
    object = websummaryRNA(
        summary_json=f"{rna_wd}/{samplename}_summary.json",
    )

    data_json = object.to_json(
        diff_table=f"{rna_wd}/step4/FindAllMarkers.xls",
        dim_table=f"{rna_wd}/step4/tsne_umi.xls",
        filtered_dir=f"{rna_wd}/step3/filtered_feature_bc_matrix",
        raw_dir=f"{rna_wd}/step3/raw_feature_bc_matrix",
        T_airr=T_airr,
        B_airr=B_airr,
    )

    if tcr_json:
        data_json.update(tcr_json)
    if bcr_json:
        data_json.update(bcr_json)

    template_dir_new = os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/report'))
    env = Environment(loader=FileSystemLoader(template_dir_new))
    template = env.get_template('base.html')

    with open(os.path.join(outdir, f'{samplename}_report.html'), 'w') as fh:
        # fh.write(template.render(websummary_json_data=data_json))
        fh.write(template.render(websummary_json_data=json.dumps(data_json).replace("5'", "5\\'").replace("3'", "3\\'")))
    shutil.copytree(f"{rna_wd}/step3/filtered_feature_bc_matrix", f"{outdir}/filtered_feature_bc_matrix", dirs_exist_ok=True)
    if bcr_wd or tcr_wd:
        ana_dir = outdir / 'Analysis'
        shutil.rmtree(ana_dir, ignore_errors=True)
