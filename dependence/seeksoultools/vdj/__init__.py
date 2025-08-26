from pathlib import Path
import click
from ..utils.helper import str2path
from .utils import get_config, get_white_list
from typing import List, Optional
from ..utils.helper import check_all

@click.group(help="assembl and annotate vdj sequence.")
@click.pass_obj
def vdj(obj):
    pass
@vdj.command()
@click.option("--samplename", required=True, type=str, help="Sample name.")
@click.option("--fq1", required=True, multiple=True, type=click.Path(exists=True), help="Path(s) to read1 file(s).")
@click.option("--fq2", required=True, multiple=True, type=click.Path(exists=True), help="Path(s) to read2 file(s).")
@click.option("--organism", type=click.Choice(["human", "mouse",]), required=False, help="Organism type.")
@click.option("--cfg", type=click.Path(exists=True), callback=str2path, required=False, help="Path to ref config.")
@click.option("--chain", type=click.Choice(["TR", "IG"]), help="Chain type, TR or IG.")
@click.option("--outdir", required=True, type=click.Path(), callback=str2path, help="Working directory.")
@click.option("--read_pair", default=False, is_flag=True, help="Use r1 and r2. default: only r2")
@click.option("--core", default=8, type=int, help="Number of threads.")
@click.option("--matrix", "matrix", type=click.Path(exists=True), help="RNA matrix directory.")
@click.option("--reads_limit", hidden=True, default=80000, type=int, show_default=True, help="Read number limit for each barcode, negative means no limit.")
@click.option("--rna_wd", type=click.Path(exists=True), callback=str2path, help="RNA working directory.")
@click.option("--keep_tmp", is_flag=True, default=False, help="Keep trust4 temporary files.")
def run(
    samplename: str,
    fq1: List[Path],
    fq2: List[Path],
    organism: str,
    chain: str,
    outdir: Path,
    read_pair: bool = False,
    cfg: Optional[Path] = None,
    matrix: Optional[Path] = None,
    rna_wd: Optional[Path] = None,
    core: int = 8,
    reads_limit: int = 80000,
    keep_tmp: bool = False
):
    from ..utils.helper import check_cores
    core = check_cores(core)
    wd = outdir / samplename 
    check_all()
    ref, fa, leader = get_config(organism, chain, cfg)
    white_list = get_white_list()
    
    from .run import run_trust4
    run_trust4(
        samplename, fq1, fq2, fa, ref, leader, organism, white_list,
        chain, wd, core, read_pair, reads_limit, matrix, rna_wd, keep_tmp
    )

@vdj.command()
@click.option("--samplename", required=True, type=str, help="Sample name.")
@click.option("--chain", type=click.Choice(["TR", "IG"]), help="Chain type, TR or IG.")
@click.option("--organism", type=click.Choice(["human", "mouse",]), required=False, help="Organism type.")
@click.option("--cfg", type=click.Path(exists=True, path_type=Path), required=False, help="Path to ref config.")
@click.option("--outdir", required=True, type=click.Path(), callback=str2path, help="Working directory.")
@click.option("--matrix", type=click.Path(exists=True), help="RNA matrix directory.")
@click.option("--rna_wd", type=click.Path(exists=True),  callback=str2path, help="RNA working directory.")
def report(
    outdir: Path,
    samplename: str,
    chain: str,
    organism: str|None,
    cfg: Path|None,
    matrix: Optional[Path] = None,
    rna_wd: Optional[Path] = None,
):
    from .report import report
    _ref, _fa, leader = get_config(organism, chain, cfg)
    report(outdir, samplename, chain, leader, matrix, rna_wd, outdir)

