import os
from ..utils.helper import logger, get_file_path
from ..utils.wrappers import cmd_execute

def do_seurat(
    matrix:str, samplename:str, outdir:str, dims:int=15, minpct:float=0.1,
    logfc:float=0.25, rscript_path:str="Rscript", **kwargs
    ):

    logger.info("seurat started!")
    outdir1 = os.path.join(outdir, "step4")
    os.makedirs(outdir1, exist_ok=True)

    Rapp = os.path.abspath(os.path.join(get_file_path(__file__),"../utils/do_seurat.R"))
    matrix = os.path.abspath(matrix)
    outdir1 = os.path.abspath(outdir1)
    args = [
        rscript_path, Rapp, "--indir", matrix, "--name", samplename, "--outdir",
        outdir1, "--dims", dims, "--minpct", minpct, "--logfc", logfc
    ]
    cmd_execute(args, check=True)

    logger.info("seurat done!")
