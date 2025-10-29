import os
import json
import numpy as np
from collections import defaultdict
from ..utils.helper import logger, get_utils_path
from ..utils.wrappers import cmd_execute
from ..utils import countUtil 

__srcdir = get_utils_path(__file__)

def count(bam, outdir, gtf, umi_correct_method, **kwargs):
    basedir = os.path.join(outdir, "step3")
    os.makedirs(basedir, exist_ok=True)
    countUtil.count(bam, basedir, gtf, umi_correct_method, **kwargs)


def cell_calling(raw_matrix, outdir, samplename, gtf, expectNum=3000, **kwargs):

    basedir = os.path.join(outdir, "step3")
    os.makedirs(basedir, exist_ok=True)

    detail_file = os.path.join(outdir, "step3", "detail.xls")
    counts_file = os.path.join(outdir, "step3", "counts.xls")
    umi_file = os.path.join(outdir, "step3", "umi.xls")
    
    filtered_matrix = os.path.join(basedir, "filtered_feature_bc_matrix")
    Rapp = os.path.join(__srcdir, "cell_identify.R")
    if kwargs["forceCell"] != None:
        args = ["Rscript", Rapp, "-i", raw_matrix, "-o", filtered_matrix, "-f", kwargs["forceCell"]]
        logger.info("force cell.")
    else:
        args = ["Rscript", Rapp, "-i", raw_matrix, "-o", filtered_matrix, "-e", expectNum, "-p", 0.01]

    logger.info("call cell started!")
    cmd_execute(args, check=False)
    logger.info("call cell done!")

    filterd_barcodes_file = os.path.join(filtered_matrix, "barcodes.tsv.gz")
    filterd_features_file = os.path.join(filtered_matrix, "features.tsv.gz")

    logger.info("calculate metrics started!")
    summary_tmp, downsample = countUtil.calculate_metrics(
        counts_file,
        detail_file,
        filterd_barcodes_file,
        filterd_features_file,
        gtf,
        basedir
    )

    with open(os.path.join(outdir, samplename+"_summary.json"), "r") as fh:
        summary = json.load(fh)
        # Total = int(summary["Sequencing"]["Number of Reads"].replace(",", ""))
        Total = summary["stat"]["total"]
    with open(os.path.join(outdir, samplename+"_summary.json"), "w") as fh:
        estimated_cell_num = summary_tmp["Estimated Number of Cells"]
        mean_reads_per_cell = Total / estimated_cell_num
        summary_tmp["Mean Reads per Cell"] = int(mean_reads_per_cell)
        del summary_tmp["Median lnc Genes per Cell"]
        summary["cells"] = summary_tmp
        downsample["Reads"] = [
            int(summary_tmp["Mean Reads per Cell"] * p)
            for p in downsample["percentage"]
        ]
        summary["downsample"] = downsample
        json.dump(
            summary,
            fh,
            indent=4,
            default=lambda o: int(o) if isinstance(o, np.int64) else o
        )

    logger.info("cell calling done!")

    return
