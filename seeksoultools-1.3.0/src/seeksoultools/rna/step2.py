import os
import json
from collections import defaultdict
from ..utils.helper import logger

def mapping_summary(STARLog, RnaSeqMetrics):
    summary = defaultdict()
    with open(STARLog, "r") as fh:
        for line in fh:
            if "Number of input reads" in line:
                summary["Number of input reads"] = int(
                    line.strip().split("\t")[-1])
            if "Uniquely mapped reads number" in line:
                summary["Uniquely mapped reads number"] = int(
                    line.strip().split("\t")[-1])
            if "Number of reads mapped to multiple loci" in line:
                summary["Number of reads mapped to multiple loci"] = int(
                    line.strip().split("\t")[-1])
            if "Number of reads mapped to too many loci" in line:
                summary["Number of reads mapped to too many loci"] = int(
                    line.strip().split("\t")[-1])

    with open(RnaSeqMetrics, "r") as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith("total alignments"):
                summary["total alignments"] = int(line.split()[-1].replace(",", ""))
            if line.startswith("reads aligned"):
                summary["reads aligned"] = int(line.split()[-1].replace(",", ""))
            if line.startswith("aligned to genes"):
                summary["aligned to genes"] = int(line.split()[-1].replace(",", ""))
            if line.startswith("no feature assigned"):
                summary["no feature assigned"] = int(line.split()[-1].replace(",", ""))
            if line.startswith("exonic"):
                summary["exonic"] = int(line.split()[-2].replace(",", ""))
            if line.startswith("intronic"):
                summary["intronic"] = int(line.split()[-2].replace(",", ""))
            if line.startswith("intergenic"):
                summary["intergenic"] = int(line.split()[-2].replace(",", ""))
                break
    return summary


def align(
    fq:list, genomeDir:str, gtf:str, samplename:str, outdir:str, region:str, sc5p:bool,
    core:int=4, star_path:str="STAR", **kwargs):

    if ("steps" not in kwargs) or (not kwargs["steps"]):
        kwargs["steps"] = ["STAR", "SortByPos", "qualimap",  "FeatureCounts", "SortByName"]

    basedir = os.path.join(outdir, "step2")
    STAR_dir = os.path.join(basedir, "STAR")
    os.makedirs(STAR_dir, exist_ok=True)
    prefix = os.path.join(STAR_dir, samplename + "_")

    if "STAR" not in kwargs["steps"]:
        logger.info("STAR skiped!")
    else:
        from ..utils.wrappers import STAR_wrapper
        logger.info("STAR started!")
        bam, STARLog = STAR_wrapper(
            fq=fq,
            core=core,
            genomeDir=genomeDir,
            prefix=prefix,
            star_path=star_path
        )
        logger.info("STAR done!")
    bam, STARLog =f"{prefix}Aligned.out.bam", f"{prefix}Log.final.out"

    # sort by pos
    if "SortByPos" not in kwargs["steps"]:
        logger.info("SortByPos skiped!")
    else:
        logger.info("SortByPos started!")
        from ..utils.wrappers import samtools_sort_wrapper
        bam = samtools_sort_wrapper(
            bam,
            f"{prefix}SortedByCoordinate.bam",
            core=core,
            clean=True,
        )
        logger.info("SortByPos done!")
    bam = f"{prefix}SortedByCoordinate.bam"

    logger.info("run_qualimap started!")
    from ..utils.wrappers import qualimap_wrapper
    RnaSeqMetrics = qualimap_wrapper(bam=bam, gtf=gtf, outdir=STAR_dir, SC5P=sc5p)
    logger.info("run_qualimap done!")

    with open(os.path.join(outdir, f"{samplename}_summary.json")) as fh:
        refpath=os.path.dirname(genomeDir.rstrip("/"))
        reffile=os.path.join(refpath, "reference.json")
        if os.path.exists(reffile):
                with open(reffile) as refjson:
                    refj=json.load(refjson)
                    genome=refj["genomes"][0]
        else:
                genome=refpath.split("/")[-1]
        summary = json.load(fh)
        summary["reference"] = genome
        Total = summary["stat"]["total"]
    
    if region=="exon":
        summary["include_introns"] = "False"
    else:
        summary["include_introns"] = "True"

    summary_tmp = defaultdict()
    tmp = mapping_summary(STARLog, RnaSeqMetrics)
    Total = tmp["Number of input reads"]
    mapped_genome_ratio = (tmp["Uniquely mapped reads number"] +
                           tmp["Number of reads mapped to multiple loci"] +
                           tmp["Number of reads mapped to too many loci"])/Total
    summary_tmp["Reads Mapped to Genome"] = mapped_genome_ratio

    mapped_confident_ratio = (tmp["aligned to genes"] + tmp["no feature assigned"]) / Total
    summary_tmp["Reads Mapped Confidently to Genome"] = mapped_confident_ratio

    mapped_intergenic_ratio = tmp["intergenic"] / Total
    summary_tmp["Reads Mapped to Intergenic Regions"] = mapped_intergenic_ratio

    mapped_intronic_ratio = tmp["intronic"] / Total
    summary_tmp["Reads Mapped to Intronic Regions"] = mapped_intronic_ratio

    mapped_exonic_ratio = tmp["exonic"] / Total
    summary_tmp["Reads Mapped to Exonic Regions"] = mapped_exonic_ratio

    with open(os.path.join(outdir, f"{samplename}_summary.json"), "w") as fh:
        summary["mapping"] = summary_tmp
        json.dump(summary, fh, indent=4)


    featureCounts_dir = os.path.join(basedir, "featureCounts")
    # FeatureCounts
    if "FeatureCounts" not in kwargs["steps"]:
        logger.info("run_featureCounts done!")
    else:
        os.makedirs(featureCounts_dir, exist_ok=True)
        logger.info("run_featureCounts started!")
        from ..utils.wrappers import featureCounts_wrapper
        bam = featureCounts_wrapper(
            bam=bam,
            samplename=samplename,
            outdir=featureCounts_dir,
            gtf=gtf,
            region=region,
            SC5P=sc5p,
            core=core
        )
        logger.info("run_featureCounts done!")
    bam = os.path.join(featureCounts_dir, f"{samplename}_SortedByCoordinate.bam.featureCounts.bam")
    # sort by name
    if "SortByName" not in kwargs["steps"]:
        logger.info("SortByName skiped!")
    else:
        logger.info("SortByName started!")
        from ..utils.wrappers import samtools_sort_wrapper
        bam = samtools_sort_wrapper(
            bam,
            os.path.join(featureCounts_dir, f"{samplename}_SortedByName.bam"),
            core=core,
            byname=True,
            clean=True,
        )
        logger.info("SortByName done!")
