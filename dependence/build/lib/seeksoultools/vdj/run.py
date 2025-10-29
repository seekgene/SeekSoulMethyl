from pathlib import Path
from ..utils.helper import logger
from .utils import (
     subset_reads, 
     count_barcode,
     readFormat2rs,
     summary_simpleqc,
     summary_enrichment_qc,
     summary_fastq_extractor,
     
)
def run_trust4(
    samplename:str, fq1:list[Path], fq2:list[Path], fa:Path, ref:Path, leader:Path,
    organism:str|None, wl:Path,  chain:str,  wd:Path, core:int=8, 
    read_pair:bool=False, reads_limit:int=80000, matrix:Path|None=None, rna_wd:Path|None=None,
    keep_tmp:bool=False
):
    summary_json = wd / "Analysis" / "summary.json"
    readFormat = "bc:0:16,um:17:28"
    # fastq qc
    # "Number of reads", "Number of valid reads", "Valid barcodes"
    # "Q30 barcodes", "Q30 UMI", "Q30 R2 read"
    qc_wd = wd / "Analysis" /"qc"
    qc_wd.mkdir(exist_ok=True, parents=True)

    from ..utils.wrappers import simpleqc_wrapper
    logger.info("simple qc")
    simpleqc_wrapper(fq1=fq1, fq2=fq2, wl=wl, rs=readFormat2rs(readFormat), outprefix=f"{qc_wd}/{samplename}_", threads=core)
    summary_simpleqc(summary_file=summary_json, wd=wd, samplename=samplename, readFormat=readFormat)
    fq1 = qc_wd/f"{samplename}_r1.fq.gz"
    fq2 = qc_wd/f"{samplename}_r2.fq.gz"

    # enrichment qc
    # "Reads mapped to any V(D)J gene" ...
    from .enrichment import enrichment_qc
    logger.info("enrichment qc")
    enrichment_qc(
        fq=[fq2,],
        chain=chain,
        fa=fa,
        samplename=samplename,
        outdir=qc_wd,
        core=core,
        star_path="STAR"
    )
    summary_enrichment_qc(summary_file=summary_json, wd=wd, samplename=samplename, chain=chain)
    # trust4 fastq extractor
    # "Number of candidate reads" ...
    turst_wd = wd /"Analysis"/"trust4"
    turst_wd.mkdir(exist_ok=True, parents=True)
    if not read_pair:
        readFormat = ",".join([e for e in readFormat.split(",") if not e.startswith("r1:")])

    from ..utils.wrappers import fastq_extractor_wrapper

    logger.info("fastq extractor")
    fastq_extractor_wrapper(samplename, fq1, fq2, fa, wl, readFormat, turst_wd, core, read_pair)
    
    # subset reads
    _candidate_reads, count_dict = count_barcode(
        turst_wd/f"{samplename}_preassemble_bc.fa",
        turst_wd/f"{samplename}_preassemble_umi.fa",
        turst_wd/f"{samplename}_bc_count.tsv",
    )
    count_dict = {k:0 for k, v in count_dict.items() if v > reads_limit}
    bc_fa = turst_wd/f"{samplename}_preassemble_bc.fa"
    umi_fa = turst_wd/f"{samplename}_preassemble_umi.fa"
    fq_fa = turst_wd/f"{samplename}_preassemble.fq"
    pfq2 = turst_wd/f"{samplename}_preassemble_2.fq"
    pfq1 = turst_wd/f"{samplename}_preassemble_1.fq"
    if len(count_dict) > 0 and reads_limit > 0:
        logger.info(f"filter the barcodes to include only those with up to {reads_limit:,} reads.")
        if not read_pair:
            pfq2 = turst_wd/f"{samplename}_preassemble.fq"
            subset_reads(
                count_dict=count_dict, bc_fa=bc_fa, umi_fa=umi_fa,
                r2_fq=pfq2, reads_limit=reads_limit
            )
        else:
            subset_reads(
                count_dict=count_dict, bc_fa=bc_fa, umi_fa=umi_fa,
                r2_fq=pfq2, r1_fq=pfq1, reads_limit=reads_limit
            )
        (turst_wd/f"{samplename}_bc_count.tsv").replace(f"{turst_wd}/_{samplename}_bc_count.tsv")
        _candidate_reads, _count_dict = count_barcode(
            turst_wd/f"{samplename}_preassemble_bc.fa",
            turst_wd/f"{samplename}_preassemble_umi.fa",
            turst_wd/f"{samplename}_bc_count.tsv",
        )
    summary_fastq_extractor(summary_file=summary_json, wd=wd, samplename=samplename)

    ## filter umis
    tmpfile = turst_wd/f"{samplename}_filter.detail.xls"
    from ..utils.wrappers import filter_umis
    logger.info("filter umis")
    bc_out = turst_wd/f"{samplename}_toassemble_bc.fa"
    umi_out = turst_wd/f"{samplename}_toassemble_umi.fa"
    fq_out = turst_wd/f"{samplename}_toassemble.fq"
    fq1_out = turst_wd/f"{samplename}_toassemble_1.fq"
    fq2_out = turst_wd/f"{samplename}_toassemble_2.fq"
    assign_f = turst_wd/f"{samplename}_assign.out"
    
    if  not read_pair:
        filter_umis(bc_fa, umi_fa, fq_fa, bc_out, umi_out, fq_out,  tmpfile, fq1=None, fq1_out=None)
    else:
        filter_umis(bc_fa, umi_fa, pfq2,   bc_out, umi_out, fq2_out, tmpfile, pfq1, fq1_out)
    
    logger.info("filter umis done")
    # trust4
    from ..utils.wrappers import trust4_wrapper
    logger.info("run trust4")
    trust4_wrapper(samplename, fq1, fq2, fa, ref, wl, readFormat, turst_wd, core, read_pair, stage=1)
    # clean tmp
    logger.info("cleaning temp files")
    raw_out = turst_wd / f"{samplename}_raw.out"
    reads = turst_wd / f"{samplename}_assembled_reads.fa"
    final_out = turst_wd / f"{samplename}_final.out"
    cleanfile = [fq_fa, pfq1, pfq2, bc_fa, umi_fa, raw_out, reads]
    for files in cleanfile:
        if files.exists():
            files.unlink()

    gzfile = [fq_out, fq1_out, fq2_out, bc_out, umi_out, assign_f]
    for files in gzfile:
        if files.exists():
            if keep_tmp:
                from ..utils.wrappers import cmd_execute
                cmd_execute(f"gzip -f {files}", check=True)
            else:
                files.unlink()
    # report 
    from .report import report
    report(wd, samplename, chain, leader, matrix, rna_wd)
