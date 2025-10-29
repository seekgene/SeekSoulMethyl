import sys
import os
from pathlib import Path
import shutil
from subprocess import run
import subprocess
import random
import dnaio
import pandas as pd
import psutil
from collections import defaultdict
from .helper import logger
from .helper import hamming_distance
from .countUtil import UMIClusterer,create_map_to_correct_umi,encode, decode
from .reassign import main
##from .barcode import mapping_report
def cmd_execute(
    args, check:bool=False, text:bool=True,
    capture_output:bool=True, env:os.environ=None
):

    if isinstance(args, list):
        args = [str(_) for _ in args]
        logger.info(" ".join(args))
        _call = run(args, check=False, text=text, capture_output=capture_output, env=env)
    elif isinstance(args, str):
        logger.info(args)
        _call = run(args, shell=True, check=False, text=text, capture_output=capture_output, env=env)
    
    if check:
        if _call.returncode!=0:
            logger.info(f"stderr: {_call.stderr}")
            logger.info(f"stdout: {_call.stdout}")
            sys.exit(1)
    return _call

def STAR_wrapper(
    fq:list, genomeDir:str, prefix:str, scoremin=None, matchnmin=None,
    core:int=4, star_path:str="STAR", readMapNumber:int=-1,
    ):
    """wrapper for STAR"""
    args = [
        star_path, '--runThreadN', core, '--limitOutSJcollapsed', 5000000, '--readMapNumber', readMapNumber,
        '--genomeDir', genomeDir, '--readFilesCommand', 'zcat', '--outFileNamePrefix',
        prefix, '--outSAMtype', 'BAM', 'Unsorted', 
    ]
    if scoremin:
        args += ['--outFilterScoreMinOverLread', scoremin,]
    if matchnmin:
        args += [ '--outFilterMatchNminOverLread', matchnmin,]

    args = args + ['--readFilesIn', ] + list(fq)
    _ = cmd_execute(args)
    STARLog = f'{prefix}Log.final.out'
    from .barcode import mapping_report
    mapping_report(STARLog)
    return f'{prefix}Aligned.out.bam', f'{prefix}Log.final.out'


def samtools_index_wrapper(inbam:str, core:int=4, samtools_path:str="samtools"):
    args = [samtools_path, "index", "-@", core, inbam]
    _call = cmd_execute(args, check=True)
    return 

def samtools_sort_wrapper(
    inbam:str, outbam:str, byname:bool=False, clean:bool=False,
    core:int=4, samtools_path:str="samtools"
)->str:
    """wrapper for samtools sort"""
    args = [
        samtools_path, "sort", "-O", "BAM", "-@",  core, "-o", outbam, inbam
    ]
    if byname:
        args.insert(2, '-n')

    _call = cmd_execute(args, check=True)

    if _call.returncode == 0:
        # index
        if not byname:
            samtools_index_wrapper(outbam,samtools_path=samtools_path)
        if clean:
            os.remove(inbam)
    return outbam

def unzip_wrapper(gzfile:str, outfile:str):
    args = [
        "gzip", "-dc", gzfile, ">", outfile
    ]
    cmd_execute(args, check=True)
    return outfile

def qualimap_wrapper(bam:str, gtf:str, outdir:str, SC5P:bool=None, qualimap_path:str="qualimap", **kwargs):
    """ -pe,--paired? -s?"""
    strand = {
        'f': 'strand-specific-forward',
        'r': 'strand-specific-reverse',
        'non': 'non-strand-specific'
    }
    s = "non"
    # 5'
    if isinstance(SC5P, bool):
        if SC5P:
            s="r"
        else:
            s="f"
    if '-pe' in kwargs :
       s = 'non'
    # gtf.gz?
    if gtf.endswith('.gz'):
        gtf_file = os.path.join(outdir, 'tmp.gtf')
        unzip_wrapper(gtf, gtf_file)
    else:
        gtf_file = gtf
    available_memory = psutil.virtual_memory().available / (1024 * 1024 * 1024)
    available_memory = int(available_memory)
    args = [
        qualimap_path, 'rnaseq', '-outformat', 'PDF', '-outdir', outdir, '-bam',
        bam, '-gtf', gtf, '-p', strand[s], f'--java-mem-size={available_memory}G'
    ]

    for k, v in kwargs.items():
        args += [f"{k}", f"{v}"]

    my_env = os.environ.copy()
    if 'DISPLAY' in my_env:
        del my_env['DISPLAY']
    
    ##cmd_execute(args, check=True, env=my_env)
    try:
        result = subprocess.run(
            args,
            env=my_env,
            check=True,
            capture_output=True,
            text=True
        )
        
        if "WARNING: out of memory" in result.stderr:
            logger.error(
                "out of memory!\n"
                "Qualimap allows to set RAM size using special argument: --java-mem-size"
            )
            sys.exit(1)
            
    except subprocess.CalledProcessError as e:
        if "WARNING: out of memory" in e.stderr:
            logger.error(
                "out of memory!\n"
                "Qualimap allows to set RAM size using special argument: --java-mem-size"
            )
            sys.exit(1)
        else:
            raise e



    return os.path.join(outdir, "rnaseq_qc_results.txt")
def check_gtf_region(gtf, region):
    region_set = set()
    logger.info("check gtf feature.")
    with open(gtf) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            tmp = line.strip().split("\t")
            region_set.add(tmp[2])
            if tmp[2] == region:
                return region
    if region == "transcript" and "gene" in region_set:
        logger.warning("No transcript feature in gtf, use gene feature instead.")
        return "gene"
    logger.warning(f"No {region} feature in gtf")
    os.exit(1)

def featureCounts_wrapper(
    bam:str, gtf:str, samplename:str, outdir:str, region:str, SC5P:bool=None,
    core:int=4, featureCounts_path="featureCounts", **kwargs
    ):
    outcounts = os.path.join(outdir, 'counts.txt')
    region = check_gtf_region(gtf, region)
    s = 0
    if isinstance(SC5P, bool):
        if SC5P:
            s = "2"
        else:
            s = "1"
    if '-p' in kwargs:
       s = 0
    args = [
        featureCounts_path, "-T", core, "-t", region , "-s", s, "-M", "-O", "-g", "gene_id",
        "--fracOverlap", 0.5, "-a", gtf, "-o", outcounts, 
    ]
    
    for k, v in kwargs.items():
        args += [f"{k}", f"{v}"]
    
    args = args + ["-R", "BAM", bam]
    try:
        cmd_execute(args, check=True)
    except Exception as e:
        logger.error('f{e}')
        sys.exit(1)
    return os.path.join(outdir, f"{outdir}/{os.path.basename(bam)}.featureCounts.bam")

def bowtie2_wrapper(
    fq:str, ref:str, bam:str, core:int=1, local_mode=True,
    bowtie2_path:str="bowtie2", samtools_path:str="samtools"
)->str:
    local_option = ""
    if local_mode:
        local_option = "--local"
    args = (
        f"{bowtie2_path} -p {core} -x {ref} -U {fq} {local_option}|"
        f"samtools view -b > {bam}"
    )
    cmd_execute(args, check=True)
    return bam

def bowtie2_build_wrapper(
    ref:str,
    outdir:str,
    bowtie2_path:str="bowtie2",
)->str:
    os.makedirs(outdir, exist_ok=True)
    if os.path.dirname(os.path.abspath(ref)) != os.path.dirname(os.path.abspath(outdir)):
        ref_l = shutil.copy(ref, outdir)
    args = [f"{bowtie2_path}-build", ref_l, ref_l]
    cmd_execute(args, check=True)
    return ref_l

def igblastn_wrapper(
    input:str,
    output:str,
    organism:str,
    chain:str,
    core:int=4,
    igblastn_path:str="igblastn",
    auxiliary_data:str="",
    internal_data:str="",
):
    outdir = os.path.dirname(output)

    cmd = f"cd {outdir}; ln -s {internal_data} ./{organism}; "
    if chain=="TR":
        cmd += (
            f"{igblastn_path} -query {input} -auxiliary_data {auxiliary_data} -outfmt 19 -ig_seqtype TCR "
            f"-germline_db_J {organism}/TR_J-REGION -germline_db_D {organism}/TR_D-REGION "
            f"-germline_db_V '{organism}/TR_L-REGION+V-REGION' -c_region_db {organism}/TR_C-REGION "
            f"-num_threads {core} -organism {organism} -out {output}"
        )
    elif chain=="IG":
        cmd += (
            f"{igblastn_path} -query {input} -auxiliary_data {auxiliary_data} -outfmt 19 -ig_seqtype Ig "
            f"-germline_db_J {organism}/IG_J-REGION -germline_db_D {organism}/IG_D-REGION "
            f"-germline_db_V {organism}/IG_L-REGION+V-REGION -c_region_db {organism}/IG_C-REGION  "
            f"-num_threads {core} -organism {organism} -out {output}"
        )
    else:
        pass

    cmd_execute(cmd, check=True)
    return output

def picard_DownsampleSam_wrapper(bam, fra, downsambam):
    cmd = "picard DownsampleSam I=%s O=%s P=%s" %(bam, downsambam, fra)
    cmd_execute(cmd)
    samtools_index_wrapper(downsambam)
    return downsambam

def geneBody_coverage_wrapper(downsambam, downdir, samplename, gtf):
    resultbed = os.path.join(downdir, samplename+".bed")
    reductionbed = os.path.join(downdir, samplename+".reduction.bed")
    pel=os.path.join(os.path.abspath(os.path.dirname(__file__)),'gtf2bed.pl')
    with open(resultbed, 'wb', 0) as bedfile:
        subprocess.call(['perl', pel, gtf], stdout=bedfile)
    with open(resultbed) as infile:
        with open(reductionbed,'w') as out:
            allbed=infile.readlines()
            lines = random.sample(allbed, min(len(allbed), 20000))
            for line in lines:
                out.write(line)
    cmd = f"cd {downdir}; geneBody_coverage.py -r {downdir}/{samplename}.reduction.bed -i {downsambam} -o {downdir}/{samplename}"
    cmd_execute(cmd)
    with open(f'{downdir}/{samplename}.geneBodyCoverage.txt','r') as infile:
        infile.readline()
        txt=infile.readline().strip().split('\t')[1:]
        numberlist=map(float,txt)
        all=sum(numberlist)
        mosttxt=txt[24:76]
        mostnumberlist=map(float,mosttxt)
        most=sum(mostnumberlist)
    return float(most)/float(all)

def fastq_extractor_wrapper(
    samplename:str, r1:str, r2:str, fa:Path, wl:Path, 
    readFormat:str, wd:Path, threads:int=10, read_pair=False,
    stage:int=0
):
    cmd = f"cd {wd}; "
    if  read_pair:
        cmd += (
            f"fastq-extractor -f {fa} -1 {r1} -2 {r2} -t {threads} "
            f"--readFormat {readFormat} --barcode {r1} --UMI {r1} "
            f"-o {samplename}_preassemble  --barcodeWhitelist {wl} "
            f"1>{wd}/_stdout.log 2>{wd}/_stderr.log "
        )
    else:
        cmd += (
            f"fastq-extractor -f {fa} -u {r2} -t {threads} "
            f"--readFormat {readFormat} --barcode {r1} --UMI {r1} "
            f"-o {samplename}_preassemble  --barcodeWhitelist {wl} "
            f"1>{wd}/_stdout.log 2>{wd}/_stderr.log "
        )
    cmd_execute(cmd)
def calculate_n50(umi_counts: dict):
    counts = sorted(umi_counts.values(), reverse=True)
    total_sum = sum(counts)
    threshold = total_sum * 0.5
    current_sum = 0
    for count in counts:
        current_sum += count
        if current_sum >= threshold:
            return count
    return 0
def filter_umis(bc_fa, umi_fa, fq_fa, bc_out, umi_out, fq_out, tmpfile, fq1=None, fq1_out=None):
    count_dict = defaultdict(lambda: defaultdict(int))
    n50_cutoff = 0.1
    with dnaio.open(bc_fa, fileformat="fasta") as fh_b, \
        dnaio.open(umi_fa, fileformat="fasta") as fh_u:
        for b1, u1 in zip(fh_b, fh_u):
            barcode = b1.name.split('_')[0]
            umi = u1.name.split('_')[1]
            #count_dict[barcode][umi] += 1
            count_dict[b1.sequence][u1.sequence] += 1
        valid_barcodes = set()
        valid_umi_counts = {}
        correction_info = {}
        umi_clusterer = UMIClusterer("adjacency")
        with open(tmpfile, 'w') as detail_fh:
            detail_fh.write("barcode\toriginal_umi\toriginal_count\tcorrected_umi\tcorrected_count\n")
            for barcode, umi_counts in count_dict.items():
                n50 = calculate_n50(umi_counts)
                threshold = n50_cutoff * n50
                filtered_counts = {encode(umi): count for umi, count in umi_counts.items() if count > threshold}
                
                if not filtered_counts:
                   continue
                #valid_barcodes.add(barcode)
                valid_umi_counts[barcode] = filtered_counts
                final_umis = umi_clusterer(filtered_counts, threshold=1)
                corrected_dict = create_map_to_correct_umi(final_umis)
                for ori_umi, new_umi in corrected_dict.items():
                    if ori_umi != new_umi:
                        filtered_counts[new_umi] += filtered_counts[ori_umi]
                        detail_fh.write(
                            f"{barcode}\t{decode(ori_umi)}\t{filtered_counts[ori_umi]}\t"
                            f"{decode(new_umi)}\t{filtered_counts[new_umi]}\n"
                            )
                        correction_info[(barcode, decode(ori_umi))] = decode(new_umi)
                        del filtered_counts[ori_umi]
                valid_umi_counts[barcode] = filtered_counts
    if not fq1 and not fq1_out:
        with dnaio.open(bc_fa, fileformat="fasta") as fh_b, \
             dnaio.open(umi_fa, fileformat="fasta") as fh_u, \
             dnaio.open(fq_fa, fileformat="fastq") as fh_fq, \
             dnaio.open(bc_out, fileformat="fasta", mode='w') as out_b, \
             dnaio.open(umi_out, fileformat="fasta", mode='w') as out_u, \
             dnaio.open(fq_out, fileformat="fastq", mode='w') as out_fq:
            for b1, u1, fq1 in zip(fh_b, fh_u, fh_fq):
                b1_name = b1.name
                barcode = b1.sequence
                original_umi = u1.sequence
                if barcode not in valid_umi_counts:
                    continue
                if encode(original_umi) not in valid_umi_counts[barcode]:
                    continue
                if (barcode, original_umi) in correction_info:
                    corrected_umi = correction_info[(barcode, original_umi)]
                    name_parts = b1.name.split('_')
                    name_parts[1] = corrected_umi
                    new_name = '_'.join(name_parts)
                    b1.name = new_name
                    u1.name = new_name
                    u1.sequence = corrected_umi
                    fq1.name = new_name
                out_b.write(b1)
                out_u.write(u1)
                out_fq.write(fq1)
    else:
        with dnaio.open(bc_fa, fileformat="fasta") as fh_b, \
             dnaio.open(umi_fa, fileformat="fasta") as fh_u, \
             dnaio.open(fq_fa, fileformat="fastq") as fh_fq, \
             dnaio.open(fq1, fileformat="fastq") as fh_fq1, \
             dnaio.open(bc_out, fileformat="fasta", mode='w') as out_b, \
             dnaio.open(umi_out, fileformat="fasta", mode='w') as out_u, \
             dnaio.open(fq_out, fileformat="fastq", mode='w') as out_fq, \
             dnaio.open(fq1_out, fileformat="fastq", mode='w') as out_fq1:
            for b1, u1, fq1, fq2 in zip(fh_b, fh_u, fh_fq, fh_fq1):
                b1_name = b1.name
                barcode = b1.sequence
                original_umi = u1.sequence
                if barcode not in valid_umi_counts:
                    continue
                if encode(original_umi) not in valid_umi_counts[barcode]:
                    continue
                if (barcode, original_umi) in correction_info:
                    corrected_umi = correction_info[(barcode, original_umi)]
                    name_parts = b1.name.split('_')
                    name_parts[1] = corrected_umi
                    new_name = '_'.join(name_parts)
                    b1.name = new_name
                    u1.name = new_name
                    u1.sequence = corrected_umi
                    fq1.name = new_name
                    fq2.name = new_name
                out_b.write(b1)
                out_u.write(u1)
                out_fq.write(fq1)
                out_fq1.write(fq2)


def trust4_wrapper2(
    samplename:str, r1:str, r2:str, fa:Path, ref:Path, wl:Path, 
    readFormat:str, wd:Path, threads:int=10, read_pair=False,
    stage:int=0
):
    cmd = f"cd {wd}; "
    if read_pair:
        cmd += (
            f"run-trust4 -f {fa} --ref {ref} -1 {r1} -2 {r2} "
            f"--readFormat {readFormat} --barcode {r1} --UMI {r1} "
            f"--od {wd} -o {samplename} --barcodeWhitelist {wl} "
            f"--outputReadAssignment -t {threads} --stage {stage} "
            f"1>>{wd}/_stdout.log 2>>{wd}/_stderr.log "
        )
    else:
        cmd += (
            f"run-trust4 -f {fa} --ref {ref} -u {r2} "
            f"--readFormat {readFormat} --barcode {r1} --UMI {r1} "
            f"--od {wd} -o {samplename}  --barcodeWhitelist {wl} "
            f"--outputReadAssignment -t {threads} --stage {stage} "
            f"1>>{wd}/_stdout.log 2>>{wd}/_stderr.log "
        )
    cmd_execute(cmd)

def trust4_wrapper(
    samplename:str, r1:str, r2:str, fa:Path, ref:Path, wl:Path,
    readFormat:str, wd:Path, threads:int=10, read_pair=False,
    stage:int=0
):
    cmd = f"cd {wd}; "
    if read_pair:
        cmd += (
            f"run-trust4 -f {fa} --ref {ref} -1 {r1} -2 {r2} "
            f"--readFormat {readFormat} --barcode {r1} --UMI {r1} "
            f"--od {wd} -o {samplename} --barcodeWhitelist {wl} "
            f"--outputReadAssignment -t {threads} --stage {stage} "
            f"1>>{wd}/_stdout.log 2>>{wd}/_stderr.log "
        )
        cmd1 = (f"trust4  -t 8 -f {fa} -o {wd}/{samplename} "
                f"-1 {wd}/{samplename}_toassemble_1.fq "
                f"-2 {wd}/{samplename}_toassemble_2.fq "
                f"--barcode {wd}/{samplename}_toassemble_bc.fa "
                f"--UMI {wd}/{samplename}_toassemble_umi.fa"
        )
    else:
        cmd1 = (f"trust4  -t 8 -f {fa} -o {wd}/{samplename} "
                f"-u {wd}/{samplename}_toassemble.fq "
                f"--barcode {wd}/{samplename}_toassemble_bc.fa "
                f"--UMI {wd}/{samplename}_toassemble_umi.fa 1>>{wd}/_stdout.log 2>>{wd}/_stderr.log")
    cmd2 = (f"annotator  -f {ref} -a {wd}/{samplename}_final.out -t 8 "
                f"-o {wd}/{samplename}  --barcode --UMI "
                f"--readAssignment {wd}/{samplename}_assign.out "
                f"-r {wd}/{samplename}_assembled_reads.fa "
                f"--airrAlignment > {wd}/{samplename}_annot.fa ") ##1>>{wd}/_stdout.log 2>>{wd}/_stderr.log")
    cmd4 = f"trust-barcoderep.pl {wd}/new.{samplename}_cdr3.out -a {wd}/{samplename}_annot.fa --chainsInBarcode 2 > {wd}/{samplename}_barcode_report.tsv "
    cmd5 = f"trust-simplerep.pl {wd}/new.{samplename}_cdr3.out  --barcodeCnt --filterBarcoderep {wd}/{samplename}_barcode_report.tsv > {wd}/{samplename}_report.tsv "
    cmd6 = f"trust-airr.pl {wd}/{samplename}_report.tsv {wd}/{samplename}_annot.fa --airr-align {wd}/{samplename}_airr_align.tsv > {wd}/{samplename}_airr.tsv"
    cmd7 = f"trust-airr.pl {wd}/{samplename}_barcode_report.tsv {wd}/{samplename}_annot.fa --format barcoderep --airr-align {wd}/{samplename}_airr_align.tsv > {wd}/{samplename}_barcode_airr.tsv"
    logger.info("run trust4")
    cmd_execute(cmd1)
    logger.info("run annotator")
    cmd_execute(cmd2)
    logger.info("run re_assign")
    main(f"{wd}/{samplename}_assign.out", f"{wd}/{samplename}_cdr3.out", f"{wd}/new.{samplename}_cdr3.out", 8)
    logger.info("run trust-report")
    cmd_execute(cmd4)
    cmd_execute(cmd5)
    cmd_execute(cmd6)
    cmd_execute(cmd7)
    logger.info("trust4 done")

def simpleqc_wrapper(fq1:list[Path], fq2:list[Path], wl:Path, rs:str, outprefix:str="./", chunk=500000, threads:int=4):
    cmd = f"simpleqc --fq1 {' '.join(fq1)} --fq2 {' '.join(fq2)} -w {wl} --rs {rs} -o {outprefix} -c {chunk} -t {threads}"
    cmd_execute(cmd)

def search_pattern_wrapper(leader_ref:Path, airr_tsv:Path, annot:Path, outfile:str="./annotations.tsv"):
    cmd = f"search-pattern --leader-ref {leader_ref} --airr {airr_tsv} --annot {annot} --outfile {outfile}"
    cmd_execute(cmd)
