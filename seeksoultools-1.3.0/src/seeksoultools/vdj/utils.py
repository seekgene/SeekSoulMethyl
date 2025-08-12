from pathlib import Path
import json
from collections import defaultdict
import dnaio
import pandas as pd
from ..utils.helper import logger 

CWD = (Path(__file__).expanduser())

WL = (CWD/"../../utils/barcode/P3CBGB/P3CB.barcode.txt.gz").resolve()

VDJ_REF = {
    "human": {
        "fa": (CWD/"../../utils/trust4_ref/Homo_sapien/bcrtcr.fa").resolve(),
        "ref": (CWD/"../../utils/trust4_ref/Homo_sapien/IMGT+C.fa").resolve(),
        "leader": {
            "IG": (CWD/"../../utils/trust4_ref/Homo_sapien/IG_L-PART1+L-PART2.fa").resolve(),
            "TR": (CWD/"../../utils/trust4_ref/Homo_sapien/TR_L-PART1+L-PART2.fa").resolve(),
        }
    },
    "mouse": {
        "fa": (CWD/"../../utils/trust4_ref/Mus_musculus/bcrtcr.fa").resolve(),
        "ref": (CWD/"../../utils/trust4_ref/Mus_musculus/IMGT+C.fa").resolve(),
        "leader": {
            "IG": (CWD/"../../utils/trust4_ref/Mus_musculus/IG_L-PART1+L-PART2.fa").resolve(),
            "TR": (CWD/"../../utils/trust4_ref/Mus_musculus/TR_L-PART1+L-PART2.fa").resolve(),
        }
    },

}

def get_config(species:str, chain:str, cfg:None|Path) -> Path:
    ref_path, fa_path, leader_path = None, None, None
    if species:
        try:
            ref_path = VDJ_REF[species].get("ref", None)
            fa_path = VDJ_REF[species].get("fa", None)
            leader_path = VDJ_REF[species]["leader"][chain]
        except KeyError as e:
            logger.error(f"Key error occurred: {e}")
    elif cfg:
        try:
            with cfg.open() as fh:
                for line in fh:
                    if line.startswith("ref"):
                        ref_path = Path(line.split(":")[1].strip())
                    if line.startswith("fa"):
                        fa_path = Path(line.split(":")[1].strip())
                    if line.startswith("leader"):
                        leader_path = Path(line.split(":")[1].strip())              
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")
    else:
        logger.error("Either --organism or --cfg must be specified.")

    if ref_path:
        if not ref_path.exists():
            logger.error(f"Ref file does not exist. File path: {ref_path}")
    else:
        logger.error("Ref file path not found.")

    if fa_path:
        if not fa_path.exists():
            logger.error(f"Fa file does not exist. File path: {fa_path}.")
    else:
        logger.error("Fa file path not found.")

    return ref_path, fa_path, leader_path

def get_white_list() -> Path:
    if WL.exists():
        return WL
    else:
        logger.error("No white list file found")

def count_barcode(bc_fa:Path, umi_fa:Path, out:Path):
    count_dict = defaultdict(lambda: defaultdict(int))
    valid_reads = 0
    with dnaio.open(bc_fa, fileformat="fasta") as fh_b, \
        dnaio.open(umi_fa, fileformat="fasta") as fh_u:
        for b, u in zip(fh_b, fh_u):
            if b.sequence != "missing_barcode":
                valid_reads += 1
                count_dict[b.sequence][u.sequence] += 1
    tmp = []
    for k, v in count_dict.items():
        umi_num = len(v.keys())
        reads_num = sum(v.values())
        tmp.append([k, umi_num, reads_num])
    del count_dict

    count_dict = {}
    with out.open("w") as fh:
        fh.write("barcode\tumi_num\treads_num\n")
        for r in sorted(tmp, key=lambda x: x[-1], reverse=True):
            k, umi_num, reads_num = r
            fh.write(f"{k}\t{umi_num}\t{reads_num}\n")
            count_dict[k] = reads_num
    return valid_reads, count_dict

def readFormat2rs(readFormat:str):
    """
    https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
    convert readFormat to Read-Structures
    readFormat: "bc:0:15,um:16:25"
    rs: "17C25M+T"
    """

    rs = ""
    tmp = sorted(
        [p.split(":") for p in readFormat.split(",")],
        key=lambda x: int(x[1])
    )

    for (p, start, end) in tmp:
        start = int(start)
        end = int(end)
        if p == "bc":
            rs += f"{end-start+1}C"
        elif p == "um":
            rs += f"{end-start+1}M"
    rs += "+T"
    return rs

def subset_reads(
    count_dict:dict, bc_fa:Path, umi_fa:Path, 
    r2_fq:Path, r1_fq:Path|None=None, reads_limit:int=80000
):
    """
    """
    _bc_fa = bc_fa.replace(f"{bc_fa.parent}/_{bc_fa.name}")
    _umi_fa = umi_fa.replace(f"{umi_fa.parent}/_{umi_fa.name}")
    _r2_fq = r2_fq.replace(f"{r2_fq.parent}/_{r2_fq.name}")
    if r1_fq is not None:
        _r1_fq = r1_fq.replace(f"{r1_fq.parent}/_{r1_fq.name}")
    bc_fh = dnaio.open(_bc_fa, fileformat="fasta")
    bc_out = dnaio.open(bc_fa, fileformat="fasta", mode="w")
    umi_fh = dnaio.open(_umi_fa, fileformat="fasta")
    umi_out = dnaio.open(umi_fa, fileformat="fasta", mode="w")
    r2_fh = dnaio.open(_r2_fq, fileformat="fastq")
    r2_out = dnaio.open(r2_fq, fileformat="fastq", mode="w")   
    if r1_fq:
        r1_fh = dnaio.open(_r1_fq, fileformat="fastq")
        r1_fh_iter = iter(r1_fh)
        r1_out = dnaio.open(r1_fq, fileformat="fastq", mode="w")
    for b, u, r2 in zip(bc_fh, umi_fh, r2_fh):
        if r1_fq:
            r1 = next(r1_fh_iter)
        if b.sequence in count_dict:
            count_dict[b.sequence] += 1
            if count_dict[b.sequence] <= reads_limit:
                bc_out.write(b)
                umi_out.write(u)
                r2_out.write(r2)
                if r1_fq: 
                    r1_out.write(r1)
        else:
            bc_out.write(b)
            umi_out.write(u)
            r2_out.write(r2)
            if r1_fq: 
                r1_out.write(r1)
    bc_fh.close()
    bc_out.close()

    umi_fh.close()
    umi_out.close()

    r2_fh.close()
    r2_out.close()

    if r1_fq:
        r1_fh.close()
        r1_out.close()

    # unlink
    _bc_fa.unlink()
    _umi_fa.unlink()
    _r2_fq.unlink()
    if r1_fq:
        _r1_fq.unlink()

def cal_valid(summary_file:Path):
    with summary_file.open() as fh:
        for line in fh:
            if line.startswith("total_num:"):
                total_num =  int(line.strip().split(":")[-1])
            if line.startswith("valid_num:"):
                valid_num = int(line.strip().split(":")[-1])
    return {
        "Number of Read Pairs": total_num,
        "Number of Valid Reads": valid_num,
        "Valid Barcodes":valid_num/total_num
    }

def cal_q30(r1_qual:Path, r2_qual:Path, readFormat:str):
    d = {}
    r1_qual_df = pd.read_csv(r1_qual).drop("pos", axis=1)
    d["Number of Read Pairs"] = int(r1_qual_df.iloc[0,:].sum())
    r1_len = r1_qual_df.shape[0]
    readFormat = [p.split(":") for p in readFormat.split(",")]
    readFormat = {x: range(int(s), int(e)+1) if int(e)>0 else range(int(s), r1_len) for (x, s, e) in readFormat}

    bc_qual_df = r1_qual_df.iloc[readFormat["bc"],:]
    d["Q30 Bases in Barcode"] = bc_qual_df.iloc[:,30: ].sum().sum()/bc_qual_df.sum().sum() if bc_qual_df.sum().sum()>0 else 0
    um_qual_df = r1_qual_df.iloc[readFormat["um"],:]
    d["Q30 Bases in UMI"] = um_qual_df.iloc[:,30: ].sum().sum()/um_qual_df.sum().sum() if um_qual_df.sum().sum()>0 else 0

    r2_qual_df = pd.read_csv(r2_qual).drop("pos", axis=1)
    d["Q30 Bases in R2 Read"] = r2_qual_df.iloc[:,30:].sum().sum()/r2_qual_df.sum().sum() if r2_qual_df.sum().sum()>0 else 0
    return d

def summary_simpleqc(summary_file:Path, wd:Path, samplename:str, readFormat:str):
    if not summary_file.exists():
        summary = {}
    else:
        with summary_file.open() as fh:
            summary = json.load(fh)
    d1 = cal_valid(wd/"Analysis"/"qc"/f"{samplename}_summary.txt")
    summary.update(d1)

    d2 = cal_q30(wd/"Analysis"/"qc"/f"{samplename}_r1_qual_stat.csv", wd/"Analysis"/"qc"/f"{samplename}_r2_qual_stat.csv", readFormat)
    summary.update(d2)

    with summary_file.open(mode="w") as fh:
        json.dump(summary, fh, indent=4)


def get_mapping_matrix(chain_summary:Path, chain:str, valid_reads:int):
    chain_map_dict = {}
    with chain_summary.open() as fh:
        for line in fh:
            tmp = line.strip().split(":")
            chain_map_dict[tmp[0]] = int(tmp[1])        
    d = {}
    if chain=="IG":
        d["Reads Mapped to IGH"] = chain_map_dict.get("IGH", 0)/valid_reads
        d["Reads Mapped to IGK"] = chain_map_dict.get("IGK", 0)/valid_reads
        d["Reads Mapped to IGL"] = chain_map_dict.get("IGL", 0)/valid_reads
        d["Reads Mapped to Any V(D)J Gene"] = (
            d["Reads Mapped to IGH"] + 
            d["Reads Mapped to IGK"] + 
            d["Reads Mapped to IGL"]
        )
    elif chain=="TR":
        d["Reads Mapped to TRA"] = chain_map_dict.get("TRA", 0)/valid_reads
        d["Reads Mapped to TRB"] = chain_map_dict.get("TRB", 0)/valid_reads
        d["Reads Mapped to Any V(D)J Gene"] = (
            d["Reads Mapped to TRA"] + 
            d["Reads Mapped to TRB"]
        )
    else:
        logger.error(f"{chain} is not supported, should be one of IG, TR.")
    return d

def summary_enrichment_qc(summary_file:Path, wd:Path, samplename:str, chain:str):
    if not summary_file.exists():
        summary = {}
    else:
        with summary_file.open() as fh:
            summary = json.load(fh)
    valid_reads = summary["Number of Valid Reads"]
    d = get_mapping_matrix(
        wd/"Analysis"/"qc"/f"{samplename}_chain_summary.txt",
        chain, valid_reads
    )
    summary.update(d)

    with summary_file.open(mode="w") as fh:
        json.dump(summary, fh, indent=4)
    val = d['Reads Mapped to Any V(D)J Gene']
    if val < 0.01:
        logger.error(f"Reads Mapped to Any V(D)J Gene: {val}, no {chain} gene in your data, please check your data and chain")
def summary_fastq_extractor(summary_file:Path, wd:Path, samplename:str):
    if not summary_file.exists():
        summary = {}
    else:
        with summary_file.open() as fh:
            summary = json.load(fh)

    df = pd.read_table(wd/"Analysis"/"trust4"/f"{samplename}_bc_count.tsv")
    candidate_reads= int(df.reads_num.sum())
    summary["Number of Candidate Reads"] = candidate_reads

    if (wd/"Analysis"/"trust4"/f"_{samplename}_bc_count.tsv").exists():
        df_full = pd.read_table(wd/"Analysis"/"trust4"/f"_{samplename}_bc_count.tsv")
        full_reads = int(df_full.reads_num.sum())
        summary["Number of Filtered Candidate Reads"] = full_reads - candidate_reads
    else:
        summary["Number of Filtered Candidate Reads"] = 0

    with summary_file.open(mode="w") as fh:
        json.dump(summary, fh, indent=4)
def anno_cr(airr, anno_file):
    airr = pd.read_table(airr, sep='\t')
    anno = pd.read_csv(anno_file, sep='\t')
    result = pd.DataFrame(columns=[
    'barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain',
    'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive',
    'fwr1', 'fwr1_nt', 'cdr1', 'cdr1_nt', 'fwr2', 'fwr2_nt', 'cdr2', 'cdr2_nt',
    'fwr3', 'fwr3_nt', 'cdr3', 'cdr3_nt', 'fwr4', 'fwr4_nt', 'reads', 'umis',
    'raw_clonotype_id', 'raw_consensus_id', 'exact_subclonotype_id'
    ])
    merged_df = pd.merge(airr, anno, on='sequence_id', suffixes=('_airr', '_anno'))
    result['barcode'] = merged_df['cell_id']
    result['is_cell'] = 'T'
    result['contig_id'] = merged_df['sequence_id']
    result['length'] = merged_df['sequence_airr'].str.len()
    result['full_length'] = merged_df['complete_vdj']
    result['chain'] = merged_df['locus']
    result['v_gene'] = merged_df['v_call']
    result['d_gene'] = merged_df['d_call']
    result['j_gene'] = merged_df['j_call']
    result['c_gene'] = merged_df['c_call']
    result['productive'] = merged_df['productive']
    result['fwr1_nt'] = merged_df['fwr1']
    result['fwr2_nt'] = merged_df['fwr2']
    result['fwr3_nt'] = merged_df['fwr3']
    result['fwr4_nt'] = merged_df['fwr4']
    result['reads'] = merged_df['consensus_count']
    result['umis'] = merged_df['duplicate_count']
    result['raw_clonotype_id'] = merged_df['clone_id']
    result['exact_subclonotype_id'] = merged_df['exact_subclonotype_id']
    result['cdr1_nt'] =  merged_df['cdr1_airr']
    result['cdr2_nt'] =  merged_df['cdr2_anno']
    result['cdr3'] = merged_df['junction_aa']
    result['cdr3_nt'] = merged_df['junction']

    result.to_csv(anno_file, sep='\t', index=False)
