from pathlib import Path
from collections import defaultdict
import pysam
import dnaio
import re
from ..utils.wrappers import cmd_execute
from ..utils.helper import logger

def STAR_build_wrapper(
    ref: str,
    outdir: str,
    core:int=4,
    star_path:str="STAR",
)->str:

    cmd = (f"{star_path} --runMode genomeGenerate --runThreadN {core}  " 
           f"--genomeFastaFiles {ref}  --genomeDir {outdir} ")
    cmd_execute(cmd, check=True)


def format_ref(ref: Path, outdir: Path) -> Path:
    regex = re.compile(r'((TR|IG)\S+)[\|\s:]')
    

    ref_l = outdir / "raw_ref.fa"


    with dnaio.open(ref, fileformat="fasta") as fh_in, \
         dnaio.open(ref_l, mode="w", fileformat="fasta") as fh_out:
        n = 0
        for record in fh_in:
            n += 1
            m = regex.search(record.name)
            if m:
                record.name = f"{m.group(1)}_{n}"
                fh_out.write(record)
    return ref_l
def count_chain(bam: Path):
    """
    """
    d = defaultdict(int)
    tmp = defaultdict(lambda: defaultdict(int))
    _default_verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam) as bamfh:
        for r in bamfh:
            chain = r.reference_name[:3]
            if r.get_tag("NH")==1:
                d[chain] += 1
            else:
                tmp[r.qname][chain] += 1
    
    for k, v in tmp.items():
        if len(v)==1:
            chain, _count = v.popitem()
            d[chain] += 1
        else:
            sorted_v = sorted(v.items(), key=lambda x: x[1], reverse=True)
            if sorted_v[0][1] > sorted_v[1][1]:
                chain, _count = sorted_v[0]
                d[chain] += 1
            else:
                d["unknown"] += 1
    
    return d
def enrichment_qc(
    fq: list,
    chain: str,
    fa: Path,
    samplename:str,
    outdir: Path,
    core: int=4,
    star_path: str="STAR"
):

    ref_dir = outdir / "fa_ref"
    ref_dir.mkdir(exist_ok=True, parents=True)
    ref_l = format_ref(fa, ref_dir)

    STAR_build_wrapper(ref_l, ref_dir, core, star_path)

    prefix = f"{outdir}/{samplename}_"
    from ..utils.wrappers import STAR_wrapper
    STAR_wrapper(
        fq=list(fq), genomeDir=ref_dir, prefix=prefix, scoremin=0.33,
        matchnmin=0.33, core = core, star_path=star_path        
    )
    bamfile = f"{prefix}Aligned.out.bam"
    d = count_chain(bamfile)
    with open(prefix + "chain_summary.txt", "w") as fh:
        for k, v in d.items():
            fh.write(f"{k}: {v}\n")
