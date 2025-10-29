import os
import shutil
import glob
import json
from functools import partial
from collections import namedtuple
import dnaio
import pysam

# Since fields with a default value must come after any fields without a default,
# the defaults are applied to the rightmost parameters. 
ColDesc= namedtuple("Colinfo_with_name", ["name", "sep", "idx"], defaults=["_", 0])
ColDesc2= namedtuple("Colinfo_no_name", ["col_idx", "sep", "idx"], defaults=["_", 0])

def refine_table(
    infile:str,
    outdir:str,
    colnames:list[ColDesc],
    map_dict:dict, 
    sep:str="\t",
) -> None:
    with open(infile) as fh, \
        open(f"{outdir}/{os.path.basename(infile)}", "w") as fh_out:
        header_line = fh.readline().strip()
        fh_out.write(f"{header_line}\n")
        header_dict = dict([(e, n) for n, e in enumerate(header_line.split(sep))])
        #print(header_dict)
        for line in fh:
            tmp = line.strip().split(sep)
            for _name, _sep, _idx in colnames:
                col_idx = header_dict[_name]
                raw = tmp[col_idx]
                if _sep:
                    _raw_tmp = raw.split(_sep)
                    _raw_tmp[_idx] = _raw_tmp[_idx].replace("-1", "")
                    _raw_tmp[_idx] = map_dict[_raw_tmp[_idx]]
                    new = _sep.join(_raw_tmp)
                else:    
                    new = map_dict[raw]
                tmp[col_idx] = new
            fh_out.write(f"{sep.join(tmp)}\n")

def refine_table_no_header(
    infile:str,
    outdir:str,
    colnames:list[ColDesc2],
    map_dict:dict,
    sep:str="\t",
)->None:
    with open(infile) as fh, \
        open(f"{outdir}/{os.path.basename(infile)}", "w") as fh_out:
        for line in fh:
            tmp = line.strip().split(sep)
            for col_idx, _sep, _idx in colnames:
                raw = tmp[col_idx]
                if _sep:
                    _raw_tmp = raw.split(_sep)
                    _raw_tmp[_idx] = _raw_tmp[_idx].replace("-1", "")
                    _raw_tmp[_idx] = map_dict[_raw_tmp[_idx]]
                    new = _raw_tmp[_idx]
                else:    
                    new = map_dict[raw]
                tmp[col_idx] = new
            fh_out.write(f"{sep.join(tmp)}\n")  

def refine_cell_barcode_json(
    infile:str,
    outdir:str,
    map_dict:dict,

)->None:
    with open(infile) as fh:
        data = json.load(fh)
    
    new_data = []
    for b in data:
        tmp = b.split("_")[0].replace("-1", "")
        tmp = map_dict[tmp]
        new_data.append(tmp)

    with open(f"{outdir}/{os.path.basename(infile)}", "w") as fh:
        json.dump(new_data, fp=fh, indent=4)

def refine_all_contig_annotations_json(
    infile:str,
    outdir:str,
    map_dict:dict,
    umi_len: int,
)->None:
    with open(infile) as fh:
        data = json.load(fh)

    for i in range(0, len(data)):
        tmp1 = data[i]["barcode"].split("_")[0].replace("-1", "")
        tmp1 = map_dict[tmp1]
        #new_barcode = "_".join(tmp1)
        data[i]["barcode"] = tmp1
        
        tmp2 = data[i]["contig_name"].split("_")
        tmp2[0] = tmp2[0].replace("-1", "")
        tmp2[0] = map_dict[tmp2[0]]
        new_contig_name = "_".join(tmp2)
        data[i]["contig_name"] = new_contig_name

        data[i]["validated_umis"] = [
            u[:umi_len] for u in data[i]["validated_umis"]
        ]

    with open(f"{outdir}/{os.path.basename(infile)}", "w") as fh:
        json.dump(data, fp=fh, indent=4)

def refine_fastx(
    infile:str,
    fastx:str,
    outdir:str,
    map_dict:dict,
    sep:str="_",
    idx:int=0,
    
):
    with dnaio.open(infile, fileformat=fastx, mode="r") as fh, \
        dnaio.open(f"{outdir}/{os.path.basename(infile)}", fileformat=fastx, mode="w") as fh_out:
        for r in fh:
            _raw_tmp = r.name.split(sep)
            _raw_tmp[idx] = _raw_tmp[idx].replace("-1", "")
            _raw_tmp[idx] = map_dict[_raw_tmp[idx]]
            r.name = sep.join(_raw_tmp)
            fh_out.write(r)

def refine_bam(
    infile:str,
    outdir:str,
    map_dict:dict,
    umi_len:int=8,
    sep:str="_"
):
    default_verbosity = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(infile, "rb")
    pysam.set_verbosity(default_verbosity)    
    
    samfile_mod = pysam.AlignmentFile(f"{outdir}/{os.path.basename(infile)}", "wb", template=samfile)
    
    for r in samfile:
        header = r.query_name
        tmp = header.split(sep)
        tmp[0] = tmp[0].replace("-1", "")
        if tmp[0] in map_dict.keys():
            tmp[0] = map_dict[tmp[0]]
        r.query_name = sep.join(tmp)
    
        raw_bc = r.get_tag("CB").replace("-1", "")
        bc = map_dict[raw_bc]
    
        if r.has_tag("CR"):
            r.set_tags([("CB", bc), 
                        ("CR", bc), 
                        ("UB", r.get_tag("UB")[:umi_len]), 
                        ("UR", r.get_tag("UR")[:umi_len])])
        else:
            r.set_tags([("CB", bc)])

        samfile_mod.write(r)
    
    samfile_mod.close() 
    samfile.close()    


def copy(
    infile:str,
    outdir:str,
)->None:
    shutil.copy(infile, f"{outdir}/{os.path.basename(infile)}")

def get_header(
    basename:str, 
    outdir:str, 
    map_file:str,
    map_dict:dict
):
    bam_file = os.path.join(outdir, basename)
    bam_tmp = os.path.join(outdir, "tmp.bam")
    header_file = os.path.join(outdir, "header.sam")
    header_mod = os.path.join(outdir, "header_mod.sam")

    cmd = "samtools view -H " + bam_file + " > " + header_file
    os.system(cmd)

    with open(header_file, "r") as f, \
            open(header_mod, "w") as fm:
        for line in f:
            tmp = line.strip().split("\t")
            if tmp[1].startswith("SN"):
                bc_tmp = tmp[1].split(":")[1].split("_")
                bc_tmp[0] = bc_tmp[0].replace("-1", "")
                bc_tmp[0] = "SN:" + map_dict[bc_tmp[0]]
                tmp[1] = "_".join(bc_tmp)
                new_hname = "\t".join(tmp)
            else:
                new_hname = "\t".join(tmp)
            fm.write(new_hname+"\n")
    cmd_reheader = f"samtools reheader {header_mod} {bam_file} > {bam_tmp} && samtools sort -@ 4 {bam_tmp} > {bam_file} && rm {bam_tmp} "
    os.system(cmd_reheader)
    return 

def get_output(
    indir:str,
    outdir:str,
    map_file:str,
    umi_len:int=8,
    clean:bool=False,
    **kwargs
)->None:

    map_dict = {}
    with open(map_file) as fh:
        for line in fh:
            tmp = line.strip().split()
            map_dict[tmp[1]] = tmp[0]
    dispatch_table = {
        
        "all_contig_annotations.csv": partial(refine_table,
            outdir=outdir, map_dict=map_dict, sep=",",
            colnames=[
                ColDesc(name="barcode"),
                ColDesc(name="contig_id")
            ],
        ),
        "filtered_contig_annotations.csv": partial(refine_table,
            outdir=outdir, map_dict=map_dict, sep=",",
            colnames=[
                ColDesc(name="barcode"),
                ColDesc(name="contig_id")
            ],
        ),
        "clonotypes.csv": partial(copy, outdir=outdir),
        "consensus_annotations.csv": partial(copy, outdir=outdir),
        "metrics_summary.csv": partial(copy, outdir=outdir),

        "all_contig.fasta": partial(refine_fastx,
            outdir=outdir, map_dict=map_dict, fastx="fasta"
        ),
        "filtered_contig.fasta": partial(refine_fastx,
            outdir=outdir, map_dict=map_dict, fastx="fasta"
        ),
        "concat_ref.fasta": partial(copy, outdir=outdir),
        "consensus.fasta": partial(copy, outdir=outdir),

        "all_contig.fastq": partial(refine_fastx,
            outdir=outdir, map_dict=map_dict, fastx="fastq"
        ),
        "filtered_contig.fastq": partial(refine_fastx,
            outdir=outdir, map_dict=map_dict, fastx="fastq"
        ),
        "airr_rearrangement.tsv": partial(refine_table,
            outdir=outdir, map_dict=map_dict, sep="\t",
            colnames=[
                ColDesc(name="cell_id"),
                ColDesc(name="sequence_id")
            ],
        ),
        "all_contig_annotations.bed": partial(refine_table_no_header,
            outdir=outdir, map_dict=map_dict, sep="\t",
            colnames=[ColDesc2(col_idx=0),],
        ),
        "all_contig_annotations.json": partial(refine_all_contig_annotations_json,
            outdir=outdir, map_dict=map_dict, umi_len=umi_len
        ),
        "cell_barcodes.json": partial(refine_cell_barcode_json,
            outdir=outdir, map_dict=map_dict
        ),
     
        "all_contig.bam": partial(refine_bam,
            outdir=outdir, map_dict=map_dict, umi_len=umi_len, sep="_"
        ),
        "consensus.bam": partial(refine_bam,
            outdir=outdir, map_dict=map_dict, umi_len=umi_len, sep="_"
        ),
        "concat_ref.bam": partial(copy, outdir=outdir)
        
    }

    os.makedirs(outdir, exist_ok=True) 
    files = glob.glob(indir + "/*", recursive=True)
    for f in files:
        basename = os.path.basename(f)
        if basename in dispatch_table:
            dispatch_table[basename](f)
            if basename == "all_contig.bam":
                get_header(basename, outdir, map_file, map_dict)
    if clean:
        shutil.rmtree(indir)
    return

def convert(
    fq:str,
    fq1:str,
    fq2:str,
    cr_barcode:str, # 737K-august-2016.txt
    map_file:str
):
    _barcodes = set()
    _barcode_invaild = ""
    _n = 0
    with open(cr_barcode) as fh:
        for l in fh:
            l = l.strip()
            if l.startswith("#"):
                continue
            if not l:
                continue

            _barcodes.add(l)
            _n += 1
            if _n==1:
                _barcode_invaild = "A"*len(l)

    barcode_map = {}
    with dnaio.open(fq, fileformat='fastq', mode='r') as fh, \
        dnaio.open(fq1, fq2, fileformat='fastq', mode='w') as outfh:
        for r in fh:
            barcode, umi, _, r2_name = r.name.split("_", 3)
            if not barcode in barcode_map:
                if len(_barcodes):
                    barcode_map[barcode] = _barcodes.pop()
                else:
                    barcode_map[barcode] = _barcode_invaild
            new_barcode = barcode_map[barcode]
            r1_name = r2_name.replace(" 1:"," 2:").replace("/2", "/1")
            # r1_sequence = barcode + umi[:10] # "length": 10,
            r1_sequence = new_barcode + (umi + "AAAAAAAA")[:12]
            r1_qualities = r.qualities[:len(r1_sequence)]
            r1 = dnaio.Sequence(
                name = r1_name,
                sequence = r1_sequence,
                qualities = r1_qualities
            )
            r.name = r2_name
            outfh.write(r1, r)

    with open(map_file, "w") as fh:
        for k, v in barcode_map.items():
            fh.write(f"{k}\t{v}\n")
