import re
from collections import defaultdict
import pysam
import click

def get_corrected_umi(umi_f):
    #read corrected umi file into dict
    with open(umi_f, "r") as f:
        umi_corrected_dict = defaultdict(lambda: defaultdict(dict))
        for l in f:
            bc, ensemble_id, umi_pair_1, umi_pair_2 = l.strip().split("\t")
            umi_1, _umi_count_1 = umi_pair_1.split(":")
            umi_2, _umi_count_2 = umi_pair_2.split(":")
            umi_corrected_dict[bc][ensemble_id][umi_1] = umi_2
    return umi_corrected_dict


def parse_alt(alt, bc, regex = re.compile(r'(\d+)([NACGT])')):
    #get raw barcode based on alteration signiture
    #alt: M means match
    #     1N3T means corrected barcode replaced N at 1st place and T at 3rd place(0 based)
    if alt == "M":
        bc_raw = bc
    else:
        groups = regex.findall(alt)
        for pos, base in groups:
            bc = bc[:int(pos)] + base + bc[int(pos)+1:]
        bc_raw = bc
    return bc_raw

def find_umi_correction(umi_raw, d):
    if umi_raw not in d:
        return umi_raw
    else:
        return d[umi_raw]

def find_umi_correction_recursive(umi_raw, d):
    if umi_raw not in d:
        return umi_raw
    else:
        return find_umi_correction_recursive(d[umi_raw], d)

#@click.command("add barcode and umi tags.")
#@click.option("--inbam", help="input bam.")
#@click.option("--outbam", help="output bam.")
#@click.option("--umifile", "umifile", help="umi.xls in step3.")
#@click.option("--recursive", "recursive", is_flag=True, default=False, show_default=True, 
#              help=(f"whether use recursive to find raw umi or not,"
#                f"set when the software version generating the umi.xls is lower than 1.1")
#            )
def add_tag(inbam:str, outbam:str, umifile:str=None, recursive:bool=False):
    """bam add CB、CR、UB、UR tag
    """
    default_verbosity = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(inbam, "rb")
    pysam.set_verbosity(default_verbosity)
    
    samfile_mod = pysam.AlignmentFile(outbam, "wb", template=samfile)
    if umifile:
        umi_corrected_dict = get_corrected_umi(umifile)
        for read in samfile:
            #read samfile
            tmp = read.query_name.split("_", 3)
            if len(tmp)==4:
                bc, umi_raw, alt, readname = tmp
                bc_raw = parse_alt(alt, bc)
            else:
                bc, umi_raw, readname = tmp
                bc_raw = bc

            if read.has_tag("XT"):
                # check if umi correction is needed
                # multi ensemble_id?
                ensemble_id = read.get_tag("XT")

                # try if ensemble_id is contained in dict
                # if contained, do umi correction
                if ensemble_id in umi_corrected_dict[bc]:
                    if recursive:
                        umi = find_umi_correction_recursive(umi_raw, umi_corrected_dict[bc][ensemble_id])
                    else:
                        umi = find_umi_correction(umi_raw, umi_corrected_dict[bc][ensemble_id])
                else:
                    umi = umi_raw
            else:
                umi = umi_raw

            read.set_tags(read.get_tags() + [("CB", bc), ("CR", bc_raw), ("UB", umi), ("UR", umi_raw)])
            # change read.query_name by removing barcode, umi and alteration
            read.query_name = readname
            samfile_mod.write(read)
    else:
        for read in samfile:
            tmp = read.query_name.split("_", 3)
            if len(tmp)==4:
                bc, umi_raw, alt, readname = tmp
                bc_raw = parse_alt(alt, bc)
            else:
                bc, umi_raw, readname = tmp
                bc_raw = bc
            umi = umi_raw
            read.set_tags(read.get_tags() + [("CB", bc), ("CR", bc_raw), ("UB", umi), ("UR", umi_raw)])
            read.query_name = readname
            samfile_mod.write(read)

    samfile_mod.close()
    samfile.close()
    return

if __name__ == "__main__":
    add_tag()
