import sys
import re
from loguru import logger
import click

def gene_gen(gtf):
    with open(gtf) as fh:
        content = []
        record_header = None

        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#"):
                yield ("__COMMENT__", line)
                continue

            tmp = line.split("\t")

            if tmp[2] in ("gene", "transcript"):   
                yield record_header, content
                record_header = tmp
                content = []
            else:
                content.append(tmp)
        
        if content:
            yield record_header, content

def prepare_content(content:list):
    tmp_str = ""
    features = {l[2] for l in content}
    for l in content:           
        if l[2] == "CDS" and ("exon" not in features):
            exon_tmp = l[:2] + ["exon", ] + l[3:]
            tmp_str += ("\t".join(exon_tmp) + "\n")
            tmp_str += ("\t".join(l) + "\n")
        else:
            tmp_str += ("\t".join(l) + "\n")
    return tmp_str 

# @click.command("gtf validator")
# @click.option("--gtf", "gtf", type=click.Path(), required=True, help="gtf path.")
# @click.option("--newgtf", "new_gtf", default=None, help="modified gtf file. if not set, only do check.")
# @click.option("--gid", "gene_id_str", default="gene_id", show_default=True, help="attribute to extract gene id.")
def validator(gtf, new_gtf:str=None, gene_id_str:str="gene_id"):
    gene_id_regex = re.compile(rf'{gene_id_str} "(.*?)";')

    gtf_gen = gene_gen(gtf)

    if new_gtf:
        outfh = open(new_gtf, "w")

    warning_line_num = 0
    gene_id_set = set()
    for (record_header, content) in gtf_gen:
        tmp_str = ""
        if not record_header:
            continue
        if record_header == '__COMMENT__':
            tmp_str += f"{content}\n"
        else:
            m = gene_id_regex.search(record_header[-1])
            if m:
                gene_id_ori = m.group(1)

                # gene_id = gene_id_ori.replace('"', '').replace("'", "").replace("-", ".").replace("_", ".")
                gene_id = gene_id_ori.replace('"', '').replace("'", "")
                if record_header[2]=="gene" and (gene_id in gene_id_set):
                    logger.error(f"duplicate gene IDs: {gene_id}")
                    sys.exit(1)
                else:
                    gene_id_set.add(gene_id)
                if gene_id_ori != gene_id:
                    record_header[-1] = record_header[-1].replace(f'"{gene_id_ori}";', f'"{gene_id}";')
                    for l in content:
                        l[-1] = record_header[-1].replace(f'"{gene_id_ori}";', f'"{gene_id}";')
                    logger.warning(f"{gene_id_ori} -> {gene_id}")
                    warning_line_num += 1
                if gene_id_str!="gene_id":
                    record_header[-1] += f' gene_id "{gene_id}";'
            else:
                logger.error(f"{gene_id_str} NOT FOUND!\n{' '.join(record_header)}")
                sys.exit(2)

            tmp_str += ("\t".join(record_header) + "\n")

            if content:
                tmp_str += prepare_content(content)

        if new_gtf:
            outfh.write(tmp_str)
    logger.info(f"total genes: {len(gene_id_set)}.")
    if new_gtf:
        outfh.close()
        logger.info(f"modified gtf file: {new_gtf}.")
    else:
        logger.info(f"total warnings: {warning_line_num}.")

# if __name__ == "__main__":
#     validator()
