import os
import re
from ..utils.helper import logger, get_file_path
from ..utils.wrappers import cmd_execute

def do_seurat(
    matrix:str, samplename:str, outdir:str, gtf:str, dims:int=15, minpct:float=0.1,
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
    logger.info('anno diff gene started')
    markergene=os.path.join(outdir1,'FindAllMarkers.xls')
    annogene=os.path.join(outdir1,'biotype_FindAllMarkers.xls')
    lncgene=os.path.join(outdir1,'lncgene_FindAllMarkers.xls')
    typedic={}
    with open(gtf, "r") as fh:
        for line in fh:
            if not line.strip() : continue
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            ids = tmp[-1]
            if "gene_id" in ids:
                    gene_id=re.search('gene_id "(.*?)";',ids).group(1)
                    genetype=re.search('(gene_type|gene_biotype) "(.*?)"',ids)
                    gene_names=re.search('gene_name "(.*?)";',ids)
                    if  genetype:
                        gene_biotype =  genetype.group(2)
                    else:
                        gene_biotype = "undefine"
                    if gene_names :
                        gene_name = gene_names.group(1)
                    else:
                        gene_name = gene_id       
                    gene_name = gene_name.replace("_","-")
                    if gene_biotype in ["lincRNA", "antisense", "lncRNA", "lnc_RNA", "lnc"]:
                          types="lncRNA"
                    else:
                          types=gene_biotype
                    typedic[gene_name] = types
    with open(markergene) as infile:
        with open(annogene,'w') as out:
            with open(lncgene,'w') as lncout:
                 head=infile.readline()
                 headlist=head.strip().split('\t')
                 headlist.insert(2,"bio_type")
                 out.write('\t'.join(headlist)+'\n')
                 lncout.write('\t'.join(headlist)+'\n')
                 for lines in infile:
                     line=lines.strip().split('\t')
                     genename=line[1]
                     if genename not in typedic:continue
                     biotype=typedic[genename]
                     line.insert(2,biotype)
                     out.write('\t'.join(line)+'\n')
                     if biotype=="lncRNA":
                         lncout.write('\t'.join(line)+'\n')
