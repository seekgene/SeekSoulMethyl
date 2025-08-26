import os
import re
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
from .wrappers import STAR_wrapper
from .helper import logger

def lowerCount(alist,b):
    return len([i for i in alist if i>b])/float(len(alist))

def actb(gtf, bam):

    default_verbosity = pysam.set_verbosity(0)
    bam_in = pysam.AlignmentFile(bam,'rb')
    pysam.set_verbosity(default_verbosity)
    
    depthlist=[]
    over_02mean = 0
    needgene=["ACTB","Actb"]
    region = {}
    with open(gtf,'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                le = line.strip().split("\t")
                _chr, ge_type, start_pos, end_pos, gene_str = le[0], le[2], int(le[3]), int(le[4]), le[-1]
                if "gene_name" not in gene_str:continue
                gene = re.findall("gene_name \"([A-Za-z0-9_\.\-\:/\ ()\[\]']+)\"",gene_str)[0]
                if gene not in needgene:continue
                if ge_type == "exon":
                    chrn = _chr
                    region.update({start_pos:end_pos})
    if len(region) == 0:
        over_02mean = 0.0

    else:
        for s in region:
            e=region[s]
            regions="%s:%s-%s" %(chrn,s,e)
            pileupcolumns = bam_in.pileup(region=regions, stepper='nofilter', truncate=True, max_depth=8000000)
            for pileupcolumn in pileupcolumns:
                depthlist.append(pileupcolumn.n)
        mean_d = np.mean(depthlist)
        if len(depthlist) == 0:
            over_02mean = 0.0
        else:
            over_02mean=lowerCount(depthlist, 0.2*mean_d)
    logger.info("ACTB:over_02mean:", over_02mean)
    return over_02mean


def rRNA_map(fq:list, rRNAgenomeDir, rRNAgtf, rRNAdir, samplename, scoremin, matchnmin, core=4, star_path='STAR'):
    rRNAprefix = os.path.join(rRNAdir, samplename+'_')
    from ..utils.wrappers import STAR_wrapper
    rRNAbam, rRNAlog = STAR_wrapper(
        fq, rRNAgenomeDir, rRNAprefix, scoremin, matchnmin, core=core, star_path=star_path, readMapNumber=800000
    )

    countfile = os.path.join(rRNAdir, 'counts.txt')
    _countout = os.path.join(rRNAdir, samplename+'_Log.final.out')

    from ..utils.wrappers import featureCounts_wrapper

    # -t: gene; region: gene!
    # "--countReadPairs": ""
    _ = featureCounts_wrapper(
        bam=rRNAbam, gtf=rRNAgtf, samplename="rRNA", outdir=rRNAdir,
        region="gene", core=core, SC5P=None, **{"-p": "", "--fraction": ""},
    )

    adic={}
    with open(rRNAgtf) as afile:
        for lines in afile:
            if "gene_id" not in lines: continue
            line=lines.strip().split('\t')
            if line[2] != "gene": continue
            ids=line[-1]
            name=re.search('gene_id "(.*?)";',ids).group(1)
            genetype=re.search('(gene_type|gene_biotype) "(.*?)"',ids)
            if  genetype:
                gene_biotype =  genetype.group(2)
            else:
                gene_biotype = "undefine"
            adic[name] = gene_biotype

    annot_df = pd.DataFrame(adic,index=['type']).T.reset_index()
    df = pd.read_table(countfile, comment='#', usecols=[0, 5, 6])
    df.columns.values[-1] = 'Count'
    allmap = df['Count'].sum()
    result = pd.merge(df, annot_df, how='left', left_on='Geneid', right_on='index').drop('index', axis=1)
    result.loc[result.Count>0,:].to_csv(os.path.join(rRNAdir, samplename+'.xls'), sep='\t', index=False)
    

    df0 = result.loc[result.Count>0,:]
    df0 = df0.groupby('type').agg({'Count':'sum'})
    df0.columns.values[-1] = samplename
    df0.loc['Assigned_Features'] = [sum(df0.loc[:,samplename])]
    if 'Mt_rRNA' not in df0.index:
        mt=0.0
    else:
        mt=df0.loc['Mt_rRNA',samplename]
    df0['type']=df0.index.values.tolist()
    if 'rRNA' not in df0.index and 'rRNA_pseudogene' not in df0.index:
        rrna=0.0
    else:
        rnslist=['rRNA','rRNA_pseudogene']
        rnadf=df0.loc[df0['type'].isin(rnslist), :].reset_index(drop=True)
        rrna=rnadf[samplename].sum()
    
    if allmap:
        Mt_rRNA= float(mt)/allmap
        rRNA = float(rrna)/allmap
    else:
        Mt_rRNA =0
        rRNA = 0
    return Mt_rRNA, rRNA

def mapping_summary(STARLog, RnaSeqMetrics):
    summary = defaultdict()
    with open(STARLog, 'r') as fh:
        for line in fh:
            if 'Number of input reads' in line:
                summary['Number of input reads'] = int(
                    line.strip().split('\t')[-1])
            if 'Uniquely mapped reads number' in line:
                summary['Uniquely mapped reads number'] = int(
                    line.strip().split('\t')[-1])
            if 'Number of reads mapped to multiple loci' in line:
                summary['Number of reads mapped to multiple loci'] = int(
                    line.strip().split('\t')[-1])
    with open(RnaSeqMetrics, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith('total alignments'):
                summary['total alignments'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('reads aligned'):
                summary['reads aligned'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('aligned to genes'):
                summary['aligned to genes'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('no feature assigned'):
                summary['no feature assigned'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('exonic'):
                summary['exonic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intronic'):
                summary['intronic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intergenic'):
                summary['intergenic'] = int(line.split()[-2].replace(',', ''))
                break
    return summary

def read_gtf(gtf):
    genetype_table = []
    mt_regex = re.compile("^(MT|mt|Mt)-")
    with open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if tmp[2] == "gene":
                if "gene_id" in tmp[-1] and 'type' in  tmp[-1]:
                    gene_id=re.findall("gene_id \"([A-Za-z0-9_\.\-\:/\ ()]+)\";",tmp[-1])[0]
                    gene_names=re.findall("gene_name \"([A-Za-z0-9_\.\-\:/\ ()]+)\"",tmp[-1])
                    genetype=re.findall("gene_type \"([A-Za-z0-9_\.\-\:/\ ()+]+)\"",tmp[-1])
                    biotype=re.findall("gene_biotype \"([A-Za-z0-9_\.\-\:/\ ()+]+)\"",tmp[-1])
                    if len(genetype) !=0 :
                        gene_biotype = genetype[0]
                    elif len(biotype) !=0 :
                        gene_biotype = biotype[0]
                    else:
                        logger.info(gene_id,gene_names)
                        continue
                    if len(gene_names)==0:
                          gene_name = gene_id
                    else:
                          gene_name = gene_names[0]
                    
                    # gene_table.append([ gene_id, gene_name ])
                    if gene_biotype == 'lincRNA' or gene_biotype == 'antisense' or gene_biotype == "lncRNA":
                        biotype="lnc"
                    else:
                        biotype="coding"
                    genetype_table.append([gene_id, biotype])
    return  genetype_table
