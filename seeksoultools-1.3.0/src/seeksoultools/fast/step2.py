import os
import json
import linecache
from collections import defaultdict
from ..utils.helper import logger

__srcdir = os.path.dirname(os.path.abspath(__file__))


def align(ofq,
          genomeDir,
          gtf,
          fra,
          samplename,
          outdir,
          sc5p,
          region,
          core=4,
          star_path="STAR",
          samtools_path="samtools",
          **kwargs):
    if ("steps" not in kwargs) or (not kwargs["steps"]):
        kwargs["steps"] = ["STAR", "SortByPos", "FeatureCounts", "SortByName"]
    basedir = os.path.join(outdir, "step2")
    STAR_dir = os.path.join(basedir, 'STAR')
    rRNAdir = os.path.join(STAR_dir, "rRNA")
    os.makedirs(rRNAdir, exist_ok=True)
    logger.info('rRNA started')
    from ..utils.fastUtil import rRNA_map
    # scoremin = kwargs['scoremin'] if  kwargs['scoremin'] else 0.66
    # matchnmin = kwargs['matchnmin'] if kwargs['matchnmin'] else 0.66
    scoremin = kwargs['scoremin']
    matchnmin = kwargs['matchnmin']
    if kwargs["rRNAgenomeDir"] != None and kwargs["rRNAgtf"] != None:
        mtrRNA, rRNA = rRNA_map(ofq, kwargs["rRNAgenomeDir"], kwargs["rRNAgtf"], rRNAdir, samplename, scoremin=scoremin, matchnmin=matchnmin, core=core, star_path=star_path)
    else:
        mtrRNA, rRNA = rRNA_map(ofq, genomeDir, gtf, rRNAdir, samplename,scoremin=scoremin, matchnmin=matchnmin, core=core, star_path=star_path)
    logger.info('rRNA done')

    prefix = os.path.join(STAR_dir, samplename + '_')
    if 'STAR' not in kwargs['steps']:
        logger.info('STAR skiped!')
    else:
        logger.info('STAR started!')
        from ..utils.wrappers import STAR_wrapper
        bam, STARLog = STAR_wrapper(fq = ofq,
                                core=core,
                                genomeDir=genomeDir,
                                prefix=prefix,
                                scoremin=scoremin,
                                matchnmin=matchnmin,
                                star_path=star_path)
        logger.info('STAR done!')
    bam, STARLog =f'{prefix}Aligned.out.bam', f'{prefix}Log.final.out'

	    # sort by pos
    if 'SortByPos' not in kwargs['steps']:
        logger.info('SortByPos skiped!')
    else:
        logger.info('SortByPos started!')
        from ..utils.wrappers import samtools_sort_wrapper
        bam = samtools_sort_wrapper(
                    inbam=bam,
                    outbam=f'{prefix}SortedByCoordinate.bam',
                    byname=False,
                    core=core,
                    clean=True,
                    samtools_path=samtools_path)
        logger.info('SortByPos done!')
	
    bam = f'{prefix}SortedByCoordinate.bam'
    downdir= os.path.join(STAR_dir,"downbam")
    os.makedirs(downdir, exist_ok=True)
    logger.info('downsample bam started!')
    from ..utils.wrappers import picard_DownsampleSam_wrapper
    downprefix=os.path.join(downdir, samplename)
    downsambam = picard_DownsampleSam_wrapper(bam, fra, f'{downprefix}.down.{fra}.bam')
    downsambam = f'{downprefix}.down.{fra}.bam'
    logger.info('downsample bam done!')

    logger.info('genebody started!')

    from ..utils.wrappers import geneBody_coverage_wrapper
    genebodycoverage = geneBody_coverage_wrapper(downsambam, downdir, samplename, gtf)
    logger.info('genebody done!')


    logger.info('ACTB coverage started!')
    from ..utils.fastUtil import actb
    actb02mean = actb(gtf, bam)
    logger.info('ACTB coverage done!')

    with open(os.path.join(outdir, f'{samplename}_summary.json')) as fh:
        refpath=os.path.dirname(genomeDir.rstrip("/"))
        reffile=os.path.join(refpath,'reference.json')
        if os.path.exists(reffile):
                with open(reffile) as refjson:
                    refj=json.load(refjson)
                    genome=refj['genomes'][0]
        else:
                #genome=genomeDir
                genome=refpath.split("/")[-1]
        summary = json.load(fh)
        summary['reference'] = genome
        if region=="transcript":
            summary['include_introns'] =  "True"
        else:
            summary['include_introns'] =  "False"
        Total = summary['stat']['total']

    # sort by name
    if "SortByName" not in kwargs["steps"]:
        logger.info("SortByName skiped!")
    else:
        logger.info("SortByName started!")
        from ..utils.wrappers import samtools_sort_wrapper
        bam = samtools_sort_wrapper(
            bam,
            f"{prefix}SortedByName.bam",
            core=core,
            byname=True,
        )
        logger.info("SortByName done!")

    logger.info('run qualimap started!')
    from ..utils.wrappers import qualimap_wrapper
    RnaSeqMetrics = qualimap_wrapper(bam=bam, gtf=gtf, outdir=STAR_dir, SC5P=sc5p, **{"-s": "", "-pe": ""})
    RnaSeqMetrics = os.path.join(STAR_dir, "rnaseq_qc_results.txt")
    logger.info('run qualimap done!')

    summary_tmp = defaultdict()
    from ..utils.fastUtil import mapping_summary
    tmp = mapping_summary(STARLog, RnaSeqMetrics)
    Total = tmp['Number of input reads']
    mapped_genome_ratio = tmp['reads aligned']/Total
    summary_tmp['Reads Mapped to Genome'] = mapped_genome_ratio
    
    if not kwargs['scoremin'] and not kwargs['matchnmin']:
        summary_tmp["rRNA% in Mapped"] = rRNA
        summary_tmp["Mt_rRNA% in Mapped"] = mtrRNA
        summary_tmp["Reads Mapped to Middle Genebody"] = genebodycoverage
        summary_tmp['Fraction over 0.2 mean coverage depth of ACTB gene'] = actb02mean
    mapped_confident_ratio = (tmp['aligned to genes'] + tmp['no feature assigned']) /2/ Total
    summary_tmp['Reads Mapped Confidently to Genome'] = mapped_confident_ratio
    mapped_intergenic_ratio = tmp['intergenic'] /2.0 / Total
    summary_tmp['Reads Mapped to Intergenic Regions'] = mapped_intergenic_ratio
    
    mapped_intronic_ratio = tmp['intronic'] / 2.0/ Total
    summary_tmp['Reads Mapped to Intronic Regions'] = mapped_intronic_ratio
    mapped_exonic_ratio = tmp['exonic'] / 2.0/ Total
    summary_tmp['Reads Mapped to Exonic Regions'] = mapped_exonic_ratio
    with open(os.path.join(outdir, f'{samplename}_summary.json'), 'w') as fh:
        summary['mapping'] = summary_tmp
        json.dump(summary, fh, indent=4)

    stringtie_dir = os.path.join(basedir, 'stringtie')
    if 'stringtie' not in kwargs['steps']:
        logger.info('skig stringtie!')
    else:
        os.makedirs(stringtie_dir, exist_ok=True)

    featureCounts_dir = os.path.join(basedir, 'featureCounts')
    # FeatureCounts
    if 'FeatureCounts' not in kwargs['steps']:
        logger.info('run_featureCounts done!')
    else:
        os.makedirs(featureCounts_dir, exist_ok=True)
        logger.info('run featureCounts started!')
        from ..utils.wrappers import featureCounts_wrapper
        # "--countReadPairs": ""
        #
        ## region:transcript, exon

        # check gtf
        bam = featureCounts_wrapper(bam=bam,
                                samplename=samplename,
                                outdir=featureCounts_dir,
                                gtf=gtf,
                                SC5P=sc5p,
                                region=region,
                                core=core,
                                **{"-p": "", })
        logger.info('run featureCounts done!')


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
    bam = os.path.join(featureCounts_dir, f"{samplename}_SortedByName.bam")
