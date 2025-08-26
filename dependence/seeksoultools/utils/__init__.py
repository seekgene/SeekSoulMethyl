import click

@click.group(help="utils.")
@click.pass_context
def utils(ctx):
    pass

@utils.command(help="gene type summary.")
@click.option('--gtf', type=click.Path(), required=True, help="gtf file.")
@click.option('--feature', default='gene', show_default=True, help="feature, e.g., gene, transcript")
@click.option('--key', default='gene_type', show_default=True, help='attribution key, e.g., gene_type')
def gtfstat(gtf, feature, key='gene_biotype'):
    from .mkref import gtfstat
    gtfstat(gtf, feature, key)

@utils.command(help="gtf filter.")
@click.option('--gtf', type=click.Path(), required=True, help="gtf file.")
@click.option('--biotype',  multiple=True, help="biotype to keep,  can specify multiple times.")
@click.option('--key', default='gene_type', show_default=True, help='attribution key, e.g., gene_type')
def gtffilter(gtf, biotype, key='gene_type'):
    from .mkref import gtffilter
    gtffilter(gtf, biotype, key)

@utils.command(help="make ref.")
@click.option('--fa', type=click.Path(), required=True, help="fa file.")
@click.option('--gtf', type=click.Path(), required=True, help="gtf file.")
@click.option('--genomeDir', 'genomeDir', type=click.Path(), required=True, help="genomeDir.")
@click.option('--runThreadN', 'runThreadN', default=8, help="runThreadN.")
@click.option('--star_path', default='STAR', help="STAR path.")
@click.option('--star_opt', default='', help="STAR opt.")
def mkref(fa, gtf, genomeDir, runThreadN=8, star_path='STAR', star_opt=''):
    from .mkref import mkref
    mkref(fa, gtf, genomeDir, runThreadN, star_path, star_opt)

@utils.command(help="addtag.")
@click.option("--inbam", type=click.Path(), required=True, help="inbam file, barcode and umi include in reads header.")
@click.option("--outbam", type=click.Path(), required=True, help="outbam file.")
@click.option("--umifile", type=click.Path(), help="umi.xls in step3.")
@click.option("--recursive", "recursive", is_flag=True, default=False, show_default=True, 
              help=(f"whether use recursive to find raw umi or not,"
                f"set when the software version generating the umi.xls is lower than 1.1")
            )
def addtag(inbam, outbam, umifile=None, recursive=False):
    from .addtag import add_tag
    add_tag(inbam, outbam, umifile, recursive)

@utils.command(help="split reference.")
@click.option("-f", "--fa", type=click.Path(), required=True, help="fasta file.")
@click.option("-g", "--gtf", type=click.Path(), required=True, help="gtf file.")
@click.option("-o", "--outdir", required=True, help="output dir.")
def splitref(fa, gtf, outdir):
    from .splitref import split_ref
    split_ref(fa, gtf, outdir)

@utils.command(help="gtf validator")
@click.option("--gtf", "gtf", type=click.Path(), required=True, help="gtf path.")
@click.option("--newgtf", "new_gtf", default=None, help="modified gtf file. if not set, only do check.")
@click.option("--gid", "gene_id_str", default="gene_id", show_default=True, help="attribute to extract gene id.")
def validate_gtf(gtf, new_gtf:str=None, gene_id_str:str="gene_id"):
    from .gtfvalidator import validator
    validator(gtf, new_gtf, gene_id_str)

@utils.command(help="multi report.")
@click.option("--outdir", type=click.Path(), required=True, help="output dir.")
@click.option("--samplename", required=True, help="samplename.")
@click.option("--rnadir", type=click.Path(), help="rna root dir.")
@click.option("--tcrdir", type=click.Path(), help="tcr root dir,  xxxx/data or xxxx/outs.")
@click.option("--bcrdir", type=click.Path(), help="bcr root dir,  xxxx/data or xxxx/outs.")
def multireport(outdir, samplename, rnadir=None, tcrdir=None, bcrdir=None):
    from .multi_report import report
    report(outdir, samplename, rnadir, tcrdir, bcrdir)
