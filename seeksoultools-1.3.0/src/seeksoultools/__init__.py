import os
import importlib
import click
from seeksoultools.rna import rna
from seeksoultools.vdj import vdj
from seeksoultools.multivdj import multivdj
from seeksoultools.fast import fast
from seeksoultools.utils import utils
from seeksoultools.utils._version import __version__
from seeksoultools.utils import helper
import psutil
from pathlib import Path 
from loguru import logger
import subprocess


def check_disk_space(path: Path) -> None:
    disk = psutil.disk_usage(str(path))
    free_gb = disk.free / (1024 * 1024 * 1024 * 1024 )
    logger.info(f"Available space: {free_gb:.1f}T")
   
    if free_gb < 0.1: 
        logger.warning(f"Warning: Available disk space is less than {free_gb:.1f}T")
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(version=__version__, message="%(version)s")
@click.option("--debug", is_flag=True, default=False, help="debug flag.")
@click.pass_context
def cli(ctx, debug):
    ctx.ensure_object(dict)
    ctx.obj["version"] = __version__
    if debug:
        os.environ["seeksoultools_debug"] = "True"
    else:
        os.environ["seeksoultools_debug"] = "False"   
    current_dir = Path.cwd()
    check_disk_space(current_dir)

cli.add_command(rna)
cli.add_command(fast)
cli.add_command(vdj)
cli.add_command(multivdj)
cli.add_command(utils)
