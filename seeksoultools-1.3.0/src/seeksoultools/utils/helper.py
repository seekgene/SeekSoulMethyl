import os
import sys
import re
import json
from pathlib import Path
from collections import defaultdict
from loguru import logger
import multiprocessing as mp
import functools
from xopen import xopen
from cutadapt.adapters import FrontAdapter, BackAdapter, Sequence, RightmostFrontAdapter
from typing import List, Optional
import importlib
import subprocess
import shutil
from ._version import __version__

# setup logger
logger.remove()
if os.environ.get("seeksoultools_debug", "False")=="True":
    logger.info("debug: True")
    logger.add(sys.stdout, level="DEBUG")
else:
    logger.info("debug: False")
    logger.add(sys.stdout, level="INFO")

logger.add(lambda _: sys.exit(1), level="ERROR")
def include_introns_callback(ctx, param:str, value:bool):
    # --include-introns
    if value:
        return "transcript"
    else:
        return "exon"

def str2path(ctx, param:str, value):
    try:
        if value:
            if isinstance(value, tuple):
                value = (Path(p).expanduser().resolve() for p in value)
            elif isinstance(value, list):
                value = [Path(p).expanduser().resolve() for p in value]
            else:
                value = Path(value).expanduser().resolve()
        else:
            return None
    except Exception as e:
        logger.error(f"{param}")
        logger.error(f"{e}")
        sys.exit(1)
    return value
def abs_path_callback(ctx, param:str, value:bool):
    return os.path.abspath(value)

def get_file_path(filename:str):
    return os.path.dirname(os.path.abspath(filename))

def get_utils_path(filename:str):
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(filename))), "utils")

def hamming_distance(s1, s2):
    return len([(i, j) for i, j in zip(s1, s2) if i != j])

def read_gtf(gtf):
    
    gene_list = []
    mt_regex = re.compile("^(MT|mt|Mt)-")
    # gene_id_regex = re.compile("gene_id \"([A-Za-z0-9_\.\-\:/\ ()]+)\"") # '\";'  ???
    # gene_name_regex = re.compile("gene_name \"([A-Za-z0-9_\.\-\:/\ ()]+)\"")
    gene_id_regex = re.compile(r'gene_id "(.*?)";')
    gene_name_regex = re.compile(r'gene_name "(.*?)";')
    id_names_dict = {}
    with xopen(gtf) as fh:
        for line in fh:
            if not line.strip(): continue
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")

            # gene，transcript and exon 
            if tmp[2] in ("gene", "transcript", "exon"):
                gene_id = gene_id_regex.search(tmp[-1])
                if not gene_id:
                    continue
                else:
                    gene_id = gene_id.group(1)

                # keep order?
                if gene_id not in id_names_dict:
                    gene_list.append(gene_id)
                else:
                    if re.sub(r"^MT-", "", id_names_dict[gene_id]) != gene_id:
                        continue
                gene_name = gene_name_regex.search(tmp[-1])

                if gene_name:
                    gene_name = gene_name.group(1)
                else:
                    gene_name = gene_id

                if tmp[0] in ("chrM", "MT", "mt"):
                    if not mt_regex.match(gene_name):
                        gene_name = "MT-" + gene_name
                id_names_dict[gene_id] = gene_name
    
    return [[g, id_names_dict[g]] for g in gene_list]

def parse_structure(string:str) -> tuple:
    regex = re.compile(r'([BLUXT])(\d+)')
    groups = regex.findall(string)
    return tuple([(_[0], int(_[1])) for _ in groups])

def read_file(file_list: list) -> dict:
    wl_dict = dict()
    for i, wl_file in enumerate(file_list):
        white_list = set()
        with xopen(wl_file, "r") as fh:
            for l in fh:
                if l.startswith("#"):
                    continue
                la = l.strip()
                if not la:
                    continue
                white_list.add(la)
        wl_dict[i] = white_list
    return wl_dict


def get_new_bc(bc:str, white_list:set, distance:int)->set:

    if distance == 1:
        BASE_LIST = ["T", "C", "G", "A"]
        mm_dict = dict()
        for i, c in enumerate(bc):
            if c == "N":
                mm_dict = { bc[:i] + base + bc[i+1:]:f"{i}{base}" for base in BASE_LIST }
                break  
            else:
                mm_dict.update({ bc[:i] + base + bc[i+1:]:f"{i}{base}" for base in BASE_LIST if base!=c })
                
        bc_set = set(mm_dict.keys()).intersection(white_list)
        # return {k: mm_dict[k] for k in bc_set}
    else:
        bc_dict = defaultdict(set)
        for bc_true in white_list:
            hmm = sum(ch1 != ch2 for ch1,ch2 in zip(bc_true,bc))
            if hmm <= distance:
                bc_dict[hmm].add(bc_true)
                
        bc_set = set()
        if len(bc_dict) != 0:
            sorted_items = sorted(bc_dict.items(), key=lambda x: x[0])
            bc_set = sorted_items[0][1]
    
    return bc_set

class AdapterFilter:
    def __init__(self, adapter1:list=[], adapter2:list=[],):
        self.adapter1 = [BackAdapter(sequence=_) if p=="3" else RightmostFrontAdapter(sequence=_) for _, p in adapter1]
        self.adapter2 = [BackAdapter(sequence=_) if p=="3" else RightmostFrontAdapter(sequence=_) for _, p in adapter2]
    
    def filter(self, r1=None, r2=None) -> tuple:
        flag = False
        if r1 and self.adapter1:
            # print("######")
            # print(len(r1))
            # print(r1.sequence)
            for _ in self.adapter1:
                m = _.match_to(r1.sequence)
                if m:
                    flag = True
                    r1 =  m.trimmed(r1)
                    # print(r1.sequence)
        if r2 and self.adapter2:
            for _ in self.adapter2:
                m = _.match_to(r2.sequence)
                if m:
                    flag = True
                    r2 =  m.trimmed(r2)
        return flag, r1, r2

class QcStat:
    def __init__(self):
        self.data = { }

    def update(self, **d):
        if not self.data:
            self.data = d
        else:
            for k, v in d.items():
                self.data[k] += v

    @staticmethod
    def _sort_gc(d):
        idx_max = max([k[0] for k in d])
        return {
            b: [d.get((i, b), 0) for i in range(idx_max+1)] for b in 'ATCGN'
        }

    @staticmethod
    def _sort_q(d, phred=33):
        idx_max = max([k[0] for k in d])
        q_max = max([ord(k[1])-phred for k in d])
        return {
            i: [d.get((i, chr(q+phred)), 0) for q in range(q_max+1)] for i in range(idx_max+1)
        }

    def save(self, path='summary.json'):
        tmp = {'__version__': __version__}
        for k in self.data:
            if k.endswith('_gc'):
                tmp[k] = self._sort_gc(self.data[k])
            elif k.endswith('_q'):
                tmp[k] = self._sort_q(self.data[k])
            else:
                tmp[k] = dict(self.data[k])
        with open(path, 'w') as fh:
            json.dump(tmp, fh, indent=4)

def check_cores(
    requested_cores: int):
    available_cores = mp.cpu_count()
    if requested_cores > available_cores:
        warning_msg = (
            f"System has {available_cores} CPU cores available, but {requested_cores} cores were requested. "
            f"Automatically adjusted to use {available_cores} cores."
        )
        logger.warning(warning_msg)
        return available_cores
    
    return requested_cores
def validate_cores(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        if 'core' in kwargs:
            kwargs['core'] = check_cores(kwargs['core'])
        return f(*args, **kwargs)
    return wrapper

def check_dependencies(
    python_packages: Optional[List[str]] = None,
    r_packages: Optional[List[str]] = None
) -> None:
    missing_py = []
    missing_r = []
    
    if python_packages:
        for package in python_packages:
            try:
                importlib.import_module(package)
            except ImportError:
                missing_py.append(package)
    
    if r_packages:
        r_check_cmd = """
        Rscript -e 'installed <- as.character(installed.packages()[,"Package"]); \
        cat(setdiff(c("%s"), installed), sep="\\n")'
        """ % '","'.join(r_packages)
        
        try:
            result = subprocess.check_output(r_check_cmd, shell=True).decode().strip().split('\n')
            missing_r = [x for x in result if x]  
        except subprocess.CalledProcessError:
            logger.error("Unable to check R package, please make sure R is installed correctly")
    
    if missing_py or missing_r:
        error_msg = []
        if missing_py:
            error_msg.append(f"loss package of Python: {', '.join(missing_py)}")
        if missing_r:
            error_msg.append(f"loss package of R: {', '.join(missing_r)}")
        raise click.ClickException('\n'.join(error_msg))
    success_msg = []
    if python_packages:
        success_msg.append(f"Python package check completed ({len(python_packages)} packages)")
    if r_packages:
        success_msg.append(f"R package check completed ({len(r_packages)} packages)")
    logger.info("✓ " + ", ".join(success_msg) + "  All dependency packages have been checked")

def check_external_tools(tools_config: dict) -> List[str]:
    missing_tools = []
    for tool, test_arg in tools_config.items():
        if not shutil.which(tool):
            missing_tools.append(tool)
            logger.error(f"Software not found: {tool}")
            continue
        try:
            result = subprocess.run(
                [tool, test_arg],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=10,
                text=True 
            )
            output = result.stdout if result.stdout else result.stderr
            output = output.strip().split('\n')
            version_line = ''
            for line in output:
                if any(v in line.lower() for v in ['version', 'v.']):
                    version_line = line
                    break
            output = version_line if version_line else output[0]
            if result.returncode == 0 or output:
                logger.info(f"{tool}: {output}")
            else:
                missing_tools.append(tool)
                logger.error(f" The software {tool}  command fails to be executed")
                
        except subprocess.TimeoutExpired:
            missing_tools.append(tool)
            logger.error(f"{tool} The software command timed out")
        except (subprocess.SubprocessError, OSError) as e:
            missing_tools.append(tool)
            logger.error(f"{tool} execution error: {str(e)}")

    return missing_tools
def check_all():
    REQUIRED_PYTHON_PACKAGES = [
    'numpy', 'pandas', 'scanpy', 'pathlib', 'json', 'dnaio', 'collections', 'shutil', 'Bio', 'collections', 'scanpy', 'dandelion', 'warnings', 'jinja2', 'loguru', 'multiprocessing', 'functools', 'xopen','cutadapt.adapters'
    ]

    REQUIRED_R_PACKAGES = [
    'Seurat', 'ggplot2', 'dplyr', 'DropletUtils', 'Matrix', 'edgeR'
    ]

    check_dependencies(
    python_packages=REQUIRED_PYTHON_PACKAGES,
    r_packages=REQUIRED_R_PACKAGES
    )
    ### check tools
    REQUIRED_TOOLS = {
    'qualimap': '--help',
    'STAR': '--version',
    'samtools': '--version',
    'run-trust4': ''
    }
    check_external_tools(REQUIRED_TOOLS)
