import os
import re
import gzip
from collections import defaultdict
from itertools import groupby, combinations
import pysam
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from intervaltree import IntervalTree
from xopen import xopen
from .helper import logger, read_gtf, hamming_distance

def count(bam, outdir, gtf, umi_correct_method, **kwargs):

    detail_file = os.path.join(outdir, "detail.xls")
    counts_file = os.path.join(outdir, "counts.xls")
    umi_file = os.path.join(outdir, "umi.xls")

    logger.info(f"umi correction method selected: {umi_correct_method}")
    gtf_tree = create_interval_tree(gtf)
    bam2table(
        bam=bam,
        detail_file=detail_file,
        counts_file=counts_file,
        umi_correct_detail=umi_file,
        umi_correct_method=umi_correct_method,
        gtf_tree=gtf_tree
    )

    raw_matrix_dir = os.path.join(outdir, "raw_feature_bc_matrix")
    os.makedirs(raw_matrix_dir, exist_ok=True)
    logger.info("write raw matrix started!")
    write_raw_matrix(counts_file, raw_matrix_dir, gtf)
    logger.info("write raw matrix done!")
    return 


def bam2table(bam, detail_file, counts_file, umi_correct_detail, umi_correct_method, gtf_tree):

    default_verbosity = pysam.set_verbosity(0)
    sam_file = pysam.AlignmentFile(bam, "rb")
    pysam.set_verbosity(default_verbosity)

    umi_clusterer = UMIClusterer(umi_correct_method)

    umi_correct_detail_fh = open(umi_correct_detail, "w")
    with open(detail_file, "w") as fh1, open(counts_file, "w") as fh2:
        fh1.write("\t".join(["cellID", "geneID", "UMI", "Num"]) + "\n")
        fh2.write("\t".join(["cellID", "geneID", "UMINum", "ReadsNum"]) + "\n")
        for barcode, g in groupby(sam_file, key=lambda x: x.qname.split("_", 1)[0]):
            counts_dict, geneid_umi_dict = umi_count(g, umi_correct_detail_fh, umi_clusterer, gtf_tree)
            for gene_id in geneid_umi_dict:
                for umi in geneid_umi_dict[gene_id]:
                    raw_umi_count = geneid_umi_dict[gene_id][umi]
                    fh1.write(f"{barcode}\t{gene_id}\t{decode(umi)}\t{raw_umi_count}\n")
            for gene_id in counts_dict:
                umi_num, reads_num = counts_dict[gene_id]
                fh2.write(f"{barcode}\t{gene_id}\t{umi_num}\t{reads_num}\n")
    umi_correct_detail_fh.close()

def multiple_alignment(reads, gene_ids, gtf_tree):

    assignment = set()
    for read in reads:
        chrm, pos = read.reference_name, read.pos
        all_matches = gtf_tree[chrm][pos]
        for match in all_matches:
            match_dict =  match.data
            if match_dict["gene_id"] in gene_ids and ("exon_id" in match_dict.keys() or match_dict["type"] == "exon"):
                assignment.add(match_dict["gene_id"])

    if len(assignment) == 1:
        return assignment.pop()
    else:
        return None

def umi_count(reads_group, umi_correct_detail_fh, umi_clusterer, gtf_tree, ispair=False):
    assigned_dict = defaultdict(lambda: defaultdict(int))
    for barcode_umi, g in groupby(reads_group, key=lambda x: x.qname):
        _, umi = barcode_umi.split("_")[:2]
        if umi == umi[0]*len(umi):   # poly
            continue
        umi = encode(umi)
        tmp_dict = defaultdict(list)
        n = 0
        for r in g:
            if ispair and r.is_read1:
                continue
            n += 1
            if r.has_tag("XT"):
                XT =r.get_tag("XT").split("XT:Z:")[-1]
                if "," in XT:
                    gene_id = multiple_alignment([r], XT.split(","), gtf_tree)
                    if not gene_id:
                        continue
                else:
                    gene_id = XT.split("XT:Z:")[-1]
                tmp_dict[gene_id].append({"mapq": r.mapping_quality,
                                          "read": r})
        if len(tmp_dict) == 1:
            gene_id = list(tmp_dict.keys())[0]
            if tmp_dict[gene_id][0]["mapq"] == 255 or n > 1:
                assigned_dict[gene_id][umi] += 1
        elif len(tmp_dict) > 1:
            gene_ids = list(tmp_dict.keys())
            reads = [v['read'] for key in tmp_dict.keys() for v in tmp_dict[key]]
            gene_id_true = multiple_alignment(reads, gene_ids, gtf_tree)
            if gene_id_true:
                assigned_dict[gene_id_true][umi] += 1

    counts_dict, geneid_umi_dict = umi_correct(assigned_dict,
                                               _,
                                              umi_correct_detail_fh,
                                              umi_clusterer)
    return counts_dict, geneid_umi_dict

def umi_correct(geneid_umi_dict, bc, umi_correct_detail_fh, umi_clusterer):
    counts_dict = defaultdict(lambda: [0, 0])
    for gene_id, umi_dict in geneid_umi_dict.items():

        final_umis = umi_clusterer(umi_dict, 1)
        corrected_dict = create_map_to_correct_umi(final_umis)

        for ori_umi, new_umi in corrected_dict.items():
            if ori_umi != new_umi:
                umi_dict[new_umi] += umi_dict[ori_umi]
                umi_correct_detail_fh.write(f"{bc}\t{gene_id}\t{decode(ori_umi)}:{umi_dict[ori_umi]}\t{decode(new_umi)}:{umi_dict[new_umi]}\n")
                del umi_dict[ori_umi]

        counts_dict[gene_id][0] = len(umi_dict.keys())
        counts_dict[gene_id][1] = sum(umi_dict.values())
    return counts_dict, geneid_umi_dict


def write_raw_matrix(counts_file, raw_matrix_dir, gtf):
    name_df = pd.DataFrame(read_gtf(gtf), columns=["geneID", "Symbol"])

    gene_dict = name_df.reset_index().set_index("geneID")["index"].to_dict()
    row, col, data = [], [], []
    barcodes = []
    n = 0
    with open(counts_file) as fh:
        fh.readline()
        for k, g in groupby(fh, lambda x: x.split("\t")[0]):
            barcodes.append(k)
            for _ in g:
                tmp = _.split("\t")
                data.append(int(tmp[2]))
                row.append(gene_dict[tmp[1]])
                col.append(n)
            n += 1
    mat = coo_matrix((data, (row, col)), shape=(name_df.shape[0], n))

    matrix_file = os.path.join(raw_matrix_dir, "matrix.mtx.gz")
    with gzip.open(matrix_file, "w") as fh:
        mmwrite(fh, mat)

    name_df["type"] = "Gene Expression"
    features_file = os.path.join(raw_matrix_dir, "features.tsv.gz")
    name_df.to_csv(features_file, sep="\t", index=False, header=False)

    barcodes_file = os.path.join(raw_matrix_dir, "barcodes.tsv.gz")
    with gzip.open(barcodes_file, "wt") as fh:
        fh.write("\n".join(barcodes))
        fh.write("\n")

def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched

def remove_umis(adj_list, cluster, nodes):
    '''removes the specified nodes from the cluster and returns
    the remaining nodes '''

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove

def get_substr_slices(umi_length, idx_size):
    '''
    Create slices to split a UMI into approximately equal size substrings
    Returns a list of tuples that can be passed to slice function
    '''
    cs, r = divmod(umi_length, idx_size)
    sub_sizes = [cs + 1] * r + [cs] * (idx_size - r)
    offset = 0
    slices = []
    for s in sub_sizes:
        slices.append((offset, offset + s))
        offset += s
    return slices

def build_substr_idx(umis, umi_length, min_edit):
    '''
    Build a dictionary of nearest neighbours using substrings, can be used
    to reduce the number of pairwise comparisons.
    '''
    substr_idx = defaultdict(
        lambda: defaultdict(set))
    slices = get_substr_slices(umi_length, min_edit + 1)
    for idx in slices:
        for u in umis:
            u_sub = u[slice(*idx)]
            substr_idx[idx][u_sub].add(u)
    return substr_idx

def iter_nearest_neighbours(umis, substr_idx):
    '''
    Added by Matt 06/05/17
    use substring dict to get (approximately) all the nearest neighbours to
    each in a set of umis.
    '''
    for i, u in enumerate(umis, 1):
        neighbours = set()
        for idx, substr_map in substr_idx.items():
            u_sub = u[slice(*idx)]
            neighbours = neighbours.union(substr_map[u_sub])
        neighbours.difference_update(umis[:i])
        for nbr in neighbours:
            yield u, nbr

class UMIClusterer:

    def __init__(self, cluster_method="adjacency"):

        if cluster_method == "adjacency":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_adjacency

        elif cluster_method == "directional":
            self.get_adj_list = self._get_adj_list_directional
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_directional

        elif cluster_method == "cluster":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_cluster

    def __call__(self, umis, threshold):
        '''umis is a dictionary that maps UMIs to their counts'''

        counts = umis
        umis = list(umis.keys())

        adj_list = self.get_adj_list(umis, counts, threshold)
        clusters = self.get_connected_components(umis, adj_list, counts)
        final_umis = [list(x) for x in
                      self.get_groups(clusters, adj_list, counts)]

        return final_umis

    def _get_best_min_account(self, cluster, adj_list, counts):
        ''' return the min UMI(s) need to account for cluster'''
        if len(cluster) == 1:
            return list(cluster)

        sorted_nodes = sorted(cluster, key=lambda x: (-counts[x], x))

        for i in range(len(sorted_nodes) - 1):
            if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
                return sorted_nodes[:i+1]

    def _get_adj_list_adjacency(self, umis, counts, threshold):
        ''' identify all umis within hamming distance threshold'''

        adj_list = {umi: [] for umi in umis}
        if len(umis) > 25:
            umi_length = len(umis[0])
            substr_idx = build_substr_idx(umis, umi_length, threshold)
            iter_umi_pairs = iter_nearest_neighbours(umis, substr_idx)
        else:
            iter_umi_pairs = combinations(umis, 2)
        for umi1, umi2 in iter_umi_pairs:
            if hamming_distance(umi1, umi2) <= threshold:
                adj_list[umi1].append(umi2)
                adj_list[umi2].append(umi1)

        return adj_list

    def _get_adj_list_directional(self, umis, counts, threshold=1):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (2 * second umi counts)-1'''

        adj_list = {umi: [] for umi in umis}
        if len(umis) > 25:
            umi_length = len(umis[0])
            substr_idx = build_substr_idx(umis, umi_length, threshold)
            iter_umi_pairs = iter_nearest_neighbours(umis, substr_idx)
        else:
            iter_umi_pairs = combinations(umis, 2)
        for umi1, umi2 in iter_umi_pairs:
            if hamming_distance(umi1, umi2) <= threshold:
                if counts[umi1] >= (counts[umi2]*2)-1:
                    adj_list[umi1].append(umi2)
                if counts[umi2] >= (counts[umi1]*2)-1:
                    adj_list[umi2].append(umi1)

        return adj_list

    def _get_connected_components_adjacency(self, umis, graph, counts):
        ''' find the connected UMIs within an adjacency dictionary'''

        # TS: TO DO: Work out why recursive function doesn't lead to same
        # final output. Then uncomment below

        # if len(graph) < 10000:
        #    self.search = breadth_first_search_recursive
        # else:
        #    self.search = breadth_first_search

        found = set()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                # component = self.search(node, graph)
                component = breadth_first_search(node, graph)
                found.update(component)
                components.append(component)
        return components

    def _group_adjacency(self, clusters, adj_list, counts):
        ''' return groups for adjacency method'''

        groups = []

        for cluster in clusters:
            if len(cluster) == 1:
                groups.append(list(cluster))

            else:
                observed = set()

                lead_umis = self._get_best_min_account(cluster, adj_list, counts)
                observed.update(lead_umis)

                for lead_umi in lead_umis:
                    connected_nodes = set(adj_list[lead_umi])
                    groups.append([lead_umi] + list(connected_nodes - observed))
                    observed.update(connected_nodes)

        return groups

    def _group_directional(self, clusters, adj_list, counts):
        ''' return groups for directional method'''

        observed = set()
        groups = []
        for cluster in clusters:
            if len(cluster) == 1:
                groups.append(list(cluster))
                observed.update(cluster)
            else:
                cluster = sorted(cluster, key=lambda x: counts[x],
                                 reverse=True)
                # need to remove any node which has already been observed
                temp_cluster = []
                for node in cluster:
                    if node not in observed:
                        temp_cluster.append(node)
                        observed.add(node)
                groups.append(temp_cluster)

        return groups

    def _group_cluster(self, clusters, adj_list, counts):
        ''' return groups for cluster or directional methods'''

        groups = []
        for cluster in clusters:
            groups.append(sorted(cluster, key=lambda x: counts[x],
                                 reverse=True))

        return groups

def create_map_to_correct_umi(cluster_list):
    """Create map to correct umi."""
    my_map = {encode(y): encode(x[0]) for x in cluster_list for y in x}
    return my_map


def decode(umi):
    try:
        new_umi = umi.decode()
    except (UnicodeDecodeError, AttributeError):
        new_umi = umi
    return new_umi


def encode(umi):
    try:
        new_umi = str.encode(umi)
    except (TypeError):
        new_umi = umi
    return new_umi

def create_interval_tree(gtf):

    gtf_tree = defaultdict(lambda: IntervalTree())
    gene_id_regex = re.compile("gene_id \"([A-Za-z0-9_\.\-\:/\ ()]+)\"")
    transcript_id_regex = re.compile("transcript_id \"([A-Za-z0-9_\.\-\:/\ ()]+)\"")
    exon_id_regex = re.compile("exon_id \"([A-Za-z0-9_\.\-\:/\ ()]+)\"")
    with xopen(gtf) as gtf_r:
        for line in gtf_r:

            if line.startswith("#"):
                continue
            if not line.strip():
                continue

            le = line.strip().split("\t")
            _chr, ge_type, start_pos, end_pos, gene_str = le[0], le[2], int(le[3]), int(le[4]), le[-1]
            if ge_type not in ("gene", "transcript", "exon"):
                continue

            if start_pos == end_pos:
                start_pos = start_pos -1

            info_dict = {}
            gene_id = gene_id_regex.search(gene_str)
            if gene_id:
                info_dict["gene_id"] = gene_id.group(1)

            transcript_id = transcript_id_regex.search(gene_str)
            if transcript_id:
                info_dict["transcript_id"] = transcript_id.group(1)

            exon_id = exon_id_regex.search(gene_str)
            if exon_id:
                info_dict["exon_id"] = exon_id.group(1)

            info_dict.update(type = ge_type)
            gtf_tree[_chr][start_pos:end_pos] = info_dict

    return gtf_tree

def calulate_chunck_size(detail_file):
    chunck_size = 5000000
    with open(detail_file) as f:
        line_count = sum(1 for _ in f)

    n_split = int(line_count/chunck_size)
    if n_split == 0:
        n_split = 1

    chunck_size = int(line_count/n_split) + 1

    return chunck_size

def get_int_type(max_value):
    if max_value <= 127:
        return np.int8
    elif max_value <= 32767:
        return np.int16
    elif max_value <= 2147483647:
        return np.int32
    else:
        return np.int64

def calculate_metrics(counts_file, detail_file, filterd_barcodes_file, filterd_features_file, gtf, basedir):
    summary = defaultdict()
    summary["Estimated Number of Cells"] = 0
    summary["Fraction Reads in Cells"] = 0
    summary["Mean Reads per Cell"] = 0
    summary["Median Genes per Cell"] = 0
    summary["Median UMI Counts per Cell"] = 0
    summary["Total Genes Detected"] = 0

    barcodes = pd.read_csv(filterd_barcodes_file, header=None, sep="\t")
    summary["Estimated Number of Cells"] = barcodes.shape[0]

    df0 = pd.read_csv(
        counts_file,
        dtype={
            "cellID": "category",
            "geneID": "category",
            "UMINum": "int32",
            "ReadsNum": "int32"
        },
        sep="\t"
    )
    summary["Sequencing Saturation"] = 1 - df0.UMINum.sum()/df0.ReadsNum.sum()
    df0 = df0.loc[df0["cellID"].isin(barcodes[0]), :].reset_index(drop=True)

    # lnc
    from ..utils.fastUtil import read_gtf
    lnclist=['lincRNA','antisenseq','lnc','lncRNA', 'lnc_RNA']
    type_df = pd.DataFrame(read_gtf(gtf), columns=["gene", "type"])
    lnc_df = type_df.loc[type_df["type"].isin(lnclist), :].reset_index(drop=True)
    if lnc_df.size == 0:
        gene_median_lnc = 0
    else:
        lncgene=lnc_df.loc[: ,'gene']       
        df0_lnc = df0.loc[df0['geneID'].isin(lncgene), :].reset_index(drop=True)
        if df0_lnc.size ==0 : gene_median_lnc = 0
        else:
            gene_median_lnc = int(df0_lnc.groupby(['cellID'], observed=True)['geneID'].nunique().median())
            del df0_lnc

    umi_median = int(df0.groupby(["cellID"], observed=True)["UMINum"].sum().median())
    summary["Median UMI Counts per Cell"] = umi_median
    gene_total = int(df0[["geneID"]].nunique())
    summary["Total Genes Detected"] = gene_total
    gene_median = int(df0.groupby(["cellID"], observed=True)["geneID"].nunique().median())
    summary["Median Genes per Cell"] = gene_median
    summary['Median lnc Genes per Cell'] = gene_median_lnc
    del df0

    chunck_size = calulate_chunck_size(detail_file)
    csv_reader = pd.read_csv(
                        detail_file,
                        dtype={
                            "cellID": "category",
                            "geneID": "category",
                            "UMI": "category",
                            "Num": "int32"
                        },
                        sep="\t",
                        chunksize=chunck_size)

    saturation_tmp = defaultdict(lambda: defaultdict(int))
    median_tmp = defaultdict(list)

    chunk_count = 0
    max_val = 0
    dtypes = np.int16

    mapped_reads_total = 0
    cell_reads_total = 0
    for df in csv_reader:
        mapped_reads_total += df["Num"].sum()
        df = df.loc[df["cellID"].isin(barcodes[0]), :].reset_index(drop=True)
        cell_reads_total += df["Num"].sum()
        rep = df["Num"]
        df = df.drop(["Num"], axis=1)
        idx = df.index.repeat(rep)
        df = df.iloc[idx].reset_index(drop=True)
        del rep, idx

        # shuffle
        df = df.sample(frac=1.0).reset_index(drop=True)

        # downsample
        n_cols_key = [str((i+1)/ 10) for i in range(0,10)]
        for n, interval in enumerate(np.array_split(np.arange(df.shape[0]), 10)):
            idx = interval[-1]
            percentage = n_cols_key[n]
            tmp_file = "tmp_" + str(chunk_count) + "_" + percentage + ".xls"
            tmp_file = os.path.join(basedir, tmp_file)

            # calculate saturation for each portion
            sampled_df = df.iloc[:idx]
            sampled_df = sampled_df.assign(**{percentage: 1})
            sampled_df = sampled_df.groupby(['cellID', 'geneID', 'UMI'], observed=True) \
                                   .sum() \
                                   .reset_index()
            np.savetxt(tmp_file, sampled_df[percentage].to_numpy(), fmt='%d')
            saturation_tmp[percentage][tmp_file] = idx

            # calculate median for each portion
            median = sampled_df.groupby([sampled_df["cellID"]],observed=True)["geneID"] \
                               .nunique() \
                               .reset_index(drop=True) \
                               .median()
            median_tmp[percentage].append(int(median))

            # refreshing int dtype
            max_curr = sampled_df[percentage].max()
            if max_val < max_curr:
                dtypes = get_int_type(max_curr)
                max_val = max_curr

        chunk_count +=1

    summary["Fraction Reads in Cells"] = cell_reads_total/mapped_reads_total

    percentage_sampling = []
    saturation_sampling = []
    median_sampling = []
    for perc, files in saturation_tmp.items():
        arr_perc = np.array([], dtype=dtypes)
        all_obs = 0
        for file_path, count in files.items():
            arr = np.loadtxt(file_path, dtype=dtypes)
            arr_perc = np.append(arr_perc, arr)
            all_obs += count
            os.remove(file_path)
        saturation = (np.sum(arr_perc[arr_perc>1] -1) - 0.0)/ all_obs * 100
        saturation_sampling.append(saturation)
    for perc, medians in median_tmp.items():
        median = int(sum(medians)/len(medians))
        median_sampling.append(median)
        percentage_sampling.append(float(perc))

    return summary, {
                      "percentage": percentage_sampling,
                      "saturation": saturation_sampling,
                      "median": median_sampling
    }
