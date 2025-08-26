from pathlib import Path
import json
import dnaio
from collections import defaultdict, OrderedDict
import pandas as pd
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from ..utils.wrappers import search_pattern_wrapper

def cal_median(s:pd.Series):
    return s.median() if s.shape[0]>0 else 0
def get_most_frequent_longest_seq(sequences):
    max_length = len(max(sequences, key=len))
    longest_seqs = [seq for seq in sequences if len(seq) == max_length]
    if len(longest_seqs) == 1:
        return longest_seqs[0]
    seq_counts = Counter(longest_seqs)
    return seq_counts.most_common(1)[0][0]

def process_airr_to_consensus(airr_file, output_fa):
    df = pd.read_csv(airr_file, sep='\t')

    exact_1_df = df[df['exact_subclonotype_id'] == 1]

    consensus_sequences = []
    for clone_id in exact_1_df['clone_id'].unique():
        clone_data = exact_1_df[exact_1_df['clone_id'] == clone_id]

        for idx, locus in enumerate(clone_data['locus'].unique(), 1):
            locus_data = clone_data[clone_data['locus'] == locus]

            sequences = locus_data['sequence'].tolist()
            consensus_seq = get_most_frequent_longest_seq(sequences)
            seq_record = SeqRecord(
                Seq(consensus_seq),
                id=f"{clone_id}_consensus_{idx}",
                description=""
            )
            consensus_sequences.append(seq_record)

    SeqIO.write(consensus_sequences, output_fa, "fasta")

def cal_clonotype(barcode_clonotype_tsv:Path, clonotype_file:Path):
    d = defaultdict(lambda: defaultdict(list))
    total = 0
    with barcode_clonotype_tsv.open() as fh:
        header_dict = {k:i for i, k in enumerate(fh.readline().strip().split("\t"))}
        if "locus_VDJ" not in header_dict:
            header_dict["locus_VDJ"] = len(header_dict)
        if "locus_VJ" not in header_dict:
            header_dict["locus_VJ"] = len(header_dict)
        if "junction_VDJ" not in header_dict:
            header_dict["junction_VDJ"] = len(header_dict)
        if "junction_VJ" not in header_dict:
            header_dict["junction_VJ"] = len(header_dict)
        if "junction_aa_VDJ" not in header_dict:
            header_dict["junction_aa_VDJ"] = len(header_dict)
        if "junction_aa_VJ" not in header_dict:
            header_dict["junction_aa_VJ"] = len(header_dict)
            
        for line in fh:
            tmp = line.strip().split("\t")
            while len(tmp) <= max(header_dict.values()):
                tmp.append("None")
                
            cdr3s_nt = ""
            cdr3s_aa = ""
            cell_id = tmp[header_dict["cell_id"]]
            total += 1
            clone_id = tmp[header_dict["clone_id"]]
            heavy_chain = tmp[header_dict["locus_VDJ"]]
            if heavy_chain!="None":
                heavy_chain_nt = tmp[header_dict["junction_VDJ"]]
                heavy_chain_aa = tmp[header_dict["junction_aa_VDJ"]]
                cdr3s_nt += f"{heavy_chain}:{heavy_chain_nt};"
                cdr3s_aa += f"{heavy_chain}:{heavy_chain_aa};"
            
            light_chain = tmp[header_dict["locus_VJ"]]
            if light_chain!="None":
                light_chain_nt = tmp[header_dict["junction_VJ"]]
                light_chain_aa = tmp[header_dict["junction_aa_VJ"]]
                cdr3s_nt += f"{light_chain}:{light_chain_nt};"
                cdr3s_aa += f"{light_chain}:{light_chain_aa};"
            cdr3s_aa = cdr3s_aa.strip(";")
            cdr3s_nt = cdr3s_nt.strip(";")
            d[clone_id][cdr3s_nt].append([cell_id, cdr3s_aa])

    clonotype_list = []
    for clone_id, clone_data in d.items():
        clone_data = sorted(clone_data.items(), key=lambda x: len(x[1]), reverse=True)
        paired_chains = []
        single_chains = []
        for cdr3s_nt, cell_data_list in clone_data:
            chains_count = len(cdr3s_nt.split(';'))
            if chains_count == 2:  
                paired_chains.append((cdr3s_nt, cell_data_list))
            else:  
                single_chains.append((cdr3s_nt, cell_data_list))
        
        if paired_chains:
            cdr3s_nt = paired_chains[0][0]
            cdr3s_aa = paired_chains[0][1][0][1]
        else:
            cdr3s_nt = clone_data[0][0]
            cdr3s_aa = clone_data[0][1][0][1]

        frequency = sum([len(x[1]) for x in clone_data])
        proportion = frequency/total if total>0 else 0
        clonotype_list.append([clone_id, frequency, proportion, cdr3s_aa, cdr3s_nt])
    with clonotype_file.open(mode="w") as fh:
        fh.write("clonotype_id\tfrequency\tproportion\tcdr3s_aa\tcdr3s_nt\n")
        clonotype_list = sorted(clonotype_list, key=lambda x: int(x[0].replace("clonotype", "")), reverse=False)     
        for (clone_id, frequency, proportion, cdr3s_aa, cdr3s_nt) in clonotype_list:
            fh.write(f"{clone_id}\t{frequency}\t{proportion}\t{cdr3s_aa}\t{cdr3s_nt}\n")

def cal_paired_clonotype_diversity(clonotypes_tsv:str):
    df = pd.read_csv(clonotypes_tsv, sep="\t")
    try:
        sum_proportion_squared = (df['proportion']**2).sum()
        if sum_proportion_squared == 0:
            return 0 
        return 1/sum_proportion_squared
    except (ZeroDivisionError, ValueError):
        return 0 

def get_clone_matrix(airr_file:Path, meta_file:Path, chain:str):
    ddl_meta = pd.read_table(meta_file)
    if "productive_VDJ" not in ddl_meta.columns:
        ddl_meta["productive_VDJ"] = "None"

    if "productive_VJ" not in ddl_meta.columns:
        ddl_meta["productive_VJ"] = "None"
        
    if "locus_VDJ" not in ddl_meta.columns:
        ddl_meta["locus_VDJ"] = "None"
        
    if "locus_VJ" not in ddl_meta.columns:
        ddl_meta["locus_VJ"] = "None"
        
    if "umi_count_B_VDJ" not in ddl_meta.columns:
        ddl_meta["umi_count_B_VDJ"] = 0
        
    if "umi_count_B_VJ" not in ddl_meta.columns:
        ddl_meta["umi_count_B_VJ"] = 0
        
    if "umi_count_abT_VDJ" not in ddl_meta.columns:
        ddl_meta["umi_count_abT_VDJ"] = 0
        
    if "umi_count_abT_VJ" not in ddl_meta.columns:
        ddl_meta["umi_count_abT_VJ"] = 0

    d = {}
    cellnum = ddl_meta.shape[0]
    d["Estimated Number of Cells"] = cellnum
    
    if chain=="IG":
        # productive_VJ vs productive_B_VJ?
        v_j_pair_df = ddl_meta.loc[(ddl_meta.productive_VDJ =="T") & (ddl_meta.productive_VJ=="T"), :]
        d["Number of Cells with Productive V-J Spanning Pair"] = v_j_pair_df.shape[0]
        d["Cells with Productive V-J Spanning Pair"] = v_j_pair_df.shape[0]/cellnum if cellnum>0 else 0

        igk_igh_pair_df = ddl_meta.loc[(ddl_meta.locus_VJ=="IGK") & (ddl_meta.locus_VDJ=="IGH"), :]
        d["Cells with Productive V-J Spanning (IGK, IGH) Pair"] = igk_igh_pair_df.shape[0]/cellnum if cellnum>0 else 0

        igl_igk_pair_df = ddl_meta.loc[(ddl_meta.locus_VJ=="IGL") & (ddl_meta.locus_VDJ=="IGH"), :]
        d["Cells with Productive V-J Spanning (IGL, IGH) Pair"] = igl_igk_pair_df.shape[0]/cellnum if cellnum>0 else 0

        igh_df = ddl_meta.loc[(ddl_meta.productive_VDJ =="T") & (ddl_meta.locus_VDJ=="IGH"), :]
        d["Cells with Productive IGH Contig"] = igh_df.shape[0]/cellnum if cellnum>0 else 0
        d["Median IGH UMIs per Cell"] = cal_median(igh_df.umi_count_B_VDJ)

        igk_df = ddl_meta.loc[(ddl_meta.productive_VJ =="T") & (ddl_meta.locus_VJ=="IGK"), :]
        d["Cells with Productive IGK Contig"] = igk_df.shape[0]/cellnum if cellnum>0 else 0
        d["Median IGK UMIs per Cell"] = cal_median(igk_df.umi_count_B_VJ)
        
        igl_df = ddl_meta.loc[(ddl_meta.productive_VJ =="T") & (ddl_meta.locus_VJ=="IGL"), :]
        d["Cells with Productive IGL Contig"] = igl_df.shape[0]/cellnum if cellnum>0 else 0
        d["Median IGL UMIs per Cell"] = cal_median(igl_df.umi_count_B_VJ)

    elif chain=="TR":
        v_j_pair_df = ddl_meta.loc[(ddl_meta.productive_VDJ =="T") & (ddl_meta.productive_VJ=="T"), :]
        d["Number of Cells with Productive V-J Spanning Pair"] = v_j_pair_df.shape[0]
        d["Cells with Productive V-J Spanning Pair"] = v_j_pair_df.shape[0]/cellnum if cellnum>0 else 0

        tra_trb_pair_df = ddl_meta.loc[(ddl_meta.locus_VDJ=="TRB") & (ddl_meta.locus_VJ=="TRA"), :]
        d["Cells with Productive V-J Spanning (TRA, TRB) Pair"] = tra_trb_pair_df.shape[0]/cellnum if cellnum>0 else 0


        tra_df = ddl_meta.loc[(ddl_meta.productive_VDJ =="T") & (ddl_meta.locus_VJ=="TRA"), :]
        d["Cells with Productive TRA Contig"] = tra_df.shape[0]/cellnum if cellnum>0 else 0
        d["Median TRA UMIs per Cell"] = cal_median(tra_df.umi_count_abT_VJ)

        trb_df = ddl_meta.loc[(ddl_meta.productive_VJ =="T") & (ddl_meta.locus_VDJ=="TRB"), :]
        d["Cells with Productive TRB Contig"] = trb_df.shape[0]/cellnum if cellnum>0 else 0
        d["Median TRB UMIs per Cell"] = cal_median(trb_df.umi_count_abT_VDJ)

    else:
        raise ValueError(f"{chain} is not supported, should be one of IG, TR.")

    return d

def get_summary(samplename:str, wd:Path, chain:str, leader_ref_file:Path|None=None, outdir:Path|None=None, pa:str="Outs"):
    """
    Get summary of the pipeline
    """
    if not outdir:
        outdir = wd
    summary_json = wd / "Analysis" / "summary.json"
    with summary_json.open() as fh:
        summary = json.load(fh)

    total_reads = summary["Number of Read Pairs"]
    valid_reads = summary["Number of Valid Reads"]
    candidate_reads = summary["Number of Candidate Reads"]
    full_candidate_reads = candidate_reads + summary["Number of Filtered Candidate Reads"]

    if (wd/"Analysis"/"trust4"/f"_{samplename}_bc_count.tsv").exists():
        df_full = pd.read_table(wd/"Analysis"/"trust4"/f"_{samplename}_bc_count.tsv")
    else:
        df_full = pd.read_table(wd/"Analysis"/"trust4"/f"{samplename}_bc_count.tsv")
    barcodes = pd.read_table(outdir/"Analysis"/"dandelion"/f"{samplename}_barcode.tsv")["cell_id"]
    cell_associated_reads = df_full.loc[df_full.barcode.isin(barcodes), "reads_num"].sum()

    summary["Fraction Reads in Cells"] = cell_associated_reads/full_candidate_reads if full_candidate_reads > 0 else 0

    d = get_clone_matrix(outdir/"Analysis"/"dandelion"/f"{samplename}_airr_rearrangement.tsv", outdir/"Analysis"/"dandelion"/f"{samplename}_clone_meta.tsv", chain)
    summary.update(d)

    metrics = OrderedDict()
    metrics['Samplename'] = samplename
    cell_number = summary['Estimated Number of Cells']
    metrics['Estimated Number of Cells'] = f"{cell_number:,}"

    metrics['Mean Read Pairs per Cell'] = f"{int(valid_reads/cell_number):,}" if cell_number > 0 else "0"
    metrics['Mean Used Read Pairs per Cell'] = f"{int(candidate_reads/cell_number):,}" if cell_number > 0 else "0"
    metrics['Fraction Reads in Cells'] = f"{(summary['Fraction Reads in Cells']):.2%}" 

    metrics['Number of Cells with Productive V-J Spanning Pair'] = f"{summary['Number of Cells with Productive V-J Spanning Pair']:,}"
    metrics['Number of Reads'] = f"{total_reads:,}"
    metrics['Valid Barcodes'] = f"{summary['Valid Barcodes']:.2%}"
    metrics['Q30 Bases in Barcode'] = f"{summary['Q30 Bases in Barcode']:.2%}"
    metrics['Q30 Bases in R2 Read'] = f"{summary['Q30 Bases in R2 Read']:.2%}"
    metrics['Q30 Bases in UMI'] = f"{summary['Q30 Bases in UMI']:.2%}"

    if chain=="IG":
        metrics['Reads Mapped to Any V(D)J Gene'] = f"{summary['Reads Mapped to Any V(D)J Gene']:.2%}"
        metrics['Reads Mapped to IGH'] = f"{summary['Reads Mapped to IGH']:.2%}"
        metrics['Reads Mapped to IGK'] = f"{summary['Reads Mapped to IGK']:.2%}"
        metrics['Reads Mapped to IGL'] = f"{summary['Reads Mapped to IGL']:.2%}"

        metrics['Median IGH UMIs per Cell'] = f"{summary['Median IGH UMIs per Cell']:,}"
        metrics['Median IGK UMIs per Cell'] = f"{summary['Median IGK UMIs per Cell']:,}"
        metrics['Median IGL UMIs per Cell'] = f"{summary['Median IGL UMIs per Cell']:,}"
        metrics['Cells with Productive V-J Spanning Pair'] = f"{summary['Cells with Productive V-J Spanning Pair']:.2%}"
        metrics['Cells with Productive V-J Spanning (IGK, IGH) Pair'] = f"{summary['Cells with Productive V-J Spanning (IGK, IGH) Pair']:.2%}"
        metrics['Cells with Productive V-J Spanning (IGL, IGH) Pair'] = f"{summary['Cells with Productive V-J Spanning (IGL, IGH) Pair']:.2%}"
        # metrics['Paired Clonotype Diversity'] = f"{1/((clonotypes['proportion']**2).sum()):,.2}"
        metrics['Cells with Productive IGH Contig'] = f"{summary['Cells with Productive IGH Contig']:.2%}"
        metrics['Cells with Productive IGK Contig'] = f"{summary['Cells with Productive IGK Contig']:.2%}"
        metrics['Cells with Productive IGL Contig'] = f"{summary['Cells with Productive IGL Contig']:.2%}"
    elif chain=="TR":
        metrics['Reads Mapped to Any V(D)J Gene'] = f"{summary['Reads Mapped to Any V(D)J Gene']:.2%}"
        metrics['Reads Mapped to TRA'] = f"{summary['Reads Mapped to TRA']:.2%}"
        metrics['Reads Mapped to TRB'] = f"{summary['Reads Mapped to TRB']:.2%}"
        metrics['Median TRA UMIs per Cell'] = f"{summary['Median TRA UMIs per Cell']:,}"
        metrics['Median TRB UMIs per Cell'] = f"{summary['Median TRB UMIs per Cell']:,}"
        metrics['Cells with Productive V-J Spanning Pair'] = f"{summary['Cells with Productive V-J Spanning Pair']:.2%}"
        metrics['Cells with Productive V-J Spanning (TRA, TRB) Pair'] = f"{summary['Cells with Productive V-J Spanning (TRA, TRB) Pair']:.2%}"
        metrics['Cells with Productive TRA Contig'] = f"{summary['Cells with Productive TRA Contig']:.2%}"
        metrics['Cells with Productive TRB Contig'] = f"{summary['Cells with Productive TRB Contig']:.2%}"

    else:
        raise ValueError(f"{chain} is not supported, should be one of IG, TR.")
    outs_dir = outdir/pa
    outs_dir.mkdir(exist_ok=True, parents=True)

    airr_file = f"{outdir}/Analysis/dandelion/{samplename}_airr_rearrangement.tsv"
    shutil.copy(airr_file, outs_dir)
    barcode_clonotype_file = outdir/"Analysis"/"dandelion"/f"{samplename}_barcode.tsv"
    #shutil.copy(barcode_clonotype_file, outs_dir)

    clonotypes_file = outdir/pa/f"{samplename}_clonotypes.tsv"
    barcode_json = outdir/pa/f"{samplename}_cell_barcodes.json"
    contig_fa = outdir/pa/f"{samplename}_filtered_contig.fasta"
    anno_file = outdir/pa/f"{samplename}_filtered_contig_annotations.tsv"
    cal_clonotype(barcode_clonotype_file, clonotypes_file)
    with open(barcode_json, 'w') as f:
        json.dump(pd.read_table(airr_file)['cell_id'].unique().tolist() , f, indent=4)
    with dnaio.open(contig_fa, fileformat="fasta", mode='w') as fa:
        for idx, row in pd.read_table(airr_file).iterrows():
           seq = dnaio.Sequence(row['sequence_id'], row['sequence'])
           fa.write(seq)
    #### consensus
    conses_fa = contig_fa = outdir/pa/f"{samplename}_consensus.fasta"
    process_airr_to_consensus(airr_file, conses_fa)

    summary['Paired Clonotype Diversity'] = cal_paired_clonotype_diversity(clonotypes_file)
    metrics['Paired Clonotype Diversity'] = f"{summary['Paired Clonotype Diversity']:.2f}"

    if leader_ref_file:
        annot_fa = wd/"Analysis"/"trust4"/f"{samplename}_annot.fa"
        search_pattern_wrapper(
            leader_ref=leader_ref_file,
            airr_tsv=airr_file,
            annot=annot_fa,
            outfile=outs_dir / f"{samplename}_filtered_contig_annotations.tsv"
        )
    from .utils import anno_cr
    from ..utils.report.websummaryVDJ import  has_valid_data
    if has_valid_data(airr_file):
        anno_cr(airr_file, anno_file)
    json.dump(
        summary, 
        open(outdir/"Analysis"/"summary.json", "w"),
        indent=4,
    )

    with (outs_dir/"metrics_summary.csv").open(mode="w") as fh: 
        keys = [f'"{e}"' if ',' in e else e for e in  metrics.keys()]
        fh.write(f"{','.join(keys)}\n")

        values =[f'"{e}"' if ',' in e else e for e in  metrics.values()]
        fh.write(f"{','.join(values)}\n")
