from pathlib import Path
from collections import defaultdict
import scanpy as sc
import dandelion as ddl
import numpy as np
import pandas as pd
import re
import os
import warnings
from ..utils.helper import logger
from ..utils.wrappers import calculate_n50
warnings.filterwarnings("ignore")

def clean_and_read_table(files):
    try:
        with open(files, 'rb') as f:
            content = f.read()
        text = content.decode('latin1') 
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
            cleaned_text = re.sub(r'[^\x00-\x7F]+', 'N', text)
            temp_file.write(cleaned_text)
            temp_path = temp_file.name

        df = pd.read_table(temp_path)
        os.unlink(temp_path)

        return df
    except Exception as e:
        print(f"the file is not correct {files}: {str(e)}")
        raise


def count_reads(final_out:Path):
    d = defaultdict(int)
    with final_out.open() as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            contig_id = header.replace(">", "").split()[0]
            _seq = fh.readline()
            list_a = fh.readline().strip().split()
            list_c = fh.readline().strip().split()
            list_g = fh.readline().strip().split()
            list_t = fh.readline().strip().split()
            pos_sum = [int(a) + int(c) + int(g) + int(t) for (a, c, g, t) in zip(list_a, list_c, list_g, list_t)]
            for i in range(len(pos_sum)):
                if i == 0:
                    total_reads = pos_sum[i]
                else:
                    n = pos_sum[i] - pos_sum[i-1]
                    if n > 0:
                        total_reads += n
            d[contig_id] = total_reads
    return d

def filter_single_chain_by_umi_threshold(vdj):
    data = vdj.data
    complete_clones = data[data['clone_id'].str.contains('VDJ.*VJ', na=False)]
    vj_only_contig = data[data['clone_id'].str.contains('^.*VJ_', na=False) &
                                ~data['clone_id'].str.contains('VDJ', na=False)]
    vdj_only_contig = data[data['clone_id'].str.contains('^.*VDJ_', na=False) &
                                 ~data['clone_id'].str.contains('_VJ_', na=False)]
    contigs_to_remove = []
    for _, contig in vj_only_contig.iterrows():
        vj_number =  contig['clone_id'].split('VJ_')[1]
        matching_complete = complete_clones[
                    complete_clones['clone_id'].str.contains(f'_VJ_{vj_number}', na=False)]
        matching_vj = vj_only_contig[
            vj_only_contig['clone_id'].str.contains(f'VJ_{vj_number}', na=False)
        ]
        matches = pd.concat([matching_complete, matching_vj])
        if len(matches) > 0:
            umi_counts = {i: count for i, count in enumerate(matches['umi_count'].values)}
            n50 = calculate_n50(umi_counts)  
            threshold = n50 * 0.01
            if contig['umi_count'] < threshold:
                contigs_to_remove.append(contig['sequence_id'])
    for _, contig in vdj_only_contig.iterrows():
        vdj_number = contig['clone_id'].split('VDJ_')[1]

        matching_complete = complete_clones[
            complete_clones['clone_id'].str.contains(f'VDJ_{vdj_number}', na=False)]
        matching_vj = vdj_only_contig[
            vdj_only_contig['clone_id'].str.contains(f'VDJ_{vdj_number}', na=False)
        ]
        matches = pd.concat([matching_complete, matching_vj])

        if len(matches) > 0:
            umi_counts = {i: count for i, count in enumerate(matches['umi_count'].values)}
            n50 = calculate_n50(umi_counts)  
            threshold = n50 * 0.01
            if contig['umi_count'] < threshold:
                contigs_to_remove.append(contig['sequence_id'])
    aa = pd.DataFrame(contigs_to_remove, columns=["contig"])
    cells_to_keep = vdj.data[~vdj.data['sequence_id'].isin(contigs_to_remove)]['cell_id'].unique()
    vdj.data = vdj.data[~vdj.data['sequence_id'].isin(contigs_to_remove)]
    vdj.metadata = vdj.metadata[vdj.metadata.index.isin(cells_to_keep)]
    #vdj.metadata = vdj.metadata[~vdj.metadata.index.isin(contigs_to_remove)]
    
    return vdj


def update_clone_id_by_chain_match(vdj, new_airr_file, meta_file,  barcode_file):
    try:
        metadata = vdj.metadata
        complete_clones = metadata[metadata['clone_id'].str.contains('VDJ.*VJ', na=False)]
        vj_only_clones = metadata[metadata['clone_id'].str.contains('^.*VJ_', na=False) &
                                ~metadata['clone_id'].str.contains('VDJ', na=False)]
        vdj_only_clones = metadata[metadata['clone_id'].str.contains('^.*VDJ_', na=False) &
                                 ~metadata['clone_id'].str.contains('_VJ_', na=False)]

        updates = []
        clone_size_map = defaultdict(int)
        for idx, row in metadata.iterrows():
            clone_size = int(row['clone_id_by_size']) if pd.notnull(row['clone_id_by_size']) else 1
            clone_size_map[row['clone_id_by_size']] = max(clone_size_map[row['clone_id_by_size']], clone_size)


        for idx, clone in vj_only_clones.iterrows():
            try:
                vj_number = clone['clone_id'].split('VJ_')[1]

                matching_complete = complete_clones[
                    complete_clones['clone_id'].str.contains(f'_VJ_{vj_number}', na=False)
                ] #.drop_duplicates(subset=['clone_id'])
                if len(set(matching_complete['clone_id'])) == 1:
                    best_match_id = matching_complete['clone_id'].iloc[0]
                    best_match = complete_clones[complete_clones['clone_id'] == best_match_id].iloc[0]
                    updates.append({
                                'index': idx,
                                'old_clone_id': clone['clone_id'],
                                'old_clone_id_by_size': clone['clone_id_by_size'],
                                'new_clone_id_by_size': best_match['clone_id_by_size'],
                                'matched_complete_clone': best_match['clone_id'],
                                'matched_clone_size': matching_complete['clone_id_by_size'].iloc[0],
                                'total_matches': 1
                                 })
                elif len(matching_complete['clone_id']) > 1:
                    clone_sizes = matching_complete['clone_id'].value_counts()
                    sorted_clone_sizes = clone_sizes.sort_values(ascending=False)
                    best_match_id = sorted_clone_sizes.index[0]
                    second_best_size = sorted_clone_sizes.iloc[1]
                    largest_size = sorted_clone_sizes.iloc[0]
                    if second_best_size / largest_size  < 0.15:
                        best_match = complete_clones[complete_clones['clone_id'] == best_match_id].iloc[0]
                        updates.append({
                                'index': idx,
                                'old_clone_id': clone['clone_id'],
                                'old_clone_id_by_size': clone['clone_id_by_size'],
                                'new_clone_id_by_size': best_match['clone_id_by_size'],
                                'matched_complete_clone': best_match['clone_id'],
                                'matched_clone_size': largest_size,
                                'total_matches': len(matching_complete['clone_id'].unique()) 
                                 })
            except Exception as e:
                continue

        for idx, clone in vdj_only_clones.iterrows():
            try:
                vdj_number = clone['clone_id'].split('VDJ_')[1]
                matching_complete = complete_clones[
                    complete_clones['clone_id'].str.contains(f'VDJ_{vdj_number}', na=False)
                ]#.drop_duplicates(subset=['clone_id'])
                if len(set(matching_complete['clone_id'])) == 1:
                    best_match_id = matching_complete['clone_id'].iloc[0]
                    best_match = complete_clones[complete_clones['clone_id'] == best_match_id].iloc[0]
                    updates.append({
                                'index': idx,
                                'old_clone_id': clone['clone_id'],
                                'old_clone_id_by_size': clone['clone_id_by_size'],
                                'new_clone_id_by_size': best_match['clone_id_by_size'],
                                'matched_complete_clone': best_match['clone_id'],
                                'matched_clone_size': matching_complete['clone_id_by_size'].iloc[0],
                                'total_matches': 1
                                 })
                elif len(matching_complete['clone_id']) > 1:
                    clone_sizes = matching_complete['clone_id'].value_counts()
                    sorted_clone_sizes = clone_sizes.sort_values(ascending=False)
                    best_match_id = sorted_clone_sizes.index[0]
                    second_best_size = sorted_clone_sizes.iloc[1]
                    largest_size = sorted_clone_sizes.iloc[0]
                    if second_best_size / largest_size  < 0.15:
                        best_match = complete_clones[complete_clones['clone_id'] == best_match_id].iloc[0]
                        updates.append({
                                'index': idx,
                                'old_clone_id': clone['clone_id'],
                                'old_clone_id_by_size': clone['clone_id_by_size'],
                                'new_clone_id_by_size': best_match['clone_id_by_size'],
                                'matched_complete_clone': best_match['clone_id'],
                                'matched_clone_size': largest_size,
                                'total_matches': len(matching_complete['clone_id'].unique())
                                })
            except Exception as e:
                continue
        for update in updates:
            metadata.loc[update['index'], 'clone_id_by_size'] = update['new_clone_id_by_size']
        vdj.metadata.update(metadata)
        return vdj

    except Exception as e:
        return vdj

def assign_exact_subclon(group):
    unique_clone_ids = group['clone_id'].unique()
    if len(unique_clone_ids) == 1:
        return pd.Series(1, index=group.index)
    clone_counts = group['clone_id'].value_counts()
    type_weights = {}
    for clone_id in unique_clone_ids:
        if 'VDJ' in clone_id and 'VJ' in clone_id:
            type_weights[clone_id] = 3
        elif 'VDJ' in clone_id:
            type_weights[clone_id] = 2
        else:
            type_weights[clone_id] = 1
    sorted_clones = sorted(
        clone_counts.items(),
        key=lambda x: (x[1], type_weights.get(x[0], 0), x[1]),
        reverse=True
    )

    id_mapping = {clone_id: i + 1 for i, (clone_id, _) in enumerate(sorted_clones)}

    return group['clone_id'].map(id_mapping)
def reassign_clone_id_by_size(vdj):
    if 'clone_id_by_size' not in vdj.metadata.columns:
        return vdj
    clone_counts = vdj.metadata['clone_id_by_size'].value_counts()
    new_ids = {old_id: new_id + 1 for new_id, (old_id, _) in 
               enumerate(clone_counts.items())}
    id_mapping = pd.DataFrame({
    'old_clone_id': list(new_ids.keys()),
    'new_clone_id': list(new_ids.values())
     })
    vdj.metadata['clone_id_by_size'] = vdj.metadata['clone_id_by_size'].map(new_ids)
    vdj.metadata['exact_subclonotype_id'] = vdj.metadata.groupby('clone_id_by_size').apply(
       lambda x: assign_exact_subclon(x)
       ).reset_index(level=0, drop=True)
    return vdj
def parse_gene_info(gene_str):
    if gene_str == '*':
        return None, None, None
    try:
        gene_parts = gene_str.split('(')[0]
        gene_name = gene_parts.split('*')[0]
        positions = gene_str.split(':')[1].strip('()')
        start, end = map(int, positions.split('-'))
        return gene_name, start, end
    except:
        return None, None, None

def add_star(df, annofile):
    airr_dict = df.set_index('sequence_id')[['v_call', 'd_call', 'j_call', 'c_call']].to_dict('index')
    """
    airr_dict = {}
    for _, row in df.iterrows():
    airr_dict[row['sequence_id']] = {
        'V': row['v_call'],
        'D': row['d_call'],
        'J': row['j_call'],
        'C': row['c_call']}
    """
    contig_positions = {}
    with annofile.open() as f:
        while True:
            header_line = f.readline().strip()
            if not header_line:
                break
            seq_line = f.readline().strip()
            if not header_line.startswith('>'):
                continue
            fields = header_line[1:].split()
            contig_id = fields[0].strip('>')
            if contig_id in airr_dict:
                v_genes = fields[3].split(',') if fields[3] != '*' else []
                d_genes = fields[4].split(',') if fields[4] != '*' else []
                j_genes = fields[5].split(',') if fields[5] != '*' else []
                c_genes = fields[6].split(',') if fields[6] != '*' else []

                positions = {
                'v_sequence_start': None, 'v_sequence_end': None,
                'd_sequence_start': None, 'd_sequence_end': None,
                'j_sequence_start': None, 'j_sequence_end': None,
                'c_sequence_start': None, 'c_sequence_end': None
                }

                for v in v_genes:
                    gene_name, start, end = parse_gene_info(v)
                    if  gene_name == df.loc[df['sequence_id'] == contig_id, 'v_call'].iloc[0]:
                        positions['v_sequence_start'] = start
                        positions['v_sequence_end'] = end
                        break
                for d in d_genes:
                    gene_name, start, end = parse_gene_info(d)
                    if  gene_name == df.loc[df['sequence_id'] == contig_id, 'd_call'].iloc[0]:
                        positions['d_sequence_start'] = start
                        positions['d_sequence_end'] = end
                        break
                for j in j_genes:
                    gene_name, start, end = parse_gene_info(j)
                    if gene_name == df.loc[df['sequence_id'] == contig_id, 'j_call'].iloc[0]:
                        positions['j_sequence_start'] = start
                        positions['j_sequence_end'] = end
                        break
                for c in c_genes:
                    gene_name, start, end = parse_gene_info(c)
                    if gene_name == df.loc[df['sequence_id'] == contig_id, 'c_call'].iloc[0]:
                        positions['c_sequence_start'] = start
                        positions['c_sequence_end'] = end
                        break
                contig_positions[contig_id] = positions
    
    new_columns = ['v_sequence_start', 'v_sequence_end', 
                   'd_sequence_start', 'd_sequence_end',
                   'j_sequence_start', 'j_sequence_end',
                   'c_sequence_start', 'c_sequence_end']
    
    for col in new_columns:
        df[col] = None
    
    for contig_id, positions in contig_positions.items():
        for col, value in positions.items():
            df.loc[df['sequence_id'] == contig_id, col] = value
    df['junction_length'] = df['junction'].str.len()
    df['junction_aa_length'] = df['junction_aa'].str.len()
    df['junction_length'] = df['junction_length'].fillna(0).astype(int)
    df['junction_aa_length'] = df['junction_aa_length'].fillna(0).astype(int)
    return df
def make_file(chain, ddl_wd, samplename):
            clonotypes_file = ddl_wd/f"{samplename}_clonotypes.tsv"
            pd.DataFrame(columns=[ 
                'clonotype_id', 'cdr3s_aa', 'cdr3s_nt', 'frequency', 'proportion'
            ]).to_csv(clonotypes_file, sep='\t', index=False)

            airr_file = ddl_wd/f"{samplename}_airr_rearrangement.tsv"
            pd.DataFrame(columns=[
            "sequence_id","sequence","rev_comp","productive","locus","v_call","d_call","j_call","c_call","sequence_alignment","germline_alignment","cdr1","cdr2","junction","junction_aa","v_cigar","d_cigar","j_cigar","c_cigar","v_identity","j_identity","cell_id","complete_vdj","consensus_count","duplicate_count","clone_id","exact_subclonotype_id","v_sequence_start","v_sequence_end","d_sequence_start","d_sequence_end","j_sequence_start","j_sequence_end","c_sequence_start","c_sequence_end","junction_length","junction_aa_length"
            ]).to_csv(airr_file, sep='\t', index=False)

            barcode_file = ddl_wd/f"{samplename}_barcode.tsv"
            clist = [
            "cell_id","clone_id","clone_id_by_size","locus_VDJ","locus_VJ","productive_VDJ","productive_VJ","v_call_VDJ","d_call_VDJ","j_call_VDJ","v_call_VJ","j_call_VJ","c_call_VDJ","c_call_VJ","junction_VDJ","junction_VJ","junction_aa_VDJ","junction_aa_VJ","v_call_abT_VDJ","d_call_abT_VDJ","j_call_abT_VDJ","v_call_abT_VJ","j_call_abT_VJ","c_call_abT_VDJ","c_call_abT_VJ","productive_abT_VDJ","productive_abT_VJ","umi_count_abT_VDJ","umi_count_abT_VJ","isotype","isotype_status","locus_status","chain_status","rearrangement_status_VDJ","rearrangement_status_VJ","exact_subclonotype_id"
            ]
            if chain == "IG":
                clist = [i.replace('abT', 'B') for i in clist]
            pd.DataFrame(columns=clist).to_csv(barcode_file, sep='\t', index=False)

            clone_meta_file = ddl_wd/f"{samplename}_clone_meta.tsv"
            mlist = ["cell_id","clone_id","clone_id_by_size","locus_VDJ","locus_VJ","productive_VDJ","productive_VJ","v_call_VDJ","d_call_VDJ",
            "j_call_VDJ","v_call_VJ","j_call_VJ","c_call_VDJ","c_call_VJ","junction_VDJ","junction_VJ","junction_aa_VDJ","junction_aa_VJ","v_call_abT_VDJ",
            "d_call_abT_VDJ","j_call_abT_VDJ","v_call_abT_VJ","j_call_abT_VJ","c_call_abT_VDJ","c_call_abT_VJ","productive_abT_VDJ","productive_abT_VJ",
            "umi_count_abT_VDJ","umi_count_abT_VJ","v_call_VDJ_main","v_call_VJ_main","d_call_VDJ_main","j_call_VDJ_main","j_call_VJ_main","c_call_VDJ_main","c_call_VJ_main","v_call_abT_VDJ_main",
            "d_call_abT_VDJ_main","j_call_abT_VDJ_main","v_call_abT_VJ_main","j_call_abT_VJ_main","isotype","isotype_status","locus_status","chain_status","rearrangement_status_VDJ","rearrangement_status_VJ","exact_subclonotype_id","_clone_id"]
            if chain == "IG":
                mlist = [i.replace('abT', 'B') for i in mlist]
            pd.DataFrame(columns=mlist).to_csv(clone_meta_file, sep='\t', index=False)

def complete_airr(
    samplename:str, raw_airr_file: Path, final_out:Path, annofile:Path, chain: str, ddl_wd:Path, matrix:Path|None, umi:int=3):

    new_airr_file_raw = ddl_wd/f"{samplename}_airr_raw.tsv" # update duplicate_count and consensus_count
    new_airr_file = ddl_wd/f"{samplename}_airr_rearrangement.tsv"
    tmpfile = ddl_wd/f"{samplename}.tmp.vdjdata.csv"
    meta_file = ddl_wd/f"{samplename}_clone_meta.tsv"
    barcode_file = ddl_wd/f"{samplename}_barcode.tsv"
    if not raw_airr_file.exists():
        raise FileNotFoundError(f"{raw_airr_file} does not exist")
    #df = pd.read_table(raw_airr_file)
    df = clean_and_read_table(raw_airr_file)

    # filter, TR,IG, other?
    if chain == 'TR':
        df = df.loc[df.locus.str.startswith('TR'), :].reset_index(drop=True)
    elif chain == 'IG':
        df = df.loc[df.locus.str.startswith('IG'), :].reset_index(drop=True)
    
    
    # add duplicate_count
    df["duplicate_count"] = df.consensus_count
    d = count_reads(final_out)
    df["consensus_count"] = df.sequence_id.map(d)

    # add cell_id
    df.cell_id = df.sequence_id.str.replace(r'_.*', '', regex=True)

    df['saturation'] = np.where(df['consensus_count'] != 0,  1 - df['duplicate_count']/df['consensus_count'], 0)
    df = df.loc[df.saturation > 0.94]
    df = df[df.groupby('cell_id')['duplicate_count'].transform('sum') >= umi]
    vdj = ddl.utl.load_data(df)
    logger.info(f"chain: {vdj.shape[0]}")
    # vdj = ddl.utl.load_data(df)
    ###print(vdj.metadata.shape if hasattr(vdj, 'metadata') else "No metadata")
    if not matrix:
        # find clones with dandelion.tl.find_clones
        
        # The default mode for calculation of junctional hamming distance is to use the CDR3 junction amino acid sequences, 
        # specified via the key option (None defaults to junction_aa). 
        # You can switch it to using CDR3 junction nucleotide sequences (key = 'junction'), 
        # or even the full V(D)J amino acid sequence (key = 'sequence_alignment_aa'), 
        # as long as the column name exists in the .data slot.
        try:
            vdj = ddl.tl.find_clones(vdj, key = "junction_aa", by_alleles=False, key_added="clone_id", identity={"ig": 0.85, "tr-ab": 1, "tr-gd": 1})
        except Exception as e:
            #logger.warning(f"Failed to find clones: {str(e)}")
            logger.warning(f"no clonotype are find, please check your data.The number of input contig {vdj.shape[0]}")
            make_file(chain, ddl_wd, samplename)
            return None
    else:
        adata = sc.read_10x_mtx(
            matrix,
            var_names='gene_symbols'
        )
        try:
            vdj, adata = ddl.pp.check_contigs(vdj, adata)
        except IndexError as e:
            logger.warning(f"No contigs passed filter.Please check whether the 5 'rna barcode and the vdj barcode of {raw_airr_file} matches")
            make_file(chain, ddl_wd, samplename)
            return vdj, None
        try: 
            vdj = ddl.tl.find_clones(vdj.data, by_alleles=False, key = "junction_aa", key_added="clone_id", identity={"ig": 0.85, "tr-ab": 1, "tr-gd": 1}, verbose=True)
        except Exception as e:
            logger.warning(f"no clonotype are find, please check your data.The number of input contig {vdj.shape[0]}")
            make_file(chain, ddl_wd, samplename)
            return None

    vdj.simplify()
    vdj = filter_single_chain_by_umi_threshold(vdj)
    vdj = update_clone_id_by_chain_match(vdj, new_airr_file, meta_file, barcode_file)
    vdj = reassign_clone_id_by_size(vdj)
    airr_data = vdj.data
    metadata = vdj.metadata
    
    required_columns = {
        "productive_VDJ": "None",
        "productive_VJ": "None",
        "locus_VDJ": "None",
        "locus_VJ": "None",
        "umi_count_B_VDJ": 0,
        "umi_count_B_VJ": 0,
        "umi_count_abT_VDJ": 0,
        "umi_count_abT_VJ": 0
    }
    
    for col, default_value in required_columns.items():
        if col not in metadata.columns:
            metadata[col] = default_value
            
    if 'clone_id_by_size' in metadata.columns:
        airr_data = airr_data.merge(metadata.loc[:,['clone_id', 'clone_id_by_size', 'exact_subclonotype_id']].drop_duplicates(), how='inner', on='clone_id')
        airr_data.clone_id = 'clonotype' + airr_data.clone_id_by_size.astype(str)
        airr_data = (
                    airr_data.drop(['rearrangement_status', 'clone_id_by_size'], axis=1)
                        .rename({'umi_count': 'duplicate_count'}, axis=1)
                )
    else:
        airr_data = airr_data.merge(metadata.loc[:,['clone_id']].drop_duplicates(), how='inner', on='clone_id')
        airr_data = (
                    airr_data.drop(['rearrangement_status'], axis=1)
                        .rename({'umi_count': 'duplicate_count'}, axis=1)
                )
    airr_data = airr_data.drop('saturation', axis=1)
    ### add gene_star_end
    airr_data = add_star(airr_data, annofile)
    airr_data.to_csv(new_airr_file, sep='\t', index=False)
    if 'clone_id_by_size' in metadata.columns:
        metadata = metadata.assign(
            _clone_id=metadata.clone_id
        ).assign(
            clone_id='clonotype' + metadata.clone_id_by_size.astype(str)
        )
    metadata.to_csv(meta_file, sep='\t', index=True, index_label="cell_id")
    columns = [c for c in metadata.columns.tolist() if not (c.startswith("_") or c.endswith("_main") or re.search(r'_[T|B]_', c))]
    metadata.loc[:,columns].to_csv(barcode_file, sep='\t', index=True, index_label="cell_id")
