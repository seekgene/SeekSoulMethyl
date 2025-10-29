import os
import json
import pandas as pd
import gzip
from scipy.io import mmread
from .websummary import websummary

basedir = os.path.dirname(os.path.abspath(__file__))

class websummaryRNA(websummary):
    def __init__(self, summary_json, logo=f"{basedir}/logo.png", description_json=f"{basedir}/rna_description.json"):
        super().__init__(logo, description_json)

        with open(summary_json, "r", encoding="utf-8") as f:
            self.summary = json.load(f)        

    def pre_sequence_data(self):
        data = {}
        data["Number of Reads"] = self.format_comma(self.summary["stat"]["total"])
        data["Valid Barcodes"] = self.format_percent(self.summary["stat"]["valid"]/self.summary["stat"]["total"])
        data["Sequencing Saturation"] = self.format_percent(self.summary["cells"]["Sequencing Saturation"])
        data["Too Short"] = self.format_comma(self.summary["stat"]["too_short"])
        data["Q30 Bases in Barcode"] = self.format_percent(self.cal_q30("barcode_q"))
        data["Q30 Bases in UMI"] = self.format_percent(self.cal_q30("umi_q"))
        self.data["Sequencing_data"] = self.format_for_web(data, comp_type="Card", title="Sequencing")

    def pre_mapping_data(self):
        data = { k: self.format_percent(v) for k, v in self.summary["mapping"].items()}
        self.data["Mapping_data"] = self.format_for_web(data, comp_type="Card", title="Mapping")
    
    def pre_cells_data(self):
        data = {k: self.format_comma(v) if isinstance(v, int) else self.format_percent(v) for k, v in self.summary["cells"].items()} 
        if "Sequencing Saturation" in data:
            del(data["Sequencing Saturation"])
        self.data["Cells_data"] = self.format_for_web(data, comp_type="Card", title="Cells")

    def pre_sample_data(self):
        data = {
            "id": self.data["id"],
            "name": "Card",
            "title": "Sample",
            "data": {
                "Name": self.summary["stat"]["samplename"],
                "Description": "",
                "Transcriptome": self.summary["reference"],
                "Chemistry": self.summary["stat"]["chemistry"],
                "Include introns": self.summary["include_introns"],
                "SeekSoul tools Version": self.summary["__version__"],
            }
        }
        data["description"] = { k: self.get_description("Sample", k) for k, v in data["data"].items()}
        self.data["Sample_data"] = data
        self.data["id"] += 1

    def pre_median_data(self):
        self.data["Median_data"] = {
            "id": self.data["id"],
            "name": "PlotlyLine",
            "title": "Median Genes per Cell",
            "data": {
                "x": [0, ] + self.summary["downsample"]["Reads"],
                "y": [0, ] + self.summary["downsample"]["median"],
                "xaxis": {
                    "title": "Mean Reads per Cell",
                },
                "yaxis": {
                    "title":"Median Genes per Cell",
                }
            },
            "description": "Performing downsample processing on the data, with Mean Reads per Cell on the x-axis and Median Genes per Cell on the y-axis"
        }
        self.data["id"] +=1

    def pre_saturation_data(self):
        self.data["Saturation_data"] = {
            "id": self.data["id"],
            "name": "PlotlyLine",
            "title": "Sequencing Saturation",
            "data": {
                "x": [0, ] + self.summary["downsample"]["Reads"],
                "y": [0, ] + self.summary["downsample"]["saturation"],
                "xaxis": {
                    "title": "Mean Reads per Cell",
                },
                "yaxis": {
                    "title":"Sequencing Saturation",
                },
            },
            "description": "Performing downsample processing on the data, with Median Genes per Cell on the x-axis and Sequencing Saturation on the y-axis"
        }
        self.data["id"] +=1

    def pre_diff_data(self, diff_table, n=50):
        df = pd.read_table(diff_table)
        tmp = df.groupby('cluster').apply(lambda x: x.sort_values(by='avg_log2FC', ascending=False).head(n))
        self.data["Diff_data"] = tmp.loc[:,['Ensembl', 'gene', 'avg_log2FC', 'p_val_adj','cluster']].to_dict('records')

    def pre_dim_data(self, dim_table, filtered_dir=None, T_airr=None, B_airr=None):
        df = pd.read_table(dim_table)
        # ordered 
        if filtered_dir:
            with gzip.open(f"{filtered_dir}/barcodes.tsv.gz", "rt") as f:
                barcodes = [l.strip().split()[0] for l in f]
            barcode_set = df.index.intersection(barcodes)
            df = df.loc[barcode_set,:]

        self.data["Dim_data"] = {
            "range":{
                "cluster": [str(i) for i in sorted(df.loc[:, "RNA_snn_res.0.8"].unique().tolist())],     
                "nCount_RNA": [int(df.loc[:,"nCount_RNA"].min()), int(df.loc[:,"nCount_RNA"].max())],
                "nFeature_RNA": [int(df.loc[:,"nFeature_RNA"].min()), int(df.loc[:,"nFeature_RNA"].max())],
                "mito": [round(float(df.loc[:,"percent.mito"].min()),2), round(float(df.loc[:,"percent.mito"].max()),2)],
            },
            "coordinate": {
                "tSNE1": df.loc[:,"tSNE_1"].tolist(),
                "tSNE2": df.loc[:,"tSNE_2"].tolist(),
            },
            "nFeature_RNA": df.loc[:,"nFeature_RNA"].tolist(),
            "nCount_RNA": df.loc[:,"nCount_RNA"].tolist(),
            "mito": [round(float(f), 2) for f in df.loc[:,"percent.mito"].tolist()],
            "cluster": [str(i) for i in df.loc[:, "RNA_snn_res.0.8"].tolist()],
        }
        self.data["id"] +=1

        df["BCR"] = False
        df["TCR"] = False
        if B_airr is not None:
            df_B = pd.read_table(B_airr)
            df["BCR"] = df.index.isin(df_B["cell_id"])
            self.data["Dim_data"]["range"]["BCR"] = sorted(df.loc[:, "BCR"].unique().tolist())
            self.data["Dim_data"]["BCR"] = df.loc[:, "BCR"].tolist()

        if T_airr is not None:
            df_T = pd.read_table(T_airr)
            df["TCR"] = df.index.isin(df_T["cell_id"])
            self.data["Dim_data"]["range"]["TCR"] = sorted(df.loc[:, "TCR"].unique().tolist())
            self.data["Dim_data"]["TCR"] = df.loc[:, "TCR"].tolist()   

    def pre_marker_data(self, marker_json):
        with open(marker_json) as fh:
            self.data["Marker_data"] = json.load(fh)

    def pre_barcode_rank_data(self, raw_dir, filtered_dir):
        with gzip.open(f"{raw_dir}/barcodes.tsv.gz") as fh:
            barcodes = [line.strip() for line in fh]

        with gzip.open(f"{filtered_dir}/barcodes.tsv.gz") as fh:
            cells = {line.strip() for line in fh if line.strip()}

        with gzip.open(f"{raw_dir}/matrix.mtx.gz") as fh:
            m = mmread(fh)

        umi_sum = m.sum(axis=0)

        df = pd.DataFrame({"barcode": barcodes, "UMI": umi_sum.A[0]})
        df["is_cell"] = False
        df.loc[df["barcode"].isin(cells), "is_cell"] = True
        df = (df.sort_values(["UMI", "is_cell"], ascending=[False, False])
                .reset_index(drop=True)
                .reset_index(names="idx"))

        '''
        df = (pd.DataFrame({"barcode": barcodes, "UMI": umi_sum.A[0]})
                .sort_values("UMI", ascending=False).reset_index(drop=True)
                .reset_index(names="idx"))
        df["is_cell"] = False
        df.loc[df["barcode"].isin(cells), "is_cell"] = True
        '''
        self.barcode_rank_data(df)

    def _barcode_rank_data(self, d):
        self.data["Barcode_rank_data"] = {
            "id": self.data["id"],
            "name": "PlotlyLineCell",
            "title": "Sequencing Saturation",
            "description": "",
            "data": {
                "data": d,
                "xaxis": {
                    "title": "Barcodes",
                },
                "yaxis": {
                    "title": "UMI counts",
                }
            }
        }  

    def barcode_rank_data(self, df):
        cell_color = "rgba(80, 80, 201, {:2})"
        bg_color = "rgba(221, 221, 221, 1.0)"
        d = []

        df_bg = df.loc[~df.is_cell,:]

        if df_bg.shape[0] == 0:
            msg = f"{100.0:.0f}% Cells<br>{df.shape[0]}/{df.shape[0]}"
            d.append({
                "x": 0,
                "y": df.UMI.tolist(),
                "text": msg,
                "color": cell_color.format(1.0),
            })
            return self._barcode_rank_data(d)
        
        df_cells = df[df.is_cell]
	
        if df_cells.shape[0] == 0:
            d.append({
                "x": 0,
                "y": df.UMI.tolist(),
                "text": "Background",
                "color": bg_color,
            })
            return self._barcode_rank_data(d)

        # BCCCC BCBBCBCCBBBBBBBBB
        # CCCCC BCBBCBCCBBBBBBBBB
        idx_first = (df_bg[df_bg.idx > df_cells.idx.iloc[0]]).idx.iloc[0]

        df_first = df[:idx_first]
        cell_ratio = df_first[df_first.is_cell].shape[0]/df_first.shape[0]
        msg = (f"{cell_ratio*100.0:.0f}% Cells<br>"
                f"{df_first[df_first.is_cell].shape[0]}/{df_first.shape[0]}")
        d.append({
            "x": 0,
            "y": df[:idx_first+1].UMI.tolist(),
            "text": msg,
            "color": cell_color.format(1.0),
        })

        # BCCCC BCBBCBCCBBBBBBBBB
        # BCCCC BCBBCBCCBBBBBBBBC
        idx_last = df_cells.idx.iloc[-1]

        # BCCCC BCBBCBCC BBBBBBBBB
        #      |        |
        #   idx-first idx-last
        n = 0
        idx_s = idx_first
        idx_e = idx_s + 1000*4**n

        while idx_e < idx_last:
            df_cells_mix = df[idx_s: idx_e]
            cell_ratio = df_cells_mix[df_cells_mix.is_cell].shape[0]/df_cells_mix.shape[0]
            msg = (f"{cell_ratio * 100.0:.0f}% Cells<br>"
                    f"{df_cells_mix[df_cells_mix.is_cell].shape[0]}/{df_cells_mix.shape[0]}")
            d.append({
                "x": int(idx_s),
                "y": df[idx_s: idx_e+1].UMI.tolist(),
                "text": msg,
                "color": cell_color.format(cell_ratio if cell_ratio>0.1 else 0.1),
            })
            idx_s = idx_e
            n += 1
            idx_e = idx_s + 1000*4**n

        if idx_s < idx_last:
            df_cells_mix = df[idx_s: idx_last+1]
            cell_ratio = df_cells_mix[df_cells_mix.is_cell].shape[0]/df_cells_mix.shape[0]
            msg = (f"{cell_ratio*100.0:.0f}% Cells<br>"
                    f"{df_cells_mix[df_cells_mix.is_cell].shape[0]}/{df_cells_mix.shape[0]}")
            d.append({
                "x": int(idx_s),
                "y": df[idx_s: idx_last+1].UMI.tolist(),
                "text": msg,
                "color": cell_color.format(cell_ratio if cell_ratio>0.1 else 0.1), 
            })

        df_last = df[idx_last:]
        if df_last.shape[0]>0:
            msg = "Background"
            d.append({
                "x":  int(idx_last),
                "y": df_last.UMI.tolist(),
                "text": msg,
                "color": "rgba(221, 221, 221, 1.0)",
            })
        return self._barcode_rank_data(d)

    def pre_summary_card1(self):
        data = {
            "id": self.data["id"],
            "name": "CardX",
            "title": "summary1",
            "data": {
                "Estimated Number of Cells": self.format_comma(self.summary["cells"]["Estimated Number of Cells"])
            },
        }
        data["description"] = { k: self.get_description("Cells", k) for k, v in data["data"].items()}
        self.data["Summary_card1"] = data
        self.data["id"] += 1

    def pre_summary_card2(self):
        data = {
            "id": self.data["id"],
            "name": "CardX",
            "title": "summary2",
            "data": {
                "Mean Reads per Cell": self.format_comma(self.summary["cells"]["Mean Reads per Cell"]),
                "Median Genes per Cell": self.format_comma(self.summary["cells"]["Median Genes per Cell"]),               
            }
        }
        data["description"] = { k: self.get_description("Cells", k) for k, v in data["data"].items()}
        self.data["Summary_card2"] = data
        self.data["id"] += 1

    def pre_mtx_data(self, filtered_dir):
        sparse_json = {
            # "mathjs": "SparseMatrix",
            "datatype": "number",
            "values": [],
            "index":[],
            "ptr":[],
            "size": [0, 0],
        }
        with gzip.open("RNA/filtered_feature_bc_matrix/matrix.mtx.gz", "rt") as fh:
            m = mmread(fh)
            csc_m = m.tocsc()
            sparse_json["values"] = csc_m.data.tolist()
            sparse_json["index"] = csc_m.indices.tolist()
            sparse_json["ptr"] = csc_m.indptr.tolist()
            sparse_json["size"] = csc_m.shape

        with gzip.open("RNA/filtered_feature_bc_matrix/barcodes.tsv.gz", "rt") as fh: 
            self.data["barcodes"] = [l.strip().split()[0] for l in fh]

        with gzip.open("RNA/filtered_feature_bc_matrix/features.tsv.gz", "rt") as fh:
            geneID = []
            symbol = []
            for l in fh:
                tmp = l.strip().split()
                geneID.append(tmp[0])
                symbol.append(tmp[1])
            self.data["geneID"] = geneID
            self.data["symbol"] = symbol

        self.data["mtx"] =  sparse_json


    def to_json(
        self, diff_table, dim_table, raw_dir, filtered_dir, 
        marker_json=f"{basedir}/marker.json", T_airr=None, B_airr=None, filename=None,
    ):
        self.pre_summary_card1()
        self.pre_summary_card2()
        self.pre_sequence_data()
        self.pre_mapping_data()
        self.pre_median_data()
        self.pre_saturation_data()
        self.pre_cells_data()
        self.pre_sample_data()
        self.pre_diff_data(diff_table)

        self.pre_dim_data(dim_table, filtered_dir, T_airr, B_airr)
        self.pre_marker_data(marker_json)
        self.pre_barcode_rank_data(raw_dir, filtered_dir)

        data_json = {
            "logo": f'data:image/png;base64, {self.data["logo"]}',
            "Rna":  [
                {
                    "name": "TwoCol",
                    "left": [ self.data["Summary_card1"],],
                    "right": [ self.data["Summary_card2"],]
                },
                {
                    "name": "TwoColX",
                    "left": [ self.data["Cells_data"],],
                    "right": [ self.data["Barcode_rank_data"],]
                },
                {
                    "name": "TwoCol",
                    "left": [ self.data["Sequencing_data"],],
                    "right":[ self.data["Mapping_data"],]
                },
                {
                    "name": "TwoCol",
                    "left": [ self.data["Median_data"],],
                    "right": [ self.data["Saturation_data"],]
                },
                {
                    "name": "TwoCol",
                    "left": [ self.data["Sample_data"],],
                    "right": [],
                }
            ],
            "diff_description": {
                "All": "The top 50 differentially expressed genes in each cluster of the dataset",
                "Marker": "Filter with some commonly used markers for human and mouse",
                "Search": "Global keyword search in the table"
            },
            "diff": self.data["Diff_data"],
            "tsne": self.data["Dim_data"],
            "marker": self.data["Marker_data"],
        }
        if filename:
            with open(filename, "w") as f:
                json.dump(data_json, f)
            return filename
        else:
            return data_json


