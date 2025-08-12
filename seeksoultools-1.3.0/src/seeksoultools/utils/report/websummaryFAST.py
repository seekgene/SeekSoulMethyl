import os
import json
import pandas as pd
from .websummaryRNA import websummaryRNA

basedir = os.path.dirname(os.path.abspath(__file__))

class websummaryFAST(websummaryRNA):
    def pre_biotype_data(self, biotype_table):
        df = pd.read_table(biotype_table)

        def group_type(t):
            if t in ("protein_coding", ):
                return "protein_coding"
            # elif t in ('lincRNA', 'lncRNA'):
            elif t in ("lncRNA", "lincRNA",
                       "non_coding", "3prime_overlapping_ncRNA", "antisense", "sense_overlapping", 
                       "retained_intron", "sense_intronic", "macro_lncRNA", "bidirectional_promoter_lncRNA",):
                return 'lncRNA'
            elif t in ('rRNA',):
                return 'rRNA'
            elif t in ('Mt_rRNA',):
                return 'Mt_rRNA'
            else:
                return 'others'

        df["new_type"] = df["type"].transform(group_type)

        tmp = df.groupby(by=['new_type'])['Count'].sum()

        self.data["Biotype_data"] = {
            "id": self.data["id"],
            "name": "PlotlyPie",
            "title": "Biotype pie",
            "data": {
                "values": tmp.tolist(),
                "labels": tmp.index.tolist(),
            },
            "description": "Proportions of reads mapped to biotypes"
        }
        self.data["id"] += 1

    def pre_genebody_data(self, genebody_file):
        df = pd.read_table(genebody_file, header=None, index_col=0)
        total = df.iloc[1,:].sum()
        self.data["Genebody_data"] = {
            "id": self.data["id"],
            "name": "PlotlyLine",
            "title": "Gene body",
            "data": {
                "x": df.iloc[0,:].tolist(),
                "y": [e/total for e in df.iloc[1,:].tolist()],
                "xaxis": {
                    "title": "Gene: 5' -> 3'",
                },
                "yaxis": {
                    "title": "Coverage",
                    "range": [0, max([e/total for e in df.iloc[1,:].tolist()])+0.02],
                }
            },
            "description": "Gene coverage distribution"
        }
        self.data["id"] += 1

    def to_json(
        self, diff_table, dim_table,  raw_dir, filtered_dir, biotype_table=None, genebody_file=None,
        marker_json=f"{basedir}/marker.json", filename=None,
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
        self.pre_dim_data(dim_table)
        self.pre_marker_data(marker_json)
        self.pre_barcode_rank_data(raw_dir, filtered_dir)
        if biotype_table:    self.pre_biotype_data(biotype_table)
        if genebody_file:    self.pre_genebody_data(genebody_file)


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
                ###{
                ###    "name": "TwoCol",
                ###    "left": [ self.data["Biotype_data"],],
                ###    "right": [ self.data["Genebody_data"]],
                ###},
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
        if biotype_table and genebody_file :
           data_json['Rna'].insert(4, {
                    "name": "TwoCol",
                    "left": [ self.data["Biotype_data"],],
                    "right": [ self.data["Genebody_data"]],})
        if filename:
            with open(filename, "w") as f:
                json.dump(data_json, f)
            return filename
        else:
            return data_json
