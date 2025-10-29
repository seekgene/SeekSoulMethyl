import os
import json
import pandas as pd
from .websummaryRNA import websummaryRNA

basedir = os.path.dirname(os.path.abspath(__file__))

def has_valid_data(file):
    if not os.path.exists(file) or os.path.getsize(file) == 0:
        return False
    df = pd.read_csv(file, sep='\t')
    return len(df) > 0

class websummaryVDJ(websummaryRNA):
    def __init__(self, summary_csv, samplename="demo", logo=f"{basedir}/logo.png", description_json=f"{basedir}/vdj_description.json"):
        self.samplename = samplename
        super(websummaryRNA, self).__init__(logo, description_json)
        df = pd.read_csv(summary_csv)
        summary_dict = df.to_dict(orient='index')[0]
        self.summary = summary_dict
        self.card1_titles = ["Estimated Number of Cells",]
        self.card2_titles = [ "Mean Read Rairs per Cell", "Number of Cells with Productive V-J Spanning Pair"]
        self.sequencing_titles = [
            "Number of Read Pairs", "Valid Barcodes", "Q30 Bases in Barcode",
            "Q30 Bases in R2 Read", "Q30 Bases in UMI"
        ]
        self.Cells_titles = ["Estimated Number of Cells", "Mean Read Pairs per Cell",
                             "Mean Used Read Pairs per Cell", "Fraction Reads in Cells"]
    def pre_data(self, report_title, data_titles, name="Card"):
        data = {
            "id": self.data["id"],
            "name": name,
            "title": report_title,
            #"data": { t: self.summary[t] for t in getattr(self, data_titles) },
            "data": { t: self.summary[t] for t in getattr(self, data_titles) if t in self.summary},
        }
        data["description"] = { k: self.get_description(report_title, k) for k, v in data["data"].items()}
        self.data["id"] += 1
        return data

    def pre_top10_clonotype_data(self, clonotypes_tsv):
        df = pd.read_table(clonotypes_tsv)
        clonotype_df = (df.head(10).rename(columns={'cdr3s_aa': 'cdr3s'}).loc[:,['clonotype_id', 'cdr3s', 'frequency', 'proportion']]
                        .assign(cdr3s=lambda x: x['cdr3s'].str.replace(';', '<br>')))
        clonotype_df["proportion"] = clonotype_df["proportion"].apply(lambda x: f"{x:.2%}")
        self.data["Top10_clonotype_data"] = clonotype_df.to_dict('records')

    def pre_barcode_rank_data(self, bc_count_tsv, airr_tsv):
        df = pd.read_table(bc_count_tsv)
        cells = pd.read_table(airr_tsv)['cell_id'].unique()
        df['is_cell'] = df.barcode.map(lambda x: x in cells)
        umi_df = (
            df.loc[:,["barcode", "umi_num", "is_cell"]]
                .rename(columns={"umi_num": "UMI"})
                .sort_values(by=['UMI','is_cell'], ascending=False)
                .reset_index().reset_index(names="idx")
                .reindex(columns=["idx", "barcode", "UMI", "is_cell"])
        )
        self.barcode_rank_data(umi_df)
        #print(json.dumps(self.data["Barcode_rank_data"], indent=2))
    def pre_clonotype_frequencies_data(self, clonotypes_tsv):
        df = pd.read_table(clonotypes_tsv)
        bar_df = df.head(n=10)
        bar_df = bar_df.assign(Clonotype_ID = bar_df['clonotype_id'].str.extract(r'(?P<id>\d+)'))
        self.data["Clonotype_frequencies_data"] = {
            "id": self.data["id"],
            "name": "PlotlyBar",
            "title": "Top 10 Clonotype Frequencies",
            "data": {
                # "x": [int(i) for i in bar_df.Clonotype_ID.tolist()],
                "x": bar_df.Clonotype_ID.tolist(),
                "y": bar_df.proportion.tolist(),
                "xaxis": {
                    "title": "Clonotype ID",
                },
                "yaxis": {
                    "title":"Fraction of Cells",
                }
            },
            "description": "The proportion of cells occupied by the top 10 most abundant clone types in the sample"           
        }
        self.data["id"] +=1

class websummaryTCR(websummaryVDJ):
    def __init__(self, summary_csv, samplename="demo", logo=f"{basedir}/logo.png", description_json=f"{basedir}/vdj_description.json"):
        super().__init__(summary_csv, samplename, logo, description_json)
        self.enrich_titles = ["Reads Mapped to Any V(D)J Gene", "Reads Mapped to TRA", "Reads Mapped to TRB"]
        self.expression_titles = ["Median TRA UMIs per Cell", "Median TRB UMIs per Cell"]
        self.annotation_titles = [
            "Cells with Productive V-J Spanning Pair", "Cells with Productive V-J Spanning (TRA, TRB) Pair",
            "Paired Clonotype Diversity", "Cells with TRA Contig", "Cells with TRB Contig",
            "Cells with Productive TRA Contig", "Cells with Productive TRB Contig"
        ]

    def to_json(self, clonotypes_tsv, bc_count_tsv, airr_tsv, filtered_dir=None, outputdir=None, filename=None):
        ###self.pre_barcode_rank_data(bc_count_tsv, airr_tsv)
        ###self.pre_top10_clonotype_data(clonotypes_tsv)
        ###self.pre_clonotype_frequencies_data(clonotypes_tsv)
        has_clonotype_data = has_valid_data(clonotypes_tsv)
        has_airr_data = has_valid_data(airr_tsv)
        if has_airr_data:
            self.pre_barcode_rank_data(bc_count_tsv, airr_tsv)
            barcode_data = self.data["Barcode_rank_data"]
        else:
            self.data["Barcode_rank_data"] = {
            "name": "BarcodeRank",
            "title": "Barcode Rank Plot",
            "description": "No data available",
            "data": []
            }
            barcode_data = self.data['Barcode_rank_data']
        tmp_data = {
         "name": "BarcodeRank",
         "title": "Barcode Rank Plot",
         "description": "No data available",
         "data": []
        }
        if has_clonotype_data:
            self.pre_top10_clonotype_data(clonotypes_tsv)
            self.pre_clonotype_frequencies_data(clonotypes_tsv)
        else:
            self.data["Clonotype_frequencies_data"] = {
            "name": "ClonalDiversity",
            "title": "Top 10 Clonotype Frequencies",
            "description": "No clonotype data available",
            "data": []
            }

            self.data["Top10_clonotype_data"] = {
            "headers": ["Clonotype ID", "CDR3s", "Frequency", "Proportion"],
            "rows": []
             }

        data_json = {
            "logo": f'data:image/png;base64, {self.data["logo"]}',
            "Tcr":  [
                {
                    "name": "TwoCol",
                    "left": [ self.pre_data("summary_card1", "card1_titles", "CardX")],
                    "right": [ self.pre_data("summary_card2", "card2_titles", "CardX")]
                },
                {
                    "name": "TwoColX",
                    "left": [ self.pre_data("Cells", "Cells_titles"),],
                    "right": [barcode_data]
                },
                {
                    "name": "TwoCol",
                    "left":[self.pre_data("Enrichment", "enrich_titles")],
                    "right":[self.pre_data("V(D)J Expression", "expression_titles")],
                },
                {
                    "name": "TwoCol",
                    "left": [self.pre_data("V(D)J Annotation", "annotation_titles")],
                    "right":[
                        self.pre_data("Sequencing", "sequencing_titles"),
                        {
                            "id": self.data["id"],
                            "name": "Card",
                            "title": "Sample",
                            "data": {
                                "Name": self.samplename,
                                # "description": ""
                            },
                        }
                    ],
                },
                {
                    "name": "OneCol",
                    "data": self.data["Clonotype_frequencies_data"],
                },
                {
                    "name": "OneCol",
                    "data": {
                        "name": "Top10Clonotype",
                        "title": "Top 10 Clonotype CDR3 Sequences",
                        "description": "The detailed information of the top 10 most abundant clone types in the sample",
                        "data": self.data["Top10_clonotype_data"],
                    }
                }
            ],
        }
        if filename:
            with open(filename, "w") as f:
                json.dump(data_json, f)
            return filename
        else:
            return data_json
            # return json.dumps(data_json).replace("5'", "5\\'").replace("3'", "3\\'")

class websummaryBCR(websummaryVDJ):
    def __init__(self, summary_csv, samplename="demo", logo=f"{basedir}/logo.png", description_json=f"{basedir}/vdj_description.json"):
        super().__init__(summary_csv, samplename, logo, description_json)
        self.enrich_titles =  ["Reads Mapped to Any V(D)J Gene", "Reads Mapped to IGH", "Reads Mapped to IGK","Reads Mapped to IGL"]
        self.expression_titles =  ["Median IGH UMIs per Cell", "Median IGK UMIs per Cell", "Median IGL UMIs per Cell"]
        self.annotation_titles =  [
            "Cells with Productive V-J Spanning Pair", "Cells with Productive V-J Spanning (IGK, IGH) Pair",
            "Cells with Productive V-J Spanning (IGL, IGH) Pair", "Paired Clonotype Diversity", 
            "Cells with Productive IGH Contig", "Cells with Productive IGK Contig",
            "Cells with Productive IGL Contig"
        ]
    def to_json(self, clonotypes_tsv, bc_count_tsv, airr_tsv, filtered_dir=None, outputdir=None, filename=None):
        ###self.pre_barcode_rank_data(bc_count_tsv, airr_tsv)
        ###self.pre_top10_clonotype_data(clonotypes_tsv)
        ###self.pre_clonotype_frequencies_data(clonotypes_tsv)
        has_clonotype_data = has_valid_data(clonotypes_tsv)
        has_airr_data = has_valid_data(airr_tsv)
        if has_airr_data:
            self.pre_barcode_rank_data(bc_count_tsv, airr_tsv)
        else:
            self.data["Barcode_rank_data"] = {
            "name": "BarcodeRank",
            "title": "Barcode Rank Plot",
            "description": "No data available",
            "data": []
            }
        if has_clonotype_data:
            self.pre_top10_clonotype_data(clonotypes_tsv)
            self.pre_clonotype_frequencies_data(clonotypes_tsv)
        else:
            self.data["Clonotype_frequencies_data"] = {
            "name": "ClonalDiversity",
            "title": "Top 10 Clonotype Frequencies",
            "description": "No clonotype data available",
            "data": []
            }
        
            self.data["Top10_clonotype_data"] = {
            "headers": ["Clonotype ID", "CDR3s", "Frequency", "Proportion"],
            "rows": []
             }
        data_json = {
            "logo": f'data:image/png;base64, {self.data["logo"]}',
            "Bcr":  [
                {
                    "name": "TwoCol",
                    "left": [ self.pre_data("summary_card1", "card1_titles", "CardX")],
                    "right": [ self.pre_data("summary_card2", "card2_titles", "CardX")]
                },
                {
                    "name": "TwoColX",
                    "left": [ self.pre_data("Cells", "Cells_titles"),],
                    "right": [ self.data["Barcode_rank_data"],]
                },
                {
                    "name": "TwoCol",
                    "left":[self.pre_data("Enrichment", "enrich_titles")],
                    "right":[self.pre_data("V(D)J Expression", "expression_titles")],
                },
                {
                    "name": "TwoCol",
                    "left": [self.pre_data("V(D)J Annotation", "annotation_titles")],
                    "right":[
                        self.pre_data("Sequencing", "sequencing_titles"),
                        {
                            "id": self.data["id"],
                            "name": "Card",
                            "title": "Sample",
                            "data": {
                                "Name": self.samplename,
                                # "description": ""
                            },
                        }
                    ],
                },
                {
                    "name": "OneCol",
                    "data": self.data["Clonotype_frequencies_data"],
                },
                {
                    "name": "OneCol",
                    "data": {
                        "name": "Top10Clonotype",
                        "title": "Top 10 Clonotype CDR3 Sequences",
                        "description": "The detailed information of the top 10 most abundant clone types in the sample",
                        "data": self.data["Top10_clonotype_data"],
                    }
                }
            ],
        }
        if filename:
            with open(filename, "w") as f:
                json.dump(data_json, f)
            return filename
        else:
            return data_json
            
