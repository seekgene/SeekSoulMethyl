import scipy.io as io
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import plotly.express as px
import operator
import base64

class plotUtils:
    @staticmethod
    def line_color(p):
        _color = ['#dbdbf4',
                '#c5c5f5',
                '#afaff6',
                '#9999f7',
                '#8383f8',
                '#6d6df9',
                '#5858fb',
                '#4141fc',
                '#2b2bfd',
                '#1616fe']
        if p==1:
            return '#0000ff'
        elif p==0:
            return '#dddddd'
        else:
            return _color[int(p/0.1)]
        
    @staticmethod
    def barcode_rank_plot(umi_df):
        # continue
        non_cell_df = umi_df.loc[~umi_df['is_cell'],]
        if non_cell_df.empty:
            plot_data = [{
                "x": list(range(1, len(umi_df) + 1)),
                "y": umi_df['umis'].tolist(),
                "text": "100% Cells",
                "hoverinfo": 'text',
                "line": dict(color='#0000ff'),
                "mode": 'lines',
                "name": 'Cells',
                "showlegend": True,
            }]
            
            config = {
                "displayModeBar": True,
                "displaylogo": False,
                'modeBarButtonsToRemove': ['lasso','zoom','pan', 'zoomin', 'zoomout', 'autoscale', 'select2d'],
                'toImageButtonOptions': {
                    'format': 'svg',
                    'filename': 'custom_image',
                    'height': 500,
                    'width': 700,
                    'scale': 1
                },
            }

            fig = go.Figure()
            for _ in plot_data:
                fig.add_trace(go.Scatter(**_))
            fig.update_layout(xaxis_type = "log", yaxis_type = "log")
            fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='#dee2e6')
            fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='#dee2e6')
            fig.layout.height = 500
            fig.layout.width = 500
            fig.update_layout(
                title=None,
                xaxis=dict(title="Barcodes"),
                yaxis=dict(title="UMI counts"),
                xaxis_type = "log",
                yaxis_type = "log",
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)',
                modebar_bgcolor='rgba(0,0,0,0)',
                modebar_color='#dee2e6'
            )
            json_dict = fig.to_plotly_json()
            json_dict['config'] = config
            return json_dict
        idx0 = umi_df.loc[~umi_df['is_cell'],].index[0]
        idx1 = umi_df[::-1].loc[umi_df['is_cell'],].index[0]
        umi_df_head = umi_df.iloc[:idx0,]
        umi_df_tail = umi_df.iloc[idx1+1:,]
        umi_df_mix = umi_df.iloc[idx0:idx1+1,]

        # head, 1 trace
        res_list = []
        tmp = []
        n = 0
        if umi_df_head.shape[0]>0:
            for row in umi_df_head.itertuples(index=False):
            # for row in demo_df.itertuples(index=False):
                if (row.count==1):
                    tmp.append((n, row.umis))
                    n += row.count
                else:
                    tmp.append((n, row.umis))
                    tmp.append((n+row.count-1, row.umis))
                    n += row.count
            res_list.append([tmp, f"100% Cells<br>{n}/{n}", 1])

        # mix
        tmp = []
        # row_num = umi_df_mix.shape[0]
        # break_num = int(np.log10(row_num)/20)
        break_num = n + 10
        step = 100
        counter = 0
        is_cell_num = 0
        barcode_num = 0
        for row in umi_df_mix.itertuples(index=False):
            if (row.count==1):
                tmp.append((n, row.umis))
                if row.is_cell:
                    is_cell_num += row.count
                n += row.count
                barcode_num += row.count
            else:
                tmp.append((n, row.umis))
                tmp.append((n+row.count-1, row.umis))
                if row.is_cell:
                    is_cell_num += row.count
                n += row.count
                barcode_num += row.count
            # counter += 1
            if  n >= break_num:
                # counter = 0
                res_list.append([tmp, f"{is_cell_num/barcode_num:.0%} Cells<br>{is_cell_num}/{barcode_num}", is_cell_num/barcode_num])
                tmp = []
                break_num = n + step
                step = step * 10
                is_cell_num = 0
                barcode_num = 0
                # print(break_num)
        else:
            if barcode_num !=0 :
                res_list.append([tmp, f"{is_cell_num/barcode_num:.0%} Cells<br>{is_cell_num}/{barcode_num}",  is_cell_num/barcode_num])

        # tail, 1 trace
        tmp = []
        for row in umi_df_tail.itertuples(index=False):
        # for row in demo_df.itertuples(index=False):
            if (row.count==1):
                tmp.append((n, row.umis, 1))
                n += row.count
            else:
                tmp.append((n, row.umis, 1))
                tmp.append((n+row.count-1, row.umis, 1))
                n += row.count
        res_list.append([tmp, "Background", 0])
        plot_data = []
        #print(res_list)
        for idx in range(len(res_list)):
            if idx > 0:
                plot_data.append({
                    "x": [plot_data[-1]['x'][-1], ]  + [_[0] for _ in res_list[idx][0]],
                    "y": [plot_data[-1]['y'][-1], ]  + [_[1] for _ in res_list[idx][0]],
                    "text": res_list[idx][1],
                    "hoverinfo": 'text',
                    "line": dict(color=plotUtils.line_color(res_list[idx][-1])),
                    "mode": 'lines',
                    "legendgroup": 'Cells',
                    "showlegend": False,
                })
            else:
                plot_data.append({
                    "x": [_[0] for _ in res_list[idx][0]],
                    "y": [_[1] for _ in res_list[idx][0]],
                    "text": res_list[idx][1],
                    "hoverinfo": 'text',
                    "line": dict(color=plotUtils.line_color(res_list[idx][-1])),
                    "mode": 'lines',
                    "legendgroup": 'Cells',
                    "showlegend": False,
                })

        plot_data[0]['name'] = 'Cells'
        plot_data[0]['showlegend'] = True
        plot_data[-1]['name'] = 'Background'
        plot_data[-1]['showlegend'] = True
        plot_data[-1]['legendgroup'] = 'Background'

        config = {
            "displayModeBar": True,
            "displaylogo": False,
            'modeBarButtonsToRemove': ['lasso','zoom','pan', 'zoomin', 'zoomout', 'autoscale', 'select2d'],
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'height': 500,
                'width': 700,
                'scale': 1
            },
        }

        fig = go.Figure()
        for _ in plot_data:
            fig.add_trace(go.Scatter(**_))
        fig.update_layout(xaxis_type = "log", yaxis_type = "log")
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='#dee2e6')
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='#dee2e6')
        fig.layout.height = 500
        fig.layout.width = 500
        fig.update_layout(
            #title='Barcode Rank',
            #title_x=0.5,
            title=None,
            xaxis=dict(title="Barcodes"),
            yaxis=dict(title="UMI counts"),
            xaxis_type = "log",
            yaxis_type = "log",
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            modebar_bgcolor='rgba(0,0,0,0)',
            modebar_color='#dee2e6'
        )
        # fig.show(config=config)
        json_dict = fig.to_plotly_json()
        json_dict['config'] = config
        return json_dict


    @staticmethod
    def barcode_rank_plot_rna(cells_gz, barcodes_gz, mtx_gz):
        cells_df = pd.read_table(cells_gz, header=None, names=['barcode'])
        all_barcode_df = pd.read_table(barcodes_gz, header=None, names=['barcode'])
        cells_idx_df = all_barcode_df.loc[all_barcode_df.barcode.isin(cells_df.barcode),:]

        mat = io.mmread(mtx_gz)
        df = pd.DataFrame({"feature": mat.row, "barcode":mat.col, "umi": mat.data})
        umi_df = df.groupby(by="barcode").agg(umis=("umi", sum)).sort_values(by='umis', ascending=False)
        umi_df = umi_df.assign(is_cell = umi_df.index.isin(cells_idx_df.index)).reset_index()
        umi_df = umi_df.groupby(by=['umis', 'is_cell']).size().to_frame(name="count").reset_index().sort_values(by=['umis', 'is_cell'], ascending=False).reset_index(drop=True)
        return plotUtils.barcode_rank_plot(umi_df)
    
    @staticmethod
    def pearson_moment_coefficient(lst):
        '''measure skewness'''
        mid_value = lst[int(len(lst)/2)]
        sigma = np.std(lst, ddof=1)
        tmp = []
        for i in lst:
                tmp.append(((i - mid_value)/sigma)**3)
        return np.mean(tmp)

    @staticmethod
    def plot_gene_body(genebody_file):
        config = { 
                "displayModeBar": True,
                "displaylogo": False,
                'modeBarButtonsToRemove': ['lasso','zoom','pan', 'zoomin', 'zoomout', 'autoscale', 'select2d'],
                'toImageButtonOptions': {
                    'format': 'svg',
                    'filename': 'custom_image',
                    'height': 500,
                    'width': 500,
                    'scale': 1
                },
            }
        dataset=[]

        for line in open(genebody_file,'r'):
                    line = line.strip()
                    if line.startswith("Percentile"):
                            continue
                    f = line.split()
                    name = f[0]
                    dat = [float(i) for i in  f[1:]]
                    skewness = plotUtils.pearson_moment_coefficient(dat)
                    #dataset.append((name, [(i -min(dat))/(max(dat) - min(dat)) for i in dat], skewness))
                    dataset.append((name,[i / sum(dat) for i in dat],skewness))
        dataset.sort(key = operator.itemgetter(2), reverse=True)
        df = pd.DataFrame({"Gene: 5' -> 3'":[ i for i in range(1,100+1)],'Coverage':dataset[0][1]})
        fig = px.line(df, x="Gene: 5' -> 3'", y="Coverage",color_discrete_sequence=px.colors.qualitative.G10)
        fig = fig.update_layout(margin = {'l': 50, 'r': 30, 't': 30,'b': 50 },width =  520, height = 330,template='plotly_white')
        fig.update_yaxes(range=[0, 0.03])
        fig_json = fig.to_plotly_json()
        fig_json['config'] = config
        for i in fig_json['data'][0].keys():
            if isinstance(fig_json['data'][0][i],np.ndarray):
                fig_json['data'][0][i] = fig_json['data'][0][i].tolist()
        return fig_json
    

    @staticmethod
    def plot_pie(biotype_file):
        config = { 
                "displayModeBar": True,
                "displaylogo": False,
                'modeBarButtonsToRemove': ['lasso','zoom','pan', 'zoomin', 'zoomout', 'autoscale', 'select2d'],
                'toImageButtonOptions': {
                    'format': 'svg',
                    'filename': 'custom_image',
                    'height': 500,
                    'width': 500,
                    'scale': 1
                },
            }
        df = pd.read_table(biotype_file)
        df0 = df.groupby('type').agg({'Count':'sum'})
        df0['type']=df0.index.values.tolist()
        needtype=['lncRNA','protein_coding','Mt_rRNA','rRNA']
        df0n=df0.loc[df0['type'].isin(needtype), :].reset_index(drop=True)
        df0o=df0.loc[df0['type'].isin(needtype)==False, :]
        df0n.loc['others']=[df0o['Count'].sum(),'others']
        fig = px.pie(df0n, values="Count", names="type",color_discrete_sequence = px.colors.qualitative.G10)
        fig = fig.update_traces(textposition='auto', textinfo='percent+label',insidetextorientation='radial')
        fig = fig.update_layout(margin = {'l': 50, 'r': 30, 't': 30,'b': 50 },width =  520, height = 330)

        fig_json = fig.to_plotly_json()
        for i in fig_json['data'][0].keys():
            if isinstance(fig_json['data'][0][i],np.ndarray):
                fig_json['data'][0][i] = fig_json['data'][0][i].tolist()
        fig_json['config'] = config
        return fig_json

    @staticmethod
    def png2base64(f):
        with open(f, 'rb') as fh:
            return base64.b64encode(fh.read()).decode()
