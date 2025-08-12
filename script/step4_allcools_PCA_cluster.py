import os
# 设置环境变量以解决scikit-learn和scipy兼容性问题
os.environ['SCIPY_ARRAY_API'] = '1'

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

from ALLCools.mcds import MCDS
from ALLCools.clustering import tsne, significant_pc_test, log_scale, lsi, binarize_matrix, filter_regions, cluster_enriched_features, ConsensusClustering, Dendrogram, get_pc_centers
from ALLCools.plot import *
import seaborn as sns
import anndata
import re
import click
   
@click.command()
@click.option(
    '--mcds_path', 
    help = 'Path to ALLCools generate-datasets output directory',
    type = click.Path(),
    required = True,
    show_default = True
    )
@click.option(
    '--samplename',
    help = 'Sample name',
    required = True,
    show_default = True
)
@click.option(
    '--var_dim',
    help = 'The name of region to redunction and cluster.',
    default = 'chrom1M',
    show_default = True
)
@click.option(
    '--outdir',
    help = 'Output directory',
    default = './',
    type=click.Path(),
    show_default = True
)
@click.option(
    '--resolution',
    help = 'Resolution for clustering',
    type = float,
    default = 0.3,
    show_default = True
)
@click.option(
    '--n_top_feature',
    help = 'Number of top features to use',
    type = int,
    default = 2500,
    show_default = True
)
@click.option(
    '--pc_cutoff',
    help = 'PC cutoff value',
    type = float,
    default = 0.1,
    show_default = True
)
@click.option(
    '--knn',
    help = 'Number of nearest neighbors',
    type = int,
    default = -1,
    show_default = True
)
@click.option(
    '--filtered_barcode_file',
    help = 'Path to filtered barcode file',
    type = click.Path(),
    default = None,
    show_default = True
)
def allcools_100k_bin_pipeline(
    mcds_path:str,
    samplename:str,
    outdir:str = './',
    var_dim:str = 'chrom1M',
    resolution:float = 0.3,
    n_top_feature:int = 2500,
    pc_cutoff:float = 0.1,
    knn:int = -1,
    filtered_barcode_file:str = None
    ):
    obs_dim = 'cell'
    use_obs = None
    os.makedirs(outdir, exist_ok = True)
    if filtered_barcode_file:
        keep_barcode = pd.read_table(
            filtered_barcode_file,
            header = None,
            sep = '\t')
        keep_barcode.columns = ['cellID']
        use_obs = keep_barcode.cellID
    
    mcds = MCDS.open(
        mcds_path, 
        obs_dim = 'cell', 
        var_dim = var_dim,
        use_obs = use_obs
    )
    total_feature = mcds.get_index(var_dim).size
    mcds.add_feature_cov_mean(var_dim = var_dim, plot = False)
    
    mcds = mcds.filter_feature_by_cov_mean(min_cov = 0)
    
    all_chrom = mcds.chrom1M.chrom1M_chrom.values
    '''
    # 暂时不做region方面的过滤吧
    keep_chroms = [
        'chr1', 'chr2', 'chr3', 'chr4',
        'chr5', 'chr6', 'chr7',
        'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13',
        'chr14', 'chr15', 'chr16',
        'chr17', 'chr18', 'chr19',
        'chr20', 'chr21', 'chr22',
        'chrX'
    ]
    exclude_chromosome = [ i  for i in all_chrom if not i in keep_chroms ]
    '''
    mcds.add_mc_frac(
        normalize_per_cell = True,  # after calculating mC frac, per cell normalize the matrix
        clip_norm_value = 10  # clip outlier values above 10 to 10
    )
    # load only the mC fraction matrix into memory so following steps is faster
    # Only load into memory when your memory size is enough to handle your dataset
    load = True
    if load and (mcds.get_index(obs_dim).size < 20000):
        mcds[f'{var_dim}_da_frac'].load()  
    mcg_hvf = mcds.calculate_hvf_svr(
        var_dim = var_dim,
        mc_type = 'CGN',
        n_top_feature = n_top_feature,
        plot = False) 
    mcg_adata = mcds.get_adata(
        mc_type='CGN',
        var_dim = var_dim,
        select_hvf = True)
    log_scale(mcg_adata)
    sc.tl.pca(mcg_adata)
    cg_n_components = significant_pc_test(mcg_adata)
    cg_pcs = mcg_adata.obsm['X_pca'][:, :cg_n_components]
    cg_pcs = cg_pcs / cg_pcs.std()
    total_pcs = cg_pcs
    adata = mcg_adata.copy()
    adata.obsm['X_pca'] = total_pcs

    del adata.uns['pca']
    del adata.varm['PCs']
    if knn == -1:
        knn = max(15, int(np.log2(adata.shape[0])*2))
    sc.pp.neighbors(adata, n_neighbors = knn)
    sc.tl.leiden(adata, resolution=resolution)
    tsne(adata,
     obsm='X_pca',
     metric='euclidean',
     exaggeration=-1,  # auto determined
     perplexity=30,
     n_jobs=-1)
    
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    _ = categorical_scatter(data = adata,
                            ax = ax,
                            coord_base = 'tsne',
                            hue = 'leiden',
                            text_anno='leiden',
                            show_legend=True)
    fig.savefig(
        f'{outdir}/{samplename}_bismark_allcools_cluster_based_on_{var_dim}_PCA_tsne.png', 
        bbox_inches='tight')
    fig.savefig(
        f'{outdir}/{samplename}_bismark_allcools_cluster_based_on_{var_dim}_PCA_tsne.pdf', 
        bbox_inches='tight')
    
    sc.tl.umap(
        adata, 
        min_dist = 0.3,
        spread = 1)
    fig, ax = plt.subplots(figsize=(4, 4), dpi = 300)
    _ = categorical_scatter(data = adata,
                            ax = ax,
                            coord_base = 'umap',
                            hue = 'leiden',
                            text_anno = 'leiden',
                            show_legend = True)
    fig.savefig(
        f'{outdir}/{samplename}_bismark_allcools_cluster_based_on_{var_dim}_PCA_umap.png', 
        bbox_inches='tight')
    fig.savefig(
        f'{outdir}/{samplename}_bismark_allcools_cluster_based_on_{var_dim}_PCA_umap.pdf', 
        bbox_inches='tight')
    adata.write_h5ad(f'{outdir}/{samplename}_bismark_allcools_based_on_{var_dim}_PCA.h5ad')

if __name__ == '__main__':
    allcools_100k_bin_pipeline()
    
       