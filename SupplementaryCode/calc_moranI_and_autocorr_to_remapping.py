import math
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import liana
sc.settings.set_figure_params(dpi=300)

meta_df = pd.read_csv('1.csv')
adata = sc.read_h5ad('adata.h5ad')
samps = meta_df.sort_values(by='Tissue', ascending=False)['Sample'].tolist()
genes = ['MMP12','CXCL10']
gene = 'CXCL10'

def moranI_calc(samp):
    pre_path = '.'
    adata = adata
    filterids = (
        (adata.obs['area'] > 10) &
        (adata.obs['area'] < 1000) &
        (adata.obs['nCounts'] > 50) &
        (adata.obs['nCounts'] < 10000)
    )
    adata = adata[filterids]
    
    # autocorr and moranI
    sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay = True)
    sq.gr.spatial_autocorr(adata, mode='moran', n_jobs=32)
    moran_df = adata.uns['moranI']
    moran_df['Sample'] = samp
    adata.uns['moranI'].to_csv(f'{samp}_MoranI.csv')

