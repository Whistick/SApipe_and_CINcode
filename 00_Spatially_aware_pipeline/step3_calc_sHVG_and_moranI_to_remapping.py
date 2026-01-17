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
samps = samp[1]
def spatial_hvg_calc(samp, p_val_cutoff=0.05, q_val_cutoff=0.05):
    """
    clac Spatial HVGs and save to csv files 
    
    samp
    p_val_cutoff
    q_val_cutoff
    """
    global adata 
    
    filterids = (
        (adata.obs['area'] > 10) &
        (adata.obs['area'] < 1000) &
        (adata.obs['nCounts'] > 50) &
        (adata.obs['nCounts'] < 10000)
    )
    temp_adata = adata[filterids].copy()

    # Spatial Neighbors
    sq.gr.spatial_neighbors(temp_adata, coord_type='generic', delaunay=True)

    # Moran's I
    # n_perms=1000 calc p-value
    print(f"Calculating Spatial Autocorr for {samp}...")
    sq.gr.spatial_autocorr(
        temp_adata, 
        mode='moran', 
        n_perms=1000, 
        n_jobs=32
    )

    # result_df :p-value  p-value adj (FDR)
    result_df = temp_adata.uns['moranI'].copy()
    
    # Moran's I > 0
    shvg_df = result_df[
        (result_df['I'] > 0) & 
        (result_df['I_p_sim_bh'] < q_val_cutoff)
    ].sort_values(by='I', ascending=False)

    shvg_df['Sample'] = samp
 
    # save
    output_name = f'{samp}_Spatial_HVGs_list.csv'
    shvg_df.to_csv(output_name)
    
    print(f"Success! Found {len(shvg_df)} spatial HVGs for {samp}. Saved to {output_name}")
    
    return shvg_df


# shvg_list = spatial_hvg_calc('Sample_01')
