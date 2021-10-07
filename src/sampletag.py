import numpy as np
import pandas as pd
import scanpy as sc
import sys
import anndata as ad

sampletag_file = sys.argv[1]
data = pd.read_csv(sampletag_file, index_col = 0)
adata = ad.AnnData(data)
adata.obs_names = adata.obs_names.astype(str)
sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.pp.neighbors(adata, n_neighbors = 30, n_pcs = 0 )
sc.tl.tsne(adata)
sc.tl.louvain(adata)
sc.pl.tsne(adata, color = ['louvain'], save = "sampletag.pdf")
#adata.write_csvs("out")
