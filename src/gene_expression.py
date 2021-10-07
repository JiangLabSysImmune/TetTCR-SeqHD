import numpy as np
import pandas as pd
import scanpy as sc
import sys
import anndata as ad

gene_expression_file = sys.argv[1]
data = pd.read_csv(gene_expression_file, comment = "#", index_col = 0)
adata = ad.AnnData(data)
adata.obs_names = adata.obs_names.astype(str)
sc.pp.filter_genes(adata, min_cells = 50)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.neighbors(adata, n_neighbors = 30, n_pcs = 0 )
sc.tl.umap(adata, min_dist = 0.01)
sc.tl.leiden(adata)
sc.pl.umap(adata, color = ['leiden'], save = "cluster.pdf")
adata.write_csvs("out")

sc.tl.rank_genes_groups(adata, 'leiden', method = 'wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
deg = pd.DataFrame({group + "_" + key:result[key][group] for group in groups for key in ['names','pvals_adj']})
deg.to_csv("out/marker_gene.csv")
sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden', n_genes=5, dendrogram=True,save = "marker.pdf")

