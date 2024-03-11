## load packages
import scanpy as sc
import numpy as np
import seaborn as sns
from scipy.sparse import save_npz
from scipy.io import mmread
from scipy.stats import median_abs_deviation
import scrublet as scr
import matplotlib.pyplot as plt
from matplotlib import colors


## read data
# read filtered data
adata=sc.read_10x_h5(pathin_filtered)
adata.var_names_make_unique()
# read raw data
adata_raw=sc.read_10x_h5(pathin_raw)
adata_raw.var_names_make_unique()


## explore QC parameters
# mark mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("mt-")

# calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=[20], log1p=True, inplace=True)

# visualize QC parameters
with plt.rc_context():
    sc.pl.violin(adata, keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True,show=False)
    plt.savefig(pathout_qc+"/qc_violin_sample"+sample+".pdf", bbox_inches="tight")


with plt.rc_context():
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',color="pct_counts_mt",show=False)
    plt.savefig(pathout_qc+"/qc_mtfrac_scatter_sample"+sample+".pdf", bbox_inches="tight")


plttmp=sns.displot(adata.obs["n_genes_by_counts"], bins=100, kde=False)
plttmp.savefig(pathout_qc+"/qc_ngene_hist_sample"+sample+".pdf")


plttmp=sns.displot(adata.obs["pct_counts_mt"], bins=100, kde=False)
plttmp.savefig(pathout_qc+"/qc_mtfrac_hist_sample"+sample+".pdf")



## perform QC 
# define QC function
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

# QC on MADs metrics
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()

# QC on mitochondria fraction
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 20
)
adata.obs.mt_outlier.value_counts()

# apply QC
print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")



## SoupX in R
# make a rough cluster for SoupX
adata_pp=adata.copy()
sc.pp.normalize_per_cell(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="soupx_groups")

# Preprocess variables for SoupX
soupx_groups = adata_pp.obs["soupx_groups"]
cells = adata.obs_names
genes = adata.var_names
data = adata.X.T
data_tod = adata_raw.X.T

# save variable for SoupX in R
cells.to_series().to_csv(pathout_data+'/cells'+sample+'.csv', header=False)
genes.to_series().to_csv(pathout_data+'/genes'+sample+'.csv', header=False)
soupx_groups.to_csv(pathout_data+'/soupx_groups'+sample+'.csv', header=False)
save_npz(pathout_data+'/data'+sample+'.npz', data)
save_npz(pathout_data+'/data_tod'+sample+'.npz', data_tod)

# read output from SoupX in R
out=mmread(pathout_data+"/out"+sample+".mtx")
out=out.toarray()

# add SoupX corrected count matrix into the layers
adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = out.T
adata.X = adata.layers["soupX_counts"]

# QC on number of genes using SoupX corrected count matrix
print(f"Total number of genes: {adata.n_vars}")
sc.pp.filter_genes(adata, min_cells=5)
print(f"Number of genes after cell filter: {adata.n_vars}")

# write QC-ed data
adata.write(pathout_data+"/finalsoupx_sample"+sample+".h5ad")
# some samples cannot run scrublet. so i save this version before applying scrublet for downstream analyses



## Doublets detection
# run scrublet
adata.obs['doublet_score']= np.zeros(adata.shape[0])
adata.obs['doublet'] = np.zeros(adata.shape[0])
scrub = scr.Scrublet(counts_matrix = adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs['doublet_score'] = doublet_scores
adata.obs['doublet'] = predicted_doublets

# save doublet plots
with plt.rc_context():
    scrub.plot_histogram();
    plt.savefig(pathout_qc+"/qc_doublets_sample"+sample+".pdf", bbox_inches="tight")

# perform doublet filter
doublet_thres=0.5 # set this value for each sample manually
print(f"before filter out doublet: {adata.n_obs}")
print(f"threshold: {doublet_thres}")
print(f"after filter out doublet: {sum(adata.obs['doublet_score']<doublet_thres)}")
adata = adata[adata.obs['doublet_score'] < doublet_thres]

# write data
adata.write(pathout_data+"/soupx_scrublet_sample"+sample+".h5ad")









