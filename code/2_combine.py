## load packages
import scanpy as sc
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


## read data 
# Note that for samples 14, 15, 23, 24 and 25, Scrublet couldn't converge. This is likely to be fine because even when Scrublet converges, we remove very few cells in each sample. As a result, we use the data wihtout doublet detection for these five samples.
# Here we show combination of condensed milk as an example: 
path_sample1 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "11" + ".h5ad")
path_sample2 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "12" + ".h5ad")
path_sample3 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "13" + ".h5ad")

# concatenate samples
adata = sample1.concatenate(sample2,sample3)
adata.obs=adata.obs.drop(["doublet","doublet_score"],axis=1)



## normalization
adata.X = adata.layers["soupX_counts"].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()

# visualize normalization
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sb.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sb.histplot(adata.layers["logcounts"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.savefig(pathout_qc+"/combine_normalization_"+title+".pdf")



## batch correction using scanorama
# define groups
batch1 = adata[adata.obs["batch"] == "0", :]
batch2 = adata[adata.obs["batch"] == "1", :]
batch3 = adata[adata.obs["batch"] == "2", :]

split = [batch1,batch2,batch3]

# define scanorama
def runScanorama(adata, hvg = None):
    import scanorama
    corrected = scanorama.correct_scanpy(split, return_dimred=True)
    corrected = corrected[0].concatenate(corrected[1:])
    return corrected

# run scanorama
adata_scanorama=runScanorama(adata,hvg=None)

# visualization
sc.pp.highly_variable_genes(adata_scanorama)
sc.tl.pca(adata_scanorama)
sc.pp.neighbors(adata_scanorama)
sc.tl.umap(adata_scanorama)

with plt.rc_context():  # Use this to set figure params like size and dpi
    sc.pl.umap(adata_scanorama,color="batch",size=15,alpha=0.5,title="batch_corrected",show=False)
    plt.savefig(pathout_scanorama+"/batchafter_"+title+".pdf", bbox_inches="tight")

# write data
adata.write(path_data+"/combine_"+title+".h5ad")


