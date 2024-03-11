## load packages
import scanpy as sc
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


## read data 
# shown as an example: for samples 15, 23, 24, 25, scrublet couldn't run, so i use the data wihtout doublet detection.
# condensed milk
path_sample1 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "11" + ".h5ad")
path_sample2 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "12" + ".h5ad")
path_sample3 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "13" + ".h5ad")
# elevated platform
path_sample4 = sc.read_h5ad(path_data + "/finalsoupx_sample" + "23" + ".h5ad")
path_sample5 = sc.read_h5ad(path_data + "/finalsoupx_sample" + "24" + ".h5ad")
path_sample6 = sc.read_h5ad(path_data + "/finalsoupx_sample" + "25" + ".h5ad")
# quinine
path_sample7 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "19" + ".h5ad")
path_sample8 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "20" + ".h5ad")
path_sample9 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "21" + ".h5ad")
# social
path_sample10 = sc.read_h5ad(path_data + "/finalsoupx_sample" + "15" + ".h5ad")
path_sample11 = sc.read_h5ad(path_data + "/soupx_scrublet_sample" + "16" + ".h5ad")

# concatenate samples
adata = sample1.concatenate(sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10,sample11)
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
batch4 = adata[adata.obs["batch"] == "3", :]
batch5 = adata[adata.obs["batch"] == "4", :]
batch6 = adata[adata.obs["batch"] == "5", :]
batch7 = adata[adata.obs["batch"] == "6", :]
batch8 = adata[adata.obs["batch"] == "7", :]
batch9 = adata[adata.obs["batch"] == "8", :]
batch10 = adata[adata.obs["batch"] == "9", :]
batch11 = adata[adata.obs["batch"] == "10", :]

split = [batch1,batch2,batch3,batch4,batch5,batch6,batch7,batch8,batch9,batch10,batch11]

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






