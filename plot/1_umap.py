## load packages
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


## read adata
adata=sc.read_h5ad(path_data + "/mapmycellannotation_combine_"+ title +".h5ad")
adata.X=adata.layers["soupX_counts"]


## UMAP
adata.uns["broadtype_colors"]=['#440154FF','#21908CFF','#FDE725FF']

with plt.rc_context({'figure.figsize':(6,6),'axes.labelsize':25,'axes.titlesize':25}):
    fig=sc.pl.umap(adata, color= ["broadtype"], size=15,alpha=0.8,legend_fontsize=25,show=False)
    plt.savefig(pathout+"/broadtype_"+title+".pdf",bbox_inches='tight')
plt.show()