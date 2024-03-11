## load packages
import pandas as pd
import scanpy as sc
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt


## read data
# read MapMyCells annotations
data=pd.read_csv(pathin_mapmycells+"/combine_"+title+".csv",skiprows=4)
keep_cols = ["cell_id", "class_name", "subclass_name", "cluster_name"]
data = data[keep_cols]

# read coresponding anndata structure
adata=sc.read_h5ad(pathin_data + "/combine_"+ title +".h5ad")


## combine MapMyCells annotations with anndata
# formatting
adata.obs.reset_index(inplace=True)
adata.obs = adata.obs.rename(columns={'index': 'cell_id'})

# Perform a left join on 'cell_id' (transfer the labels)
adata.obs = adata.obs.merge(data, on='cell_id', how='left')

# write output
adata.obs.to_csv(path_or_buf=pathout+"/mapmycells_"+title+".csv",sep="\t", index=False, header=True,quoting=None)

# make a broader cell type annotation
adata.obs["broadtype"]=adata.obs["class_name"]
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('Glut'),"Glut",adata.obs["broadtype"])
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('GABA'),"GABA",adata.obs["broadtype"])
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('Dopa'),"Others",adata.obs["broadtype"])
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('Sero'),"Others",adata.obs["broadtype"])
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('Astro-Epen'),"Others",adata.obs["broadtype"])
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('OPC-Oligo'),"Others",adata.obs["broadtype"])
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('OEC'),"Others",adata.obs["broadtype"])
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('Vascular'),"Others",adata.obs["broadtype"])
adata.obs["broadtype"]=np.where(adata.obs["class_name"].str.contains('Immune'),"Others",adata.obs["broadtype"])


## visualization
# UMAP for detailed cell types
with plt.rc_context({'figure.figsize': (3, 3)}):
    fig=sc.pl.umap(adata, color= ["class_name"], size=15,show=False) # or the "broadtype" key to show broad cell types
    plt.savefig(pathout+"/mapmycells_umap_"+title+".pdf",bbox_inches='tight')
plt.show()

# percentage barplot
# Count the number of rows for each cell type
cell_type_counts = adata.obs['class_name'].value_counts() # or the "broadtype" key to show broad cell types
# Create a bar plot
plt.rcParams["figure.figsize"] = (3,6)
fig, ax = plt.subplots()
ax.barh(cell_type_counts.index,cell_type_counts.values)
ax.invert_yaxis()
# Set X and Y axis labels and a title
plt.xlabel('n_cells')
plt.title('n_cells per cell type',fontsize=10, fontweight="bold")
plt.savefig(pathout+"/mapmycells_ncells_"+title+".pdf", bbox_inches='tight')
plt.show()


## comparison across all samples
# read all data
container=[]
for i in ["all","allcases", "con", "ctrl_con", "ctrl_ep", "ctrl_q", "ctrl_soc", "ctrl", "ep", "q", "soc"]:
    tmp=pd.read_csv(pathout+"/mapmycells_"+str(i)+".csv",sep="\t")
    tmp=tmp.sort_values(by="class_name")
    tmp=tmp['class_name'].value_counts().index.tolist()
    tmp=sorted(tmp)
    tmp=[i]+tmp
    
    container.append(tmp)

# combine all samples
maxlength=max(len(sublist) for sublist in container)
merge=[sublist + [None] * (maxlength - len(sublist)) for sublist in container]
merge=pd.DataFrame(merge)
merge=merge.transpose()
merge.columns=merge.iloc[0]
merge=merge[1:]
merge.to_csv(pathout+"/celltype_table_combine.csv",sep="\t", index=False, header=True,quoting=None)


