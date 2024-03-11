## load packages
import scanpy as sc
import decoupler as dp
import matplotlib.pyplot as plt
import numpy as np
# Import DESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


## read data
adata=sc.read_h5ad(path_data + "/mapmycellannotation_combine_"+ title +".h5ad")
print(title)
print(adata.shape)
print(adata.obs["status"].cat.categories)
print(adata.obs["broadtype"].cat.categories)



## make pseudobulks
# run decoupler
adata_dp = dp.get_pseudobulk(adata,sample_col="batch", 
                             groups_col="broadtype",
                             layer="soupX_counts",
                             mode="sum",min_cells=0,min_counts=0)
# visualize pseudobulks
dp.plot_psbulk_samples(adata_dp, groupby=['status', 'broadtype'], 
                       figsize=(9, 3),
                      save=pathout+"/pseudobulk_"+title+".pdf")



## check pseudobulk profiles
# calcualte visualization parameters
pp_pdata = adata_dp.copy()
sc.pp.normalize_total(pp_pdata, target_sum=1e6)
sc.pp.log1p(pp_pdata)
sc.pp.scale(pp_pdata, max_value=10)
sc.tl.pca(pp_pdata, n_comps=10)
# plot first 2 PCs
with plt.rc_context():
    sc.pl.pca(pp_pdata, color=['status', 'broadtype'], ncols=2, show=False, size=300)
    plt.savefig(pathout+"/top2PCs_"+title+".pdf")



## run DEG analyses
# define case group and control group
strcontrol = adata.obs["status"].cat.categories[0]
strcontrast = adata.obs["status"].cat.categories[1]
title=strcontrast+"_compareto_"+strcontrol

# run DESeq2
for i in adata.obs["broadtype"].cat.categories:
    print(i)
    
    try:
        # select cell profiles
	    tmp = adata_dp[adata_dp.obs['broadtype'] == i].copy()
	    # plot filter
	    dp.plot_filter_by_expr(tmp, group='status', 
	                       min_count=10, min_total_count=15,
	                       save=pathout+"/featureqc_"+title+"_"+i+".pdf")
	    
	    # Build DESeq2 object
	    dds = DeseqDataSet(
	        adata=tmp,
	        design_factors='status',
	        ref_level=["status",strcontrol],
	        refit_cooks=True,
	        n_cpus=8
	    )	

	    # Compute LFCs
	    dds.deseq2()
	    
	    # Extract contrast between two conditions
	    stat_res = DeseqStats(dds, 
	                          contrast=['status',strcontrast,strcontrol], 
	                          n_cpus=8)

	    # Compute Wald test
	    stat_res.summary()	

	    # Shrink LFCs
	    stat_res.lfc_shrink(coeff='status_'+strcontrast+'_vs_'+strcontrol)	

	    # Extract results
	    results_df = stat_res.results_df
        
        # global pvalue correction
	    results_df["pglobadj"]=results_df["padj"]*10
	    results_df=results_df.sort_values(by="pglobadj")
        
        # write csv
	    results_df.to_csv(pathout+"/DEGlist_"+title+"_"+i+".csv",index=True)
    
        # plot volcano
	    dp.plot_volcano_df(results_df, x='log2FoldChange', y='padj',top=0,
                          save=pathout+"/volcanoplot_"+title+"_"+i+".pdf")

    except ValueError:
        continue



