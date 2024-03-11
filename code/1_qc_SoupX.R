## load packages
library(data.table)
library(reticulate)
library(SoupX)
library(Matrix)

scipy_sparse = import("scipy.sparse")
np = import("numpy")


## read data
# read data_tod
data_tod=scipy_sparse$load_npz(paste0(pathin,"/data_tod",sample,".npz"))

# read soupx_groups
soupx_groups = fread(paste0(pathin,"/soupx_groups",sample,".csv"))
soupx_groups=data.frame(soupx_groups)
soupx_groups=py_run_string("import pandas as pd", convert = TRUE)$pd$Series(soupx_groups$V2)

# read and format other needed data
data=scipy_sparse$load_npz(paste0(pathin,"/data",sample,".npz"))

cells = fread(paste0(pathin,"/cells",sample,".csv"),header=FALSE)
cells=data.frame(cells)
cells=cells[,c(1)]

genes = fread(paste0(pathin,"/genes",sample,".csv"),header=FALSE)
genes=data.frame(genes)
genes=genes[,c(1)]

rownames(data)=genes
colnames(data)=cells

data=as(data,"sparseMatrix")


## run SoupX
# Generate SoupChannel Object
sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

# add extra meta data to the SoupChannel
soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
sc = setSoupProfile(sc, soupProf)

# Set cluster information in SoupChannel
sc = setClusters(sc, soupx_groups)

pdf(paste0(pathout_plot,"/soupxrho_sample",sample,".pdf"))
# Estimate contamination fraction
sc  = autoEstCont(sc, doPlot=TRUE)
dev.off()

# Infer corrected table of counts and rount to integer
out = adjustCounts(sc, roundToInt = TRUE)


## write output
writeMM(out,paste0(pathout,"/out",sample,".mtx"))




