## load packages
library(data.table)
library(tidyr)
library(ggplot2)


## define function that does rename for broad cell type
rename_cell <- function(value) {
  if (value %like% "Glut") {
    return("Glut")
  } else if (value %like% "GABA") {
    return("GABA")
  } else if (value %like% "Dopa") {
    return("Others")
  } else if (value %like% "Sero") {
    return("Others")
  } else if (value %like% "Astro-Epen") {
    return("Others")
  } else if (value %like% "OPC-Oligo") {
    return("Others")
  } else if (value %like% "OEC") {
    return("Others")
  } else if (value %like% "Vascular") {
    return("Others")
  } else if (value %like% "Immune") {
    return("Others")
  } else {
    return("NA")
  }
}


## read all data
container=NULL
for (i in seq(11,26)) {
	data=fread(paste0(pathin,"/mapmycells_",i,".csv"))
	data=data.frame(data)

	keep_cols=c("cell_id", "class_name", "subclass_name", "cluster_name")
	data=data[,keep_cols]

	data = as.data.frame(lapply(data, function(col) sapply(col, rename_cell)))
	print(i)

	tmp=data.frame(table(data$class_name))
	tmp$treatment=i
	container=rbind(container,tmp)
}


## format data
colnames(container)=c("celltype","n_cell","treatment")
combine=spread(container,key="treatment","n_cell")
combine$ctrl=rowMeans(combine[,c("14","17","18","22","26")])
combine$con=rowMeans(combine[,c("11","12","13")])
combine$ep=rowMeans(combine[,c("23","24","25")])
combine$q=rowMeans(combine[,c("19","20","21")])
combine$soc=rowMeans(combine[,c("15","16")])
# save a copy
write.table(combine,paste0(pathout,"/ncell_eachsample.txt"),row.names=F,col.names=T,sep="\t",quote=F)


## further format
useful=combine[,c("celltype","ctrl", "con", "ep", "q", "soc")]
useful=gather(useful,"treatment","n_cell",2:6)



## plot
# Specify the desired order for the x-axis
treatment_order = c("ctrl", "con", "ep", "q", "soc")
container$treatment = factor(container$treatment, levels = treatment_order)
# customize color
custom_colors = c("GABA" = "#8B0000", "Glut" = "#00008B", "Others" = "#FF00FF")
# plot
pdf(paste0(pathout,"/stack_celltype.pdf"),width=5,height=5)
ggplot(useful, aes(fill=celltype, y=n_cell, x=treatment)) + 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  scale_fill_manual(values = custom_colors)+
  labs(title = "cell type in each treatment", x = " ", y = "Percentage")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
