
library(data.table)
library(ggplot2)


pathin="/Users/lianyun.huang/OneDrive/PhD_program/results2023format/singlecell/final_v4/DEGlist"
pathout="/Users/lianyun.huang/OneDrive/PhD_program/results2023format/singlecell/final_v4/DEGvolcano"


# ep_compareto_con, q_compareto_con, soc_compareto_con, 
# q_compareto_ep, soc_compareto_ep, soc_compareto_q
# con_compareto_rest, ep_compareto_rest, q_compareto_rest, soc_compareto_rest

comparison="soc_compareto_rest"
celltype="Glut"
#mycolor=c("#00AFBB", "grey", "#bb0c00") # "#00AFBB", "grey", "#bb0c00"

if (celltype=="GABA") {
  if (comparison %in% c("q_compareto_ep","q_compareto_con")) {
    mycolor=c("#00AFBB", "grey", "#bb0c00") # "#00AFBB", "grey", "#bb0c00"
  } else if (comparison %in% c("q_compareto_rest")) {
    mycolor=c("#00AFBB", "grey")
  } else if (comparison %in% c("ep_compareto_con","soc_compareto_con","soc_compareto_q","ep_compareto_rest")) {
    mycolor=c("grey","#bb0c00")
  } else if (comparison %in% c("soc_compareto_ep","con_compareto_rest","soc_compareto_rest")) {
    mycolor=c("grey")
  }
}


if (celltype=="Glut") {
  if (comparison %in% c("ep_compareto_con","q_compareto_con","q_compareto_ep","soc_compareto_q","ep_compareto_rest","q_compareto_rest")) {
    mycolor=c("#00AFBB", "grey", "#bb0c00") # "#00AFBB", "grey", "#bb0c00"
  } else if (comparison %in% c("soc_compareto_con","soc_compareto_ep")) {
    mycolor=c("#00AFBB", "grey")
  } else if (comparison %in% c("con_compareto_rest")) {
    mycolor=c("grey","#bb0c00")
  } else if (comparison %in% c("soc_compareto_rest")) {
    mycolor=c("grey")
  }
}



df1=fread(paste0(pathin,"/DEGlist_",comparison,"_",celltype,".csv"))
df1=data.frame(df1)

df1[,"diffexpressed"]="NO"

if (dim(df1[which(df1$log2FoldChange < -0.5 & df1$pglobadj<0.05),])[1]==0) {
  print("no value")
} else if (dim(df1[which(df1$log2FoldChange < -0.5 & df1$pglobadj<0.05),])[1]!=0) {
  df1[which(df1$log2FoldChange < -0.5 & df1$pglobadj<0.05),]$diffexpressed="DOWN"
} 

if (dim(df1[which(df1$log2FoldChange > 0.5 & df1$pglobadj<0.05),])[1]==0) {
} else if (dim(df1[which(df1$log2FoldChange > 0.5 & df1$pglobadj<0.05),])[1]!=0) {
  df1[which(df1$log2FoldChange > 0.5 & df1$pglobadj<0.05),]$diffexpressed="UP"
}

df1[,"delabel"]=df1$V1
df1[which(df1$diffexpressed == "NO"),]$delabel=NA

table(df1$diffexpressed)


# "#00AFBB", "grey", "#bb0c00"
pdf(paste0(pathout,"/volcano_",celltype,"_",comparison,".pdf"),width=5,height=5)
ggplot(data = df1, aes(x = log2FoldChange, y = -log10(pglobadj), col = diffexpressed, label = delabel))+
geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed')+
geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed')+
geom_point(size = 2)+
scale_color_manual(values = mycolor)+
geom_text(hjust = 0.8, vjust = 1, check_overlap = TRUE, na.rm=TRUE)+
ggtitle(paste0(celltype,": ",comparison))+
labs(x="log2FC",y="-log10(padj)")+
coord_cartesian(ylim = c(min(df1$pglobadj)-1, max(df1$pglobadj)+1), xlim = c(min(df1$log2FoldChange)-1, max(df1$log2FoldChange)+1))+
theme_classic()+
theme(legend.position="none",
  axis.text.x = element_text(size = 22),
  axis.title.x = element_text(size = 22),
  axis.text.y = element_text(size = 22),
  axis.title.y = element_text(size = 22),
  plot.title=element_text(size=22)
)
dev.off()








