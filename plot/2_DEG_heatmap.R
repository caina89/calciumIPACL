## load packages
library(data.table)
library(tidyr)
library(ggplot2)


## read data
tmp_con=fread(pathin_con)
tmp_con=data.frame(tmp_con)
tmp_con=tmp_con[,c("names","logfoldchanges","pvals_adj")]
colnames(tmp_con)[which(colnames(tmp_con)=="logfoldchanges")]="logfoldchanges_con"
colnames(tmp_con)[which(colnames(tmp_con)=="pvals_adj")]="pvals_adj_con"


tmp_ep=fread(pathin_ep)
tmp_ep=data.frame(tmp_ep)
tmp_ep=tmp_ep[,c("names","logfoldchanges","pvals_adj")]
colnames(tmp_ep)[which(colnames(tmp_ep)=="logfoldchanges")]="logfoldchanges_ep"
colnames(tmp_ep)[which(colnames(tmp_ep)=="pvals_adj")]="pvals_adj_ep"

tmp_q=fread(pathin_q)
tmp_q=data.frame(tmp_q)
tmp_q=tmp_q[,c("names","logfoldchanges","pvals_adj")]
colnames(tmp_q)[which(colnames(tmp_q)=="logfoldchanges")]="logfoldchanges_q"
colnames(tmp_q)[which(colnames(tmp_q)=="pvals_adj")]="pvals_adj_q"

tmp_soc=fread(pathin_soc)
tmp_soc=data.frame(tmp_soc)
tmp_soc=tmp_soc[,c("names","logfoldchanges","pvals_adj")]
colnames(tmp_soc)[which(colnames(tmp_soc)=="logfoldchanges")]="logfoldchanges_soc"
colnames(tmp_soc)[which(colnames(tmp_soc)=="pvals_adj")]="pvals_adj_soc"



## format
data=merge(tmp_con,tmp_ep,by=c("names"),all=TRUE)
data=merge(data,tmp_q,by=c("names"),all=TRUE)
data=merge(data,tmp_soc,by=c("names"),all=TRUE)

data=data[,c("names","logfoldchanges_con","logfoldchanges_ep","logfoldchanges_q","logfoldchanges_soc")]
colnames(data)=c("names","con","ep","q","soc")

data_long=tidyr::gather(data, key = "column", value = "value", -names)


## plot
pdf(paste0(pathout,"/deg_heatmap.pdf"),width=5,height=0.15*nrow(data))
ggplot(data_long, aes(x = column, y = names, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#709AE1", high = "#FD7446",midpoint=0,na.value="black") +  # Adjust the color gradient
  theme_classic()+
  labs(title = "DEG Heatmap", x = "log fold changes", y = " ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()