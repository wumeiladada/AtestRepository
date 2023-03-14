#给蛋白质互作文件去重 已经确定了全是2倍重复的 重复的边打分也一致 这里用unique行不行？
w<-read.table("D:\\下载\\9606.protein.links.v11.5.txt",sep=" ",header=T)
nn<-apply(w,1,function(x){c(sort(x[1:2]),x[3])})
library("dplyr")
bb<-as.data.frame(t(nn))
aa<-distinct(bb)
saveRDS(aa,"D:\\去重后的蛋白质互作文件.rds")

#把icp和蛋白质互作文件中的id匹配起来
aa<-read.table("D:\\学校电脑\\9606.protein.aliases.v11.5.txt",sep = "\t",fill=T,quote="")
ll<-read.table("D:\\学校电脑\\匹配好的icp的ensembl-ncbi-symbol.txt",sep="\t",header = T)
a<-unique(aa[aa[,2] %in% ll[,1],1:2])
b<-unique(aa[aa[,2] %in% ll[,2],1:2])
ccc<-merge(a,b,by="V1")
saveRDS(ccc,"D:\\学校电脑\\ENSP-ENSG-ncbi.rds")

#挑出与icp互作的基因  可以叫“关联有罪”策略
ppi<-readRDS("D:\\学校电脑\\去重后的蛋白质互作文件.rds")
EEN<-readRDS("D:\\学校电脑\\ENSP-ENSG-ncbi.rds")
b<-apply(ppi,1,function(x){
  if(x[1] %in% EEN[,1] | x[2] %in% EEN[,1]){
    return(x)
  }
})
all_edges<-matrix(unlist(b),ncol=3,byrow = T)
all_nodes<-unique(c(all_edges[,1],all_edges[,2]))


#去掉lnc-gene网中类型不对的基因
lnc_gene<-read.table("D:\\学校电脑\\Net_LPI(GY).txt",sep="\t",header = T)
zhus<-read.table("D:\\最新注释文件.txt",sep="\t")
nn<-sapply(zhus[,1],function(x){unlist(strsplit(x,"[.]"))[1]})
lnc_in_gc<-zhus[which(zhus[,2]=="lncRNA"),]
s1<-lnc_gene[lnc_gene[,1] %in% nn[as.numeric(rownames(lnc_in_gc))],]
library("clusterProfiler")
library("org.Hs.eg.db")
gene1<-bitr(s1[,2],fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb="org.Hs.eg.db")
s2<-merge(s1,gene1,by.y="SYMBOL",by.x="gene")
s2<-s2[,-3]

#提取基因、边
eg<-aa[grep("Ensembl_gene",aa[,3]),]
P_E<-unique(eg[,1:2])
pp<-match(all_nodes,P_E[,1])
d<-match(P_E[pp,2],s2[,3])
dd<-s2[s2[,3] %in% P_E[pp,2],]
saveRDS(dd,file="D:\\学校电脑\\提取的lnc-gene互作总网.rds")


























