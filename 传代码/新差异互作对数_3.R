#用TCGA数据计算差异表达、相关系数

#去掉样本小于10的阶段，查看分布！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
library("rio")

clinlic_files<-list.files(path="D:\\学校电脑\\临床tsv\\",full.names =T,recursive=T)
clin<-import_list(clinlic_files)
#筛选有stage的样本   每个癌型的stage列不同，这里匹配定义一下
col_of_stage<-lapply(clin,function(x){grep("Stage",x)})
col_of_stage[[8]]<-28
col_of_stage[[10]]<-28
col_of_stage[[13]]<-28
col_of_stage[[28]]<-28
col_of_stage[[33]]<-28

#暂时去掉stage不明晰的数据集
jsd<-setdiff(1:33,c(9,14,15,22,23,25))

sample_And_Stage<-lapply(jsd,function(aa){
  x<-clin[[aa]]
  m1<-unique(x[,c(2,col_of_stage[[aa]])])
  m2<-m1[grep("Stage",m1[,2]),]
  m2
})
names(sample_And_Stage)<-names(clin)[jsd]

cancers<-read.table("D:\\cancers.txt")

for(i in cancers[jsd,1]){               
  lujing1<-paste0("D:\\学校电脑\\整合完的表达谱\\",i,'.txt')
  aa<-read.table(lujing1,header = T,row.names = 1,check=F)
  w<-grep(i,names(sample_And_Stage))
  
  
  #临床信息里有stage的样本和表达谱里的样本有一些对不上，这里重新取交集
  bothIn<-intersect(sample_And_Stage[[w]][,1],colnames(aa))
  ins<-match(bothIn,sample_And_Stage[[w]][,1])
  inaa<-match(bothIn,colnames(aa))
  profile_with_stage<-aa[,inaa]
  rm(aa)
  ss<-sample_And_Stage[[w]][ins,2]
  st1<-gsub("[A-Ea-e]","",ss)
  stagess<-gsub("[0-9]","",st1)
  stagess<-gsub("Stg IS","Stg I",stagess)
  
  als<-table(stagess)
  remain_stage<-names(als)[als>=10]
  if(length(remain_stage)>1){
    rema<-stagess %in% remain_stage
    new_stage<-stagess[rema]
    new_profile<-profile_with_stage[,rema]
    rm(profile_with_stage)
    number_of_stage<-length(unique(new_stage))
    calc_fangcha<-apply(new_profile,1,function(x){
      h<-log2(x)
      zhengtai<-tapply(h,new_stage,function(xiaozu){shapiro.test(xiaozu)$p.value})
      if(length(which(zhengtai>0.05))==number_of_stage){
        qici<-bartlett.test(h~new_stage)$p.value
        if(qici>0.05){
          return(1)
        } else{return(0)}
      } else{return(0)}
    })
    #大于70%的基因符合正态+方差齐就方差分析
    if(sum(calc_fangcha)>=dim(new_profile)[1]*0.7){
      pvals<-apply(new_profile,1,function(n){
        h<-log2(n)
        e<-aov(h~new_stage)
        p_of_fangcha<-summary(e)[[1]][1,5]
        p_of_fangcha
      })
    } else{
      pvals<-apply(new_profile,1,function(n){
        kruskal.test(n~new_stage)$p.value
      })
    }
  
    b<-sapply(rownames(new_profile),function(x){unlist(strsplit(x,"[.]"))[1]})
    a<-data.frame(genename=b,pval=pvals)
    lujing2<-paste0("D:\\学校电脑\\所有基因阶段差异检验结果\\",i,".txt")
    write.table(a,file=lujing2,row.names = T,col.names = T,sep="\t",quote=F)
  }
}
    
    
    
#计算相关性
library("qvalue")
okNet<-readRDS("D:\\学校电脑\\提取的lnc-gene互作总网.rds")
all_genes<-unique(c(okNet[,2],okNet[,3]))
cg<-list.files(path="D:\\学校电脑\\所有基因阶段差异检验结果\\",full.names =T,recursive=T)

jsd<-setdiff(1:33,c(9,14,15,22,23,25))


for(i in cancers[jsd,1]){               
  lujing1<-paste0("D:\\学校电脑\\整合完的表达谱\\",i,'.txt')
  aa<-read.table(lujing1,header = T,row.names = 1,check=F)
  w<-grep(i,names(sample_And_Stage))
  
  
  #临床信息里有stage的样本和表达谱里的样本有一些对不上，这里重新取交集
  bothIn<-intersect(sample_And_Stage[[w]][,1],colnames(aa))
  ins<-match(bothIn,sample_And_Stage[[w]][,1])
  inaa<-match(bothIn,colnames(aa))
  profile_with_stage<-aa[,inaa]
  rm(aa)
  ss<-sample_And_Stage[[w]][ins,2]
  st1<-gsub("[A-Ea-e]","",ss)
  stagess<-gsub("[0-9]","",st1)
  stagess<-gsub("Stg IS","Stg I",stagess)
  
  als<-table(stagess)
  remain_stage<-names(als)[als>=10]
  if(length(remain_stage)>1){
    rema<-stagess %in% remain_stage
    new_stage<-stagess[rema]
    new_profile<-profile_with_stage[,rema]
    nn<-sapply(rownames(new_profile),function(x){unlist(strsplit(x,"[.]"))[1]})
    rm(profile_with_stage)

    cyf<-cg[grep(i,cg)]
    a<-read.table(file=cyf,sep="\t",header=T)
    b<-a[a[,1] %in% all_genes,]
jiaozheng<-qvalue(b[,2])
final_pval<-jiaozheng$qvalues
asd<-which(final_pval<0.05)
length(asd)

s1<-okNet[,3] %in% b[asd,1]
s2<-okNet[,2] %in% b[asd,1]
final_net<-okNet[s1 | s2,2:3]


#计算差异对的spearman
lujing2<-paste0("D:\\学校电脑\\新互作边\\",i,".txt")
stageee<-tapply(1:length(new_stage),new_stage,function(stg){stg})
cat("lnc","gene",paste(rep(names(stageee),each=2),c("$estimate","$p.value")),sep="\t",file=lujing2)
cat("\n",file=lujing2,append = T)
mm<-apply(final_net,1,function(bian){
  a<-which(nn==bian[1])
  b<-which(nn==bian[2])
  if(length(a)!=0 & length(b)!=0){
    cat(bian[1],bian[2],sep="\t",file=lujing2,append = T)
    for(kk in 1:length(stageee)){
      xx<-stageee[[kk]]
      d<-cor.test(as.numeric(new_profile[a,xx]),as.numeric(new_profile[b,xx]),method="pearson")
      cat("\t",file=lujing2,append = T)  
      cat(as.numeric(d$estimate),d$p.value,file=lujing2,sep="\t",append = T)
      
    }
    cat("\n",file=lujing2,append = T)
  }
}) 
  }}    
    
    
 
  
#统计各个网络的显著边信息  
library("qvalue")  
lujing3<-paste0("D:\\学校电脑\\筛选显著\\",i,'.txt')
one_file<-"D:\\学校电脑\\癌症各个阶段的统计信息.txt"

for(i in cancers[jsd,1]){               
  lujing1<-paste0("D:\\学校电脑\\新互作边\\",i,'.txt')
  if(file.exists(lujing1)){

a<-read.table(lujing1,sep="\t",header=T)
p_in_lie<-grep("p.value",colnames(a))
jpvs<-apply(a[,p_in_lie],2,function(ph){
    a1<-qvalue(ph)$qvalues
    xz<-which(a1<0.05)
    pair_num<-length(xz)
    lnc_num<-length(unique(a[xz,1]))
    gene_num<-length(unique(a[xz,2]))
    return(c(pair_num,lnc_num,gene_num))
  })
cat("\n",i,"\n",file=one_file,append=T)
rownames(jpvs)<-c("pair_number","lnc_number","gene_number")
colnames(jpvs)<-gsub('..p.value',"",colnames(jpvs))
write.table(jpvs,file=one_file,append=T,quote=F,col.names = T,row.names = T,sep="\t")

}}





#功能富集
library("qvalue")  
library("clusterProfiler")
library("org.Hs.eg.db")

okNet<-readRDS("D:\\学校电脑\\提取的lnc-gene互作总网.rds")

for(i in cancers[jsd,1]){               
  lujing1<-paste0("D:\\学校电脑\\新互作边\\",i,'.txt')
  if(file.exists(lujing1)){
    
    a<-read.table(lujing1,sep="\t",header=T)
    p_in_lie<-grep("p.value",colnames(a))
    huzuogene<-apply(a[,p_in_lie],2,function(ph){
      a1<-qvalue(ph)$qvalues
      xz<-which(a1<0.05)
      ge<-unique(a[xz,2])
      intergenes<-unique(okNet[okNet[,3] %in% ge,4])
      go1<-enrichGO(intergenes,OrgDb='org.Hs.eg.db',ont="BP",pAdjustMethod="fdr",pvalueCutoff=0.05,keyType="ENTREZID")
      kegg1<-enrichKEGG(intergenes,pAdjustMethod="fdr",pvalueCutoff=0.05,keyType="kegg")
      gogs<-go1$Description
      kgs<-kegg1$Description
      n1<-paste0("************",i,"--","GO-BP","************")
      n2<-paste0("************",i,"--","kegg","************")
      return(c(n1,gogs,n2,kgs))
    })
    
    
uio<-do.call(cbind,lapply(huzuogene,'length<-',max(lengths(huzuogene))))    
# uio<-do.call(cbind,lapply(huzuogene,function(x){
#   length(x)<-max(lengths(huzuogene))
#   x
#   }))    
#  这个和上面那行等同，只不过function的写法不同


    colnames(uio)<-gsub('..p.value',"",colnames(uio))
    lujing3<-paste0("D:\\学校电脑\\功能富集\\",i,'.txt')
    write.table(uio,file=lujing3,quote=F,col.names = T,row.names = F,sep="\t")
    
  }}







  
  
  













































