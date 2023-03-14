#循环一千次

s<-list.files("D:\\新互作边\\")
ns<-gsub(".txt","",s)
one_file<-"D:\\随机后统计对.txt"

for(i in ns){ 
  
  lujing4<-paste0("D:\\新互作边\\",i,".txt")
  bia<-read.table(file=lujing4,sep="\t",header=T)
  #计算真实网络的评分
  p_in_lie<-grep("p.value",colnames(bia))
  pvals<-bia[,p_in_lie]
  cors<-bia[,p_in_lie-1]
  lie<-dim(pvals)[2]
  nns<-gsub("..estimate","",colnames(cors))
  nns<-gsub("[.]"," ",nns)
  okdata<-sapply(1:dim(pvals)[1],function(z1){
    val<-rep(0,lie)
    xz<-which(pvals[z1,]<0.05)
    if(length(xz)!=0){
      val<-cors[z1,]
      val[setdiff(1:lie,xz)]<-0
    }
    return(val)
  })

  zong_sc<-apply(okdata,2,function(val){
    zong_score<-0
    val<-abs(unlist(val))
    if(length(which(val!=0))!=0){
      aa<-val/max(val)
      specific_score<-(lie-sum(aa))/(lie-1)
      strong_score<-max(val)/0.75
      zong_score<-specific_score*strong_score
    }
    return(zong_score)
  })

  nl<-apply(okdata,2,function(val){
    lx<-"none"
    if(length(which(val!=0))!=0){
      aa<-abs(unlist(val))
      lx<-nns[which.max(aa)]
    }
    return(lx)
  })

  #开始与随机结果比较
  lujing_s<-paste0("D:\\随机矩阵66\\",i,".rds")
  sdfg<-readRDS(lujing_s)
  sjpf<-do.call(cbind,sdfg[1,])
  sjsl<-do.call(cbind,sdfg[2,])
  bj<-sapply(1:length(nl),function(x){
    if(nl[x]!="none"){
      w<-which(sjsl[x,]==nl[x])
    return(length(which(sjpf[x,w]>zong_sc[x]))/1000)
    }
  })
  
  dd<-sapply(1:length(bj),function(s){
    if(!is.null(bj[[s]])){
      if(bj[[s]]<0.05 & zong_sc[s]>0.5){
        return(s)
      }
    }
  })
  
  ud<-unlist(dd)
  nets<-bia[ud,1:2]
  lb<-nl[ud]
  test1<-tapply(1:length(lb), lb, function(gs){
    pair_num<-length(gs)
    lnc_num<-length(unique(nets[gs,1]))
    gene_num<-length(unique(nets[gs,2]))
    return(c(pair_num,lnc_num,gene_num))
  })
  
  tj<-as.data.frame.list(test1)
  rownames(tj)<-c("edge_number","lnc_number","gene_number")
  cat("\n",i,"\n",file=one_file,append=T)
  write.table(tj,file=one_file,append=T,quote=F,col.names = T,row.names = T,sep="\t")
}


