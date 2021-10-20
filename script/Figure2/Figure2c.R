
plot.ko.fpkm <- function(x){
  pfc <- c()
  kk <- c()
  ti <- c()
  n=c()
  require(stringr)
  require(reshape2)
  require(ggplot2)
  ko_fpkm <- read.table("ko_fpkm.prl.txt",sep="\t",header = T,row.names = 1,check.names = F)
  for(i in x){
    fpkm.tmp <- ko_fpkm[grep(i,rownames(ko_fpkm)),] #||K14506|K21604
    fpkm.tmp$ko_id <- str_sub(row.names(fpkm.tmp),1,6)
    tmp.for.plot <- melt(fpkm.tmp)
    tmp.for.plot$tissue <- str_sub(tmp.for.plot$variable,-1,-1)
    tmp.for.plot$region <- str_sub(tmp.for.plot$variable,1,1)
    tmp.for.plot$region2 <- ifelse(tmp.for.plot$region=="M"|tmp.for.plot$region=="W"|tmp.for.plot$region=="T","Xinhui","Others")
    tmp.for.plot$region_t <- paste(tmp.for.plot$region2,tmp.for.plot$tissue,sep = "_")
    tmp.for.plot$region_t <- factor(tmp.for.plot$region_t,levels = c("Others_F","Xinhui_F","Others_L","Xinhui_L","Others_R","Xinhui_R"))
    pdf(paste0("ko-fpkm-",i,".pdf")  ,height = 3,width = 3.5)#_K21604
    p <- ggplot(tmp.for.plot,aes(region2,value,fill=region2))+geom_boxplot(outlier.alpha=0)+ylab("FPKM")+xlab("")+theme_bw()+theme(axis.text.x=element_text(angle = 45,hjust=1,vjust=1))+facet_grid(~ko_id,scales = "free")+geom_point(position=position_jitter(width = 0.2),aes(shape=region2))+scale_shape_manual(values = c(1,1))+theme_classic()+facet_grid(.~tissue,scales="free")
    print (p)
    dev.off()
    xinhui.f <- grep("^[M|T|W].+F$",colnames(fpkm.tmp),perl=T)
    others.f<- grep("^[D|X|L|W].+F$",colnames(fpkm.tmp),perl=T)
    xinhui.l <- grep("^[M|T|W].+L$",colnames(fpkm.tmp),perl=T)
    others.l <- grep("^[D|X|L|W].+L$",colnames(fpkm.tmp),perl=T)
    xinhui.r <- grep("^[M|T|W].+R$",colnames(fpkm.tmp),perl=T)
    others.r <- grep("^[D|X|L|W].+R$",colnames(fpkm.tmp),perl=T)
    pf <- wilcox.test(as.numeric(fpkm.tmp[1,xinhui.f]),as.numeric(fpkm.tmp[1,others.f]),exact = F)[3] %>% as.numeric()
    pl <-  wilcox.test(as.numeric(fpkm.tmp[1,xinhui.l]),as.numeric(fpkm.tmp[1,others.l]),exact = F)[3] %>% as.numeric()
    pr <- wilcox.test(as.numeric(fpkm.tmp[1,xinhui.r]),as.numeric(fpkm.tmp[1,others.r]),exact = F)[3] %>% as.numeric()
    pfc <- c(pfc,pf,pl,pr)
    k <- rep(i,3)
    kk <- c(kk,k)
    t <- c("F","L","R")
    ti <- c(ti,t)
    n=paste(n,i,sep = "_")
  }
  print(data.frame(pfc=pfc,kk=kk,ti=ti))
  write.table(data.frame(pfc=pfc,kk=kk,ti=ti),file=paste0(n,".wilcox.test.result.txt"),sep="\t",quote=F,row.names = F)
}

plot.ko.fpkm(c("K02183","K18834","K20725","K13459"))
