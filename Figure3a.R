setwd("./")
library(reshape2)

soil <- read.table("soil.txt",header = T)
soil.rd <- soil[soil$Soiltype=="Rhizosphere",]
soil.rd$region <-str_sub(soil.rd$Number,1,1)
soil.rd <- soil.rd[,5:17]
m.soil.rd <- melt(soil.rd)
m.soil.rd$region <- factor(m.soil.rd$region,levels = c("D","X","L","S","M","T","W"))
m.soil.rd$region2 <- ifelse(m.soil.rd$region=="D"|m.soil.rd$region=="X"|m.soil.rd$region=="L"|m.soil.rd$region=="S","Non-core","Core")
pdf("soil.pdf",width =5,height = 5)
ggplot(m.soil.rd,aes(region2,value,fill=region2))+geom_boxplot(outlier.alpha=0)+facet_wrap(variable~.,scales = "free_y",nrow=3)+geom_point(position=position_jitter(width = 0.2),aes(shape=region2))+scale_shape_manual(values = c(1,1))+theme_classic()
dev.off()

core <- grep("^[M|T|W]",soil.rd$region,perl=T)
non-core<- grep("^[D|X|L|S]",soil.rd$region,perl=T)

pfc <- c()
ti <- c()
for(i in colnames(soil.rd)[1:12]){
  pf <- mean(as.numeric(soil.rd[core,i]))/mean(as.numeric(soil.rd[non-core,i]))
  pfc <- c(pfc,pf)
  ti <- c(ti,i)
}

write.table(data.frame(pfc,ti),"wilcox.test.result.soil.rd.v0609.txt",sep="\t",quote=FALSE,row.names = F)
write.table(data.frame(pfc,ti),"fold.change.of.soil.rd.txt",sep="\t",quote=FALSE,row.names = F)


