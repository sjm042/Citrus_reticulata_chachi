
#### Figure 2b and wilcox test
setwd("./")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(stringr)
library(reshape2)
oil <- read.table("Essential.oils.csv",sep=",",header = T)
oil[is.na(oil)] <- 0

oil$region <- str_sub(oil$sample,1,1)
oil$site <- ifelse(oil$region=="D"|oil$region=="X"|oil$region=="L"|oil$region=="S","Non-core","Core")
a <- oil[,!colnames(oil)%in%c("sample","region")]
a <- melt(a)
data.plot <- a[!a$variable%in%c("D.Limonene","γ.Terpinene"),] 
data.plot <- data.plot[data.plot$value<5,]
data.plot$variable <- reorder(data.plot$variable,data.plot$value,FUN=mean)
data.plot$variable <- factor(data.plot$variable,levels = rev(levels(data.plot$variable)))
pdf("essential.oils.11.pdf",width =7,height = 3.5)
ggplot(data.plot,aes(site,value,fill=site))+geom_boxplot(outlier.alpha=0)+facet_grid(.~variable,scales = "free_y")+geom_point(position=position_jitter(width = 0.2),aes(shape=site))+scale_shape_manual(values = c(1,1))+theme_classic()
dev.off()
pdf("essential.oils.for.D.Limonene.pdf",width =2,height = 3.5)
ggplot(a[a$variable=="D.Limonene",],aes(site,value,fill=site))+geom_boxplot(outlier.alpha=0)+facet_grid(.~variable,scales = "free_y")+geom_point(position=position_jitter(width = 0.2),aes(shape=site))+scale_shape_manual(values = c(1,1))+theme_classic()
dev.off()

pdf("essential.oils.for.γ.Terpinene.pdf",width =2,height = 3.5)
ggplot(a[a$variable=="γ.Terpinene",],aes(site,value,fill=site))+geom_boxplot(outlier.alpha=0)+facet_grid(.~variable,scales = "free_y")+geom_point(position=position_jitter(width = 0.2),aes(shape=site))+scale_shape_manual(values = c(1,1))+theme_classic()
dev.off()

#### wilcox test ####
oilname <- c()
pvalue <- c()
for(i in colnames(oil)[2:14]){
  p <-wilcox.test(oil.pca[,i]~site,data=oil,exact = F)[[3]]
  pvalue <- c(pvalue,p)
  oilname <- c(oilname,i)
}
oil.wilcox <- data.frame(oilname=oilname,pvalue=pvalue)
write.table(oil.wilcox,"essential.oil.wilcox.txt",quote = F,sep="\t",row.names = F)
