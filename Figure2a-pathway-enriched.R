setwd("./")
library(Hmisc)
library(ggplot2)
library(dplyr)


enrich.f<-read.table("Peel.diff.gene.list.path.xls",header = T,sep="\t")
enrich.f <- filter(enrich.f,Pvalue<=0.05) %>% .[,1:5] 
colnames(enrich.f) <- c("Pathway","Diff.genes","All.genes","Pvalue","Qvalue")
enrich.f <- enrich.f[!enrich.f$Pathway%in%c("Metabolic pathways","Biosynthesis of secondary metabolites"),]
enrich.f$Richfactor <- enrich.f$Diff.genes/enrich.f$All.genes
enrich.f$Pathway <- droplevels(enrich.f$Pathway)
enrich.f$Pathway <- reorder(enrich.f$Pathway,enrich.f$Pvalue,FUN=mean)
enrich.f$tissue <- "peel"
enrich.l<-read.table("Leaf.diff.gene.list.path.xls",header = T,sep="\t")
enrich.l <- filter(enrich.l,Pvalue<=0.05) %>% .[,1:5] 
colnames(enrich.l) <- c("Pathway","Diff.genes","All.genes","Pvalue","Qvalue")
enrich.l <- enrich.l[!enrich.l$Pathway%in%c("Metabolic pathways","Biosynthesis of secondary metabolites"),]
enrich.l$Richfactor <- enrich.l$Diff.genes/enrich.l$All.genes
enrich.l$Pathway <- droplevels(enrich.l$Pathway)
enrich.l$Pathway <- reorder(enrich.l$Pathway,enrich.l$Pvalue,FUN=mean)
enrich.l$tissue <- "leaf"
enrich.r<-read.table("Root.diff.gene.list.path.xls",header = T,sep="\t")
enrich.r <- filter(enrich.r,Pvalue<=0.05) %>% .[,1:5] 
colnames(enrich.r) <- c("Pathway","Diff.genes","All.genes","Pvalue","Qvalue")
enrich.r <- enrich.r[!enrich.r$Pathway%in%c("Metabolic pathways","Biosynthesis of secondary metabolites"),]
enrich.r$Richfactor <- enrich.r$Diff.genes/enrich.r$All.genes
enrich.r$Pathway <- droplevels(enrich.r$Pathway)
enrich.r$Pathway <- reorder(enrich.r$Pathway,enrich.r$Pvalue,FUN=mean)
enrich.r$tissue <- "root"
enrich.rfl <- rbind(enrich.f,enrich.l,enrich.r)
enrich.rfl$tissue <- factor(enrich.rfl$tissue,levels=c("peel","leaf","root"))
pdf("pathway.enriched.rfl.pdf",width = 7.5,height = 3.5)
ggplot(enrich.rfl,aes(Richfactor,Pathway))+
  geom_point(aes(size=Diff.genes,color=-1*log10(Pvalue)))+
  scale_color_gradient(low="blue",high="red")+
  facet_grid(~tissue)
dev.off()
