setwd("./")
ko.compare <- read.table("Core-vs-non-core.wilcox.reporterscore.list.stat.txt",sep="\t",row.names = 1,header = T)
colnames(ko.compare) <- c("pathway","zscore")

ko.compare.filter <- ko.compare[abs(ko.compare$zscore)>=1.7,]
ko.compare.filter <-ko.compare.filter[order(ko.compare.filter$zscore,decreasing = TRUE),,drop=FALSE]
ko.compare.filter$pathway <- factor(ko.compare.filter$pathway,levels =ko.compare.filter$pathway )
ko.compare.filter$enrich <- ifelse(ko.compare.filter$zscore>0,"Core","Non-core")

mydata <- ko.compare.filter %>%mutate(x_t=ifelse(zscore>=0,-2*1, 2*1))

mytheme <- theme(panel.background = element_blank(),
                 axis.title = element_text(size=rel(1.5)),
                 text=element_text(colour = "black"),
                 axis.text = element_text(colour = "black",size=rel(1.0)))

pdf("meta.ko.enrich.core.non-core.pdf",width = 5.5,height =5)
ggplot(mydata,aes(zscore,pathway,fill=enrich))+geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#B66795","#5F9769"))+mytheme+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  geom_text(aes(x=x_t,label=pathway), color="black",size=3)+
  scale_x_continuous(breaks = c(-6,-5,-4,-3,-1.7,0,1.7,3,4,5,6))
dev.off()
