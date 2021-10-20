
setwd("./")
library(ade4)
library(ggplot2)
library(RColorBrewer)
library(ggord)
library(stringr)
#PCA
oil.pca <- read.table("Essential.oils.csv",sep=",",header = T,row.names = 1)
oil.pca[is.na(oil.pca)] <- 0
oil.pca <- oil.pca[,c("α.Pinene","α.Thujene","β.Pinene","α.Terpinene","γ.Terpinene","o.Cymene","δ.Carene")]
pca<- dudi.pca(oil.pca, center = TRUE,scale =TRUE,scannf = FALSE)
pca_eig <- (pca$eig)[1:2] / sum(pca$eig)
sample_site <- data.frame({pca$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCA1', 'PCA2')
sample_site$region <- str_sub(sample_site$names,1,1) 
sample_site$region <- factor(sample_site$region,levels=c("D","X","L","S","M","T","W"))
sample_site$region2<-ifelse(sample_site$region=="D"|sample_site$region=="X"|sample_site$region=="L"|sample_site$region=="S","Non-core","Core") 
color4pca <- c("#c51b8a","#f03b20","#c994c7","#fec44f","#92C953","#356094","#94B3D7")
pdf("Figure1a.pdf",width = 4.5,height = 3.5)
ggplot(sample_site, aes(PCA1, PCA2)) +geom_point(aes(shape=region2,fill=region2),size = 3)+
  scale_shape_manual(values = c(25,21,17))+
  theme_classic()+
  geom_vline(xintercept = 0, color = 'black', size = 0.4,lty=2) + 
  geom_hline(yintercept = 0, color = 'black', size = 0.4,lty=2) +
  scale_color_manual(values = color4pca) + 
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank())+
 stat_ellipse(aes(fill=region2),type="norm",geom="polygon",level=0.8,alpha=0.5)+ 
  labs(x = paste('PC1 (', round(100 * pca_eig[1], 2), '%)'), y = paste('PC2 (', round(100 * pca_eig[2], 2), '%)')) 
dev.off()





