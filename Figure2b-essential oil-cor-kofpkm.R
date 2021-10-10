
library(dplyr)
pathogen.gene.r <- read.table("pathogen-up-core.txt",header = T,sep="\t")
hormone.gene.r <- read.table("hormone-up-core.txt",header = T,sep="\t")
mapk.gene.r <- read.table("mapk-up-core.txt",header = T,sep="\t")




gene2ko <- read.table("gene2ko.txt",sep="\t",header = T)

three.kos <- gene2ko[gene2ko$geneid%in%unique(c(as.character(hormone.gene.r$hormone),as.character(mapk.gene.r$mapk),as.character(pathogen.gene.r$pathogen))),] 
k3 <- three.kos$koid %>% unique() %>% as.character()
k3 <- c(k3,"K07478")  ## add putative ATPase--K07478
samples48 <- read.table("samples48.txt",header = T) %>% .[,1] %>% as.character()
kos.fpkm <- ko_fpkm[rownames(ko_fpkm)%in%k3,]
kos.fpkm.f <- dplyr::select(kos.fpkm,ends_with("F"))
colnames(kos.fpkm.f) <- gsub("F$","",colnames(kos.fpkm.f))
kos.fpkm.f <- kos.fpkm.f[,samples48]
kos.fpkm.l <- dplyr::select(kos.fpkm,ends_with("L"))
colnames(kos.fpkm.l) <- gsub("L$","",colnames(kos.fpkm.l))
kos.fpkm.l <- kos.fpkm.l[,samples48]

kos.fpkm.r <- dplyr::select(kos.fpkm,ends_with("R"))
colnames(kos.fpkm.r) <- gsub("R$","",colnames(kos.fpkm.r))
kos.fpkm.r <- kos.fpkm.r[,samples48]

oil.rel <- read.table("Essential.oils.csv",sep=",",header = TRUE,row.names = 1)
oil.rel[is.na(oil.rel)] <- 0
oil.rel <- oil.rel[samples48,]
oil.rel <- oil.rel[,c("α.Pinene","α.Thujene","β.Pinene","β.Myrcene","α.Terpinene","o.Cymene","δ.Carene")]


#### loop for cor.test for oil and ko fpkm ####
r.cor.result <- data.frame(tmp=1:(ncol(oil.rel)*nrow(kos.fpkm.r)))
k=1
for(i in 1:ncol(oil.rel)){
  for(j in 1:nrow(kos.fpkm.r)){
    p <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.r[j,]),method = "spearman",exact = FALSE)$p.value
    r.cor.result[k,1] <- colnames(oil.rel)[i]
    r.cor.result[k,2] <- rownames(kos.fpkm.r)[j]
    r.cor.result[k,3] <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.r[j,]),method = "spearman",exact = FALSE)$estimate %>% as.numeric()
    r.cor.result[k,4] <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.r[j,]),method = "spearman",exact = FALSE)$p.value
    k=k+1
  }
}
##### end for loop ####
colnames(r.cor.result) <- c("oil","kos","cor","pvalue")
r.cor.result <- na.omit(r.cor.result)
r.cor.result$p.adjust <- p.adjust(r.cor.result$pvalue,method = "BH")
#### filter result #####
r.cor.result.filter <- r.cor.result[r.cor.result$p.adjust<=0.05,]
r.cor.result.filter$abs_cor <- abs(r.cor.result.filter$cor)
r.cor.result.filter$cor2 <- ifelse(r.cor.result.filter$cor>0,1,-1)

################################ leaf #######

l.cor.result <- data.frame(tmp=1:(ncol(oil.rel)*nrow(kos.fpkm.l)))
k=1

for(i in 1:ncol(oil.rel)){
  for(j in 1:nrow(kos.fpkm.l)){
    p <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.l[j,]),method = "spearman",exact = FALSE)$p.value
    l.cor.result[k,1] <- colnames(oil.rel)[i]
    l.cor.result[k,2] <- rownames(kos.fpkm.l)[j]
    l.cor.result[k,3] <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.l[j,]),method = "spearman",exact = FALSE)$estimate %>% as.numeric()
    l.cor.result[k,4] <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.l[j,]),method = "spearman",exact = FALSE)$p.value
    k=k+1
  }
}
colnames(l.cor.result) <- c("oil","kos","cor","pvalue")
l.cor.result <- na.omit(l.cor.result)
l.cor.result$p.adjust <- p.adjust(l.cor.result$pvalue,method = "BH")

l.cor.result.filter <- l.cor.result[l.cor.result$p.adjust<=0.05,]
l.cor.result.filter$abs_cor <- abs(l.cor.result.filter$cor)
l.cor.result.filter$cor2 <- ifelse(l.cor.result.filter$cor>0,1,-1)

################################ peels #######

p.cor.result <- data.frame(tmp=1:(ncol(oil.rel)*nrow(kos.fpkm.f)))
k=1
for(i in 1:ncol(oil.rel)){
  for(j in 1:nrow(kos.fpkm.f)){
    p <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.f[j,]),method = "spearman",exact = FALSE)$p.value
    
    p.cor.result[k,1] <- colnames(oil.rel)[i]
    p.cor.result[k,2] <- rownames(kos.fpkm.l)[j]
    p.cor.result[k,3] <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.f[j,]),method = "spearman",exact = FALSE)$estimate %>% as.numeric()
    p.cor.result[k,4] <- cor.test(as.numeric(oil.rel[,i]),as.numeric(kos.fpkm.f[j,]),method = "spearman",exact = FALSE)$p.value
    k=k+1
  }
}
colnames(p.cor.result) <- c("oil","kos","cor","pvalue")
p.cor.result <- na.omit(p.cor.result)
p.cor.result$p.adjust <- p.adjust(p.cor.result$pvalue,method = "BH")

p.cor.result.filter <- p.cor.result[p.cor.result$p.adjust<=0.05,]
p.cor.result.filter$abs_cor <- abs(p.cor.result.filter$cor)
p.cor.result.filter$cor2 <- ifelse(p.cor.result.filter$cor>0,1,-1)



###  
