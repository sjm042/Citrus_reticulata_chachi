library(amplicon)
setwd("./")
genus.abd <- read.table('bacterial.genus.abdunce.txt',header = T, row.names = 1,sep="\t")
metadata.group <- data.frame(region=str_sub(colnames(genus.abd),1,1))
metadata.group$region2 <- ifelse(metadata.group$region=="D"|metadata.group$region=="L"|metadata.group$region=="S"|metadata.group$region=="X","Non-core","Core")
rownames(metadata.group) <- colnames(genus.abd)
soil <- read.table("soil.txt",header = T,sep="\t")
rownames(soil) <- paste0(soil$Number,soil$Soiltype)
soil.env <- soil[intersect(colnames(genus.abd),rownames(soil)),]
soil.env <- soil.env[,5:16]
result = RDA_CCA(genus.abd,metadata.group,ps=NULL,env=soil.env,group="region2")
p=result[[1]]
p
ggsave("RDA.genus.abundance.region2.pdf", p, width=180, height=150, units="mm")
