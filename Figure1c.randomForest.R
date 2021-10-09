setwd("./")
library(MASS)
library(randomForest)
library(Boruta)
library(RColorBrewer)
library(Cairo) #install.packages("Cairo")  用于保存PDF图片
library(pheatmap)
library(stringr)
library(dplyr)
oil <- read.csv("Essential.oils.csv",header = TRUE,sep = ",",row.names = 1)
oil[is.na(oil)] <- 0
oil$region <- str_sub(rownames(oil),1,1)
oil$region2 <- ifelse(oil$region=="D"|oil$region=="X"|oil$region=="L"|oil$region=="S","Non-core","Core")


### leave one out cv ######################
oil.forest <- oil[,colnames(oil)!=c("region")] 
oil.forest$region2 <- as.factor(oil.forest$region2)
loocv_tmp <- data.frame()
fs <- list()
rf.imp <- list()
rf.imp.name <- list()
set.seed(1205)
for (k in 1:56) {
  train <- oil.forest[-k, ]
  test <- oil.forest[k, ]
  #feature selection with only training set
  boruta.train <- Boruta(region2~., data = train, doTrace = 0,getImp=getImpRfZ,maxRuns=999)
  final.boruta <- TentativeRoughFix(boruta.train)
  boruta.df <- attStats(final.boruta)
  a <- boruta.df[which(boruta.df$decision =="Confirmed"),]
  #the selected OTUs ordered by importance
  a <- a[order(a$meanImp,decreasing=TRUE),] 
  select <- rownames(a)
  fs[[k]] <- select
  #the train set with only selected OTUs
  fitted_models <- randomForest(train[,select],train[,'region2'],importance=TRUE, ntree=1000) 
  imp <- randomForest::importance(fitted_models)
  imp <- imp[order(-imp[,3]),]
  rf.imp[[k]] <- imp 
  rf.imp.name[[k]] <- rownames(imp)
  pred <- predict(fitted_models, test[,-16])
  if (test[,14]== as.character(pred)){  #oil.forest[,14]
    loocv_tmp[k,1] <- 1
  }else{
    loocv_tmp[k,1] <- 0
  }
  loocv_tmp[k,2] <- names(pred)
}
(loocv <- mean(loocv_tmp[,1]))
write.table(loocv_tmp,"loocv.ratio.xls",sep="\t",quote=FALSE,row.names = FALSE)
indx <- sapply(fs, length)
fs.df <- as.data.frame(do.call(rbind,lapply(fs, `length<-`,max(indx))))
write.csv(fs.df,"loocv.feature_selection.seed1205.oil.csv")
indx <- sapply(rf.imp.name, length)
rf.imp.name.df <- as.data.frame(do.call(rbind,lapply(rf.imp.name, `length<-`,max(indx))))
write.csv(rf.imp.name.df,"loocv.feature_selection.rf_imp_name.seed12.5.oil.csv")
rf.imp.df <- as.data.frame(do.call(rbind, rf.imp))
write.csv(rf.imp.df,"loocv.feature_selection.rf_imp.seed1205.oil.csv")

oil.all <- unlist(rf.imp.name, recursive = FALSE)
oil.select <- as.matrix(table(oil.all))
oil.sort<-oil.select[rev(order(oil.select[,1])),]
oil.sort
#pick the first n OTUs
(pick <- names(oil.select[oil.select[,1] >= 56*0.9,]))

#### "o.Cymene","α.Pinene","α.Terpinene","α.Thujene","β.Myrcene","β.Pinene","δ.Carene" 
#### end of leave one out ####


################ Feature Selection in R with the Boruta R Package ###

set.seed(1205) 
Boruta.oil <- Boruta(region2 ~ ., data = oil.forest, doTrace =0, maxRuns = 999,getImp = getImpRfZ) 

final.boruta <- TentativeRoughFix(Boruta.oil)
boruta.df <- attStats(final.boruta)
boruta.Confirmed <- boruta.df[which(boruta.df$decision =="Confirmed"),]
boruta.Confirmed <- boruta.Confirmed[order(boruta.Confirmed$meanImp,decreasing=TRUE),] 
select <- rownames(boruta.Confirmed)
###

fitted_models <- randomForest(oil.forest[,select],oil.forest[,'region2'],importance=TRUE, ntree=1000) 
imp <- randomForest::importance(fitted_models)
imp <- imp[order(-imp[,3]),]


##### plot importance ####
final.boruta <- TentativeRoughFix(Boruta.oil)
imp.h <- final.boruta$ImpHistory
imp.m <- reshape2::melt(imp.h)
colnames(imp.m) <- c("number","volatile.oil","importance")

imp.m2 <-imp.m %>% .[.$importance!="-Inf",] %>% group_by(volatile.oil) %>% dplyr::summarize(importance=mean(importance))  %>% .[order(.$importance,decreasing = FALSE),]
imp.m2$volatile.oil
imp.m$volatile.oil <- factor(imp.m$volatile.oil,levels = c("shadowMin","β.Phellandrene","α.terpineol","shadowMean","D.Limonene","α.Sinensal","Methyl.methanthranilate","shadowMax","γ.Terpinene","β.Myrcene","δ.Carene","α.Pinene","α.Thujene","α.Terpinene","o.Cymene","β.Pinene"))

imp.m$volatile.oil <- reorder(imp.m$volatile.oil,imp.m$importance,FUN=mean)
mycolor2=c(rep("#dd1c77",1),rep("#756bb1",2),rep("#dd1c77",1),rep("#756bb1",3),"#dd1c77","#FFED6F",rep("#8DD3C7",7))
library(ggplot2)

ggplot(imp.m,aes(x=volatile.oil,y=importance,fill=volatile.oil))+geom_boxplot()+scale_fill_manual(values=mycolor2)+coord_flip()+guides(fill=FALSE)+xlab("Essential oils")+ylab("Importance")+theme_classic()+geom_hline(yintercept = 2.669011, color = 'gray', size = 0.4,lty=2)
dev.off()


####  2021.1.12  #####
set.seed(1205)
train <- sample(nrow(oil), 0.8*nrow(oil))
oil.train <- oil.forest[train,]
oil.validate <- oil.forest[-train,]
table(oil.train$region2)
table(oil.validate$region2)
fit.forest <- randomForest(region2~α.Pinene + α.Thujene + β.Pinene + α.Terpinene + δ.Carene + β.Myrcene + o.Cymene, data=oil.train, na.action=na.roughfix,importance=TRUE)  
fit.forest$predicted
importance(fit.forest, type=2)                          
forest.pred <- predict(fit.forest, oil.validate)
forest.perf <- table(oil.validate$region2, forest.pred, dnn=c("Actual", "Predicted"))
forest.perf

yes_ratio <-paste(round((forest.perf[1]/(forest.perf[1]+forest.perf[2])),4)*100,"%",sep="")
no_ratio <- paste(round(forest.perf[4]/(forest.perf[3]+forest.perf[4]),4)*100,"%",sep="")
data.frame(yes_ratio,no_ratio)



###########







