require(tidyr)
require(dplyr)
require(gplots)
#require(xlsx)

FileName <- "C:/Metabolomics raw data/working_folder/PoolSize"

DataInput <- read.csv(file=paste(FileName,".csv",sep=""),check.names = FALSE)

CleanInput <- DataInput[,c(1,6:ncol(DataInput))]

# for (i in ncol(CleanInput):2){
#   Column_name <- unlist(strsplit(colnames(CleanInput)[i], "_"))
#   colnames(CleanInput)[i] <- paste(Column_name[c(1,3,2)],collapse = "_")
#   
# }
#   
# CleanInput <- CleanInput[,sort(names(CleanInput))]
# CleanInput <- CleanInput[,c(1,80:97,37:55,56:61,68:79,14:31)]

DataMatrix <- data.matrix(CleanInput[,2:ncol(CleanInput)])
row.names(DataMatrix) <- CleanInput[,1]

DataMatrix2 <- log(DataMatrix+1000,2)

RowAvg <- t(apply(DataMatrix2,MARGIN=1,FUN=function(x) x-mean(x)))

par(cex.main=0.8)
heatmap.2(RowAvg,trace="none",Colv=FALSE,margins=c(10,12),dendrogram = "row",col=bluered,cexCol=0.4,
          cexRow = 0.4,srtCol=45,main="Pool Size in Liver after Infusion (Ordered)")  

heatmap.2(RowAvg,trace="none",Colv=TRUE,margins=c(10,12),dendrogram = "both",col=bluered,cexCol=0.4,
          cexRow = 0.4,srtCol=45,main="Pool Size in Liver after Infusion (Clustered)")  
