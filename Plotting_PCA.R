require(tidyr)
require(dplyr)
require(gplots)
require(ggplot2)

FileName <- "C:/Metabolomics raw data/working_folder/PoolSize"

DataInput <- read.csv(file=paste(FileName,".csv",sep=""),check.names = FALSE)

CleanInput <- DataInput[,c(1,18:26,6:17,42:44,46:53,27:41)]

DataMatrix <- data.matrix(CleanInput[,2:ncol(CleanInput)])
row.names(DataMatrix) <- CleanInput[,1]


PC<-prcomp(t(DataMatrix),scale=T)
eigs <- PC$sdev^2
PC1v <- round(eigs[1] / sum(eigs),3) *100
PC2v <- round(eigs[2] / sum(eigs),3) *100

PCi<-data.frame(PC$x,Genotype=c(rep("GFP-RD",9),rep("GFP-HFD",12),rep("PKA-RD",11),rep("PKA-HFD",15)))
PCi$Genotype <- factor(PCi$Genotype,levels=c("GFP-RD","GFP-HFD","PKA-RD","PKA-HFD"))

ggplot(PCi,aes(x=PC1,y=PC2,col=Genotype,label=row.names(PCi))) +
  geom_point(size=5,alpha=0.7) + 
  # geom_text(size=3) +
  theme_classic() +
  ggtitle("Principal Component Analysis - GFP vs. PKA & RD vs. HFD") +
  scale_x_continuous(paste("PC1 (",PC1v,"%)"))+
  scale_y_continuous(paste("PC2 (",PC2v,"%)"))+
  theme(axis.line =element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        legend.position = c(0.95, 0.8),legend.key.size = unit(6,units = "mm"))


FileName2 <- "C:/Metabolomics raw data/working_folder/Cleaninput"
Namelist <- read.csv(file=paste(FileName2,".csv",sep=""),check.names = FALSE)
PCr<-data.frame(PC$rotation,Group = Namelist$Group)
PCr$Group <- factor(PCr$Group, levels = c("AA","TCA"))


ggplot(PCr,aes(x=PC1,y=PC2,col=Group,label=row.names(PCr))) +
  geom_point(size=3,alpha=0.7) +
  # geom_text(data=subset(PCr,Group == "AA"),size=5) +
  scale_x_continuous(paste("PC1 (",PC1v,"%)"))+
  scale_y_continuous(paste("PC2 (",PC2v,"%)"))+
  theme_classic() +
  ggtitle("Principal Component Analysis - Loading Plot") +
  theme(axis.line =element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        legend.position = c(0.9, 0.9),legend.key.size = unit(6,units = "mm"))

