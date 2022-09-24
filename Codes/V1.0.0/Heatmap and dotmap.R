library(ggplot2)
library(tidyverse)
library(patchwork)
library(enrichplot)
library(Scillus)
library(Seurat)
        Treg<-subset(SP9w,idents = c("Treg"))
        Treg<-subset(PP.combined,idents = c("Treg"))
        Treg<-subset(SP10d,idents = c("Treg"))
        Treg<-subset(Tregcccc,idents = c("Treg"))
        Treg<-subset(Tregcccc2,idents = c("Treg"))
        Idents(Treg)<-"orig.ident"
Data<-Treg@assays$RNA@scale.data %>% data.frame
gene=c("Foxp3","Gata3","Ctla4","Icos","Il2ra","Nrp1","Tnfrsf4","Tnfrsf9","Tnfrsf18","Tgfb1")
gene=c("Foxp3","Gata3","Ctla4","Icos","Il2ra","Nrp1","Tnfrsf4","Tnfrsf9","Tnfrsf18","Ccr2","Cxcr3","Tgfb1")
gene=c("Foxp3","Gata3","Ctla4","Icos","Il2ra","Nrp1","Tnfrsf4","Tnfrsf9","Tnfrsf18","Tgfb1","Ccr7","Sell","S100a6","S100a4")
gene=c("Tbx21","Cxcr3","Ifng")
gene=c("Rorc","Il17a")
gene=c("Lag3","Nt5e","Pdcd1","Cxcr5")
gene=c("Ccr2","Cxcr3")
gene=c("Klf2","Ccr7","S1pr1","S1pr4","Sell")
gene=c("Cd40lg","Il4","Il21","Cxcl13")
gene=c("Dapl1")
gene=c("NFATC1","FOXO1","BLIMP1","TBX21","TOX","TOX2")
gene=c("Foxo1","Ets1","Klf2","Il7r","Ccr7","Sell","S1pr1","S1pr4")
gene=c("Rorc","Il17a","Il17f","Il22")
gene=c("Pdcd1","Cxcr5","Icos","Tox","Tox2","Cd40lg","Il4","Il21")
#Data<-markers[gene,]
Idents(Treg)<-"orig.ident"
Treg<-RenameIdents(object = Treg,"Nt9wPP"="Nt9wPP","At9wPP"="At9wPP","Nt9wSP"="Nt9wSP","At9wSP"="At9wSP")
VlnPlot(Treg,features = c("Foxp3","Gata3","Ccr7","Sell","Ctla4","Icos","Il2ra","Nrp1","Tnfrsf4","Tnfrsf9","Tnfrsf18","Tgfb1"),pt.size = 0,stacked = T,group.by = "orig.ident")
Data<-Data[gene,]
colnames(Data)<-gsub("\\.","-",colnames(Data))
ClusterID<- lapply(split(Treg@meta.data, list(Treg@meta.data$orig.ident)), function(x)rownames(x))
#ClusterID<- lapply(split(Treg@meta.data, list(Treg@meta.data$old_Ident)), function(x)rownames(x))
#ClusterID1<-ClusterID$Nt12DSP
#ClusterID2<-ClusterID$Nt9wSP
#ClusterID1<-ClusterID$Nt12DSP
#ClusterID2<-ClusterID$Nt9wSP

ClusterID <- lapply(ClusterID, function(x) ifelse(!is.na(x), paste0("X", x), x))
#ClusterID$Nt12DSP<-ClusterID1
#ClusterID$Nt9wSP<-ClusterID2
#ClusterMean<-aggregate()
#ClusterMean<-vag_exp$RNA
#ClusterMean<-scale(ClusterMean)
#ClusterMean<-as.matrix(ClusterMean)
ClusterMean <- t(apply(Data,1,function(x){
        lapply(ClusterID,function(ID){
                mean(x[ID]) 
        }) %>% unlist
}))
ClusterMean1<-ClusterMean[,1:2]
ClusterMean2<-ClusterMean[,3:4]

#使用默认渐变色画热图；
pheatmap::pheatmap(ClusterMean2, cluster_row = FALSE,cluster_col=FALSE )
mycolor <- c( 'lightgrey', 'darkorange1', 'red')
DotPlot(Treg,features = c("Ccr2","Cxcr3","Itgb1","Ccr4","Gpr15","Itgal","Itga4","Ccr5","Ccr6","Ccr8","Itgb2","Itgb7"),group.by = "orig.ident")+coord_flip()
DotPlot(Treg,features = c("Cxcr5","Pdcd1","Nt5e","Lag3","Ccr2","Cxcr3"),cols = mycolor)+coord_flip()
DotPlot(Treg,features = c("Tbx21","Ifng","Cxcr3"),cols = mycolor,group.by = "orig.ident")+coord_flip()
DotPlot(Treg,features = c("Pdcd1","Cxcr5","Icos","Tox","Tox2","Cd40lg","Il4","Il21"),cols = mycolor,group.by = "orig.ident",dot.scale = 30)
DotPlot(Treg,features = c("S1pr4","S1pr1","Ccr7","Sell","Klf2"),cols = mycolor,group.by = "orig.ident",dot.scale = 20)+coord_flip()
DotPlot(Tregcccc,features = rev(gene),cols = mycolor,group.by = "orig.ident",dot.scale = 20)+coord_flip()
DotPlot(Treg,features = rev(gene),cols = mycolor,group.by = "orig.ident")+coord_flip()
DotPlot(Treg,features = c("Batf","Rora","Il22"),cols = mycolor,dot.scale =30 )+coord_flip()
