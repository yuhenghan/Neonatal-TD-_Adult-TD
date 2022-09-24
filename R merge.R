library(ggplot2)
library(showtext) 
library(clusterProfiler)
library(org.Mm.eg.db) ##加载小鼠
library(org.Hs.eg.db) ##加载人类
library(fgsea)
library(dplyr)
library(stringr)
library(msigdbr) #提供MSigdb数据库基因集
library(fgsea)
library(dplyr)
library(tibble)
library(Seurat)
library(dplyr)
library(Seurat)
library(DOSE)
library(patchwork)
  SP9w<-RenameIdents(SP9w, "0"="Nt9wSP","1"="At9wSP")
  PP.combined<-RenameIdents(object = PP.combined,"0"="Nt9wPP","1"="At9wPP")
SP10d<-RenameIdents(object = SP10d,"0"="Nt12DSP","1"="At12DSP")
SP9w@meta.data$old_Ident<-SP9w@active.ident
SP10d@meta.data$old_Ident<-SP10d@active.ident
PP.combined@meta.data$old_Ident<-PP.combined@active.ident
Tregcccc <- merge(SP9w, y = c(PP.combined),add.cell.ids = c("neonatal12d", "adult12d"))
Tregcccc@meta.data$orig.ident<-Tregcccc@active.ident
Tregcccc <- NormalizeData(object = Tregcccc, normalization.method = "LogNormalize", scale.factor = 10000)
Tregcccc <- FindVariableFeatures(Tregcccc, selection.method = "vst", nfeatures = 5000)
Tregcccc <- ScaleData(Tregcccc, verbose = FALSE)
Tregcccc <- ScaleData(object = Tregcccc, vars.to.regress = 'percent.mt')
Tregcccc <- RunPCA(Tregcccc, verbose = FALSE)
Tregcccc <- Tregcccc %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
#Treg2 <- Treg2 %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE,sigma=10,theta=20,lamba=20)
   Tregcccc <- Tregcccc %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution =1) %>% 
  identity()#3.4
   #sp+pp
  Tregcccc <- Tregcccc %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution =2) %>% 
    identity()#3.4
  #sp9w+12d
  Tregcccc2 <- merge(SP9w, y = c(SP10d),add.cell.ids = c("neonatal12d", "adult12d"))
  Tregcccc2@meta.data$orig.ident<-Tregcccc2@active.ident
  Tregcccc2 <- NormalizeData(object = Tregcccc2, normalization.method = "LogNormalize", scale.factor = 10000)
  Tregcccc2 <- FindVariableFeatures(Tregcccc2, selection.method = "vst", nfeatures = 5000)
  Tregcccc2 <- ScaleData(Tregcccc2, verbose = FALSE)
  Tregcccc2 <- ScaleData(object = Tregcccc2, vars.to.regress = 'percent.mt')
  Tregcccc2 <- RunPCA(Tregcccc2, verbose = FALSE)
  Tregcccc2 <- Tregcccc2 %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
  Tregcccc2 <- Tregcccc2 %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution =2.2) %>% 
    identity()#3.4 
   Tregcccc<-RenameIdents(object = Tregcccc,'0'  ='Tn' ,'1'  ='Tn' , '2' = 'Tn','3'= 'Tn','4'='Tn','5'='Tm','6' ='Tm' , '7'= 'Tn' ,'8'='Tn' ,'9'='Unknown'  ,'10'='Tn'  ,'11'='Tn'  ,'12'='Th1'  ,'13' ='Treg','14'='Th1','15'='Th17' ,'16'='Tfh', '17'='Treg'  ,'18'  ='Unknown'  , '19'='Treg' , '20'='Unknown' , '21'='Tn', '22'='Treg', '23'='Unknown', '24'='Unknown', '25'='Th1', '26'='Unknown', '27'='Unknown', '28'='Unknown')
  Tregcccc<-RenameIdents(object = Tregcccc,'Tn' = 'Tn','Treg' = 'Treg', 'Tm' = 'Tm', 'Th1' = 'Th1','Th17' = 'Th17', 'Tfh' = 'Tfh','Unknown' = 'Unknown')
  SP9wcol<-c('#000EFF', '#53A85F', '#F1BB72','#BD956A',  "#476D87",  "#E95C59",'#AB3282')
Tregcccc<-RenameIdents(object = Tregcccc, '0' = 'Tn', '1' = 'Tn', '2' = 'Tn', '3' = 'Tn','6' = 'Treg', '4' = 'Tm', '10' = 'Th1','5' = 'Th17',  '7' = 'Tn', '8' = 'Treg', '12' = 'Th1', '14' = 'Tfh', '15' = 'Unknown',  '9' = 'Unknown', '11' = 'Unknown',  '13' = 'Unknown')
Tregcccc2<-RenameIdents(object = Tregcccc2, '0' = 'Tn', '1' = 'Tn', '2' = 'Tfh', '3' = 'Tn', '4' = 'Tn','5' = 'Tfh', '6' = 'Tn', '7' = 'Tfh','8' = 'Treg',  '9' = 'Tm', '10' = 'Treg', '11' = 'Tfh', '12' = 'Th1', '13' = 'Th17', '14' = 'Unknown','15' = 'Unknown', '16' = 'Th1','17' = 'Tfh', '18' = 'Tn','19' = 'Unknown','20' = 'Th17')
Tregcccc2<-RenameIdents(object = Tregcccc2,'Tn' = 'Tn','Treg' = 'Treg', 'Tm' = 'Tm', 'Th1' = 'Th1','Th17' = 'Th17', 'Tfh' = 'Tfh','Unknown' = 'Unknown')
DimPlot(Tregcccc, reduction = "umap", pt.size = .1,label = T,repel =T)
DimPlot(Tregcccc, reduction = "umap", pt.size = .1,label = T,repel = T)+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) +guides(color=F)

Tregcccc@meta.data$new_Ident<-Tregcccc@active.ident
Idents(object = Tregcccc) <- 'orig.ident'
table(Tregcccc@active.ident)
Tregcccc_1<-subset(Tregcccc, downsample=1759)
Tregcccc_2<-subset(Tregcccc, idents=c("Nt12DSP"))
Tregcccc_3<-subset(Tregcccc, idents=c("At12DSP"))
Tregcccc_4<-subset(Tregcccc, idents=c("Nt9wSP"))
Tregcccc_5<-subset(Tregcccc, idents=c("At9wSP"))
Idents(object = Tregcccc_2) <- 'old_Ident'
DimPlot(Tregcccc_2, reduction = "umap", pt.size = .1,label = F,repel = F,cols = SP9wcol)+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) +guides(color=F)
p1
