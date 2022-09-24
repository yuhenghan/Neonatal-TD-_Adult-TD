library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(ggplot2)
library(Rmisc)
library(car)
library(purrr)
library(ggpubr)
library(colourpicker)
library(paletteer)
library(scPred)
library(sqldf)
Sample4.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix") 
Sample3.data<-Read10X(data.dir = "~/mm10")
# Add meta data of group
AdultSP<- CreateSeuratObject(counts =  Sample4.data, project = "ASP9W", min.cells = 3, min.features = 200)
NeonatalSP <- CreateSeuratObject(counts =  Sample3.data, project = "NSP9W", min.cells = 3, min.features = 200)
SP9w <- merge(NeonatalSP, y = AdultSP, add.cell.ids = c("NSP9W", "ASP9W"), project = "PBMC12K")
#reduce batch effect，not change results
#counts <- GetAssayData(PP.combined, assay = "RNA")
#counts <- counts[-(which(rownames(counts) %in% c('Gnb2l1','Gltscr2',"Gm8730","Rpl23a-ps3","Gm9843","Gm10036","Sepw1","Rpl13-ps3","Rps12-ps3","Tceb2","Hmha1","Atp5o","Shfm1","Myeov2","2700060E02Rik","Selk","Tonsl","Rps18-ps3","Wdr89","Nhp2l1","Sep15","Fam101b","Sepp1","Gm10073","Whsc1l1","Zfos1","Gm10076","Gm8797","Gm26917"))),]
#SP9w <- subset(SP9w, features = rownames(counts))
SP9w[["percent.mt"]] <- PercentageFeatureSet(SP9w, pattern = "^mt-")
SP9w[["percent.B"]] <- PercentageFeatureSet(SP9w, pattern = "Cd79")
SP9w[["percent.Rps"]] <- PercentageFeatureSet(SP9w, pattern = "^Rps|^Rpl|^Gm")
SP9w[["percent.CD8"]] <- PercentageFeatureSet(SP9w, pattern = "Cd8a|Cd8b1")
VlnPlot(SP9w, features = c("nFeature_RNA", "percent.Rps", "percent.mt","percent.CD8"), ncol = 4, pt.size = 0.0001)
 FeatureScatter(SP9w, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
plot2 <- FeatureScatter(SP9w, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
SP9w <- subset(SP9w, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 8 & percent.B < 0.02& percent.CD8 < 0.02)
SP9w <- NormalizeData(object = SP9w, normalization.method = "LogNormalize", scale.factor = 10000)
SP9w <- FindVariableFeatures(SP9w, selection.method = "vst", nfeatures = 2000)
#SP9w <- FindVariableFeatures(SP9w, selection.method = "vst", nfeatures = 14319)
SP9w <- ScaleData(SP9w, verbose = FALSE)
SP9w <- ScaleData(object = SP9w, vars.to.regress = 'percent.mt')
SP9w <- RunPCA(SP9w, verbose = FALSE)
#options(repr.plot.height = 2.5, repr.plot.width = 6)
SP9w <- SP9w %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
SP9w <- SP9w %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.9) %>% 
  identity()
SP9w<-subset(x = SP9w, idents =c("0","1","2","3","4","5","6","7","8","9","10","16","17"))
SP9w <- SP9w %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 2) %>% 
  identity()
DimPlot(SP9w,reduction = "umap")
SP9w<-RenameIdents(object = SP9w, '0' = 'Klf3L Tn', '1' = 'Klf3L Tn', '2' = 'Klf3L Tn', '3' = 'Klf3L Tn', '4' = 'Klf3L Tn', '5' = 'Tm', '6' = 'Tm', '7' = 'Klf3H Tn', '8' = 'Klf3L Tn', '4' = 'Treg', '9' = 'Ikzf2H Treg', '10' = 'Th1', '11' = 'Ikzf2L Treg', '12' = 'Tfh', '13' = 'Th17', '14' = 'Isg+ T', '15' = 'Th1', '16' = 'Ccr2+ Treg', '17' = 'Nme+ T', '18' = 'Pro T')
SP9w<-RenameIdents(object = SP9w,'Pro T' = 'Unknown','Isg+ T' = 'Unknown','Nme+ T' = 'Unknown')
SP9w<-RenameIdents(object = SP9w,'Tm' = 'Tm','Th1' = 'Th1','Th1' = 'Th1','Th17' = 'Th17','Tfh' = 'Tfh')
SP9w<-RenameIdents(object = SP9w, 'Klf3L Tn' = 'Tn', 'Klf3H Tn' = 'Tn','Ikzf2H Treg' = 'Treg','Ikzf2L Treg' = 'Treg','Ccr2+ Treg' = 'Treg')
SP9w<-RenameIdents(object = SP9w, 'Tn' = 'Tn','Treg' = 'Treg', 'Tm' = 'Tm', 'Th1' = 'Th1','Th1' = 'Th1',  'Th17' = 'Th17', 'Tfh' = 'Tfh', 'Unknown' = 'Unknown')
SP9wcol<-c('#000EFF', '#53A85F', '#F1BB72','#BD956A',  "#476D87",  "#E95C59",'#AB3282')
DimPlot(SP9w, reduction = "umap", pt.size = .1,label = F,repel = T,cols = SP9wcol)+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) +guides(color=F)
DimPlot(SP9w, reduction = "umap", pt.size = .1,label = F,cols = myidentcol,group.by = "orig.ident")+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) +guides(color=F)+ labs(title ="")#无标题，加框，去注释框

