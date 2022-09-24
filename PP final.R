library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(cowplot)
library(ggplot2)
library(MySeuratWrappers)
Sample2.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix") 
Sample1.data<-Read10X(data.dir = "~/mm10")
# Add meta data of group
AdultPP<- CreateSeuratObject(counts =  Sample2.data, project = "APP9W", min.cells = 3, min.features = 200)
NeonatalPP <- CreateSeuratObject(counts =  Sample1.data, project = "NPP9W", min.cells = 3, min.features = 200)
PP.combined <- merge(NeonatalPP, y = AdultPP, add.cell.ids = c("NPP9W", "APP9W"), project = "PBMC12K")
#reduce batch effect£¬not change results
counts <- GetAssayData(PP.combined , assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('Gnb2l1','Gltscr2',"Gm8730","Rpl23a-ps3","Gm9843","Gm10036","Sepw1","Rpl13-ps3","Rps12-ps3","Tceb2","Hmha1","Atp5o","Shfm1","Myeov2","2700060E02Rik","Selk","Tonsl","Rps18-ps3","Wdr89","Nhp2l1","Sep15","Fam101b","Sepp1","Gm10073","Whsc1l1","Zfos1","Gm10076","Gm8797","Gm26917"))),]
PP.combined  <- subset(PP.combined , features = rownames(counts))
PP.combined[["percent.mt"]] <- PercentageFeatureSet(PP.combined, pattern = "^mt-")
PP.combined[["percent.B"]] <- PercentageFeatureSet(PP.combined, pattern = "Cd79")
PP.combined[["percent.Rps"]] <- PercentageFeatureSet(PP.combined, pattern = "^Rps|^Rpl|^Gm")
PP.combined[["percent.N"]] <- PercentageFeatureSet(PP.combined, pattern = "Gnb2l1|Gltscr2|Sepw1|Tceb2|Hmha1|Atp5o|Shfm1|Hn1|Myeov2|2310036022Rik|Tonsl|2700060E02Rik|Selk|Wdr89|Nhp2l1|Sep15|Sepp1")
PP.combined[["percent.CD8"]] <- PercentageFeatureSet(PP.combined, pattern = "Cd8a|Cd8b1")
VlnPlot(PP.combined1, features = c("nFeature_RNA", "percent.Rps", "percent.mt","percent.CD8"), ncol = 4, pt.size = 0.0001)
plot1 <- FeatureScatter(PP.combined1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PP.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
PP.combined <- subset(PP.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.B < 0.02& percent.CD8 < 0.02)
PP.combined <- NormalizeData(object = PP.combined, normalization.method = "LogNormalize", scale.factor = 10000)
PP.combined <- FindVariableFeatures(PP.combined, selection.method = "vst", nfeatures = 2000)
PP.combined <- ScaleData(PP.combined, verbose = FALSE)
PP.combined <- ScaleData(object = PP.combined, vars.to.regress = 'percent.mt')
PP.combined <- RunPCA(PP.combined, verbose = FALSE)
PP.combined <- PP.combined %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
PP.combined <- PP.combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution =3) %>% 
  identity()
PP.combined<-subset(x = PP.combined, idents =c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","25","26","27","28","29"))
PP.combined <- PP.combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution =3.4) %>% 
  identity()#3.4
PP.combined<-subset(x = PP.combined, idents =c("0","1","2","3","4","5","6","7","9","10","11","12","13","14","16","17","18","19","20","21","22","23","24","25","26","27","28","29"))

PP.combined<-RenameIdents(object = PP.combined, '0' = 'Tfh', '1' = 'Tfh', '2' = 'Tn', '3' = 'T act', '4' = 'Tn', '5' = 'Tfh', '6' = 'Tfh','7' = 'Tfh', '9' = 'Treg','10' = 'T act', '11' = 'Egr+ T', '12' = 'Tfh', '13' = 'T act', '14' = 'Ikzf2L Treg',  '16' = 'Nme1+ T','17' = 'Tfh', '18' = 'T act','19' = 'Isg+ T', '20' = 'T act', '21' = 'Th1', '22' = 'Tm17', '23' = 'Isg+ T', '24' = 'Nme1+ T', '25' = 'Tm17', '26' = 'Tm', '27' = 'Th17', '28' = 'Pro T', '29' = 'Ccr2+ Treg')
PP.combined<-RenameIdents(object = PP.combined, 'Tn' = 'Tn','Treg' = 'Treg', 'Tm' = 'Tm', 'Th1' = 'Th1',  'Th17' = 'Th17', 'Tfh' = 'Tfh','Pre-Tfh' = 'Tfh', 'Pro T' = 'Unknown')
PP.combined<-RenameIdents(object = PP.combined, 'Tn' = 'Tn','Treg' = 'Treg','Tm' = 'Tm','Th1' = 'Th1', "Th17"="Th17","Tfh"="Tfh","Unknown"="Unknown")
SP9wcol<-c('#000EFF', '#53A85F', '#F1BB72','#BD956A',  "#476D87",  "#E95C59",'#AB3282')

markers <- c( 'Cd3g', 'Sell', 'Lef1', 'Foxp3', 'Ikzf2', 'Ccr2',"Tnfsf8", 'Egr1',"Egr2","Tmem176a","Tmem176b", 'S100a6', 'S100a4',"Ccl5","Nkg7", 'Birc5')
DimPlot(PP.combined, reduction = "umap", pt.size = .1,label = F,cols = PP9wcol)+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) +guides(color=F)
DimPlot(PP.combined, reduction = "umap", pt.size = .1,label = T)


PP.combined@meta.data$new_Ident<-PP.combined@active.ident
Idents(object = PP.combined1) <- 'orig.ident'
PP.combined<-subset(PP.combined, downsample=2912)
PP.combined1<-subset(PP.combined, ident="0")
Idents(object = PP.combined1) <- 'new_Ident'
PP.combined1<-subset(PP.combined, ident="1")
Idents(object = PP.combined1) <- 'new_Ident'
FeaturePlot(object=PP.combined1,cols = mycolor,features = "Bcl6",pt.size = 0.1)&annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)&annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) &guides(color=F)&labs(title ="")
FeaturePlot(object=PP.combined1,cols = mycolor,features = "Nrn1",pt.size = 0.1)&annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)&annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) &guides(color=F)&labs(title ="")










