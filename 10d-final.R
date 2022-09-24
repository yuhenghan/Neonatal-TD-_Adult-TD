#10d spleen
library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(ggplot2)
library(Rmisc)
library(car)
library(purrr)
library(ggpubr)
Sample2.data <- Read10X(data.dir = "~s/filtered_feature_bc_matrix") 
Sample1.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix") 

NeonatalSP <- CreateSeuratObject(counts =  Sample1.data, project = "0", min.cells = 3, min.features = 200)
AdultSP<- CreateSeuratObject(counts =  Sample2.data, project = "1", min.cells = 3, min.features = 200)
SP10d <- merge(NeonatalSP, y = AdultSP, add.cell.ids = c("neonatal12d", "adult12d"), project = "PBMC12K")
SP10d[["percent.mt"]] <- PercentageFeatureSet(SP10d, pattern = "^mt-")
SP10d[["percent.B"]] <- PercentageFeatureSet(SP10d, pattern = "Cd79")
SP10d[["percent.Rps"]] <- PercentageFeatureSet(SP10d, pattern = "^Rps|^Rpl|^Gm")
SP10d[["percent.N"]] <- PercentageFeatureSet(SP10d, pattern = "Gnb2l1|Gltscr2|Sepw1|Tceb2|Hmha1|Atp5o|Shfm1|Hn1|Myeov2|2310036022Rik|Tonsl|2700060E02Rik|Selk|Wdr89|Nhp2l1|Sep15|Sepp1")
SP10d[["percent.CD8"]] <- PercentageFeatureSet(SP10d, pattern = "Cd8a|Cd8b1")
#B.genes <- grep(pattern = "^Cd79", x = rownames(x = SP10d), value = TRUE)
#percent.B <- Matrix::colSums(x = GetAssayData(object = SP10d, slot = 'counts')[B.genes, ]) / Matrix::colSums(x = GetAssayData(object = SP10d, slot = 'counts'))
#SP10d[['percent.RPS']] <- percent.B
VlnPlot(PP.combined, features = c("nFeature_RNA", "percent.Rps", "percent.mt"), ncol = 4, pt.size = 0.0001)
plot1 <- FeatureScatter(SP10d, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SP10d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
SP10d <- subset(SP10d, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5 )
SP10d <- NormalizeData(object = SP10d, normalization.method = "LogNormalize", scale.factor = 10000)
SP10d <- FindVariableFeatures(SP10d, selection.method = "vst", nfeatures = 2000)
SP10d <- FindVariableFeatures(SP10d, selection.method = "vst", nfeatures = 14564)
SP10d <- ScaleData(SP10d, verbose = FALSE)
SP10d <- ScaleData(object = SP10d, vars.to.regress = 'percent.mt')
SP10d <- RunPCA(SP10d, verbose = FALSE)
#options(repr.plot.height = 2.5, repr.plot.width = 6)
SP10d <- SP10d %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
SP10d <- SP10d %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 2) %>% identity() 
SP10d<-subset(x = SP10d, idents =c("0","1","2","3","4","5","6","7","8","9","10","11","13","14","15","16","17","18","19","20","21"))
SP10d <- SP10d %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1.5) %>% identity() 

SP10d <- SP10d %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 3) %>% identity() 
SP10d<-subset(x = SP10d, idents =c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","16","17","18","19","20","21","22","23","24","25","26","27"))
  SP10d<-RenameIdents(object = SP10d,'0'  ='Tn' ,'1'  ='Tn' , '2' = 'Tn','3'= 'Tn','4'='Unknown','5'='Tn','6' ='Tn' , '7'= 'Tn' ,'8'='Unknown' ,'9'='Tn'  ,'10'='Treg'  ,'11'='Tn'  ,'12'='Tn'  ,'13' ='Tn','14'='Tn' ,'16'='Tn', '17'='Th1'  ,'18'  ='Unknown'  , '19'='Unknown' , '20'='Treg' , '21'='Unknown', '22'='Tm', '23'='Th17', '24'='Unknown', '25'='Unknown', '26'='Th1', '27'='Tn')
SP10d<-RenameIdents(object = SP10d,'0'  ='Klf3L Tn' ,'1'  ='Klf3L Tn' , '2' = 'Klf3L Tn','3'= 'Klf3H Tn','4'='Klf3H Tn','5'='Myb+ T','6' ='Klf3L Tn' , '7'= 'Ikzf2H Treg' ,'8'='Klf3H Tn' ,'9'='Klf3H Tn'  ,'10'='Klf3H Tn'  ,'11'='Klf3H Tn'  ,'12'='Klf3H Tn'  ,'13' ='Isg+ T','14'='Klf3H Tn' ,'15'='Klf3H Tn','16'='Nme1+ T', '17'='mtH T'  ,'18'  ='Klf3L Tn'  , '19'='Th1' , '20'='Ikzf2L Treg' , '21'='Egr+ T', '22'='Isg+ T', '23'='Pro T', '24'='Tm', '25'='Pro T', '26'='Pro T', '27'='Klf3L Tn')
SP10d<-RenameIdents(object = SP10d,'Pro T' = 'Unknown','Isg+ T' = 'Unknown','Nme1+ T' = 'Unknown','Egr+ T' = 'Unknown','Myb+ T' = 'Unknown')
SP10d<-RenameIdents(object = SP10d, 'Klf3L Tn' = 'Tn', 'Klf3H Tn' = 'Tn','Ikzf2H Treg' = 'Treg','IkzfL Treg' = 'Treg','Tm' = 'Tm','Th1' = 'Th1')
SP10d<-subset(x = SP10d, idents =c("Tn","Treg","Tm","Th1","Unknown"))
SP10dcol<-c('#000EFF', '#53A85F', '#F1BB72','#768A00', '#E59CC4', '#AB3282')
DimPlot(SP10d, reduction = "umap", pt.size = .1,label = F,cols = SP10dcol)#+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) +guides(color=F)
mycolor <- c( 'lightgrey', 'darkorange1', 'red') #设置颜色 大小6:6
FeaturePlot(object=SP10d,cols = mycolor,features = "Cd3g",pt.size = 0.5)&annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)&annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) &guides(color=F)&labs(title ="")
FeaturePlot(object=SP9w,cols = mycolor,features = "Gzmb",pt.size = 0.1)&annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)&annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) &guides(color=F)&labs(title ="")
FeaturePlot(object=SP10d,cols = mycolor,features = "Tbx21",pt.size = 0.1)&annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)&annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) &guides(color=F)&labs(title ="")
FeaturePlot(object=SP10d,cols = mycolor,features = "Tnfrsf9",pt.size = 0.1)&annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)&annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) &guides(color=F)&labs(title ="")
FeaturePlot(object=SP10d,cols = mycolor,features = "Mif",pt.size = 0.1)&annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)&annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) &guides(color=F)&labs(title ="")
FeaturePlot(object=SP10d,cols = mycolor,features = "Ifng",pt.size = 0.1)&annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)&annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) &guides(color=F)&labs(title ="")
DimPlot(SP10d, reduction = "umap", pt.size = .1,label = F,cols = myidentcol,group.by = "orig.ident")+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) +guides(color=F)

DimPlot(SP10d, reduction = "umap", pt.size = .1,label = T)
+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) +guides(color=F)






