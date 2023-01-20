library(org.Mm.eg.db) ##¼ÓÔØÐ¡Êó
library(org.Hs.eg.db) ##¼ÓÔØÈËÀà
library(ggplot2)
library(clusterProfiler)
library(dplyr)
library(msigdbr)
library(fgsea)
library(tibble)
library(Seurat)
library(dplyr)
library(Matrix)
markers <- FindMarkers(Treg, ident.1 = "0", ident.2 = "1", logfc.threshold = 0.2,min.pct =0.02)
markers1 <- FindMarkers(Treg1, ident.1 = "0", ident.2 = "1", logfc.threshold = 0.2,min.pct =0.02)
up <-rownames(markers[intersect(which(markers [,1]<0.05),which(markers [,2]>=0.2)),])
down <-rownames(markers[intersect(which(markers [,1]<0.05),which(markers [,2]<=(-0.2))),])
up1 <-rownames(markers[intersect(which(markers1 [,1]<0.05),which(markers [,2]>=0.2)),])
down1 <-rownames(markers[intersect(which(markers1 [,1]<0.05),which(markers [,2]<=(-0.2))),])
gs1 = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
gs2 = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
gs3 = bitr(up1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
gs4 = bitr(down1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ego.bp1 = enrichGO(gene=gs1$ENTREZID, OrgDb = org.Mm.eg.db,ont= "ALL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
ego.bp2 = enrichGO(gene=gs2$ENTREZID, OrgDb = org.Mm.eg.db,ont= "ALL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
ego.bp3 = enrichGO(gene=gs1$ENTREZID, OrgDb = org.Mm.eg.db,ont= "ALL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
ego.bp4 = enrichGO(gene=gs2$ENTREZID, OrgDb = org.Mm.eg.db,ont= "ALL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
a1<-ego.bp@result
a3<-ego.bp3@result
generatio<-c(52/605,20/605,4/605,39/605,6/605)
generatio1<-c(34/363,11/363,9/363,5/363,28/363,5/363,4/363)
a1$GeneRatio<-generatio
a3$GeneRatio<-generatio1
write.csv(a1,file = "Neonatal-TD 9W PP.csv")
write.csv(a3,file = "Adult-TD 9W PP.csv")
#CandidateGOTerms <- c("GO:0042110", "GO:0042098", "GO:0007179", "GO:0002507", "GO:0007159","GO:0045066")
CandidateGOTerms <- c("GO:0042110", "GO:0042098", "GO:0007179", "GO:0002507", "GO:0007159","GO:0045066","GO:0071604")
CandidateGOTerms <- c("GO:0042110", "GO:0042098",  "GO:0002507", "GO:0007159","GO:0045066")
a1<-a1[CandidateGOTerms,]
a3<-a3[CandidateGOTerms,]
p1<-ggplot(a1,aes(y=Count,x=reorder(Description,-Count)))+
  geom_col(aes(fill=-log10(p.adjust)))+theme(axis.line.x=element_line(linetype=1,color="black",size=1),axis.line.y=element_line(linetype=1,color="black",size=1))+theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))                    
p1 <- p1 + geom_line(aes(colour=p.adjust)
p1 <- p1 + scale_colour_continuous(low="blue",high="red")
p1<-ggplot(a1,aes(y=-log10(pvalue),x=reorder(Description,log10(p.adjust))))+
  geom_col(aes(fill=GeneRatio))+theme(axis.line.x=element_line(linetype=1,color="black",size=1),axis.line.y=element_line(linetype=1,color="black",size=1))+theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradient2(low = "#FFFFCC", high = "brown",mid = "orange")                    


                   
