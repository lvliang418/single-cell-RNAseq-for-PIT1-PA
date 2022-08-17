library(dplyr)
library(Seurat)
library(patchwork)
#step0 scCancer QC
#step1 further check

B1 <- readRDS("/my_path/scCancerQC/B1/scAnnotation/expr.RDS")
B2 <- readRDS("my_path/scCancerQC/B2/scAnnotation/expr.RDS")
B3 <- readRDS("my_path/scCancerQC/B3/scAnnotation/expr.RDS")
B4 <- readRDS("my_path/scCancerQC/B4/scAnnotation/expr.RDS")

B1<-subset(B1, nFeature_RNA <6000 & nFeature_RNA > 200)
dim(B1)
B1<-subset(B1, nCount_RNA <20000)
dim(B1)
B1<-subset(B1, mito.percent<0.25)
dim(B1)

B2<-subset(B2, nFeature_RNA <5600 & nFeature_RNA > 200)
dim(B2)
B2<-subset(B2, nCount_RNA <20000)
dim(B2)
B2<-subset(B2, mito.percent<0.35)
dim(B2)


B3<-subset(B3, nFeature_RNA <6000 & nFeature_RNA > 200)
dim(B3)
B3<-subset(B3, nCount_RNA <20000)
dim(B3)

B4<-subset(B4, nFeature_RNA <5000 & nFeature_RNA > 200)
dim(B4)
B4<-subset(B4, nCount_RNA <17187)
dim(B4)
B4<-subset(B4, mito.percent<0.25)
dim(B4)

#integration
seuratlist<-list(B1,B2,B3,B4)
for (i in 1:length(seuratlist)) {
  seuratlist[[i]]<-NormalizeData(seuratlist[[i]])
  seuratlist[[i]]<-FindVariableFeatures(seuratlist[[i]],nfeatures = 2000,selection.method = "vst")
}

#CCA-MMN batch effect
PA.anchors <- FindIntegrationAnchors(object.list = seuratlist, anchor.features= 2000,dims = 1:30)
rm(seuratlist)
gc()

PA.integrated <-IntegrateData(anchorset = PA.anchors, dims = 1:30)
dim(PA.integrated)
save(PA.integrated,file = "PA.integrated2.RData")

library(ggplot2)
ggplot(PA.integrated@meta.data, aes(doublet.score))+
  geom_histogram(binwidth = 0.001)
ggplot(PA.integrated@meta.data, aes(mito.percent))+
  geom_histogram(binwidth = 0.001)
ggplot(PA.integrated@meta.data, aes(ribo.percent))+
  geom_histogram(binwidth = 0.001)
ggplot(PA.integrated@meta.data, aes(nCount_RNA))+
  geom_histogram(binwidth = 100)
ggplot(PA.integrated@meta.data, aes(nFeature_RNA))+
  geom_histogram(binwidth = 10)
PA.subset<-subset(PA.integrated, subset= doublet.score <0.25)

#RBC removal
hemoglobin<-c("HBB","HBA1","HBG1","HBG2","HBA2","HBD","HBE1")
PA.subset<-AddModuleScore(PA.subset, features = list(hemoglobin), assay = "RNA", name = "RBC_score", search = T)
library(ggplot2)
ggplot(PA.subset@meta.data, aes(RBC_score1))+
  geom_histogram(binwidth = 0.1)+geom_vline(xintercept = 0.055)

data<-data.frame(RBC_score=PA.subset@meta.data$RBC_score1,
                 row.names = rownames(PA.subset@meta.data))

data$RBC<-data$RBC_score > 0.055
highCells<-rownames(subset(data, RBC==T))
cell<-setdiff(Cells(PA.subset), highCells)
save(PA.subset,file="PA.subset.RBC.RData")

PA.subset<-subset(PA.subset, cells = cell)

save(PA.subset,file="PA.subset2.RData")
dim(PA.subset)

rm(B1)
rm(B2)
rm(B3)
rm(B4)
rm(PA.anchors)
gc()

DefaultAssay(PA.subset)<-"integrated"
PA.subset <- ScaleData(PA.subset)

#PCA
PA.subset <- RunPCA(PA.subset, features = VariableFeatures(object = PA.subset))
ElbowPlot(PA.subset) 
DimHeatmap(PA.subset, dims = 1:15, cells = 100, balanced = TRUE)
PA.subset <- FindNeighbors(PA.subset, dims = 1:9)
PA.subset <- FindClusters(PA.subset, resolution = 2)

#tsne umap
PA.subset <- RunUMAP(PA.subset, dims = 1:9)
DimPlot(PA.subset, reduction = "umap", label = T)

PA.subset <- RunTSNE(PA.subset, dims = 1:9)
DimPlot(PA.subset, reduction = "tsne", label = T)
save(PA.subset,file="PA.subset2.RData")

#clinical information
barcode<-data.frame(barcode=rownames(PA.subset@meta.data),
                    ID=substring(rownames(PA.subset@meta.data), 18,18),
                    Sample_ID=substring(rownames(PA.subset@meta.data), 18,18),
                    Pathologic_diagnosis=substring(rownames(PA.subset@meta.data), 18,18),
                    PRL_staining=substring(rownames(PA.subset@meta.data), 18,18),
                    TSH_staining=substring(rownames(PA.subset@meta.data), 18,18),
                    Ki67_index=substring(rownames(PA.subset@meta.data), 18,18),
                    Gonad_Function=substring(rownames(PA.subset@meta.data), 18,18),
                    Thyroid_Function=substring(rownames(PA.subset@meta.data), 18,18),
                    PRL_level=substring(rownames(PA.subset@meta.data), 18,18),
                    OGTT_Test=substring(rownames(PA.subset@meta.data), 18,18))
Sample_ID<-c("B1","B2","B3","B4")
Pathologic_diagnosis<-c("Mixed somatotroph-lactotroph adenoma",
                        "Densely granulated somatotroph adenoma",
                        "Sparsely granulated somatotroph adenoma",
                        "Plurihormonal adenoma")
PRL_staining<-c("Positive","Negative","Negative","Positive")
TSH_staining<-c("Negative","Negative","Negative","Positive")
Ki67_index<-c(2.5,1.5,1.5,1)
Gonad_Function<-c("Hypogonadism","Normal","Normal","Normal")
Thyroid_Function<-c("Hypothyroidism","Normal","Normal","Hyperthyroidism")
PRL_level<-c("Hyperprolactinemia","Normal","Normal","Normal")
OGTT_Test<-c("Positve","Positve","Positve","Negative")
ID<-c(1,2,3,4)

barcode$Sample_ID<-plyr::mapvalues(barcode$ID, from = ID, to = Sample_ID)
barcode$Pathologic_diagnosis<-plyr::mapvalues(barcode$ID, from = ID, to = Pathologic_diagnosis)
barcode$PRL_staining<-plyr::mapvalues(barcode$ID, from = ID, to = PRL_staining)
barcode$TSH_staining<-plyr::mapvalues(barcode$ID, from = ID, to = TSH_staining)
barcode$Ki67_index<-plyr::mapvalues(barcode$ID, from = ID, to = Ki67_index)
barcode$Gonad_Function<-plyr::mapvalues(barcode$ID, from = ID, to = Gonad_Function)
barcode$Thyroid_Function<-plyr::mapvalues(barcode$ID, from = ID, to = Thyroid_Function)
barcode$PRL_level<-plyr::mapvalues(barcode$ID, from = ID, to = PRL_level)
barcode$OGTT_Test<-plyr::mapvalues(barcode$ID, from = ID, to = OGTT_Test)
rownames(barcode)<-barcode$barcode
barcode<-barcode[,-c(1,2)]

PA.subset<-AddMetaData(PA.subset,barcode)
save(PA.subset,file="PA.subset2.RData")

PA.markers <- FindAllMarkers(PA.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(PA.markers, "cluster_markers.csv")
top3<-PA.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
genes_to_check <- unique(top3$gene)
tab1 <- cbind(as.data.frame(PA.subset@meta.data$Sample_ID),as.data.frame(PA.subset@active.ident))
colnames(tab1) <- c("Sample", "Cluster")
library(ggplot2)
pdf("DEgenes_Dotplot.pdf", width = 9, height = 12)
DotPlot(PA.subset, assay= "RNA", features = genes_to_check) + coord_flip()+
  theme(panel.border = element_rect(colour="black",fill=NA))
ggplot(tab1) +
  aes(x = Cluster, fill = factor(Sample)) +
  geom_bar(position = "fill")+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

pdf("DimPlot.pdf", onefile = T, width = 10, height = 8)
DimPlot(PA.subset, reduction = "tsne", group.by= "Sample_ID",label = F, pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
DimPlot(PA.subset, reduction = "tsne", group.by= "Malign.type",label = F, pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
DimPlot(PA.subset, reduction = "tsne", group.by= "PRL_staining",label = F, pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
DimPlot(PA.subset, reduction = "tsne", group.by= "TSH_staining",label = F, pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
DimPlot(PA.subset, reduction = "tsne", label = T,  pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
DimPlot(PA.subset, reduction = "umap", label = T, pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
DimPlot(PA.subset, reduction = "tsne", group.by= "Cell.Type",label = F, pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

pdf("FeaturePlot.pdf", onefile = T, width = 9, height = 8)
FeaturePlot(PA.subset, reduction = "tsne", features = "doublet.score", pt.size = 1.5, cols = c("#FFFAF6", "#50802F"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.subset, reduction = "tsne", features = "Malign.score", pt.size = 1.5,cols = c("#FFFAF6", "#F57E87"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.subset, reduction = "tsne", features = "CellCycle.score", cols = c("#FFFAF6", "#139D49"), pt.size = 1.5,)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.subset, reduction = "tsne", features = "Stemness.score", cols = c("#FFFAF6", "#FF920E"), pt.size = 1.5)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

pdf("QCPlot.pdf", onefile = T, width = 9, height = 8)
FeaturePlot(PA.subset, reduction = "tsne", features = "nCount_RNA", pt.size = 0.6,cols = c("#FFFAF6", "#50802F"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.subset, reduction = "tsne", features = "nFeature_RNA", pt.size = 0.6, cols = c("#FFFAF6", "#F57E87"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.subset, reduction = "tsne", features = "mito.percent", pt.size = 0.6,cols = c("#FFFAF6", "#139D49"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.subset, reduction = "tsne", features = "ribo.percent", pt.size = 0.6,cols = c("#FFFAF6", "#FF920E"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.subset, reduction = "tsne", features = "diss.percent", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()


library(ggplot2)
DefaultAssay(PA.subset)<-"RNA" 
#tsne
p1<- FeaturePlot(PA.subset, reduction = "tsne", features = "EPCAM", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<- FeaturePlot(PA.subset, reduction = "tsne", features = "PLVAP", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<- FeaturePlot(PA.subset, reduction = "tsne", features = "COL1A2", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p4<-FeaturePlot(PA.subset, reduction = "tsne", features = "PTPRC", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA)) 
p5<-FeaturePlot(PA.subset, reduction = "tsne", features = "CD4", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA)) 
p6<-FeaturePlot(PA.subset, reduction = "tsne", features = "CD8A", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA)) 
p7<-FeaturePlot(PA.subset, reduction = "tsne", features = "IGKC", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA)) 
p8<-FeaturePlot(PA.subset, reduction = "tsne", features = "LYZ", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA)) 
pdf("FourMainClusters_combine.pdf", onefile = T, width = 18, height = 16)
CombinePlots(list(p1,p2,p3,p4,p5,p6,p7,p8), ncol = 3, legend = "right")
dev.off()

FeaturePlot(PA.subset, reduction = "tsne", features = "SOX2", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))


library(RColorBrewer)
library(ggplot2)
#violin
p1<-VlnPlot(PA.subset, features = "EPCAM", pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<-VlnPlot(PA.subset, features = "PLVAP",pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-VlnPlot(PA.subset, features = "COL1A2", pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p4<-VlnPlot(PA.subset, features = "PTPRC", pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p5<-VlnPlot(PA.subset, features = "CD4", pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p6<-VlnPlot(PA.subset, features = "CD8A", pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p7<-VlnPlot(PA.subset, features = "IGKC", pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p8<-VlnPlot(PA.subset, features = "LYZ", pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
pdf("FourMainClusters_combine_violin.pdf", width = 20, height = 20)
CombinePlots(list(p1,p2,p3,p4,p5,p6,p7,p8), ncol = 2, legend = "right")
dev.off()

VlnPlot(PA.subset, features = "PTTG1", pt.size = 0)+
  geom_boxplot(width=0.1,alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))

pdf("Pit-1.pdf", onefile = T, width = 10, height = 8)
FeaturePlot(PA.subset, reduction = "tsne", features = "POU1F1", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA)) 
FeaturePlot(PA.subset, reduction = "tsne", features = "GH1", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA)) 
FeaturePlot(PA.subset, reduction = "tsne", features = "PRL", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.subset, reduction = "tsne", features = "TSHB", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

#clustering
main.annotation<- c("Epithelial","Epithelial","Epithelial","Epithelial",
                    "Epithelial","Epithelial","Epithelial","Immune",
                    "Epithelial", "Immune","Epithelial","Epithelial",
                    "Epithelial","Epithelial","Immune","Epithelial", 
                    "Epithelial","Immune","Immune", "Immune",
                    "Immune","Immune","Epithelial","Fibroblast",
                    "Immune","Immune","Epithelial","Epithelial",
                    "Epithelial","Doublet","Endothelial")
cluster.ids <- sort(as.numeric(unique(as.character(PA.subset@meta.data$seurat_clusters))))
PA.subset@meta.data$main.annotation <- PA.subset@meta.data$seurat_clusters
PA.subset@meta.data$main.annotation <- plyr::mapvalues(x = PA.subset@meta.data$main.annotation, from = cluster.ids, to = main.annotation)
PA.subset@meta.data$main.annotation<-factor(PA.subset@meta.data$main.annotation,
                                            levels = c("Epithelial","Immune","Fibroblast","Endothelial","Doublet"))

color<-c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0")
pdf("Main.annotation.pdf", onefile = T, width = 10, height = 8)
DimPlot(PA.subset, reduction = "tsne", label = FALSE, group.by = "main.annotation", cols = color, pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()
save(PA.subset, file = "PA.subset2.RData")

p1<-VlnPlot(PA.subset, features = c("EPCAM"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
p2<-VlnPlot(PA.subset, features = c("KRT8"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
p3<-VlnPlot(PA.subset, features = c("KRT18"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)

p4<-VlnPlot(PA.subset, features = c("PTPRC"), pt.size = 0,group.by = "main.annotation",  
        cols =c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
p5<-VlnPlot(PA.subset, features = c("NKG7"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
p6<-VlnPlot(PA.subset, features = c("CD3D"), pt.size = 0,group.by = "main.annotation",  
        cols =c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)

p7<-VlnPlot(PA.subset, features = c("COL1A2"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
p8<-VlnPlot(PA.subset, features = c("ACTA2"), pt.size = 0,group.by = "main.annotation",  
        cols =c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
p9<-VlnPlot(PA.subset, features = c("FN1"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)

p10<-VlnPlot(PA.subset, features = c("PLVAP"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
p11<-VlnPlot(PA.subset, features = c("VWF"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
p12<-VlnPlot(PA.subset, features = c("CD34"), pt.size = 0,group.by = "main.annotation",  
        cols = c("#F1664E", "#388DFD", "#1FC678","#FF59AF","#90E0E0"), ncol = 1)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+labs(y="EXpression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)
pdf("violinplot_main_annotation.pdf", onefile = T, width = 9, height = 11)
CombinePlots(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12), ncol = 3, legend = "right")
dev.off()

#step2 
library(ggplot2)
p1<-ggplot(PA.subset@meta.data, aes(main.annotation)) + 
  geom_bar(fill=brewer.pal(n = 5, name = "Blues")[5:1], colour="black") +
  coord_flip()+ 
  theme(legend.position = "bottom",panel.border = element_rect(colour="black",fill=NA)) + 
  xlab("main.annotation") +
  ylab("Count") 
ggsave("annotation_count.pdf",p1, width = 6,height = 4)
Number.type <- table(PA.subset@meta.data$main.annotation) %>% as.data.frame()
write.csv(Number.type, "Number.type.csv", row.names = F)
#????nFeature_RNA vs nCount_RNA
p2<- ggplot(PA.subset@meta.data, aes(x = nFeature_RNA, y = nCount_RNA)) + 
  geom_point() + theme_bw() + scale_y_log10()
ggsave("annotation_count vs feature.pdf",p2, width = 8.27,height = 4.11)


Epithelial <- row.names(PA.subset@meta.data)[which(PA.subset@meta.data[["main.annotation"]]=='Epithelial')]
table(PA.subset@meta.data[["main.annotation"]])
length(Epithelial)
PA.Epithelial<- subset(PA.subset, cells=Epithelial)

Immune <- row.names(PA.subset@meta.data)[which(PA.subset@meta.data[["main.annotation"]]=='Immune')]
table(PA.subset@meta.data[["main.annotation"]])
length(Immune) 
PA.Immune<- subset(PA.subset, cells=Immune)

Fibroblast <- row.names(PA.subset@meta.data)[which(PA.subset@meta.data[["main.annotation"]]=='Fibroblast')]
table(PA.subset@meta.data[["main.annotation"]])
length(Fibroblast) 
PA.Fibroblast<- subset(PA.subset, cells=Fibroblast)

Endothelial <- row.names(PA.subset@meta.data)[which(PA.subset@meta.data[["main.annotation"]]=='Endothelial')]
table(PA.subset@meta.data[["main.annotation"]])
length(Endothelial) 
PA.Endothelial<- subset(PA.subset, cells=Endothelial)


dir.create("../Step2 subset subcluster")
save(PA.Epithelial,file = "../Step2 subset subcluster/PA.Epithelial.RData")
save(PA.Immune,file = "../Step2 subset subcluster/PA.Immune.RData")
save(PA.Fibroblast,file = "../Step2 subset subcluster/PA.Fibroblast.RData")
save(PA.Endothelial,file = "../Step2 subset subcluster/PA.Endothelial.RData")
rm(list= ls()) 
gc() 


load("~/scRNA/PA_combination_mito0.2/Double.score0.4/Step0-1/PA.subset2.RData")
Idents(PA.subset)<-"main.annotation"
subtype.markers <- FindAllMarkers(PA.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(subtype.markers, "subtype_markers.csv")

top20 <- subtype.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
matrix <- GetAssayData(PA.subset[["RNA"]], slot = "counts")
new_cluster <- sort(PA.subset@active.ident)
matrix_top20 <- as.matrix(matrix[top20$gene, names(new_cluster)])
Meta <- data.frame(Cluster.ID=PA.subset@meta.data[,"seurat_clusters"],
                   Sample.ID=PA.subset@meta.data[["Sample_ID"]],
                   rownames = rownames(PA.subset@meta.data))
annotation_col <- data.frame(main.annotation=new_cluster)
annotation_col$rownames<-rownames(annotation_col)
annotation_col <- inner_join(annotation_col, Meta, by="rownames")
rownames(annotation_col) <- annotation_col$rownames
annotation_col <- annotation_col[,c("main.annotation","Cluster.ID","Sample.ID")]

library(pheatmap)
matrix_top20 <- NormalizeData(matrix_top20, scale.factor = 10000)
matrix_top20 <- ScaleData(matrix_top20, scale.max = 2.5) 
plot<-pheatmap(matrix_top20,show_colnames =F,show_rownames = T, cellheight = 3,
               cluster_rows = T,cluster_cols = F, fontsize = 3,
               annotation_col=annotation_col)
ggsave("Heatmap.pdf",plot, width = 6,height = 4)

#prepare velocyto filtriting barcode
barcode<-rownames(PA.subset@meta.data)
barcode<-as.data.frame(gsub("_", "-", barcode))
write.csv(barcode, "filtered_barcode.csv")
