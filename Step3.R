library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#my functions
GSVA_barplot<-function(GSVA_result, p, n){
  rownames(GSVA_result)<-gsub("HALLMARK_","",rownames(GSVA_result))
  GSVA_result <- subset(GSVA_result, adj.P.Val < p) %>% arrange(logFC)
  GSVA_result$Hallmark <- rownames(GSVA_result)
  GSVA_result$Hallmark <- factor(GSVA_result$Hallmark, levels= GSVA_result$Hallmark)
  GSVA_result$Diff[GSVA_result$logFC>0] = "Up-regulated"
  GSVA_result$Diff[GSVA_result$logFC<0] = "Down-regulated"
  GSVA_result <- GSVA_result %>% group_by(Diff) %>% top_n(n = n, wt = abs(logFC)) %>% data.frame()
  rownames(GSVA_result) <- GSVA_result$Hallmark
  plot <- ggplot2::ggplot(GSVA_result, aes(Hallmark, logFC)) + geom_bar(stat = 'identity') + 
    coord_flip() + 
    geom_col(aes(fill= adj.P.Val))+
    labs(x="Hallmark Pathway", y="logFC", title="Hallmark Pathways from GSVA") +
    theme_minimal()+
    theme(axis.text.y = element_text(size = 7))
  return(plot)
}

HeatmapforFindAllMarkers<-function(Markers, seuratobj,category, annotation_col, number=10, scalefactor = 3, pch=3, pf=3){
  library(Seurat)
  library(pheatmap)
  top <- Markers %>% group_by(cluster) %>% top_n(n = number, wt = avg_log2FC)
  matrix <- GetAssayData(seuratobj[["RNA"]], slot = "counts")
  Idents(seuratobj)<-category
  new_cluster<-sort(seuratobj@active.ident)
  matrix_top <- as.matrix(matrix[top$gene, names(new_cluster)])
  matrix_top <- NormalizeData(matrix_top, scale.factor = 10000)
  matrix_top <- ScaleData(matrix_top, scale.max = scalefactor) 
  plot<-pheatmap(matrix_top,show_colnames =F,show_rownames = T, cellheight = pch,
                 cluster_rows = T,cluster_cols = F, fontsize = pf,
                 annotation_col=annotation_col)
  return(plot)
}

GSVA_heatmap<-function(gsva_matrix, seuratobj, category, 
                       annotation_col, scalefactor = 4.5, pfr=5){
  library(pheatmap)
  es.heatmap <- as.data.frame(gsva_matrix)
  rownames(es.heatmap) <-gsub("HALLMARK_", "", rownames(es.heatmap))
  Idents(seuratobj)<-category
  new_cluster <- sort(seuratobj@active.ident)
  es.heatmap <- as.matrix(es.heatmap[rownames(es.heatmap), names(new_cluster)])
  es.heatmap <- ScaleData(es.heatmap, scale.max = scalefactor)
  heatmap <- pheatmap(es.heatmap, show_rownames=T,
                      show_colnames=F, 
                      annotation_col=annotation_col,
                      cluster_cols=F, cluster_rows = T,
                      fontsize_row=pfr, width=10, height=5)
  return(heatmap)
}

setwd("/my_path/Step3 reannotation of Epithelial")

load("../Step2 subset subcluster/PA.Epithelial.RData")


DefaultAssay(PA.Epithelial)<-"integrated"
PA.Epithelial <- ScaleData(PA.Epithelial)

#scores calculations for doublet detection
immunlist<-c("PTPRC","IGKC","IGHG1","IGHG3","IGLC1","LYZ","IGHM","JCHAIN","IGHA1","CD79A",
             "CD3D","CD3G","CD3E","CD4","CD8A","CD8B","NKG7")
TAFgenes<-c("ACTA2", "CCN2","COL1A1","COL4A1","COL1A2","COL4A2","COL5A1","COL5A2","COL5A3",
            "COL6A1","COL6A2","COL10A1","COL15A1","COL18A1","TAGLN","MYL9","TPM1",
            "MMP11","HOPX","POSTN","TPM2","PDGFRA","CFD","DPT","LMNA","AGTR1","HAS1")
Endothelialgenes<-c("VWF","PECAM1","ACE","CD34","ENG","ICAM1","NCAM1","SELE","SELP",
                    "VEGFR1","VEGFR2","TEK","PLVAP","CD248","FN1","TFRC","VTN") 
PA.Epithelial<-AddModuleScore(PA.Epithelial, features = list(immunlist), assay = "RNA", name = "Immune_score", search = F)
PA.Epithelial<-AddModuleScore(PA.Epithelial, features = list(TAFgenes), assay = "RNA", name = "TAF_score", search = F)
PA.Epithelial<-AddModuleScore(PA.Epithelial, features = list(Endothelialgenes), assay = "RNA", name = "Endo_score", search = F)

library(ggplot2)
p1<-ggplot(PA.Epithelial@meta.data, aes(Immune_score1))+
  geom_histogram(binwidth = 0.001)+geom_vline(xintercept = 0.10)
p2<-ggplot(PA.Epithelial@meta.data, aes(TAF_score1))+
  geom_histogram(binwidth = 0.001)+geom_vline(xintercept = 0.15)
p3<-ggplot(PA.Epithelial@meta.data, aes(Endo_score1))+
  geom_histogram(binwidth = 0.001)+geom_vline(xintercept = 0.167)
p<-p1|p2|p3
ggsave(filename = "QC_score.pdf",p,width = 8,height = 3)

data<-data.frame(Immune_score=PA.Epithelial@meta.data$Immune_score1,
                 TAF_score=PA.Epithelial@meta.data$TAF_score1,
                 Endo_score=PA.Epithelial@meta.data$Endo_score1,
                 row.names = rownames(PA.Epithelial@meta.data))

data$Doublet1<-data$Immune_score > 0.10
data$Doublet2<-data$TAF_score > 0.15
data$Doublet3<-data$Endo_score > 0.167
highCells1<-rownames(subset(data, Doublet1==T))
highCells2<-rownames(subset(data, Doublet2==T))
highCells3<-rownames(subset(data, Doublet3==T))

highCells<-unique(c(highCells1,highCells2,highCells3))
cell<-setdiff(Cells(PA.Epithelial), highCells)
PA.Epithelial<-SetIdent(PA.Epithelial, cells=highCells, value = "Doublet")
PA.Epithelial<-SetIdent(PA.Epithelial, cells=cell, value = "Single")
DimPlot(PA.Epithelial, reduction = "tsne", label=FALSE)
save(PA.Epithelial,file="PA.Epithelial_predrop.RData")

PA.Epithelial<-subset(PA.Epithelial, cells = cell)
dim(PA.Epithelial)
save(PA.Epithelial,file="PA.Epithelial.RData")


DefaultAssay(PA.Epithelial)<-"RNA"
pdf("Cluster_marker_post-check.pdf", width = 6, height = 4.5)
FeaturePlot(PA.Epithelial, reduction = "tsne", features = "COL1A2", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.Epithelial, reduction = "tsne", features = "PLVAP", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.Epithelial, reduction = "tsne", features = "PTPRC", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.Epithelial, reduction = "tsne", features = "CD4", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.Epithelial, reduction = "tsne", features = "CD8A", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.Epithelial, reduction = "tsne", features = "IGKC", pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()


DefaultAssay(PA.Epithelial)<-"integrated"
PA.Epithelial <- ScaleData(PA.Epithelial)
##PCA
PA.Epithelial <- RunPCA(PA.Epithelial, features = VariableFeatures(object = PA.Epithelial))
ElbowPlot(PA.Epithelial) 
#Neighbor graph
PA.Epithelial <- FindNeighbors(PA.Epithelial, dims = 1:5)

PA.Epithelial <- FindClusters(PA.Epithelial, resolution = 2)
DimPlot(PA.Epithelial, reduction = "pca",label = TRUE)

PA.Epithelial <- RunUMAP(PA.Epithelial, dims = 1:5)
DimPlot(PA.Epithelial, reduction = "umap",label = TRUE)
PA.Epithelial <- RunTSNE(PA.Epithelial, dims = 1:5)
DimPlot(PA.Epithelial, reduction = "tsne", label=TRUE)
save(PA.Epithelial,file="PA.Epithelial.RData")

load("./PA.Epithelial.RData")

library(ggplot2)
pdf("Dimplot.pdf", width = 6, height = 4.5)
DimPlot(PA.Epithelial, reduction = "umap", label = T, pt.size = 0.6)+
        theme(panel.border = element_rect(colour="black",fill=NA))
DimPlot(PA.Epithelial, reduction = "tsne", label = T, pt.size = 0.6)+
        theme(panel.border = element_rect(colour="black",fill=NA))
DimPlot(PA.Epithelial, reduction = "tsne", group.by = "Sample_ID", pt.size = 0.6)+
        theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

#cluster biomarkers
Epithelial.markers <- FindAllMarkers(PA.Epithelial, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Epithelial.markers, "Epithelial.markers.csv")

top3<-Epithelial.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
genes_to_check <- unique(top3$gene)
tab1 <- cbind(as.data.frame(PA.Epithelial@meta.data$Sample_ID),as.data.frame(PA.Epithelial@active.ident))
colnames(tab1) <- c("Sample", "subcluster")

pdf("DEgenes_Dotplot.pdf", width = 10, height = 12)
DotPlot(PA.Epithelial, assay= "RNA", features = genes_to_check) + coord_flip()+
        theme(panel.border = element_rect(colour="black",fill=NA))
ggplot(tab1) +
        aes(x = subcluster, fill = factor(Sample)) +
        geom_bar(position = "fill")+
        theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

top5<-Epithelial.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
pdf("Cluster_Doheatmap.pdf", width = 22, height = 15)
DoHeatmap(PA.Epithelial, features = top5$gene, assay = "integrated")
dev.off()

DefaultAssay(PA.Epithelial)<- "RNA"
#violin
p1<-VlnPlot(PA.Epithelial, features = "POU1F1", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p2<-VlnPlot(PA.Epithelial, features = "SF1", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p3<-VlnPlot(PA.Epithelial, features = "GH1",pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p4<- VlnPlot(PA.Epithelial, features = "PRL", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p5<-VlnPlot(PA.Epithelial, features = "TSHB", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p6<-VlnPlot(PA.Epithelial, features = "VIM", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p7<-VlnPlot(PA.Epithelial, features = "PTTG1", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p8<-VlnPlot(PA.Epithelial, features = "SOX2", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p9<-VlnPlot(PA.Epithelial, features = "SOX9", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
pdf("Epithelial_cluster_identification_violin.pdf", width = 20, height = 10)
CombinePlots(list(p1,p2,p3,p4,p5,p6,p7,p8,p9), ncol = 3, legend = "right")
dev.off()
#Feutureplot
cols <- brewer.pal(n=4, name = "OrRd")
p1<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "POU1F1",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p2<- FeaturePlot(PA.Epithelial, reduction = "umap", features = "SF1",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p3<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "GH1",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p4<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "PRL",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p5<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "TSHB",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p6<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "VIM",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p7<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "PTTG1",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p8<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "SOX2",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p9<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "SOX9",cols = cols)+
        theme(panel.border = element_rect(colour="black",fill=NA))
pdf("Epithelial_cluster_identification_featureplot.pdf", width = 8, height = 7)
CombinePlots(list(p1,p2,p3,p4,p5,p6,p7,p8,p9), ncol = 3, legend = "right")
dev.off()

p1<-VlnPlot(PA.Epithelial, features = "doublet.score", pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p2<-VlnPlot(PA.Epithelial, features = "Malign.score",pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p3<-VlnPlot(PA.Epithelial, features = "CellCycle.score",pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
p4<-VlnPlot(PA.Epithelial, features = "Stemness.score",pt.size = 0)+
        geom_boxplot(width=0.3,alpha=1,outlier.colour = NA)+
        stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
        theme(panel.border = element_rect(colour="black",fill=NA))
pdf("Epithelial_cluster_identification_metadata_violin.pdf", width = 20, height = 10)
CombinePlots(list(p1,p2,p3,p4), ncol = 2, legend = "right")
dev.off()
#Feutureplot
p1<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "doublet.score")+
        theme(panel.border = element_rect(colour="black",fill=NA))
p2<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "Malign.score")+
        theme(panel.border = element_rect(colour="black",fill=NA))
p3<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "CellCycle.score")+
        theme(panel.border = element_rect(colour="black",fill=NA))
p4<-FeaturePlot(PA.Epithelial, reduction = "umap", features = "Stemness.score")+
        theme(panel.border = element_rect(colour="black",fill=NA))
pdf("Epithelial_cluster_identification_metadata_featureplot.pdf", width = 10, height = 9)
CombinePlots(list(p1,p2,p3,p4), ncol = 2, legend = "right")
dev.off()

genescore<-read.csv("./41467_2020_19012_MOESM4_ESM.csv")
E_score<-subset(genescore, clusters=="epi")$genes
M_score<-subset(genescore, clusters=="mesen")$genes
S_score<-setdiff(genescore$genes, c(E_score, M_score))
PA.Epithelial<-AddModuleScore(PA.Epithelial, features = list(E_score), assay = "RNA", name = "E_score", search = F)
PA.Epithelial<-AddModuleScore(PA.Epithelial, features = list(M_score), assay = "RNA", name = "M_score", search = F)
PA.Epithelial<-AddModuleScore(PA.Epithelial, features = list(S_score), assay = "RNA", name = "S_score", search = F)


p1<-RidgePlot(PA.Epithelial, features = "E_score1",cols = color)
p2<-RidgePlot(PA.Epithelial, features = "M_score1",cols = color)
p3<-RidgePlot(PA.Epithelial, features = "S_score1",cols = color)
p<-CombinePlots(list(p1,p2,p3), ncol = 1, legend = "right")
ggsave("./Cluster_score_ridgeplot.pdf", width = 5, height = 10)


Epithelial.markers<-subset(Epithelial.markers, p_val_adj<0.05)
GH<-subset(Epithelial.markers,gene=="GH1")
PRL<-subset(Epithelial.markers,gene=="PRL")
TSH<-subset(Epithelial.markers,gene=="TSHB")

Epithelial.annotation<- c("GH","GH_PRL","Triple_P","GH",
                          "GH_PRL","GH","GH","Triple_N",
                          "GH", "GH","GH","Triple_N",
                          "Triple_N","EMT_like_A","Triple_N","GH", 
                          "Triple_P","GH_PRL","PRL","GH",
                            "GH","EMT_like_A","EMT_like_A","EMT_like_B",
                          "GH","Triple_P","Stem_cell_B","GH",
                          "PRL","EMT_like_B","Stem_cell_C","Stem_cell_A")
PA.Epithelial@meta.data$Epithelial.annotation <- PA.Epithelial@meta.data[["seurat_clusters"]]
cluster.ids <- sort(as.numeric(unique(as.character(PA.Epithelial@meta.data[["seurat_clusters"]]))))
PA.Epithelial@meta.data$Epithelial.annotation <- plyr::mapvalues(x = PA.Epithelial@meta.data$Epithelial.annotation, from = cluster.ids, to = Epithelial.annotation)
table(PA.Epithelial@meta.data$Epithelial.annotation)
PA.Epithelial@meta.data[["Epithelial.annotation"]]<-factor(PA.Epithelial@meta.data[["Epithelial.annotation"]], levels = c("GH","PRL","GH_PRL","Triple_P","Triple_N","EMT_like_A","EMT_like_B","Stem_cell_A","Stem_cell_B","Stem_cell_C"))
save(PA.Epithelial,file="PA.Epithelial.RData")
pdf("Epithelial.annotation.pdf", onefile = T, width = 7, height = 5)
DimPlot(PA.Epithelial, reduction = "umap", label = FALSE, group.by = "Epithelial.annotation", 
        cols =  c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35")
        ,pt.size = 0.6)+
        theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

dat <- data.frame(row.names = rownames(PA.Epithelial@meta.data),
                  Subtype = PA.Epithelial@meta.data$Epithelial.annotation,
                  Sample = PA.Epithelial@meta.data$Sample_ID)
p<-ggplot(dat) +
  aes(x = Sample, fill = Subtype) +
  geom_bar(position = "fill")+
  labs(y="Proportion")+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 8))+
  scale_fill_manual(values=c("#E41A1C","#4376AC","#49A75A","#87638F","#D77F32","#F58879","#D690C6","#B17A7D","#847A74","#4385BF"))
ggsave(filename = "Epithelial_Subtype_contribution.pdf", p, width = 4, height = 4)

#violin
Idents(PA.Epithelial)<-"Epithelial.annotation"
color= c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35")
p1<-VlnPlot(PA.Epithelial, features = "GH1",pt.size = 0, cols = color)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  labs(y="Expression")+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<- VlnPlot(PA.Epithelial, features = "PRL", pt.size = 0, cols = color)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  labs(y="Expression")+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-VlnPlot(PA.Epithelial, features = "TSHB", pt.size = 0, cols = color)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  labs(y="Expression")+
  theme(panel.border = element_rect(colour="black",fill=NA))
p4<-VlnPlot(PA.Epithelial, features = "VIM", pt.size = 0, cols = color)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  labs(y="Expression")+
  theme(panel.border = element_rect(colour="black",fill=NA))
p5<-VlnPlot(PA.Epithelial, features = "SOX2", pt.size = 0, cols = color)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  labs(y="Expression")+
  theme(panel.border = element_rect(colour="black",fill=NA))
p6<-VlnPlot(PA.Epithelial, features = "SOX9", pt.size = 0, cols = color)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  labs(y="Expression")+
  theme(panel.border = element_rect(colour="black",fill=NA))
pdf("Epithelial_markers_violin.pdf", width = 8, height = 5.5)
CombinePlots(list(p1,p2,p3,p5), ncol = 2, legend = "right")
dev.off()

p1<-RidgePlot(PA.Epithelial, features = "M_score1",cols = color)
p2<-VlnPlot(PA.Epithelial, features = "M_score1", pt.size = 0, cols = color)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  labs(y="Score",title ="M_score")
ggsave(filename = "Epithelial_mscore_ridge.pdf", p1,width = 6.5,height = 4)
ggsave(filename = "Epithelial_mscore_vln.pdf", p2,width = 6,height = 3)

Idents(PA.Epithelial)<-"Epithelial.annotation" #????active.idents
Subtype.markers <- FindAllMarkers(PA.Epithelial, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Subtype.markers, "Subtype.markers.csv")

top10 <- Subtype.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
matrix <- GetAssayData(PA.Epithelial[["RNA"]], slot = "counts")
new_cluster <- rownames(arrange(PA.Epithelial@meta.data, Epithelial.annotation, malignType))
matrix_top10 <- as.matrix(matrix[top10$gene, new_cluster])
annotation_col<-data.frame(Epithelial.annotation=PA.Epithelial@meta.data[["Epithelial.annotation"]],
                           Cluster_ID=PA.Epithelial@meta.data[["seurat_clusters"]],
                            malignType=PA.Epithelial@meta.data[["malignType"]],
                           Sample_ID=PA.Epithelial@meta.data[["Sample_ID"]],
                           row.names = rownames(PA.Epithelial@meta.data))
library(pheatmap)
matrix_top10 <- NormalizeData(matrix_top10, scale.factor = 10000)
matrix_top10 <- ScaleData(matrix_top10, scale.max = 6.5) 

plot<-pheatmap(matrix_top10,show_colnames =F,show_rownames = T, cellheight = 3,
         cluster_rows = T,cluster_cols = F, fontsize = 3,
         annotation_col=annotation_col)
ggsave("Heatmap.pdf", plot, width = 5, height = 8)


#GSVA
library(GSVA)

library(msigdbr)

m_df<- msigdbr(species = "Homo sapiens", category = "H")
gset.list<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

DefaultAssay(PA.Epithelial) <- "RNA" 
expr <- as.matrix(GetAssayData(PA.Epithelial, slot = 'data')) 
save(expr,gset.list, file="./my_path/ubutu.es.Epithelial.RData")
es <- gsva(expr, gset.list, mx.diff=TRUE, verbose=TRUE, parallel.sz=12)
save(es, file = "es_Epithelial.RData")


load("~/scRNA/PA_combination_mito0.2/Double.score0.4/Step3 reannotation of Epithelial/es_Epithelial.RData")
load("~/scRNA/PA_combination_mito0.2/Double.score0.4/Step3 reannotation of Epithelial/PA.Epithelial.RData")

es.heatmap<-as.data.frame(es)
rownames(es.heatmap) <-gsub("HALLMARK_", "", rownames(es.heatmap))
PA.Epithelial<-AddMetaData(PA.Epithelial, t(es.heatmap), col.name =rownames(es.heatmap))
save(PA.Epithelial, file = "PA.Epithelial.RData")

Idents(PA.Epithelial)<-"Epithelial.annotation"
meta<- data.frame(row.names =rownames(PA.Epithelial@meta.data),
                  Epithelial.annotation=PA.Epithelial@meta.data[["Epithelial.annotation"]],
                  malignType=PA.Epithelial@meta.data[["malignType"]],
                  Cluster_ID=PA.Epithelial@meta.data[["seurat_clusters"]])
new_cluster <- rownames(arrange(PA.Epithelial@meta.data, Epithelial.annotation, malignType)) #对表达矩阵排�?
es.heatmap <- as.matrix(es.heatmap[rownames(es.heatmap), new_cluster])

es.heatmap <- ScaleData(es.heatmap, scale.max = 5)
heatmap <- pheatmap(es.heatmap, show_rownames=1,
                    show_colnames=0, 
                    annotation_col=meta, fontsize = 3.5,
                    cluster_cols=F, cluster_rows = T,
                    fontsize_row=4.5, width=10, height=5)
ggsave("GSVA_heatmap.pdf",heatmap, width = 10, height = 5.7)


#infercnv
load("/my_path/Step3 reannotation of Epithelial/PA.Epithelial.RData")
load("/my_path/Step4 reannotation of Fibroblast/PA.Fibroblast.RData") #control

dfcount <- as.data.frame(PA.Epithelial@assays$RNA@counts)
TAFdf<-as.data.frame(PA.Fibroblast@assays$RNA@counts)
dfcount<-cbind(dfcount, TAFdf)

groupinfo<-data.frame(CellID=rownames(PA.Epithelial@meta.data),
                      annotation=PA.Epithelial@meta.data$Epithelial.annotation)
TAFinfo<-data.frame(CellID=rownames(PA.Fibroblast@meta.data),
                    annotation=PA.Fibroblast@meta.data$Fibroblast.annotation)
groupinfo<-rbind(groupinfo, TAFinfo)
table(groupinfo$annotation)

library(AnnoProbe)
geneInfor<-annoGene(rownames(dfcount),"SYMBOL",'human')
geneInfor<-geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor<-geneInfor[!duplicated(geneInfor[,1]),]
geneInfor$chr<-factor(geneInfor$chr, levels = c("chr1","chr2","chr3","chr4","chr5",
                                                "chr6","chr7","chr8","chr9","chr10",
                                                "chr11","chr12","chr13","chr14","chr15",
                                                "chr16","chr17","chr18","chr19","chr20",
                                                "chr21","chr22","chrX","chrY"))
geneInfor<-arrange(geneInfor, chr)

dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 

dir.create("./infercnv")
dir.create("./infercnv/input")
expFile='./infercnv/input/Epithelial_counts_matrix.txt'
write.table(dfcount ,file = expFile,sep = '\t',quote = F)
groupFiles='./infercnv/input/Epithelial_groupFiles.txt'
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
geneFile='./infercnv/input/Epithelial_geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)


#run infercnv in screen
rm(list=ls())
gc()
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./infercnv/input/Epithelial_counts_matrix.txt",
                                    annotations_file="./infercnv/input/Epithelial_groupFiles.txt",
                                    delim="\t",
                                    gene_order_file= "./infercnv/input/Epithelial_geneFile.txt",
                                    ref_group_names=c('mTAF_C',"mTAF_B","mTAF_A","iTAF")) 

save(infercnv_obj, file="./infercnv/infercnv_obj.RData")

# Run inferCNV
infercnv_all = infercnv::run(infercnv_obj,
                             cutoff=0.1,  
                             out_dir= "./infercnv/Epithelial", 
                             cluster_by_groups=T,   
                             output_format = "pdf",
                             num_threads=16,
                             denoise=T,
                             HMM=T)
save(infercnv_all, file="./infercnv/infercnv_all.RData")

Epithelial_cnv<-read.table("./infercnv/Epithelial/infercnv.observations.txt", header=T)
refer_cnv <-read.table("./infercnv/Epithelial/infercnv.references.txt", header=T)
malignScore <- colSums((Epithelial_cnv - 1)^2)*10
malignScore <- malignScore / dim(Epithelial_cnv)[1]
referScore <- colSums((refer_cnv - 1)^2)*10
referScore <- referScore / dim(refer_cnv)[1]

library(scCancer)
#quantile(malignScore)

getBimodalThres <- function(scores){
  x.density <- density(scores)
  d.x.density <- diff(x.density$y)
  d.sign <- (d.x.density > 0) + 0
  
  ext.pos <- which(d.sign[2:length(d.sign)] - d.sign[1:(length(d.sign)-1)] != 0)
  ext.density <- x.density$y[ext.pos]
  y.max <- max(ext.density)
  if(length(ext.pos) >= 3){
    del.ix <- c()
    for(ei in 2:length(ext.density)){
      if(abs(ext.density[ei] - ext.density[ei - 1]) < y.max * 0.001){
        del.ix <- c(del.ix, ei - 1, ei)
      }
    }
    sel.ix <- !(1:length(ext.density) %in% unique(del.ix))
    ext.density <- ext.density[sel.ix]
    ext.pos <- ext.pos[sel.ix]
  }
  
  if(length(ext.pos) >= 3){
    t.ext.density <- c(0, ext.density, 0)
    ext.height <- sapply(2:(length(ext.pos) + 1), FUN = function(x){
      return(min(abs(t.ext.density[x] - t.ext.density[x-1]), abs(t.ext.density[x] - t.ext.density[(x+1)])))
    })
    ext <- data.frame(x = ext.pos, y = ext.density, height = ext.height)
    max.ix <- order(ext.density, decreasing = T)
    if(ext.height[max.ix[2]] / ext.height[max.ix[1]] > 0.01){
      cut.df <- ext[c(max.ix[2]:max.ix[1]), ]
      threshold <- x.density$x[cut.df[which.min(cut.df$y), ]$x]
    }else{
      threshold <- NULL
    }
  }else{
    threshold <- NULL
  }
  
  return(threshold)
}

malignPlot <- function(obserScore, referScore, malign.thres = NULL){
  scoreDF <- data.frame(malignScore = c(obserScore, referScore),
                        sets = c(rep("Observation", length(obserScore)),
                                 rep("Reference", length(referScore))))
  p <- ggplot() +
    geom_histogram(data = subset(scoreDF, sets == "Observation"),
                   mapping = aes(x = malignScore, fill = "Observation"),
                   bins = 150, alpha = 0.6) +
    geom_histogram(data = subset(scoreDF, sets == "Reference"),
                   mapping = aes(x = malignScore, fill = "Reference"),
                   bins = 150, alpha = 0.6) +
    labs(x = "Malignancy score", y = "Droplets count") +
    scale_fill_manual(name = "Cells sets", guide = "legend",
                      values = c("Observation"="#2e68b7", "Reference"="grey")) +
    theme_classic() +
    ggplot_config(base.size = 7) +
    theme(legend.justification = c(1.12,1.12), legend.position = c(1,1))
  if(!is.null(malign.thres)){
    p <- p + geom_vline(xintercept = malign.thres, colour = "red", linetype = "dashed")
  }
  return(p)
}



cutoff<-getBimodalThres(malignScore)#cutoff
#cutoff<-0.006

p<-malignPlot(malignScore, referScore, cutoff)
ggsave("./infercnv/malignPlot.pdf",p, width = 8, height = 4)

data<-as.data.frame(malignScore)
data$malignType<-data$malignScore > cutoff
data$malignType <- plyr::mapvalues(x = data$malignType, from = c(T,F), to = c("highCNV","lowCNV"))
table(data$malignType)
PA.Epithelial<-AddMetaData(PA.Epithelial, data)
table(PA.Epithelial@meta.data$malignType)
save(PA.Epithelial, file = "PA.Epithelial.RData")

VlnPlot(PA.Epithelial, features = "malignScore", pt.size = 0, 
        group.by = "Epithelial.annotation",
        split.by = "malignType",split.plot = F)
DimPlot(PA.Epithelial, reduction = "umap", group.by = "malignType")

plot<-VlnPlot(PA.Epithelial, features = "malignScore", pt.size = 0, 
        group.by = "Epithelial.annotation",
        cols = c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  theme(panel.border = element_rect(colour="black",fill=NA))
ggsave("./infercnv/Epithelial_malignScore.pdf",plot, width = 6, height = 3.5)

library(RColorBrewer)
brewer.pal.info
display.brewer.pal(8,"Accent")
cols<-brewer.pal(8,"Accent")[c(1,3,4,6)]
p<-VlnPlot(PA.Epithelial, features = "malignScore", pt.size = 0, 
        group.by = "Sample_ID",cols = cols)+labs(y="Score")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(angle=0,hjust = 0.5))
ggsave("./infercnv/Sample_malignScore.pdf",p, width = 4, height = 2.5)


color<-c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E")
p<- VlnPlot(PA.Epithelial, features = "malignScore", pt.size = 0, 
        group.by = "Epithelial.annotation",cols = color)+labs(y="Score")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(angle=45,hjust = 1, size = 8))
ggsave("./infercnv/Epithelial_malignScore.pdf",p, width = 6, height = 2.8)

tab <- cbind(as.data.frame(PA.Epithelial@meta.data$malignType),
             as.data.frame(PA.Epithelial@meta.data$Epithelial.annotation),
             as.data.frame(PA.Epithelial@meta.data$Sample_ID))
colnames(tab) <- c("malignType", "Epithelial.annotation","Sample_ID")

p1<-ggplot(tab) +
  aes(x = Epithelial.annotation, fill = malignType) +
  geom_bar(position = "fill")+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1))+
  scale_fill_manual(values=c("#807DBA","#6BAED6"))
p2<-ggplot(tab) +
  aes(x = Sample_ID, fill = malignType) +
  geom_bar(position = "fill")+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1))+
  scale_fill_manual(values=c("#807DBA","#6BAED6"))
ggsave("./infercnv/Epithelial_distribution_barplot.pdf",p1, width = 4, height = 4.5)
ggsave("./infercnv/Sample_distribution_barplot.pdf",p2, width = 5, height = 5)

Idents(PA.Epithelial)<-"malignType"
malignType.markers <- FindMarkers(PA.Epithelial, assay = "RNA", ident.1 = "highCNV", ident.2 = "lowCNV", logfc.threshold = 0.25)
write.csv(malignType.markers, "./infercnv/malignType.markers.csv")

malignType.markers<-subset(malignType.markers, abs(avg_log2FC) > 1 & p_val_adj < 0.05)

matrix <- GetAssayData(PA.Epithelial[["RNA"]], slot = "counts")
new_cluster <- sort(PA.Epithelial@active.ident) 
matrix <- as.matrix(matrix[rownames(malignType.markers), names(new_cluster)])
annotation_col<-data.frame(malignType=PA.Epithelial@meta.data[["malignType"]],
                           Epithelial.annotation=PA.Epithelial@meta.data[["Epithelial.annotation"]],
                           row.names = rownames(PA.Epithelial@meta.data))

library(pheatmap)
matrix <- NormalizeData(matrix, scale.factor = 10000)
matrix <- ScaleData(matrix, scale.max = 6.5) 
plot<-pheatmap::pheatmap(matrix,show_colnames =F,show_rownames = T,
               cluster_rows = T,cluster_cols = F, fontsize = 3,
               annotation_col=annotation_col)
ggsave("./infercnv/Heatmap_malignType.pdf", plot, width = 6,height = 4) #热图还是不是很好�?

library(ComplexHeatmap)
meta<-PA.Epithelial@meta.data[names(new_cluster),] %>% select(c("malignType"))
group<-HeatmapAnnotation(df = meta,which = "column",simple_anno_size = unit(3, "mm"),
                         col = list(malignType=c("highCNV"="#807DBA","lowCNV"="#6BAED6")),
                         annotation_name_gp=gpar(fontsize=6),annotation_label="malignType")
#topmarkergenes
malignType.markers$Significance <- ifelse(malignType.markers$avg_log2FC>0, "Up", "Down")
malignType.markers$gene<-rownames(malignType.markers)
genelist<-malignType.markers%>% group_by(Significance)%>% top_n(n=5, wt=avg_log2FC)
genelist<-genelist$gene
genelist
markgene<-c("COL4A4", "GNAS","HLA−A", "SERF2","METRN", "AMIGO2", "THBS2", "HES6","CD79B", 
            "GH1") %>% union(genelist) %>% unique()
index<-which(rownames(matrix) %in% markgene)
labels <- rownames(matrix)[index]

lab = rowAnnotation(ano = anno_mark(at = index,
                                    labels = labels,
                                    labels_gp = gpar(fontsize = 4.5),
                                    lines_gp = gpar()))

library(RColorBrewer)
brewer.pal.info
display.brewer.pal(11,"RdBu")
cols<-c(brewer.pal(11,"RdBu")[9],"white",brewer.pal(11,"RdBu")[2])

pdf(file = "./infercnv/Heatmap_malignType_chm.pdf",width = 6,height = 4)
Heatmap(matrix,
        col = colorRampPalette(colors = cols)(100),
        column_split = meta$malignType, column_gap = unit(0.7,"mm"),
        cluster_rows = T,cluster_columns = F, show_column_names = F,show_row_dend = F,
        show_row_names = FALSE,right_annotation = lab,
        top_annotation = group, use_raster = F,
        name = "Exp")
dev.off()

#(FC>1, DAVID)
dat<-read.csv("./infercnv/DAVID/all_up.csv", header = T) %>% subset(FDR<0.05) %>% arrange(Count)
rownames(dat)<-dat$Term

#GO
go<-subset(dat, Category %in% c("GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT"))
top<-go%>% group_by(Category) %>% top_n(n=3, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,12)
top$Term<- factor(top$Term, levels = top$Term)

library(RColorBrewer)
p1<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= PValue))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  #scale_fill_continuous(low=brewer.pal(9,"BuPu")[7], high=brewer.pal(9,"BuPu")[2], limits=c(0,1.360e-06))
  labs(x="GO Term", y="Count", title="GO Term Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 7))

#kegg
kegg<-subset(dat, Category =="KEGG_PATHWAY") 
top<-kegg%>% group_by(Category) %>% top_n(n=9, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,10)
top$Term<- factor(top$Term, levels = top$Term)

p2<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= PValue))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])
  labs(x="KEGG Pathway", y="Count", title="KEGG Pathway Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 7))
ggsave("./infercnv/Epithelial_malignType_up_enrichment.pdf", p1, width = 4.5, height = 2)


dat<-read.csv("./infercnv/DAVID/all_down.csv", header = T) %>% subset(FDR<0.05) %>% arrange(desc(Count))
rownames(dat)<-dat$Term

#GO
go<-subset(dat, Category %in% c("GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT"))
top<-go%>% group_by(Category) %>% top_n(n=3, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,12)
top$Term<- factor(top$Term, levels = top$Term)

p1<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= PValue))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+#颜色变化不明�?
  #scale_fill_continuous(low=brewer.pal(9,"BuPu")[7], high=brewer.pal(9,"BuPu")[2], limits=c(0,1.360e-06))
  labs(x="GO Term", y="Count", title="GO Term Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 7))

#kegg
kegg<-subset(dat, Category =="KEGG_PATHWAY")

ggsave("./infercnv/Epithelial_malignType_down_enrichment.pdf", p1, width = 6, height = 2.5)

#GSVA
load("~/scRNA/PA_combination_mito0.2/Double.score0.4/Step3 reannotation of Epithelial/es_Epithelial.RData")
library(limma)

meta <- PA.Epithelial@meta.data[,c("malignType")]
group <- factor(meta,levels = c("lowCNV","highCNV"),ordered = F)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(es,design)
contrast.matrix <- makeContrasts(highCNV-lowCNV,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff.pathway <- topTable(fit2,adjust="fdr", number=Inf)
rownames(diff.pathway)<-gsub("HALLMARK_", "", rownames(diff.pathway))
write.csv(diff.pathway,"./infercnv/Epithelial_highCNV_lowCNV_GSVA.csv")

plot <- GSVA_barplot(diff.pathway, p = 0.05, n=5)
ggsave("./infercnv/Epithelial_highCNV_lowCNV_GSVA.pdf",plot, width = 7.5, height = 3.5)


View(cc.genes)
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(PA.Epithelial)) 
s_genes <- cc.genes$s.genes   
s_genes <- CaseMatch(search=s_genes, match=rownames(PA.Epithelial)) 

PA.Epithelial <- CellCycleScoring(PA.Epithelial, g2m.features=g2m_genes, s.features=s_genes)
save(PA.Epithelial, file = "PA.Epithelial.RData")

Idents(PA.Epithelial)<-"Epithelial.annotation"
cols<-c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35")
DimPlot(PA.Epithelial, reduction = "umap", label = F, group.by = "Phase",
        cols = c("#B452CD","#826DAE","#37B8BD"),pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))

p1<-VlnPlot(PA.Epithelial, features = c("G2M.Score"), pt.size = 0,slot = "data", group.by = "Epithelial.annotation",cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Scores", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size=10))
p2<-VlnPlot(PA.Epithelial, features = c("S.Score"), pt.size = 0,slot = "data", group.by = "Epithelial.annotation",cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Scores", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA), 
        axis.text.x = element_text(size=10))
p<-CombinePlots(list(p1,p2),ncol = 1,legend = "right")
ggsave("./Epithelial_CellCycleScores.pdf", p, width = 4, height = 4.8)

ggplot(PA.Epithelial@meta.data, aes(S.Score))+
  geom_histogram(binwidth = 0.001)+geom_vline(xintercept = 0.15)

ggplot(PA.Epithelial@meta.data, aes(G2M.Score))+
  geom_histogram(binwidth = 0.001)+geom_vline(xintercept = 0.15)


phenodata<-data.frame(S.Score=PA.Epithelial@meta.data[["S.Score"]],
                      G2M.Score=PA.Epithelial@meta.data[["G2M.Score"]],
                      Phase=PA.Epithelial@meta.data$Phase,
                      Cell.Type=PA.Epithelial@meta.data$Epithelial.annotation,
                      Sample_ID=PA.Epithelial@meta.data$Sample_ID,
                      row.names = rownames(PA.Epithelial@meta.data))
phenodata$Cycling.Type[phenodata$S.Score>0.15 | phenodata$G2M.Score>0.15]="Cycling_Cell"
phenodata$Cycling.Type[phenodata$S.Score<0.15 & phenodata$G2M.Score<0.15]="Resting_Cell"

p1<-ggplot(phenodata)+
  geom_point(aes(x=S.Score,y=G2M.Score, colour=Cycling.Type), size=0.6)+
  scale_color_manual(values =c("#FF59AF", "#388DFD"))+
  theme_bw()

p2<-ggplot(phenodata) +
  geom_bar(aes(x = Phase, fill = Cell.Type), position = "fill")+
  labs(y="Proportion")+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text = element_text(size = 8))+
  scale_fill_manual(values=cols)

p3<-ggplot(phenodata) +
  aes(x = Cycling.Type, fill = Cell.Type) +
  geom_bar(position = "fill")+
  labs(y="Proportion")+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))+
  scale_fill_manual(values=cols)
p4<-ggplot(phenodata) +
  aes(x = Sample_ID, fill = Cycling.Type) +
  geom_bar(position = "fill")+
  labs(y="Proportion")+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))+
  scale_fill_manual(values=c("#FF59AF", "#388DFD"))
p5<-p2|p3|p4
ggsave("Cycling.point.pdf", p1, width = 3.5, height = 2)
ggsave("Cycling.bar.pdf", p5, width = 8.5, height = 3)

meta<-select(phenodata, "Cycling.Type")
PA.Epithelial<-AddMetaData(PA.Epithelial, metadata = meta)
save(PA.Epithelial, file = "PA.Epithelial.RData")

Idents(PA.Epithelial)<-"Cycling.Type"
Cycling.markers <- FindAllMarkers(PA.Epithelial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


Idents(PA.Epithelial)<-"Epithelial.annotation"
DefaultAssay(PA.Epithelial)<-"RNA"
malignType_Epithelial.markers<-list()
subset.item<-c("GH","PRL","GH_PRL","Triple_P","Triple_N","EMT_like_A","EMT_like_B")
for (i in subset.item) {
  malignType_Epithelial.markers[[i]]<-FindMarkers(PA.Epithelial,assay = "RNA", ident.1 = "highCNV", group.by = "malignType", subset.ident=i)
}

write.csv(malignType_Epithelial.markers[["GH"]], "./infercnv/malign_GH_markers.csv")
write.csv(malignType_Epithelial.markers[["PRL"]], "./infercnv/malign_PRL_markers.csv")
write.csv(malignType_Epithelial.markers[["GH_PRL"]], "./infercnv/malign_GH_PRL_markers.csv")
write.csv(malignType_Epithelial.markers[["Triple_P"]], "./infercnv/malign_Triple_P_markers.csv")
write.csv(malignType_Epithelial.markers[["Triple_N"]], "./infercnv/malign_Triple_N_markers.csv")
write.csv(malignType_Epithelial.markers[["EMT_like_A"]], "./infercnv/malign_EMT_like_A_markers.csv")
write.csv(malignType_Epithelial.markers[["EMT_like_B"]], "./infercnv/malign_EMT_like_B_markers.csv")


#top10
for (i in subset.item) {
  malignType_Epithelial.markers[[i]]$Significance <-as.factor(ifelse(malignType_Epithelial.markers[[i]]$avg_log2FC>0 ,'Up','Down'))
  malignType_Epithelial.markers[[i]]$gene <- rownames( malignType_Epithelial.markers[[i]])
}
save(malignType_Epithelial.markers, file = "./infercnv/malignType_Epithelial.RData")

genelist<-NULL
for (i in subset.item) {
  genelist[[i]]<-subset(malignType_Epithelial.markers[[i]], abs(avg_log2FC) >0.5 & p_val_adj <0.05)%>% group_by(Significance)%>% top_n(n=10,wt=abs(avg_log2FC))
}

genelist<-unique(c(genelist[[1]]$gene,genelist[[2]]$gene,genelist[[3]]$gene,genelist[[4]]$gene,
                   genelist[[5]]$gene,genelist[[6]]$gene,genelist[[7]]$gene))


metadata <- arrange(PA.Epithelial@meta.data, Epithelial.annotation, malignType)
cell <- row.names(metadata)[which(metadata$Epithelial.annotation %in% subset.item)]

matrix<-data.frame(PA.Epithelial@assays$RNA@counts[genelist, cell])

annotation_col<-data.frame(malignType=PA.Epithelial@meta.data[["malignType"]],
                          Epithelial.annotation=PA.Epithelial@meta.data[["Epithelial.annotation"]],
                          row.names = rownames(PA.Epithelial@meta.data))[cell,]

library(pheatmap)
matrix.heatmap <- NormalizeData(matrix, scale.factor = 10000)
matrix.heatmap <- ScaleData(matrix.heatmap, scale.max = 6.5) 
plot<-pheatmap(matrix.heatmap,show_colnames =F,show_rownames = T,
               cluster_rows = T,cluster_cols = F, fontsize = 3,
               annotation_col=annotation_col)
ggsave("./infercnv/Heatmap_malignType_Epithelial.pdf", plot, width = 4.5,height = 3) 

annotation_col$Type<-paste0(annotation_col$Epithelial.annotation, "-",annotation_col$malignType)
matrix.mean<-t(data.frame(PA.Epithelial@assays$RNA@counts[genelist, cell])) %>% cbind(select(annotation_col,Type))

matrix.mean <-aggregate(x = select(matrix.mean, -"Type"), by = list(Type=matrix.mean$Type), FUN = mean)
matrix.mean$Type<-factor(matrix.mean$Type,levels = c("GH-highCNV","GH-lowCNV","PRL-highCNV","PRL-lowCNV","GH_PRL-highCNV","GH_PRL-lowCNV",
                                                     "Triple_P-highCNV","Triple_P-lowCNV","Triple_N-highCNV","Triple_N-lowCNV",
                                                     "EMT_like_A-highCNV","EMT_like_A-lowCNV","EMT_like_B-highCNV","EMT_like_B-lowCNV"))
matrix.mean<-arrange(matrix.mean,Type)
rownames(matrix.mean)<-matrix.mean$Type

annotation_col<-select(matrix.mean,"Type")
annotation_col$Type<-as.character(annotation_col$Type)

list_col<-strsplit(annotation_col$Type, "-")
for (i in 1:length(list_col)) {
  annotation_col[i,"malignType"]=list_col[[i]][2]
  annotation_col[i,"Epithelial.annotation"]=list_col[[i]][1]
}
annotation_col<-select(annotation_col, c("malignType","Epithelial.annotation"))
annotation_col$Epithelial.annotation<-factor(annotation_col$Epithelial.annotation, levels = c("GH","PRL","GH_PRL", "Triple_P","Triple_N","EMT_like_A","EMT_like_B"))

matrix.mean<-matrix.mean%>% select(-"Type")%>%t()

plot<-pheatmap::pheatmap(matrix.mean,show_colnames =F,show_rownames = T,scale = "row",
               cluster_rows = T,cluster_cols = F,fontsize = 3.4,
               cellwidth = 18, cellheight = 3,border_color = NA,
               annotation_col = annotation_col)
ggsave("./infercnv/Heatmap_malignType_Epithelial_mean.pdf", plot, width = 5.5,height = 3.5) 


genelist1<-rownames(subset(malignType_Epithelial.markers[["GH"]], abs(avg_log2FC) >1 &p_val_adj <0.05)) 
genelist2<-rownames(subset(malignType_Epithelial.markers[["PRL"]], abs(avg_log2FC) >1 & p_val_adj <0.05))
genelist3<-rownames(subset(malignType_Epithelial.markers[["GH_PRL"]], abs(avg_log2FC) >1 & p_val_adj <0.05))
genelist4<-rownames(subset(malignType_Epithelial.markers[["Triple_P"]], abs(avg_log2FC) >1 & p_val_adj <0.05)) 
genelist5<-rownames(subset(malignType_Epithelial.markers[["Triple_N"]], abs(avg_log2FC) >1 & p_val_adj <0.05))
genelist6<-rownames(subset(malignType_Epithelial.markers[["EMT_like_A"]], abs(avg_log2FC) >1 & p_val_adj <0.05)) 
genelist7<-rownames(subset(malignType_Epithelial.markers[["EMT_like_B"]], abs(avg_log2FC) >1 & p_val_adj <0.05)) 

genelist<-list(GH=genelist1,
               PRL=genelist2,
               GH_PRL=genelist3,
               Triple_P=genelist4,
               Triple_N=genelist5,
               EMT_like_A=genelist6,
               EMT_like_B=genelist7)


library(UpSetR)
pdf("./infercnv/Upset_Epithelial.pdf", width = 5.5, height = 3)
upset(fromList(genelist), sets= c("GH","PRL","GH_PRL","Triple_P","Triple_N","EMT_like_A","EMT_like_B"),
      keep.order = T,mb.ratio = c(0.55, 0.45), order.by = "freq",number.angles = 0, line.size = 1, 
      mainbar.y.label = "Intersection Size", sets.x.label = "Gene number",
      queries = list(list(query = intersects, params = list("GH", "PRL","GH_PRL", "Triple_P"), color = "orange", active = T), 
                     list(query = intersects, params = list("GH", "PRL","GH_PRL", "Triple_P","Triple_N"), active = T)),
      query.legend = "top")
dev.off()


MHC_1<-c("HLA.A","HLA.B","HLA.C")
counts<-GetAssayData(PA.Epithelial, assay = "RNA",slot = "counts")
data<-GetAssayData(PA.Epithelial, assay = "RNA",slot = "data")
rownames(counts)<-gsub("-",".",rownames(counts))
rownames(data)<-gsub("-",".",rownames(data))
PA.Epithelial[["RNA"]] <- CreateAssayObject(counts = counts)
PA.Epithelial[["RNA"]] <- CreateAssayObject(data = data)
save(PA.Epithelial, file = "./infercnv/PA.Epithelial.genesymbol.adj.RData")


Idents(PA.Epithelial)<-"Epithelial.annotation"
cols<-c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35")
p1<-VlnPlot(PA.Epithelial, features = "HLA.A", pt.size = 0,slot = "data", group.by = "Epithelial.annotation", cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Expression", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size=10))
p2<-VlnPlot(PA.Epithelial, features = "HLA.B", pt.size = 0,slot = "data", group.by = "Epithelial.annotation", cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Expression", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size=10))
p3<-VlnPlot(PA.Epithelial, features = "HLA.C", pt.size = 0,slot = "data", group.by = "Epithelial.annotation", cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Expression", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size=10))
p<-CombinePlots(list(p1,p2,p3), ncol=1, legend = "right")
ggsave("./infercnv/Sample_MHC_vln.pdf", p, width = 6,height = 8) 


DefaultAssay(PA.Epithelial)<-"RNA"
library(RColorBrewer)
brewer.pal.info
display.brewer.pal(8,"Spectral")
cols<-c(brewer.pal(9,"Purples")[6],brewer.pal(9,"Blues")[5])
p1<-VlnPlot(PA.Epithelial, features = "HLA.A", pt.size = 0, split.by = "malignType", 
           idents =subset.item, 
        group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p2<-VlnPlot(PA.Epithelial, features = "HLA.B", pt.size = 0, split.by = "malignType", 
           idents =subset.item, 
           group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p3<-VlnPlot(PA.Epithelial, features = "HLA.C", pt.size = 0, split.by = "malignType", 
           idents =subset.item, 
           group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p<-CombinePlots(list(p1,p2,p3), ncol=1, legend = "right")
ggsave("./infercnv/Triple_P_MHC_vln.pdf", p, width = 6,height = 8) 

FeaturePlot(PA.Epithelial, features = MHC_1, reduction = "umap")+
  theme(panel.border = element_rect(colour="black",fill=NA))

#vlnplot
MHC_1<-c("HLA.A","HLA.B","HLA.C")
p1<-RidgePlot(PA.Epithelial, features = MHC_1, group.by = "malignType", idents = "GH",cols = cols)
p2<-RidgePlot(PA.Epithelial, features = MHC_1, group.by = "malignType", idents = "PRL",cols = cols)
p3<-RidgePlot(PA.Epithelial, features = MHC_1, group.by = "malignType", idents = "GH_PRL",cols = cols)
p4<-RidgePlot(PA.Epithelial, features = MHC_1, group.by = "malignType", idents = "Triple_P",cols = cols)
p5<-RidgePlot(PA.Epithelial, features = MHC_1, group.by = "malignType", idents = "Triple_N",cols = cols)
p6<-RidgePlot(PA.Epithelial, features = MHC_1, group.by = "malignType", idents = "EMT_like_A",cols = cols)
p7<-RidgePlot(PA.Epithelial, features = MHC_1, group.by = "malignType", idents = "EMT_like_B",cols = cols)
p<-CombinePlots(list(p1,p2,p3,p4,p5,p6,p7), ncol=2, legend = "right")
ggsave("./infercnv/Triple_P4_MHC.pdf", p, width = 15,height = 9) 
rm(PA.Epithelial)
load("./PA.Epithelial.RData")

#FOLDCHANGE>1
genelist<-unique(c(genelist1,genelist2,genelist3,genelist4,
                   genelist5,genelist6,genelist7))
neuro_oncology<-c("GADD45G", "CDKN2A", "CCND1","AMIGO2","SERF2") 
gene_check<-neuro_oncology[which(neuro_oncology %in% genelist)]

library(RColorBrewer)
brewer.pal.info
display.brewer.pal(8,"Spectral")
cols<-c(brewer.pal(9,"Purples")[6],brewer.pal(9,"Blues")[5])

p1<-VlnPlot(PA.Epithelial, features = "AMIGO2", pt.size = 0, split.by = "malignType",
        idents = subset.item,
        group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p2<-VlnPlot(PA.Epithelial, features = "SERF2", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p3<-VlnPlot(PA.Epithelial, features = "CCND1", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p4<-VlnPlot(PA.Epithelial, features = "CDKN2A", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
#oncogenes
oncogenes<-c("GNAS", "USP9X", "AIP","USP8","PTTG1","MEG3","DLK1")
gene_check<-oncogenes[which(oncogenes %in% genelist)]
p5<-VlnPlot(PA.Epithelial, features = "GNAS", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p6<-VlnPlot(PA.Epithelial, features = "MEG3", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p7<-VlnPlot(PA.Epithelial, features = "DLK1", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p<-CombinePlots(list(p1,p3,p2,p4,p5,p6), ncol=2, legend = "right")
ggsave("./infercnv/Epithelial_malignType_oncogenes.pdf", p, width = 6.5,height = 8) 


genes_5g<-genelist1%>% intersect(genelist2) %>% intersect(genelist3) %>% intersect(genelist4) %>% intersect(genelist5)
p8<-VlnPlot(PA.Epithelial, features = genes_5g, pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")

genes_4g<-genelist1%>% intersect(genelist2) %>% intersect(genelist3)%>% intersect(genelist4)
write.csv(data.frame(setdiff(genes_4g,genes_5g)),file = "./infercnv/genes_4g.csv")

p<-VlnPlot(PA.Epithelial, features = genes_4g, pt.size = 0, split.by = "malignType",idents = subset.item,
        group.by = "Epithelial.annotation",cols = cols, ncol = 3)+labs(y="Expression")
ggsave("./infercnv/genes_4g.pdf", p, width = 9,height = 10.5) 

p9<-VlnPlot(PA.Epithelial, features = "TIMP2", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p10<-VlnPlot(PA.Epithelial, features = "HTRA1", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p11<-VlnPlot(PA.Epithelial, features = "THBS2", pt.size = 0, split.by = "malignType",
            idents = subset.item,
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+
  theme(panel.border = element_rect(colour="black",fill=NA))+labs(y="Expression")
p<-CombinePlots(list(p8,p10,p9,p11), ncol=2, legend = "right")
ggsave("./infercnv/Epithelial_malignType_4g_GOrelated.pdf", p, width = 6.5,height = 5.5) 


p<-CombinePlots(list(p1,p2,p5,p6,p8,p10,p9,p11), ncol=2, legend = "right")
ggsave("./infercnv/Epithelial_malignType_geneselect_combine.pdf", p, width = 6.5,height = 10) 

#GSVA
library(GSVA)
library(msigdbr)
m_df<- msigdbr(species = "Homo sapiens", category = "H")
gset.list<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

DefaultAssay(PA.Epithelial) <- "RNA" 
expr <- as.matrix(GetAssayData(subset(PA.Epithelial, Epithelial.annotation %in% subset.item), 
                               slot = 'data')) 
es.DC <- gsva(expr, gset.list, mx.diff=TRUE, verbose=TRUE, parallel.sz=16)
save(es.DC, file = "es.DC.RData")

#GSVA
library(limma)
load("~/my_path/Step3 reannotation of Epithelial/es.DC.RData")
malignType_Epithelial.GSVA <-list()
subset.item<-c("GH","PRL","GH_PRL","Triple_P","Triple_N","EMT_like_A","EMT_like_B")
for (i in subset.item) {
  seuratobj<- subset(PA.Epithelial, Epithelial.annotation==i)
  meta <- seuratobj@meta.data[,c("malignType")]
  group <- factor(meta,levels = c("lowCNV","highCNV"),ordered = F)## levels???�??Ѷ?????????ǰ??
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  es.sub<-es.DC[,Cells(seuratobj)]
  rownames(es.sub) <-gsub("HALLMARK_", "", rownames(es.sub))
  fit <- lmFit(es.sub,design)
  contrast.matrix <- makeContrasts(highCNV-lowCNV,levels=design)#???????Ǻ??ߣ?
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  malignType_Epithelial.GSVA[[i]]<-topTable(fit2,adjust="fdr", number=Inf)
} 
save(malignType_Epithelial.GSVA, file = "./infercnv/malignType_Epithelial.GSVA.RData")

pathway1<-rownames(subset(malignType_Epithelial.GSVA[["GH"]], adj.P.Val <0.05)) 
pathway2<-rownames(subset(malignType_Epithelial.GSVA[["PRL"]],  adj.P.Val <0.05))
pathway3<-rownames(subset(malignType_Epithelial.GSVA[["GH_PRL"]],  adj.P.Val <0.05))
pathway4<-rownames(subset(malignType_Epithelial.GSVA[["Triple_P"]], adj.P.Val <0.05)) 
pathway5<-rownames(subset(malignType_Epithelial.GSVA[["Triple_N"]], adj.P.Val <0.05)) 
pathway6<-rownames(subset(malignType_Epithelial.GSVA[["EMT_like_A"]], adj.P.Val <0.05))
pathway7<-rownames(subset(malignType_Epithelial.GSVA[["EMT_like_B"]], adj.P.Val <0.05))

pathwaylist<-list(GH=pathway1,
                  PRL=pathway2,
                  GH_PRL=pathway3,
                  Triple_P=pathway4,
                  Triple_N=pathway5,
                  EMT_like_A=pathway6,
                  EMT_like_B=pathway7)

library(UpSetR)
pdf("./infercnv/Upset_Epithelial_GSVA_H50.pdf", width = 7, height = 3.5)
upset(fromList(pathwaylist), sets= c("GH","PRL","GH_PRL","Triple_P","Triple_N","EMT_like_A","EMT_like_B"),
      keep.order = T,mb.ratio = c(0.55, 0.45), order.by = "freq",number.angles = 0, line.size = 1, 
      mainbar.y.label = "Intersection Size", sets.x.label = "Pathway number",
      queries = list(list(query = intersects, params = list("GH","PRL","GH_PRL","Triple_P","Triple_N","EMT_like_A","EMT_like_B"), color = "orange", active = T), 
                     
                     list(query = intersects, params = list("GH", "PRL","GH_PRL", "Triple_P","Triple_N"), active = T)),
      query.legend = "top")
dev.off()

Reduce(intersect,list(pathway1,pathway2,pathway3,pathway4,pathway5)) %>% setdiff(pathway6) %>% setdiff(pathway7)
#"REACTIVE_OXYGEN_SPECIES_PATHWAY"

PA.sub<-subset(PA.Epithelial,Epithelial.annotation %in% subset.item)
es.data<-data.frame(t(es.DC))
PA.sub<-AddMetaData(PA.sub, es.data,col.name = colnames(es.data))
save(PA.sub,file = "./infercnv/PA.sub.RData")

p1<-VlnPlot(PA.sub, features = "HALLMARK_ANGIOGENESIS", pt.size = 0, split.by = "malignType",
        group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+labs(y="GSVA Score")+
        theme(panel.border = element_rect(colour="black",fill=NA))
p2<-VlnPlot(PA.sub, features = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", pt.size = 0, split.by = "malignType",
              group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+labs(y="GSVA Score")+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-VlnPlot(PA.sub, features = "HALLMARK_INFLAMMATORY_RESPONSE", pt.size = 0, split.by = "malignType",
            group.by = "Epithelial.annotation",  cols = cols, ncol = 1)+labs(y="GSVA Score")+
  theme(panel.border = element_rect(colour="black",fill=NA))
p<-CombinePlots(list(p1,p2,p3),ncol = 1,legend = "right")
ggsave("./infercnv/GSVA_select.pdf", p, width = 5,height = 8)


dir.create("./SCENIC/")
write.csv(t(as.matrix(PA.Epithelial@assays$RNA@counts)), file = "./SCENIC/Epithelial_pyscenic.csv")
#pyscenic

PA.Epithelial@meta.data$BinaryType <- PA.Epithelial@meta.data[["Epithelial.annotation"]]
PA.Epithelial@meta.data$BinaryType<-ifelse(PA.Epithelial@meta.data$BinaryType %in% c("GH","PRL","GH_PRL","Triple_P","Triple_N"), "DC",
                                           ifelse(PA.Epithelial@meta.data$BinaryType %in% c("EMT_like_A","EMT_like_B"),"EMT", "PC"))
table(PA.Epithelial@meta.data$BinaryType)
PA.Epithelial@meta.data[["BinaryType"]]<-factor(PA.Epithelial@meta.data[["BinaryType"]], levels = c("DC","EMT","PC"))
save(PA.Epithelial,file="PA.Epithelial.RData")

dir.create("./Stem")
Idents(PA.Epithelial)<-"BinaryType"
BinaryType.markers <- FindAllMarkers(PA.Epithelial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(BinaryType.markers, "./Stem/BinaryType.markers.csv")
BinaryType.markers<- subset(BinaryType.markers, avg_log2FC >1 & p_val_adj<0.05)

#top25
top25 <- BinaryType.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
matrix <- GetAssayData(PA.Epithelial[["RNA"]], slot = "counts")
new_cluster <- arrange(PA.Epithelial@meta.data, BinaryType, Epithelial.annotation) %>% select(BinaryType, Epithelial.annotation)


matrix_top <- as.matrix(matrix[top25$gene, rownames(new_cluster)])
annotation_col<-data.frame(CellType=PA.Epithelial@meta.data[["BinaryType"]],
                           Epithelial.annotation=PA.Epithelial@meta.data[["Epithelial.annotation"]],
                           row.names = rownames(PA.Epithelial@meta.data))

library(pheatmap)
matrix_top <- NormalizeData(matrix_top, scale.factor = 10000)
matrix_top <- ScaleData(matrix_top, scale.max = 2.5) 
plot<-pheatmap(matrix_top,show_colnames =F,show_rownames = T, cellheight = 3.5,
               cluster_rows = T,cluster_cols = F, fontsize = 3.5,
               annotation_col=annotation_col)
ggsave("./Stem/Epithelial_BinaryType_heatmap.pdf", plot, width = 5,height = 5)


#top4
top4<- top_n(subset(BinaryType.markers, cluster=="PC"), n=4, wt = avg_log2FC)$gene
VlnPlot(PA.Epithelial, features = top4, pt.size = 0,
        group.by = "BinaryType", ncol = 2)
plot<-FeaturePlot(PA.Epithelial, features = top4, reduction = "umap")+
  theme(panel.border = element_rect(colour="black",fill=NA))
ggsave("./Stem/Epithelial_BinaryType_top4.pdf", plot, width = 6.7,height = 5.5)


#GSVA差异分析
library(limma)
load("~/my_path/Step3 reannotation of Epithelial/es_Epithelial.RData")
rownames(es) <-gsub("HALLMARK_", "", rownames(es))

meta <- PA.Epithelial@meta.data[,c("BinaryType")]
group <- factor(meta,levels = c("DC","EMT","PC"),ordered = F)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(es,design)
contrast.matrix <- makeContrasts(PC-DC,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
Epithelial_BinaryType_GSVA <-topTable(fit2,adjust="fdr", number=Inf)
write.csv(Epithelial_BinaryType_GSVA,"./Stem/Epithelial_BinaryType_GSVA.csv")

#barplot
plot<-GSVA_barplot(GSVA_result = Epithelial_BinaryType_GSVA, p = 0.05, n = 5)
ggsave("./Stem/Epithelial_BinaryType_GSVA.pdf", plot, width = 8,height = 3)


cols<-c("#E05B52","#B455CD","#52B79A")
p1<-RidgePlot(PA.Epithelial, features = "E_score1",cols = cols)
p2<-RidgePlot(PA.Epithelial, features = "M_score1",cols = cols)
p3<-RidgePlot(PA.Epithelial, features = "S_score1",cols = cols)
p<-CombinePlots(list(p1,p2,p3), ncol = 1, legend = "right")
ggsave("./Stem/Epithelial_score_ridgeplot.pdf", width = 3, height = 6)

library(ggplot2)
ggplot(PA.Epithelial@meta.data,aes(x=E_score1,y=malignScore))+
  stat_binhex()+
  scale_fill_gradient(low="lightblue",high="red",breaks=c(0,250,500,1000,2000,4000,6000),limits=c(0,6000))

PC <- row.names(PA.Epithelial@meta.data)[which(PA.Epithelial@meta.data[["Epithelial.annotation"]] %in% c("Stem_cell_A","Stem_cell_B","Stem_cell_C","EMT_like_A","EMT_like_B"))]
table(PA.Epithelial@meta.data[["Epithelial.annotation"]])
length(PC)
PA.PC<- subset(PA.Epithelial, cells=PC)
save(PA.PC,file="./Stem/PA.PC.RData")

Idents(PA.PC)<-"Epithelial.annotation"
cols<-c("#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E")
p1<-RidgePlot(PA.PC, features = "S_score1",cols = cols)
p2<-RidgePlot(PA.PC, features = "E_score1",cols = cols)
p3<-RidgePlot(PA.PC, features = "M_score1",cols = cols)
p<-CombinePlots(list(p1,p2,p3), ncol = 1, legend = "right")
ggsave("./Stem/PC_score_ridgeplot.pdf", width = 4, height = 6)

Idents(PA.PC)<-"Epithelial.annotation"
PC.markers <- FindAllMarkers(PA.PC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(PC.markers, "./Stem/PC.markers.csv")
PC.markers<- subset(PC.markers, avg_log2FC >1 & p_val_adj<0.05)

#top25
annotation_col<-data.frame(Stem_ID=PA.PC@meta.data[["Epithelial.annotation"]],
                           malignType=PA.PC@meta.data[["malignType"]],
                           row.names = rownames(PA.PC@meta.data))

p<-HeatmapforFindAllMarkers(Markers = PC.markers, seuratobj = PA.PC,
                         annotation_col = annotation_col,category = "Epithelial.annotation",
                         number = 25)

ggsave("./Stem/Heatmap_PC.pdf", p, width = 6,height = 6)



#GSVA
library(GSVA)

library(msigdbr)

m_df<- msigdbr(species = "Homo sapiens", category = "H")
gset.list<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

DefaultAssay(PA.PC) <- "RNA" 
expr <- as.matrix(GetAssayData(PA.PC, slot = 'data')) 
es <- gsva(expr, gset.list, mx.diff=TRUE, verbose=TRUE, parallel.sz=4)
save(es, file = "./Stem/es_PC.RData")

annotation_col <- data.frame(row.names =rownames(PA.PC@meta.data),
                             BinaryType=PA.PC@meta.data[["BinaryType"]],
                             Stem_ID=PA.PC@meta.data[["Epithelial.annotation"]])
plot<-GSVA_heatmap(es, PA.PC,annotation_col,category = "Epithelial.annotation")

ggsave("./Stem/PC_GSVA_heatmap.pdf",plot, width = 10, height = 5.7)

library(limma) 

meta <- PA.PC@meta.data[,c("Epithelial.annotation")]
group <- factor(meta, levels = unique(PA.PC@meta.data[["Epithelial.annotation"]]),ordered = F)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(es,design)
contrast.matrix <- makeContrasts(Stem_cell_B-Stem_cell_A, Stem_cell_C-Stem_cell_B, EMT_like_A-EMT_like_B,Stem_cell_C-(EMT_like_A+EMT_like_B),levels=design)#???????Ǻ??ߣ?
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff.all<-topTable(fit2, number = 243, resort.by = "logFC", adjust.method = "fdr") %>% subset(adj.P.Val<0.05)

StemB_A <- topTable(fit2, coef=1, resort.by = "logFC", adjust="fdr", number=Inf)
StemC_B <- topTable(fit2, coef=2, resort.by = "logFC", adjust="fdr", number=Inf)
EMTA_EMTB <- topTable(fit2, coef=3, resort.by = "logFC", adjust="fdr", number=Inf)
StemC_EMT <- topTable(fit2, coef=4, resort.by = "logFC", adjust="fdr", number=Inf)

library(ggplot2)
datalist<-list(StemB_A,StemC_B,EMTA_EMTB,StemC_EMT)
plotlist<-NULL 
for (i in 1:length(datalist)) {
  plotlist[[i]]=GSVA_barplot(datalist[[i]],p=0.05,n=5)
}
plot <- CombinePlots(list(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]]), ncol = 1, legend = "right")
ggsave("./Stem/GSVA_barplot_combine.pdf",plot, width = 8, height = 8)

ggsave("./Stem/GSVA_barplot_B-A.pdf",plotlist[[1]], width = 6, height = 2.4)
ggsave("./Stem/GSVA_barplot_EMTA-B.pdf",plotlist[[3]], width = 6, height = 2.4)

p1<-VlnPlot(PA.Epithelial, features = c("G2M.Score"), pt.size = 0,slot = "data", group.by = "Epithelial.annotation",cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Scores", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size=10))
p2<-VlnPlot(PA.Epithelial, features = c("S.Score"), pt.size = 0,slot = "data", group.by = "Epithelial.annotation",cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Scores", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA), 
        axis.text.x = element_text(size=10))
p<-CombinePlots(list(p1,p2),ncol = 1,legend = "right")
ggsave("./Epithelial_CellCycleScores.pdf", p, width = 4, height = 4.8)


ggplot(PA.PC@meta.data, aes(S.Score))+
  geom_histogram(binwidth = 0.001)+geom_vline(xintercept = 0.1)

ggplot(PA.PC@meta.data, aes(G2M.Score))+
  geom_histogram(binwidth = 0.001)+geom_vline(xintercept = 0.1)


phenodata<-data.frame(S.Score=PA.PC@meta.data[["S.Score"]],
                      G2M.Score=PA.PC@meta.data[["G2M.Score"]],
                      Phase=PA.PC@meta.data$Phase,
                      Cell.Type=PA.PC@meta.data$Epithelial.annotation,
                      row.names = rownames(PA.PC@meta.data))
phenodata$Cycling.Type[phenodata$S.Score>0.1 | phenodata$G2M.Score>0.1]="Cycling_Cell"
phenodata$Cycling.Type[phenodata$S.Score<0.1 & phenodata$G2M.Score<0.1]="Resting_Cell"

cols<-c("#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E")
p1<-ggplot(phenodata)+
  geom_point(aes(x=S.Score,y=G2M.Score, colour=Cycling.Type), size=0.6)+
  scale_color_manual(values =c("#FF59AF", "#388DFD"))+
  theme_bw()

p2<-ggplot(phenodata) +
  geom_bar(aes(x = Cell.Type, fill = Phase), position = "fill")+
  labs(y="Proportion")+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))+
  scale_fill_manual(values=cols)

p3<-ggplot(phenodata) +
  aes(x = Cell.Type, fill = Cycling.Type) +
  geom_bar(position = "fill")+
  labs(y="Proportion")+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))+
  scale_fill_manual(values=c("#37B8BD","#B452CD"))
p4<-p2|p3
ggsave("./Stem/Cycling.point.pdf", p1, width = 3.5, height = 2)
ggsave("./Stem/Cycling.bar.pdf", p4, width = 6, height = 3)

p1<-VlnPlot(PA.PC, features = c("G2M.Score"), pt.size = 0,slot = "data", group.by = "Epithelial.annotation",cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Scores", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size=8))
p2<-VlnPlot(PA.PC, features = c("S.Score"), pt.size = 0,slot = "data", group.by = "Epithelial.annotation",cols = cols)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(y="Scores", x="")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=1.5)+
  theme(panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size=8))
p<-p1/p2
ggsave("./Stem/G2Mscore.pdf", p, width = 5, height = 5)

rest<-subset(phenodata, Cycling.Type == "Resting_Cell")
Cycling<-subset(phenodata, Cycling.Type == "Cycling_Cell")
dat<-data.frame(Cycling_Cell=table(Cycling$Cell.Type),
                Resting_Cell=table(rest$Cell.Type))
rownames(dat)<-dat$Cycling_Cell.Var1
dat<-dat[6:10,c(2,4)] %>% t()
chisq.test(dat)
#Pearson's Chi-squared test

#data:  dat
#X-squared = 11.942, df = 4, p-value = 0.01779
fisher.test(dat)

#data:  dat
#p-value = 0.009771
#alternative hypothesis: two.sided

dat<- t(dat) %>% data.frame()
dat$CyclingP<-dat$Cycling_Cell.Freq/sum(dat$Cycling_Cell.Freq)
dat$RestP<-dat$Resting_Cell.Freq/sum(dat$Resting_Cell.Freq)
write.csv(dat, file = "Stem/Cycling.Type_statics.csv")

G1<-subset(phenodata, Phase == "G1")
G2M<-subset(phenodata, Phase == "G2M")
S<-subset(phenodata, Phase == "S")
dat<-data.frame(G1_Cell=table(G1$Cell.Type),
                G2M_Cell=table(G2M$Cell.Type),
                S_Cell=table(S$Cell.Type))
rownames(dat)<-dat$G1_Cell.Var1
dat<-dat[6:10,c(2,4,6)] %>%t()

chisq.test(dat)
#X-squared = 115.26, df = 8, p-value < 2.2e-16


library(Seurat)
library(dplyr)
expr <- GetAssayData(PA.PC, assay = "RNA", slot = 'data') 
dir.create("./Stem/GSEA")   
dir.create("./Stem/GSEA/input")

expr <- data.frame(NAME=rownames(expr), Description=rep('na', nrow(expr)), expr, stringsAsFactors=F)
write('#1.2', "./Stem/GSEA/input/expr.gct", ncolumns=1)
write(c(nrow(expr),(ncol(expr)-2)), "./Stem/GSEA/input/expr.gct", ncolumns=2, append=T, sep='\t')
write.table(expr, "./Stem/GSEA/input/expr.gct", row.names=F, sep='\t', append=T, quote=F)

tmp <- table(PA.PC@meta.data$Epithelial.annotation)
line.1 <- c((ncol(expr)-2), length(unique(PA.PC@meta.data[["Epithelial.annotation"]])), 1) 
line.3 <- as.character(PA.PC@meta.data$Epithelial.annotation)
line.2 <- c("#", unique(line.3))#line.2 ?? line.3?е?????˳??????һ??
write(line.1, './Stem/GSEA/input/group.cls', ncolumns=length(line.1), append=T, sep='\t')
write(line.2, './Stem/GSEA/input/group.cls', ncolumns=length(line.2), append=T, sep='\t')
write(line.3, './Stem/GSEA/input/group.cls', ncolumns=length(line.3), append=T, sep='\t')

write.csv(PA.PC@meta.data,file="./Stem/PC.Metadata.csv")

#GSEA
library(RColorBrewer)
#StemC vs. EMT
StemC_go <- read.table(file = "./Stem//GSEA/output/GO/StemC_EMT.Gsea.1639149130288/gsea_report_for_Stem_cell_C_1639149130288.tsv",
                    sep = '\t', header = T)
EMT_go <- read.table(file = "./Stem//GSEA/output/GO/StemC_EMT.Gsea.1639149130288/gsea_report_for_EMT_1639149130288.tsv",
                    sep = '\t', header = T)
StemC_go$Category<-substr(StemC_go$NAME, 1,4)
StemC_go$NAME<-substring(StemC_go$NAME,6)
StemC_go<-subset(StemC_go, FDR.q.val<0.05)

EMT_go$Category<-substr(EMT_go$NAME, 1,4)
EMT_go$NAME<-substring(EMT_go$NAME,6)
EMT_go<-subset(EMT_go, FDR.q.val<0.05)

top5P<-StemC_go %>% group_by(Category) %>% top_n(n = 3, wt = abs(NES))
top5R<-EMT_go %>% group_by(Category) %>% top_n(n = 3, wt = abs(NES))
top5<-rbind(top5P, top5R)
top5<-arrange(top5, Category, NES)
top5$NAME<-factor(top5$NAME,levels =top5$NAME)
p1 <- ggplot(top5, aes(NAME, NES)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR.q.val))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO term", y="NES", title="Stem_cell_C vs EMT:GO term NES from GSEA") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 7))
ggsave("./Stem/GSEA_barplot_StemC_EMT.pdf",p1, width = 13.5, height = 7)

#KEGG
StemC_kegg <- read.table(file = "./Stem//GSEA/output/KEGG/StemC_EMT.Gsea.1639147826243/gsea_report_for_Stem_cell_C_1639147826243.tsv",
                       sep = '\t', header = T)
EMT_kegg <- read.table(file = "./Stem//GSEA/output/KEGG/StemC_EMT.Gsea.1639147826243/gsea_report_for_EMT_1639147826243.tsv",
                      sep = '\t', header = T)
StemC_kegg$Category<-substr(StemC_kegg$NAME, 1,4)
StemC_kegg$NAME<-substring(StemC_kegg$NAME,6)
StemC_kegg<-subset(StemC_kegg, FDR.q.val<0.05)

EMT_kegg$Category<-substr(EMT_kegg$NAME, 1,4)
EMT_kegg$NAME<-substring(EMT_kegg$NAME,6)
EMT_kegg<-subset(EMT_kegg, FDR.q.val<0.05)

top5P<-StemC_kegg %>% group_by(Category) %>% top_n(n = 5, wt = abs(NES))
top5R<-EMT_kegg %>% group_by(Category) %>% top_n(n = 5, wt = abs(NES))
top5<-rbind(top5P, top5R)
top5<-arrange(top5, Category, NES)
top5$NAME<-factor(top5$NAME,levels =top5$NAME)
p2 <- ggplot(top5, aes(NAME, NES)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR.q.val))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="KEGG Pathway", y="NES", title="StemC:KEGG Pathway NES from GSEA") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 7))
ggsave("./Stem/GSEA_barplot_kegg_Stem_cell_C vs EMT.pdf",p2, width = 8.5, height = 3.5)

#hybrid epithelial/mesenchymal state and transition
load("~/scRNA/PA_combination_mito0.2/Double.score0.4/Step3 reannotation of Epithelial/PA.Epithelial.RData")
Idents(PA.Epithelial)<-"BinaryType"
p1<-VlnPlot(PA.Epithelial, features = c("KRT8","EPCAM","KRT18","CDH2","VIM","E_score1","M_score1","S_score1"), 
        pt.size = 0,cols = c("#E05B52","#B455CD","#52B79A"))
p2<-RidgePlot(PA.Epithelial, features = c("KRT8","EPCAM","KRT18","CDH2","VIM","E_score1","M_score1","S_score1"),
          group.by = "BinaryType",cols = c("#E05B52","#B455CD","#52B79A"))
plot<-CombinePlots(list(p1,p2),ncol = 2, legend = "right")
ggsave("./Stem/PC-DC_Stem_markers_vln.pdf", plot, width = 15, height = 7)

p1<-VlnPlot(PA.Epithelial, features = c("KRT8"), 
            pt.size = 0,cols = c("#E05B52","#B455CD","#52B79A"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(axis.text.x = element_text(angle = 0))
p2<-VlnPlot(PA.Epithelial, features = c("KRT18"), 
            pt.size = 0,cols = c("#E05B52","#B455CD","#52B79A"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(axis.text.x = element_text(angle = 0))
p3<-VlnPlot(PA.Epithelial, features = c("CDH2"), 
            pt.size = 0,cols = c("#E05B52","#B455CD","#52B79A"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(axis.text.x = element_text(angle = 0))
p4<-VlnPlot(PA.Epithelial, features = c("VIM"), 
            pt.size = 0,cols = c("#E05B52","#B455CD","#52B79A"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(axis.text.x = element_text(angle = 0))

plot<-CombinePlots(list(p1,p2,p3,p4),ncol = 2, legend = "right")
ggsave("./Stem/PC-DC_Stem_markers_vln_select.pdf", plot, width = 4.5, height = 4.8)

Idents(PA.PC)<-"Epithelial.annotation"
p1<-VlnPlot(PA.PC, features = c("KRT8","EPCAM","KRT18","CDH2","VIM","E_score1","M_score1","S_score1"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"), pt.size = 0,ncol = 3, cols = c("#B452CD","#826DAE","#37B8BD"))

p2<-RidgePlot(PA.PC, features = c("KRT8","EPCAM","KRT18","CDH2","VIM","E_score1","M_score1","S_score1"),
              idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"), group.by = "Epithelial.annotation",ncol = 3, cols = c("#B452CD","#826DAE","#37B8BD"))
plot<-CombinePlots(list(p1,p2),ncol = 2, legend = "right")
ggsave("./Stem/PC_subtype_markers_vln.pdf", plot, width = 18, height = 9)

p1<-VlnPlot(PA.PC, features = c("KRT8"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"), pt.size = 0,ncol = 1, cols = c("#B452CD","#826DAE","#37B8BD"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(x="",y="Expression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<-VlnPlot(PA.PC, features = c("KRT18"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"), pt.size = 0,ncol = 1, cols = c("#B452CD","#826DAE","#37B8BD"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(x="",y="Expression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-VlnPlot(PA.PC, features = c("CDH2"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"), pt.size = 0,ncol = 1, cols = c("#B452CD","#826DAE","#37B8BD"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(x="",y="Expression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p4<-VlnPlot(PA.PC, features = c("VIM"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"), pt.size = 0,ncol = 1, cols = c("#B452CD","#826DAE","#37B8BD"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  labs(x="",y="Expression")+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot<-CombinePlots(list(p1,p2,p3,p4),ncol = 2, legend = "right")
ggsave("./Stem/PC_subtype_markers_vln_select.pdf", plot, width = 5.4, height = 5.6)

DimPlot(PA.Epithelial, reduction = "pca", label = F, group.by = "Epithelial.annotation",
        cols =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35")
        ,pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))

#RNA velocity
Epithelial_obs<-read.csv("~/my_path/Python/scvelo/PC/PC_obs.csv", row.names = 1)
Epithelial_obs<-select(Epithelial_obs,c("nCount_spliced","nFeature_spliced","nCount_unspliced","nFeature_unspliced","nCount_ambiguous","nFeature_ambiguous",
                                        "initial_size_spliced","initial_size_unspliced","initial_size","n_counts",
                                        "velocity_self_transition","root_cells","end_points","velocity_pseudotime","latent_time"))
PA.PC<-AddMetaData(PA.PC,Epithelial_obs)
save(PA.PC, file="./Stem/PA.PC.RData")

p<-RidgePlot(PA.PC,features = "latent_time",assay = "RNA",group.by = "Epithelial.annotation",
          cols = c("#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E"))
ggsave(filename = "./Stem/latent_time_ridge.pdf", width = 6,height = 3)

Epithelial_obs<-read.csv("~/scRNA/Python/scvelo/PC/PC_obs.csv", row.names = 1)
a<-data.frame(latent_time=Epithelial_obs$latent_time,
              Scores=Epithelial_obs$E_score1,
              Type="E_score")
b<-data.frame(latent_time=Epithelial_obs$latent_time,
              Scores=Epithelial_obs$M_score1,
              Type="M_score")
c<-data.frame(latent_time=Epithelial_obs$latent_time,
              Scores=Epithelial_obs$S_score1,
              Type="S_score")

dat<-rbind(a,b,c)
p1<-ggplot(dat) +
  geom_smooth(aes(x=latent_time,y=Scores,color=Type),method = "gam",se = T, size=0.8, show.legend = T)+
  ylab(label="Scores")+theme_bw()

DefaultAssay(PA.PC)<-"RNA"
data <-PA.PC %>% GetAssayData(slot = "data")%>%data.frame()%>%t()
matrix<-data[,c("KRT8","EPCAM","KRT18","CDH2","VIM","COL11A2","SOX9","SOX2","LHX3")]%>%data.frame()

d<-data.frame(latent_time=Epithelial_obs$latent_time,
              Expression=matrix$EPCAM,
              Type="EPCAM")
e<-data.frame(latent_time=Epithelial_obs$latent_time,
              Expression=matrix$KRT8,
              Type="KRT8")
f<-data.frame(latent_time=Epithelial_obs$latent_time,
              Expression=matrix$KRT18,
              Type="KRT18")
dat<-rbind(d,e,f)
p2<-ggplot(dat) +
  geom_smooth(aes(x=latent_time,y=Expression,color=Type),method = "gam",se = T, size=0.8, show.legend = T)+theme_bw()


g<-data.frame(latent_time=Epithelial_obs$latent_time,
              Expression=matrix$SOX9,
              Type="SOX9")
h<-data.frame(latent_time=Epithelial_obs$latent_time,
              Expression=matrix$SOX2,
              Type="SOX2")
i<-data.frame(latent_time=Epithelial_obs$latent_time,
              Expression=matrix$LHX3,
              Type="LHX3")
dat<-rbind(g,h,i)
p3<-ggplot(dat) +
  geom_smooth(aes(x=latent_time,y=Expression,color=Type),method = "gam",se = T, size=0.8, show.legend = T)+theme_bw()


j<-data.frame(latent_time=Epithelial_obs$latent_time,
              Expression=matrix$CDH2,
              Type="CDH2")
k<-data.frame(latent_time=Epithelial_obs$latent_time,
              Expression=matrix$VIM,
              Type="VIM")
l<-data.frame(latent_time=Epithelial_obs$latent_time,
                 Expression=matrix$COL11A2,
                 Type="COL11A2")
dat<-rbind(j,k,l)
p4<-ggplot(dat) +
  geom_smooth(aes(x=latent_time,y=Expression,color=Type),method = "gam",se = T, size=0.8, show.legend = T)+theme_bw()

plot<-CombinePlots(list(p1,p2,p3,p4),ncol = 2)
ggsave("./Stem/PC_markers_smoothline.pdf", plot, width = 6.5, height = 3.5)

#GSEA

#StemB vs. StemA
StemA <- read.table(file = "./Stem//GSEA/output/GO/StemB_StemA.Gsea.1639151109954/gsea_report_for_Stem_cell_A_1639151109954.tsv",
                    sep = '\t', header = T)
StemB <- read.table(file = "./Stem//GSEA/output/GO/StemB_StemA.Gsea.1639151109954/gsea_report_for_Stem_cell_B_1639151109954.tsv",
                    sep = '\t', header = T)

StemA$Category<-substr(StemA$NAME, 1,4)
StemA$NAME<-substring(StemA$NAME,6)
StemA<-subset(StemA, FDR.q.val<0.05)

StemB$Category<-substr(StemB$NAME, 1,4)
StemB$NAME<-substring(StemB$NAME,6)
StemB<-subset(StemB, FDR.q.val<0.05)

top5A<-StemA %>% group_by(Category) %>% top_n(n = 5, wt = abs(NES))
top5B<-StemB %>% group_by(Category) %>% top_n(n = 5, wt = abs(NES))
top5<-rbind(top5A, top5B)
top5<-arrange(top5, Category, NES)
top5$NAME<-factor(top5$NAME,levels =top5$NAME)
p1 <- ggplot(top5, aes(NAME, NES)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR.q.val))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO term", y="NES", title=paste0(deparse(substitute(StemB)),":GO term NES from GSEA")) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 5))

library(stringr)
StemA<-StemA[which(str_detect(StemA$NAME, "STEM|DIFFERENTIATION")==TRUE & str_detect(StemA$NAME, "SYSTEM")==FALSE),]
StemB<-StemB[which(str_detect(StemB$NAME, "STEM|DIFFERENTIATION")==TRUE & str_detect(StemB$NAME, "SYSTEM")==FALSE),]
topA<-StemA %>% group_by(Category) %>% top_n(n = 10, wt = abs(NES))
topB<-StemB %>% group_by(Category) %>% top_n(n = 10, wt = abs(NES))
top10<-rbind(topA, topB)
top10<-arrange(top10, Category, NES)
top10$NAME<-factor(top10$NAME,levels =top10$NAME)
p2 <- ggplot(top10, aes(NAME, NES)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR.q.val))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO term", y="NES", title=paste0(deparse(substitute(StemB)),":GO term NES from GSEA")) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 5))
plot <- CombinePlots(list(p1,p2), ncol = 1, legend = "right")
ggsave("./Stem/GSEA_barplot_B_A.pdf",plot, width = 7, height = 9)

#StemC vs. StemB
StemB <- read.table(file = "./Stem//GSEA/output/GO/StemC_StemB.Gsea.1639151649708/gsea_report_for_Stem_cell_B_1639151649708.tsv",
                    sep = '\t', header = T)
StemC <- read.table(file = "./Stem//GSEA/output/GO/StemC_StemB.Gsea.1639151649708/gsea_report_for_Stem_cell_C_1639151649708.tsv",
                    sep = '\t', header = T)
StemC$Category<-substr(StemC$NAME, 1,4)
StemC$NAME<-substring(StemC$NAME,6)
StemC<-subset(StemC, FDR.q.val<0.05)

StemB$Category<-substr(StemB$NAME, 1,4)
StemB$NAME<-substring(StemB$NAME,6)
StemB<-subset(StemB, FDR.q.val<0.05)
top5C<-StemC %>% group_by(Category) %>% top_n(n = 3, wt = abs(NES))
top5A<-StemB %>% group_by(Category) %>% top_n(n = 3, wt = abs(NES))
top5<-rbind(top5C, top5A)
top5<-arrange(top5, Category, NES)
top5$NAME<-factor(top5$NAME,levels =top5$NAME)
top5$NES=-top5$NES
p3 <- ggplot(top5, aes(NAME, NES)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR.q.val))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO term", y="NES", title="B.vs.C:GO term NES from GSEA") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6))
ggsave("./Stem/GSEA_barplot_B_C.pdf",p3, width = 15, height = 6)

library(stringr)
StemB<-StemB[which(str_detect(StemB$NAME, "STEM|DIFFERENTIATION")==TRUE & str_detect(StemB$NAME, "SYSTEM")==FALSE),]
StemC<-StemC[which(str_detect(StemC$NAME, "STEM|DIFFERENTIATION")==TRUE & str_detect(StemC$NAME, "SYSTEM")==FALSE),]
topA<-StemB %>% group_by(Category) %>% top_n(n = 10, wt = abs(NES))
topC<-StemC %>% group_by(Category) %>% top_n(n = 10, wt = abs(NES))
top10<-rbind(topA, topC)
top10<-arrange(top10, Category, NES)
top10$NAME<-factor(top10$NAME,levels =top10$NAME)
top10$NES=-top10$NES
p4 <- ggplot(top10, aes(NAME, NES)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR.q.val))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO term", y="NES", title=paste0(deparse(substitute(StemB)),":GO term NES from GSEA")) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 7))
ggsave("./Stem/GSEA_barplot_B_C_DIFF.pdf",p4, width = 13, height = 5)

Idents(PA.PC)<-"Epithelial.annotation"
p1<-VlnPlot(PA.PC, features = c("KRT8","EPCAM","KRT18","SOX9","SOX2","LHX3"),cols = c("#B452CD","#826DAE","#37B8BD"),
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),pt.size = 0, ncol = 3)#ABC3群表达的干细胞marker谱不尽相�?
p2<-VlnPlot(PA.PC, features = c("COL1A2","PTPRC","IGKC","CD4","CD8A"), cols = c("#B452CD","#826DAE","#37B8BD"),
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),pt.size = 0, ncol = 3)
p3<-VlnPlot(PA.PC, features = c("CD74","HLA-A","HLA-B","HLA-C","HLA-DRA","SLPI"), cols = c("#B452CD","#826DAE","#37B8BD"),
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),pt.size = 0, ncol = 3) #StemC具有最强的免疫抑制作用（可能在免疫细胞分泌的细胞因子的作用下干细胞产生的亚群）
plot<-CombinePlots(list(p1,p3,p2),ncol = 2, legend = "right")
ggsave("./Stem/Stem_markers_vln.pdf", plot, width = 14, height = 10)

p1<-VlnPlot(PA.PC, features = c("CD74"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
            pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<-VlnPlot(PA.PC, features = c("HLA-A"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
           pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-VlnPlot(PA.PC, features = c("HLA-B"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
            pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p4<-VlnPlot(PA.PC, features = c("HLA-C"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
            pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p5<-VlnPlot(PA.PC, features = c("HLA-DRA"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
            pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p6<-VlnPlot(PA.PC, features = c("SLPI"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
            pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot<-CombinePlots(list(p1,p2,p3,p4,p5,p6),ncol = 3, legend = "right")
ggsave("./Stem/Stem_markers_vln_select.pdf", plot, width = 6.4, height = 5.5)


cells<-Cells(subset(PA.PC, Epithelial.annotation %in% c("Stem_cell_A","Stem_cell_B","Stem_cell_C")))
FeaturePlot(PA.PC, reduction = "umap", features = c("KRT8","EPCAM","KRT18"), 
            cells = cells,pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.PC, reduction = "umap", features = c("SOX9","SOX2","LHX3"), 
            cells = cells,pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
FeaturePlot(PA.PC, reduction = "umap", features = c("CD74","HLA-A","HLA-B","HLA-C","HLA-DRA","SLPI"), 
            cells = cells, pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))

#violin
p1<-VlnPlot(PA.PC, features = c("SOX2"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
            pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<-VlnPlot(PA.PC, features = c("SOX9"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
            pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-VlnPlot(PA.PC, features = c("LHX3"), 
            idents = c("Stem_cell_A","Stem_cell_B","Stem_cell_C"),
            cols = c("#B452CD","#826DAE","#37B8BD"),
            pt.size = 0, ncol = 1)+labs(y="Expression")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot<-CombinePlots(list(p1,p2,p3),ncol = 3, legend = "right")
ggsave("./Stem/Stem_markers_vln_Stem_select.pdf", plot, width = 7.5, height = 2.8)
#featureplot
library(RColorBrewer)
cols <- brewer.pal(n=6, name = "Oranges")
cells<-Cells(subset(PA.PC, Epithelial.annotation %in% c("Stem_cell_A","Stem_cell_B","Stem_cell_C")))
p1<-FeaturePlot(PA.PC, reduction = "umap", features = "SOX2", pt.size = 0.6, cols = cols, cells = cells)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<-FeaturePlot(PA.PC, reduction = "umap", features = "SOX9", pt.size = 0.6, cols = cols, cells = cells)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-FeaturePlot(PA.PC, reduction = "umap", features = "LHX3", pt.size = 0.6, cols = cols, cells = cells)+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot<-CombinePlots(list(p1,p2,p3),ncol = 3, legend = "right")
ggsave("./Stem/Stem_markers_featureplot_Stem_select.pdf", plot, width = 6.4, height = 2.2)

#SCENIC + pyscenic
dir.create("./Stem/SCENIC")
write.csv(t(as.matrix(PA.PC@assays$RNA@counts)), file = "./Stem/SCENIC/PC_pyscenic.csv")
library(tidyverse)
library(Seurat)
library(data.table)
library(SCENIC)
library(AUCell)
library(SCopeLoomR)
dir.create("./Stem/SCENIC")

loom<-open_loom("/home/liang/scRNA/Python/pyscenic/Epithelial/Stem/PC.SCENIC.loom")
regulons_incidMat<-get_regulons(loom, column.attr.name = "Regulons")
regulons<-regulonsToGeneLists(regulons_incidMat)
regulonAUC<-get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

TFtoGene<-data.frame(gene=unlist(regulons))
for (i in 1:dim(TFtoGene)[1]) {
  TFtoGene[i,"TF"] <- unlist(strsplit(rownames(TFtoGene)[i], "[(]"))[1]
}
write.csv(TFtoGene,"./Stem/SCENIC/regulonTargetsInfo.csv")

load("./Stem/PA.PC.RData")
library(pheatmap)
TF_matrix<-regulonAUC[,] %>% getAUC() 

PA.PC@assays$pyscenic<-CreateAssayObject(counts = TF_matrix)


Idents(PA.PC)<-"BinaryType"
EMT.TFs<-FindAllMarkers(PA.PC,assay = "pyscenic",slot = "data",logfc.threshold = 0)
write.csv(EMT.TFs, "./Stem/SCENIC/EMT.TFs.csv")

Idents(PA.PC)<-"Epithelial.annotation"
PC.TFs<-FindAllMarkers(PA.PC,assay = "pyscenic",slot = "data",logfc.threshold = 0)
write.csv(PC.TFs, "./Stem/SCENIC/PC.TFs.csv")
PC.TFs <- subset(PC.TFs, p_val_adj<0.05)

EMT.TFs <- subset(EMT.TFs, p_val_adj<0.05)
EMT_top <- EMT.TFs%>% group_by(cluster)%>% top_n(n=25, wt=avg_log2FC)

top<-EMT_top$gene
TF_matrix<-GetAssayData(PA.PC,assay = "pyscenic",slot = "counts")
new_cluster <- sort(PA.PC@active.ident)
AUC.matrix.top<-TF_matrix[top,names(new_cluster)]

AUC.matrix.top <- NormalizeData(AUC.matrix.top, scale.factor = 10)
AUC.matrix.top <- ScaleData(AUC.matrix.top, scale.max = 4) 

annotation_col<-data.frame(Stem_ID=PA.PC@meta.data$Epithelial.annotation,
                           CellType=PA.PC@meta.data[["BinaryType"]],
                           row.names = rownames(PA.PC@meta.data))

library(RColorBrewer)
brewer.pal.info
display.brewer.pal(11,"RdBu")
cols<-c(brewer.pal(11,"RdBu")[9],"white",brewer.pal(11,"RdBu")[2])

plot<-pheatmap(AUC.matrix.top,
               show_colnames =F,show_rownames = T, 
               cluster_rows = T,cluster_cols = F, 
               fontsize = 6,
               color = colorRampPalette(colors = cols)(100),
               annotation_col = annotation_col)
ggsave("./Stem/SCENIC/TFs_Heatmap_EMT_top.pdf",plot, width = 6.5,height = 4.5)

Idents(PA.PC)<-"Epithelial.annotation"
data<-FindMarkers(PA.PC,assay = "pyscenic",ident.1 = "EMT_like_A",ident.2 = "EMT_like_B",logfc.threshold = 0)
write.csv(data,"./Stem/SCENIC/EMTA-B.csv")
dat<-data%>%subset(p_val_adj<0.05)

genes<-subset(TFtoGene,TF %in% c("FOXO1","SP4"))
write.csv(genes,"./Stem/SCENIC/EMTA-B_topTFtarget.csv")
rm(PA.PC)


load("./Stem/PA.PC.RData")
PC <- row.names(PA.PC@meta.data)[which(PA.PC@meta.data[["BinaryType"]]=='PC')]
seuratobj<-subset(PA.Epithelial, cells  = PC)

seuratobj@assays$pyscenic<-CreateAssayObject(counts = TF_matrix[,PC])
save(seuratobj, file="./Stem/PA.Stem.TF.RData")

Idents(seuratobj)<-"Epithelial.annotation"
Stem.TFs<-FindAllMarkers(seuratobj,assay = "pyscenic",slot = "data",logfc.threshold = 0)
write.csv(Stem.TFs, "./Stem/SCENIC/Stem.TFs.csv")

library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

Stem.TFs <- subset(Stem.TFs, p_val_adj<0.05)
Stem.TFs <- Stem.TFs%>% group_by(cluster)%>% top_n(n=10, wt=avg_log2FC)

top<-Stem.TFs$gene
TF_matrix<-GetAssayData(seuratobj,assay = "pyscenic",slot = "counts")
new_cluster <- sort(seuratobj@active.ident)
AUC.matrix.top<-TF_matrix[top,names(new_cluster)]

AUC.matrix.top <- NormalizeData(AUC.matrix.top, scale.factor = 10)
AUC.matrix.top <- ScaleData(AUC.matrix.top, scale.max = 3.5) 

annotation_col<-data.frame(Stem_ID=seuratobj@meta.data$Epithelial.annotation,
                           malignType=seuratobj@meta.data[["malignType"]],
                           row.names = rownames(seuratobj@meta.data))

library(RColorBrewer)
brewer.pal.info
display.brewer.pal(9,"Oranges")
display.brewer.pal(9,"Blues")
cols<-c(brewer.pal(9,"Blues")[5],"white",brewer.pal(9,"Oranges")[6])

plot<-pheatmap(AUC.matrix.top, cellheight = 3,
               show_colnames =F,show_rownames = T, 
               cluster_rows = T,cluster_cols = F, 
               fontsize = 3,
               color = colorRampPalette(colors = cols)(100),
               annotation_col = annotation_col)
ggsave("./Stem/SCENIC/TFs_Heatmap_Stem_top.pdf",plot, width = 4,height = 5)

StemB_A<-FindMarkers(seuratobj, assay = "pyscenic", ident.1 = "Stem_cell_B",ident.2 = "Stem_cell_A",logfc.threshold = 0)
StemC_B<-FindMarkers(seuratobj, assay = "pyscenic", ident.1 = "Stem_cell_C",ident.2 = "Stem_cell_B",logfc.threshold = 0)

StemB_A$TFs<-rownames(StemB_A)
StemB_A$logpvalue= -log10(StemB_A$p_val_adj)
StemB_A$Significance = as.factor(ifelse(StemB_A$p_val_adj < 0.05, ifelse(StemB_A$avg_log2FC > 0 ,'Up','Down'),'No-diff'))

StemC_B$TFs<-rownames(StemC_B)
StemC_B$logpvalue= -log10(StemC_B$p_val_adj)
StemC_B$Significance = as.factor(ifelse(StemC_B$p_val_adj < 0.05, ifelse(StemC_B$avg_log2FC > 0 ,'Up','Down'),'No-diff'))

library(ggplot2)
library(ggrepel)
p1<-ggplot(StemB_A,aes(x=avg_log2FC,y=logpvalue,colour=Significance))+
  labs(x="log2(fold change)",y="-log10 (FDR)",title="StemB.vs.StemA:Differential TFs")+
  geom_point(size=2.5,alpha=0.6)+
  scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+xlim(c(-0.35, 0.35))+
  geom_text_repel(data = top_n(group_by(subset(StemB_A, Significance!="No-diff"), Significance),n=7,wt = abs(avg_log2FC)),
    aes(label = TFs),
    size = 3, max.overlaps = 20,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)

p2<-ggplot(StemC_B,aes(x=avg_log2FC,y=logpvalue,colour=Significance))+
  labs(x="log2(fold change)",y="-log10 (FDR)",title="StemC.vs.StemB:Differential TFs")+
  geom_point(size=2.5,alpha=0.6)+
  scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+xlim(c(-0.3, 0.3))+
  geom_text_repel(data = top_n(group_by(subset(StemC_B, Significance!="No-diff"), Significance),n=7,wt = abs(avg_log2FC)),
                  aes(label = TFs),
                  size = 3, max.overlaps = 25,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)

p<-CombinePlots(list(p1,p2),ncol = 2, legend = "right")
ggsave("./Stem/SCENIC/PC_subtype_TFs_volcano.pdf", p, width = 8, height = 4.5)

TF_matrix<-regulonAUC[,] %>% getAUC() 
AUC.data<-as.data.frame(t(TF_matrix))
PA.PC<- AddMetaData(PA.PC, AUC.data)
save(PA.PC,file = "./Stem/PA.PC.RData")

DimPlot(PA.PC, reduction = "umap", label = F, group.by = "Epithelial.annotation",
        cols = c("#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E"),pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))

library(RColorBrewer)
DefaultAssay(PA.PC)<-"RNA"
#StemB_A
col1 <- brewer.pal(n=6, name = "Blues")
#StemB
top_n(group_by(Stem.TFs,cluster),n=5, wt=avg_log2FC) #top5TFs
p1<-FeaturePlot(PA.PC, reduction = "umap", features = "SOX2...", pt.size = 0.6, cols = col1)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<-FeaturePlot(PA.PC, reduction = "umap", features = "SOX9...", pt.size = 0.6, cols = col1)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-FeaturePlot(PA.PC, reduction = "umap", features = "NFIX...", pt.size = 0.6, cols = col1)+
  theme(panel.border = element_rect(colour="black",fill=NA))
#StemC_A
col2 <- brewer.pal(n=6, name = "Purples")
#StemC
p4<-FeaturePlot(PA.PC, reduction = "umap", features = "SP4...", pt.size = 0.6, cols = col2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p5<-FeaturePlot(PA.PC, reduction = "umap", features = "ZFP90...", pt.size = 0.6, cols = col2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p6<-FeaturePlot(PA.PC, reduction = "umap", features = "EGR4...", pt.size = 0.6, cols = col2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
#StemA
col4 <- brewer.pal(n=6, name = "OrRd")
p7<-FeaturePlot(PA.PC, reduction = "umap", features = "FOXA2...", pt.size = 0.6, cols = col4)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p8<-FeaturePlot(PA.PC, reduction = "umap", features = "MESP1...", pt.size = 0.6, cols = col4)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p9<-FeaturePlot(PA.PC, reduction = "umap", features = "EHF...", pt.size = 0.6, cols = col4)+
  theme(panel.border = element_rect(colour="black",fill=NA))
#EMT
top_n(group_by(EMT.TFs,cluster),n=5, wt=avg_log2FC) 
col5 <- brewer.pal(n=6, name = "Greens")
p10<-FeaturePlot(PA.PC, reduction = "umap", features = "POU1F1...", pt.size = 0.6, cols = col5)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p11<-FeaturePlot(PA.PC, reduction = "umap", features = "HES6...", pt.size = 0.6, cols = col5)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p12<-FeaturePlot(PA.PC, reduction = "umap", features = "FOXO1...", pt.size = 0.6, cols = col5)+
  theme(panel.border = element_rect(colour="black",fill=NA))

plot<-(p7|p8|p9)/(p1|p2|p3)/(p10|p11|p12)/(p4|p5|p6)
ggsave("./Stem/SCENIC/PC_TFs_Featureplot.pdf",plot, width = 8.7,height = 9.3)

regulon.A <-c("FOXA2","MESP1","EHF")
regulon.B <-c("SOX2","SOX9","NFIX")
regulon.EMT <-c("POU1F1","HES6","FOXO1")

tableSubset.A<- subset(TFtoGene, TF %in% regulon.A)
tableSubset.B <- subset(TFtoGene, TF %in% regulon.B)
tableSubset.EMT <- subset(TFtoGene, TF %in% regulon.EMT)

MEyellow<-read.csv("./scWGCNA/module_yellow.csv")
MEgreen<-read.csv("./scWGCNA/module_green.csv")
MEyellow<-MEyellow$x
MEgreen<-MEgreen$x
targetGenesA<- subset(tableSubset.A, gene %in% MEyellow)
targetGenesB<- subset(tableSubset.B, gene %in% MEgreen)
targetGenesA$ME<-"MEyellow"
targetGenesB$ME<-"MEgreen"
write.csv(rbind(targetGenesA,targetGenesB),file = "./Stem/SCENIC/targetgenesAB.csv")

library(dplyr)
library(igraph)
library(ggraph)

net <- igraph::graph_from_data_frame(d=targetGenesA,directed = F)

igraph::V(net)$deg <- igraph::degree(net) 
igraph::V(net)$size <- igraph::degree(net)/5 
p1<-ggraph(net,layout = "stress")+ 
  geom_edge_fan(color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="#E7542B", alpha=0.7)+
  geom_node_text(aes(label=name), size = 3, repel = T, max.overlaps = 100)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10))+
  guides(size=F)+
  theme_graph()
ggsave("./Stem/SCENIC/StemA-GENE_network.pdf", p1, width = 5.5, height = 4.5)

net <- igraph::graph_from_data_frame(d=targetGenesB,directed = F)

igraph::V(net)$deg <- igraph::degree(net) 
igraph::V(net)$size <- igraph::degree(net)/5 
p2 <- ggraph(net,layout = "stress")+ 
  geom_edge_fan(color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="#A060ED", alpha=0.7)+
  geom_node_text(aes(label=name), size = 3, repel = T, max.overlaps = 100)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10))+
  guides(size=F)+
  theme_graph()
ggsave("./Stem/SCENIC/StemB-GENE_network.pdf", p2, width = 7.5, height = 4.5)


#monocle
library(monocle)
#B4
load("./PA.B4.RData")

pdf("./monocle/B4.Epithelial.annotation.pdf", onefile = T, width = 7, height = 5)
DimPlot(PA.B4, reduction = "umap", label = FALSE, group.by = "Epithelial.annotation", 
        cols =  c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35")
        ,pt.size = 0.6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

plot<-VlnPlot(PA.B4, features = "malignScore", pt.size = 0, 
              group.by = "Epithelial.annotation",
              cols = c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
ggsave("./monocle/B4_Epithelial_malignScore.pdf",plot, width = 6, height = 3.5)

data <- PA.B4@assays$RNA@counts
pd <- new('AnnotatedDataFrame', data = PA.B4@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
Epithelial.B4.cds <- newCellDataSet(data, phenoData = pd,
                                 featureData = fd,
                                 lowerDetectionLimit = 0.3,
                                 expressionFamily=negbinomial.size())
Epithelial.B4.cds <- estimateSizeFactors(Epithelial.B4.cds)
Epithelial.B4.cds <- estimateDispersions(Epithelial.B4.cds)
dim(Epithelial.B4.cds)

#QC
Epithelial.B4.cds <- detectGenes(Epithelial.B4.cds) 

expressed_genes <- row.names(subset(fData(Epithelial.B4.cds),num_cells_expressed >= 100))
Epithelial.B4.cds <- Epithelial.B4.cds[expressed_genes,]
dim(Epithelial.B4.cds)

valid_cells <- row.names(subset(pData(Epithelial.B4.cds),num_genes_expressed >= 200))
Epithelial.B4.cds <- Epithelial.B4.cds[,valid_cells] 
dim(Epithelial.B4.cds)

expdt <- exprs(Epithelial.B4.cds)

pData(Epithelial.B4.cds)$Total_mRNAs <- Matrix::colSums(expdt)
head(pData(Epithelial.B4.cds))

valid_cells2 <- pData(Epithelial.B4.cds)$Total_mRNAs < 1e6
Epithelial.B4.cds <- Epithelial.B4.cds[,valid_cells2]

dim(Epithelial.B4.cds)

disp_table <- dispersionTable(Epithelial.B4.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
Epithelial.B4.cds <- setOrderingFilter(Epithelial.B4.cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(Epithelial.B4.cds)
plot_pc_variance_explained(Epithelial.B4.cds, return_all = F) 

Epithelial.B4.cds <- reduceDimension(Epithelial.B4.cds, max_components = 2, num_dim = 5, 
                                  reduction_method = 'tSNE', verbose = T)
Epithelial.B4.cds <- clusterCells(Epithelial.B4.cds, num_clusters = 16)
plot_cell_clusters(Epithelial.B4.cds, color = "Cluster", show_cell_names = T)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dir.create("./monocle")
save(Epithelial.B4.cds, file = "./monocle/Epithelial.B4.cds.RData")

#marker
plot_cell_clusters(Epithelial.B4.cds, 1, 2, color = "Cluster", markers = "SOX2")+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_clusters(Epithelial.B4.cds, 1, 2, color = "Cluster", markers = "SOX9")+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_clusters(Epithelial.B4.cds, 1, 2, color = "Cluster", markers = "LHX3")+
  theme(panel.border = element_rect(colour="black",fill=NA))

pdf("./monocle/Monocle_Cell_cluster.pdf", width = 4, height = 4)
plot_cell_clusters(Epithelial.B4.cds, color = "Cluster", show_cell_names = T)+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_clusters(Epithelial.B4.cds, color = "seurat_clusters", show_cell_names = T)+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_clusters(Epithelial.B4.cds, color = "Epithelial.annotation", show_cell_names = T,
                   cols = c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

expressed_genes <- row.names(subset(fData(Epithelial.B4.cds),num_cells_expressed >= 100))
diff_test_res <- differentialGeneTest(Epithelial.B4.cds[expressed_genes,], fullModelFormulaStr = "~Epithelial.annotation")
diff_test_res <- diff_test_res %>% arrange(qval)
write.csv(diff_test_res, "./monocle/Epithelial_B4_annotation_DEGgenes.csv") 

library(ggplot2)
top20 <- diff_test_res[c(1:20),]$gene_short_name
top20 <- plot_cell_clusters(Epithelial.B4.cds, 1, 2, color = "Epithelial.annotation", markers = top20)+
  theme(panel.border = element_rect(colour="black",fill=NA))
ggsave("./monocle/B4_top20_cell_cluster.pdf", top20, width = 10,height = 8)

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
Epithelial.B4.cds <- setOrderingFilter(Epithelial.B4.cds, ordering_genes)
plot_ordering_genes(Epithelial.B4.cds)

#pseudotime trajectory building
Epithelial.B4.cds <- reduceDimension(Epithelial.B4.cds, max_components = 2,
                                  method = 'DDRTree') 
Epithelial.B4.cds <- orderCells(Epithelial.B4.cds)
save(Epithelial.B4.cds,file="./monocle/Epithelial.B4.cds.RData")

pdf("./monocle/B4_pseudotime_trajectory_preroot.pdf", width = 5,height = 5)
plot_cell_trajectory(Epithelial.B4.cds, color_by = "State",cell_size=0.8)+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_trajectory(Epithelial.B4.cds, color_by = "Epithelial.annotation",cell_size=0.8)+
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_trajectory(Epithelial.B4.cds, color_by = "malignType",cell_size=0.8)+
  scale_color_manual(values =c("orchid3","dodgerblue"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_trajectory(Epithelial.B4.cds, color_by = "Pseudotime",cell_size=0.8)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

pdf("./monocle/B4_pseudotime_trajectory_seperated.pdf", width = 5.8,height = 6)
plot_cell_trajectory(Epithelial.B4.cds, color_by = "Epithelial.annotation") +
  facet_wrap(~Epithelial.annotation, nrow = 3)+
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

Epithelial.B4.cds <- orderCells(Epithelial.B4.cds, root_state = 2, reverse = T)
pdf("./monocle/B4_pseudotime_trajectory_postroot.pdf", width = 5,height = 5)
plot_cell_trajectory(Epithelial.B4.cds, color_by = "State",cell_size=0.8)+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_trajectory(Epithelial.B4.cds, color_by = "Epithelial.annotation",cell_size=0.8)+
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  theme(panel.border = element_rect(colour="black",fill=NA))
plot_cell_trajectory(Epithelial.B4.cds, color_by = "Pseudotime",cell_size=0.8)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

pdf("./monocle/B4_pseudotime_trajectory_seurat_clusters.pdf", width = 8,height = 13)
plot_cell_trajectory(Epithelial.B4.cds, color_by = "seurat_clusters") +
  facet_wrap(~seurat_clusters, nrow = 6)+
  theme(panel.border = element_rect(colour="black",fill=NA))
dev.off()

BEAM_1 <- BEAM(Epithelial.B4.cds, branch_point = 1, cores = 16)
BEAM_1 <- BEAM_1[order(BEAM_1$qval),]
write.csv(BEAM_1, "./monocle/BEAM_1.csv")
BEAM_1 <- BEAM_1[,c("gene_short_name", "pval", "qval")]

library(RColorBrewer)
heatmap_BEAM_1 <-plot_genes_branched_heatmap(Epithelial.B4.cds[row.names(subset(BEAM_1,qval<1e-4)),],
                                             branch_point = 1,
                                             num_clusters = 5, 
                                             cores = 32,
                                             branch_labels = c("Triple_N branch", "Triple_P branch"),
                                             hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                             branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2
                                             use_gene_short_name = T,
                                             show_rownames = F,
                                             return_heatmap = T)
save(heatmap_BEAM_1, file="./monocle/heatmap_BEAM_1.RData")

p<-heatmap_BEAM_1$ph_res
ggsave("./monocle/branched_heatmap.pdf", p, width = 8.5,height = 7.5)
gene_group=heatmap_BEAM_1$annotation_row
gene_group$gene=rownames(gene_group)
table(gene_group$Cluster)
write.csv(gene_group, "./monocle/gene_group.csv")

#Cluster genes DAVID
#Triple_P branch (cluster1 2)
#cluster1
dat<-read.csv("./monocle/DAVID/cluster1.csv", header = T) %>% subset(FDR <0.05) %>% arrange(Count)
rownames(dat)<-dat$Term

#GO
go<-subset(dat, Category %in% c("GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT"))
top<-go%>% group_by(Category) %>% top_n(n=3, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,12)
top$Term<- factor(top$Term, levels = top$Term)

library(RColorBrewer)
p1<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO Term", y="Count", title="GO Term Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))

#kegg
kegg<-subset(dat, Category =="KEGG_PATHWAY")
top<-kegg%>% group_by(Category) %>% top_n(n=9, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,10)
top$Term<- factor(top$Term, levels = top$Term)

p2<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[2])+
  labs(x="KEGG Pathway", y="Count", title="KEGG Pathway Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))

p=p1/p2
ggsave("./monocle/Cluster1_enrichment.pdf", p, width = 6, height = 5)

#cluster2
dat<-read.csv("./monocle/DAVID/cluster2.csv", header = T) %>% subset(FDR <0.05) %>% arrange(Count)
rownames(dat)<-dat$Term

#GO
go<-subset(dat, Category %in% c("GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT"))
top<-go%>% group_by(Category) %>% top_n(n=3, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,12)
top$Term<- factor(top$Term, levels = top$Term)

library(RColorBrewer)
p3<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+#颜色变化不明�?
  #scale_fill_continuous(low=brewer.pal(9,"BuPu")[7], high=brewer.pal(9,"BuPu")[2], limits=c(0,1.360e-06))
  labs(x="GO Term", y="Count", title="GO Term Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))

#kegg
kegg<-subset(dat, Category =="KEGG_PATHWAY")
top<-kegg%>% group_by(Category) %>% top_n(n=9, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,10)
top$Term<- factor(top$Term, levels = top$Term)

p4<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[2])+
  labs(x="KEGG Pathway", y="Count", title="KEGG Pathway Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))

p=p3/p4
ggsave("./monocle/Cluster2_enrichment.pdf", p, width = 6, height = 5)

#Triple_N branch(cluster3 5)
#cluster3
dat<-read.csv("./monocle/DAVID/cluster3.csv", header = T) %>% subset(FDR <0.05) %>% arrange(Count)
rownames(dat)<-dat$Term

#GO
go<-subset(dat, Category %in% c("GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT"))
top<-go%>% group_by(Category) %>% top_n(n=3, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,12)
top$Term<- factor(top$Term, levels = top$Term)

library(RColorBrewer)
p5<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+#颜色变化不明�?
  #scale_fill_continuous(low=brewer.pal(9,"BuPu")[7], high=brewer.pal(9,"BuPu")[2], limits=c(0,1.360e-06))
  labs(x="GO Term", y="Count", title="GO Term Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))

#kegg
kegg<-subset(dat, Category =="KEGG_PATHWAY")
top<-kegg%>% group_by(Category) %>% top_n(n=9, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,10)
top$Term<- factor(top$Term, levels = top$Term)

p6<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[2])+
  labs(x="KEGG Pathway", y="Count", title="KEGG Pathway Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))

p=p5/p6
ggsave("./monocle/Cluster3_enrichment.pdf", p, width = 6, height = 5)

#cluster5
dat<-read.csv("./monocle/DAVID/cluster5.csv", header = T) %>% subset(FDR <0.05) %>% arrange(Count)
rownames(dat)<-dat$Term

#GO
go<-subset(dat, Category %in% c("GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT"))
top<-go%>% group_by(Category) %>% top_n(n=3, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,12)
top$Term<- factor(top$Term, levels = top$Term)

library(RColorBrewer)
p7<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+#颜色变化不明�?
  #scale_fill_continuous(low=brewer.pal(9,"BuPu")[7], high=brewer.pal(9,"BuPu")[2], limits=c(0,1.360e-06))
  labs(x="GO Term", y="Count", title="GO Term Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))

#kegg
kegg<-subset(dat, Category =="KEGG_PATHWAY")
top<-kegg%>% group_by(Category) %>% top_n(n=9, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,10)
top$Term<- factor(top$Term, levels = top$Term)

p8<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[2])+
  labs(x="KEGG Pathway", y="Count", title="KEGG Pathway Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))

p=p7/p8
ggsave("./monocle/Cluster5_enrichment.pdf", p, width = 6, height = 5)

p<-CombinePlots(list(p1,p3,p5,p7),ncol=2,legend = "right")
ggsave("./monocle/Cluster_enrichment_go_combine.pdf", p, width = 11, height = 5)

p<-CombinePlots(list(p2,p6,p4,p8),ncol=2,legend = "right")
ggsave("./monocle/Cluster_enrichment_kegg_combine.pdf", p, width = 11, height = 5)


#GSVA
library(GSVA)
library(msigdbr)
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
m_df<-m_df[which(str_detect(m_df$gs_name, "METABOLISM")==TRUE),]
gset.list<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
DefaultAssay(PA.B4) <- "RNA"
expr <- as.matrix(GetAssayData(subset(PA.B4, Epithelial.annotation %in% c("Triple_P","Triple_N")), slot = 'data'))
es.PN <- gsva(expr, gset.list, mx.diff=TRUE, verbose=TRUE, parallel.sz=16)
save(es.PN, file = "./monocle/es.PN.RData")

library(limma)
seuratobj<- subset(PA.B4, Epithelial.annotation %in% c("Triple_P","Triple_N"))
meta <- seuratobj@meta.data[,c("Epithelial.annotation")]
group <- factor(meta,levels = c("Triple_P","Triple_N"),ordered = F)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
es.sub<-es.PN
rownames(es.sub) <-gsub("HALLMARK_", "", rownames(es.sub))
fit <- lmFit(es.sub,design)
contrast.matrix <- makeContrasts(Triple_P-Triple_N,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
METABOLISM<-topTable(fit2,adjust="fdr", number=Inf)

METABOLISM<-subset(METABOLISM,adj.P.Val<0.05)

library(pheatmap)
es.heatmap <- as.data.frame(es.PN)[rownames(METABOLISM),]
rownames(es.heatmap) <-gsub("KEGG_", "", rownames(es.heatmap))
Idents(seuratobj)<-"Epithelial.annotation"
new_cluster <- sort(seuratobj@active.ident) 
es.heatmap <- as.matrix(es.heatmap[rownames(es.heatmap), names(new_cluster)])

es.heatmap <- ScaleData(es.heatmap, scale.max = 10)

annotation_col<- data.frame(row.names =rownames(seuratobj@meta.data),
                            Epithelial.annotation=seuratobj@meta.data[["Epithelial.annotation"]])

library(RColorBrewer)
brewer.pal.info
display.brewer.pal(11,"RdBu")
cols<-c(brewer.pal(11,"RdBu")[9],"white",brewer.pal(11,"RdBu")[2])

heatmap <- pheatmap::pheatmap(es.heatmap, show_rownames=T,
                    show_colnames=F, cellheight = 3.1,
                    color = colorRampPalette(colors = cols)(100),
                    annotation_col=annotation_col,
                    cluster_cols=F, cluster_rows = T,
                    fontsize_row=3)
ggsave("./monocle/GSVA_Heatmap_PN.pdf",heatmap, width = 6.5,height = 4)

library(ComplexHeatmap)
meta<-seuratobj@meta.data[names(new_cluster),] %>% select(c("Epithelial.annotation"))
group<-HeatmapAnnotation(df = meta,which = "column",simple_anno_size = unit(3, "mm"),
                         col = list(Epithelial.annotation=c("Triple_P"="#E868A1","Triple_N"="#3C81C4")),
                         annotation_name_gp=gpar(fontsize=6),annotation_label="Epithelial.annotation")

library(RColorBrewer)
brewer.pal.info
display.brewer.pal(11,"RdBu")
cols<-c(brewer.pal(11,"RdBu")[9],"white",brewer.pal(11,"RdBu")[2])

pdf(file = "./monocle/GSVA_Heatmap_PN_chm.pdf",width = 8.5,height = 4.5)
Heatmap(es.heatmap,
        col = colorRampPalette(colors = cols)(100),
        column_split = meta$Epithelial.annotation, column_gap = unit(0.7,"mm"),row_names_gp = gpar(fontsize = 6),
        cluster_rows = T,cluster_columns = F, show_column_names = F,show_row_dend = F,
        show_row_names = TRUE, use_raster = F,
        top_annotation = group,
        name = "Exp")
dev.off()

seuratobj<-AddMetaData(seuratobj, t(es.heatmap), col.name =rownames(es.heatmap))
save(seuratobj, file = "./monocle/PA.PN.GSVA.RData")

p1<-VlnPlot(seuratobj, features = c("PYRUVATE_METABOLISM"), 
            cols = c("#E868A1","#3C81C4"),
            pt.size = 0, ncol = 1)+labs(y="GSVA_score")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA), title =element_text(size = 6))
p2<-VlnPlot(seuratobj, features = c("GLUTATHIONE_METABOLISM"), 
            cols = c("#E868A1","#3C81C4"),
            pt.size = 0, ncol = 1)+labs(y="GSVA_score")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA), title =element_text(size = 6))
p3<-VlnPlot(seuratobj, features = c("NITROGEN_METABOLISM"), 
            cols = c("#E868A1","#3C81C4"),
            pt.size = 0, ncol = 1)+labs(y="GSVA_score")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA), title =element_text(size = 6))
p4<-VlnPlot(seuratobj, features = c("TAURINE_AND_HYPOTAURINE_METABOLISM"), 
            cols = c("#E868A1","#3C81C4"),
            pt.size = 0, ncol = 1)+labs(y="GSVA_score")+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.colour = NA)+
  stat_summary(fun=mean,geom="point",fill="white",shape=21,size=2)+
  theme(panel.border = element_rect(colour="black",fill=NA), title =element_text(size = 6))
plot<-CombinePlots(list(p1,p2,p3,p4),ncol = 2, legend = "right")
ggsave("./monocle//PN_GSVA_select_vln.pdf", plot, width = 4, height = 4.5)



load("~/my_path/Step3 reannotation of Epithelial/PA.B4.RData")
dir.create("./monocle/SCENIC")
write.csv(t(as.matrix(PA.B4@assays$RNA@counts)), file = "./monocle/SCENIC/B4_pyscenic.csv")

library(tidyverse)
library(data.table)
library(SCENIC)
library(AUCell)
library(SCopeLoomR)

loom<-open_loom("/my_path/Python/pyscenic/Epithelial/B4.SCENIC.loom")
regulons_incidMat<-get_regulons(loom, column.attr.name = "Regulons")
regulons<-regulonsToGeneLists(regulons_incidMat)
regulonAUC<-get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
TF_matrix<-regulonAUC[,] %>% getAUC() 

TFtoGene<-data.frame(gene=unlist(regulons))#转换regulons格式
for (i in 1:dim(TFtoGene)[1]) {
  TFtoGene[i,"TF"] <- unlist(strsplit(rownames(TFtoGene)[i], "[(]"))[1]
}
write.csv(TFtoGene,"./monocle/SCENIC/regulonTargetsInfo.csv")

PA.B4[["pyscenic"]] <- CreateAssayObject(counts = TF_matrix)
save(PA.B4, file="./monocle/SCENIC/PA.B4.RData")

diff.TFs<-FindAllMarkers(PA.B4,assay = "pyscenic",slot = "data",logfc.threshold = 0)
diff.TFs<-subset(diff.TFs,p_val_adj<0.05)

library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
topTFs <- diff.TFs %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)
Idents(PA.B4)<-"Epithelial.annotation"
new_cluster <- sort(PA.B4@active.ident)
AUC.matrix <- GetAssayData(PA.B4[["pyscenic"]], slot = "counts")

AUC.matrix <- as.matrix(AUC.matrix[topTFs$gene, names(new_cluster)])

AUC.matrix <- NormalizeData(AUC.matrix, scale.factor = 10)
AUC.matrix <- ScaleData(AUC.matrix, scale.max = 3) 
annotation_col<-data.frame(Epithelial.annotation=PA.B4@meta.data$Epithelial.annotation,
                           CellType=PA.B4@meta.data[["BinaryType"]],
                           row.names = rownames(PA.B4@meta.data))

library(RColorBrewer)
brewer.pal.info
display.brewer.pal(11,"RdBu")
cols<-c(brewer.pal(11,"RdBu")[9],"white",brewer.pal(11,"RdBu")[2])

plot<-pheatmap(AUC.matrix,cellheight = 3,
               show_colnames =F,show_rownames = T, 
               cluster_rows = T,cluster_cols = F,
               color = colorRampPalette(colors = cols)(100),
               fontsize = 3, 
               annotation_col = annotation_col)
ggsave("./monocle/SCENIC/TFs_Heatmap_alltop7.pdf",plot, width = 4.5,height = 6)

Idents(PA.B4)<-"Epithelial.annotation"
EMTA_EMTB<-FindMarkers(PA.B4,assay = "pyscenic",slot = "data",logfc.threshold = 0, ident.1 = "EMT_like_A",ident.2 = "EMT_like_B")
P_N<-FindMarkers(PA.B4,assay = "pyscenic",slot = "data",logfc.threshold = 0, ident.1 = "Triple_P",ident.2 = "Triple_N")


P_N$TFs<-rownames(P_N)
P_N$logpvalue= -log10(P_N$p_val_adj)
P_N$Significance = as.factor(ifelse(P_N$p_val_adj < 0.05, ifelse(P_N$avg_log2FC > 0 ,'Up','Down'),'No-diff'))

library(ggplot2)
library(ggrepel)
p1<-ggplot(P_N,aes(x=avg_log2FC,y=logpvalue,colour=Significance))+
  labs(x="log2(fold change)",y="-log10 (FDR)",title="Triple_P.vs.Triple_N:Differential TFs")+
  geom_point(size=2.5,alpha=0.6)+
  scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+xlim(c(-0.15, 0.15))+
  geom_text_repel(data = P_N %>% group_by(Significance) %>% top_n(n=5, wt=abs(avg_log2FC))%>% subset(Significance!="No-diff"),
                  aes(label = TFs),
                  size = 3, max.overlaps = 20,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
ggsave("./monocle/SCENIC/P-N_TFs_volcano.pdf", p1, width = 4.5, height = 4.5)

library(RColorBrewer)
DefaultAssay(PA.B4)<-"pyscenic"
#Triple_P
col1 <- brewer.pal(n=8, name = "Blues")
top_n(group_by(P_N,Significance),n=5, wt=abs(avg_log2FC)) #top5TFs
p1<-FeaturePlot(PA.B4, reduction = "umap", features = "CREB3(+)", pt.size = 0.6, cols = col1)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p2<-FeaturePlot(PA.B4, reduction = "umap", features = "NELFE(+)", pt.size = 0.6, cols = col1)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p3<-FeaturePlot(PA.B4, reduction = "umap", features = "TAF7(+)", pt.size = 0.6, cols = col1)+
  theme(panel.border = element_rect(colour="black",fill=NA))
#Triple_N
col2 <- brewer.pal(n=8, name = "OrRd")
p4<-FeaturePlot(PA.B4, reduction = "umap", features = "SOX7(+)", pt.size = 0.6, cols = col2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p5<-FeaturePlot(PA.B4, reduction = "umap", features = "ZNF91(+)", pt.size = 0.6, cols = col2)+
  theme(panel.border = element_rect(colour="black",fill=NA))
p6<-FeaturePlot(PA.B4, reduction = "umap", features = "TAF1(+)", pt.size = 0.6, cols = col2)+
  theme(panel.border = element_rect(colour="black",fill=NA))

plot<-(p1|p2|p3)/(p4|p5|p6)
ggsave("./monocle/SCENIC/P-N_TFs_Featureplot.pdf",plot, width = 8.7,height = 4.8)

PA.B4@meta.data[["Monocle.Pseudotime"]]<-Epithelial.B4.cds@phenoData@data[["Pseudotime"]]
PA.B4@meta.data[["Monocle.State"]]<-Epithelial.B4.cds@phenoData@data[["State"]]
VlnPlot(PA.B4, features = "Monocle.Pseudotime", sort = F,
        group.by = "Monocle.State",pt.size = 0)
VlnPlot(PA.B4, features = "Monocle.Pseudotime", sort = F,
        group.by = "Epithelial.annotation",pt.size = 0)

plot<-RidgePlot(PA.B4, features = "Monocle.Pseudotime",cols = c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))
ggsave("./monocle/B4_Pseudotime_ridge.pdf", plot, width = 6.5,height = 4)
save(PA.B4,file="./monocle/PA.B4.RData")

phenodata <- select(PA.B4@meta.data, c("Epithelial.annotation","seurat_clusters","malignScore",
                                 "Monocle.Pseudotime","Monocle.State"))

AUC.matrix <- GetAssayData(PA.B4[["pyscenic"]], slot = "counts")
AUC.data<-AUC.matrix[c("CREB3(+)","NELFE(+)","TAF7(+)",
                       "SOX7(+)","ZNF91(+)","TAF1(+)"),]%>%t()%>%data.frame()
phenodata<-cbind(phenodata,AUC.data)
colnames(phenodata) <- c("Epithelial.annotation","seurat_clusters","malignScore",
                         "Monocle.Pseudotime","Monocle.State",
                         "CREB3","NELFE","TAF7",
                         "SOX7","ZNF91","TAF1")
library(ggplot2)
TFs_selected<-c("CREB3","NELFE","TAF7",
                "SOX7","ZNF91","TAF1","malignScore")                        

plotlist<-list()
for (i in 1:length(TFs_selected)) {
  plotlist[[i]]<-ggplot(phenodata)+
    geom_point(aes_string(x="Monocle.Pseudotime",y=TFs_selected[i], colour="Epithelial.annotation"), size=0.6)+
    scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
    geom_smooth(aes_string(x="Monocle.Pseudotime",y=TFs_selected[i]), method = "gam",se = F, colour= brewer.pal(9,"Blues")[8], size=0.5)+
    theme_bw()
}
plot<-CombinePlots(list(plotlist[[1]],plotlist[[2]],plotlist[[3]],
                        plotlist[[4]],plotlist[[5]],plotlist[[6]],
                        plotlist[[7]]), ncol = 3,legend = "right")
ggsave("./monocle/SCENIC/B4_SCENIC_TFs.pdf",plot, width = 9.5,height = 7)

gene_group<-read.csv("./monocle/gene_group.csv", row.names = 1)
BEAM_1<-read.csv("./monocle/BEAM_1.csv", row.names = 1)

PA_TFs<-c("TBX19","NR5A1","GATA2","SF1","POU1F1")
top_genes <- row.names(subset(fData(Epithelial.B4.cds), gene_short_name %in% PA_TFs))
p1<-plot_genes_branched_pseudotime(Epithelial.B4.cds[top_genes,],
                                   branch_point = 1,
                                   branch_labels = c("Triple_N branch", "Triple_P branch"),
                                   color_by = "Epithelial.annotation",
                                   ncol = 3)+ 
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))
p2<-plot_genes_in_pseudotime(Epithelial.B4.cds[top_genes,], 
                             color_by = "Epithelial.annotation",
                             ncol = 3)+ 
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))
p<-CombinePlots(list(p1,p2), ncol = 1, legend = "right")
ggsave("./monocle/monocle_PA_TFs.pdf", p, width = 8, height = 3.5)


TFs_N<-c("SOX7","ZNF91","TAF1")

targetGeneN<-subset(TFtoGene, TF %in% TFs_N)

genelistN<-intersect(targetGeneN$gene, subset(gene_group, Cluster%in%c(3,5))$gene)

#Triple_P Triple_N
Idents(PA.B4)<- "Epithelial.annotation"
DEGs <- FindMarkers(PA.B4, assay = "RNA", ident.1 = "Triple_P", ident.2 = "Triple_N", logfc.threshold = 0.25)
write.csv(DEGs, "./monocle/P-N_degs.csv")

DEGs$gene<-rownames(DEGs)
DEGs$Significance = as.factor(ifelse(DEGs$p_val_adj < 0.05& abs(DEGs$avg_log2FC)> 1, ifelse(DEGs$avg_log2FC> 1 ,'Up','Down'),'No-diff'))

genelistN<-intersect(genelistN,subset(DEGs, Significance=="Down")$gene)#19�?

BEAM_N<-BEAM_1[genelistN,] %>% arrange(desc(num_cells_expressed),qval)
genelistN<-head(BEAM_N,10)$gene

m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
m_df<-m_df[which(str_detect(m_df$gs_name, "METABOLISM")==TRUE),]
genelistP<-intersect(targetGeneP$gene, subset(gene_group, Cluster%in%c(1,2))$gene)
genelistP<-intersect(genelistP, m_df$human_gene_symbol)

BEAM_P<-BEAM_1[genelistP,] %>% arrange(desc(num_cells_expressed),qval)
genelistP<-head(BEAM_P,10)$gene

top_genes <- row.names(subset(fData(Epithelial.B4.cds), gene_short_name %in% unique(c(genelistP,genelistN))))
p1<-plot_genes_branched_pseudotime(Epithelial.B4.cds[top_genes,],
                                   branch_point = 1,
                                   branch_labels = c("Triple_N branch", "Triple_P branch"),
                                   color_by = "Epithelial.annotation",
                                   ncol = 4)+ 
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))
p2<-plot_genes_in_pseudotime(Epithelial.B4.cds[top_genes,], 
                             color_by = "Epithelial.annotation",
                             ncol = 4)+ 
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))
p<-CombinePlots(list(p1,p2), ncol = 1, legend = "right")
ggsave("./monocle/monocle_DEGs.pdf", p, width = 12, height = 12)
  
P<-c("GNAS","ASAH1","MGST3")
N<-c("MEG3","N4BP2L2","MALAT1")
top_genes <- row.names(subset(fData(Epithelial.B4.cds), gene_short_name %in% P))
p1<-plot_genes_branched_pseudotime(Epithelial.B4.cds[top_genes,],
                                   branch_point = 1,
                                   branch_labels = c("Triple_N branch", "Triple_P branch"),
                                   color_by = "Epithelial.annotation",
                                   ncol = 3)+ 
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  labs(y="Expression")
p2<-plot_genes_in_pseudotime(Epithelial.B4.cds[top_genes,], 
                             color_by = "Epithelial.annotation",
                             ncol = 3)+ 
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  labs(y="Expression")
#
top_genes <- row.names(subset(fData(Epithelial.B4.cds), gene_short_name %in% N))
p3<-plot_genes_branched_pseudotime(Epithelial.B4.cds[top_genes,],
                                   branch_point = 1,
                                   branch_labels = c("Triple_N branch", "Triple_P branch"),
                                   color_by = "Epithelial.annotation",
                                   ncol = 3)+ 
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  labs(y="Expression")
p4<-plot_genes_in_pseudotime(Epithelial.B4.cds[top_genes,], 
                             color_by = "Epithelial.annotation",
                             ncol = 3)+ 
  scale_color_manual(values =c("#DE302E","#E05B41","#F39F16","#E868A1","#3C81C4","#B452CD","#826DAE","#37B8BD","#53B77A","#C6BC3E","#D37C35"))+
  labs(y="Expression")

p<-CombinePlots(list(p1,p3,p2,p4), ncol = 2, legend = "bottom")
ggsave("./monocle/monocle_DEGs_select.pdf", p, width = 9.5, height = 3.6)

library(ggplot2)
library(ggrepel)
p1<-ggplot(DEGs,aes(x=avg_log2FC,y=logpvalue,colour=Significance))+
  labs(x="log2(fold change)",y="-log10 (p.adjust)",title="Triple_P.vs.Triple_N:DEGs")+
  geom_point(size=2,alpha=0.6)+
  scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_text_repel(data = DEGs %>% group_by(Significance) %>% top_n(n=5, wt=abs(avg_log2FC)) %>% subset(Significance !="No-diff"),
                  aes(label = gene),
                  size = 3, max.overlaps = 40,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
ggsave("./monocle/P-N_DEGs_vocalno.pdf",p1, width = 4, height = 4)


#WGCNA
rm(list = ls())
gc()
options(stringsAsFactors = F)

library(WGCNA)
library(Seurat)
library(tidyverse)
library(reshape2)
library(stringr)
library(flashClust)
enableWGCNAThreads(nThreads = 16)

load("~/my_path/Step3 reannotation of Epithelial/PA.Epithelial.RData")
datadf <- as.matrix(PA.Epithelial@assays$RNA@data)
meta<-PA.Epithelial@meta.data[,c("Epithelial.annotation","seurat_clusters")]

#top5000
PA.Epithelial <- FindVariableFeatures(PA.Epithelial, selection.method = "vst", nfeatures = 5000)
datadf<-datadf[Seurat::VariableFeatures(PA.Epithelial),]
datExpr<-t(datadf)

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Choosing a soft-threshold to fit a scale-free topology to the network
powers = c(c(1:10), seq(from = 12, to=30, by=2));
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,
                      powerVector = powers,corFnc = cor,
                      corOptions = list(use = 'p'),
                      networkType = "unsigned")

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

power = sft$powerEstimate
softPower  = power
softPower

type="unsigned"
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}
softPower  = power
softPower=1
#calclute the adjacency matrix
adj= adjacency(datExpr,type = "unsigned", power = softPower);

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower);

colnames(TOM) =rownames(TOM) =colnames(datExpr)
dissTOM=1-TOM

#Module detection
#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree 
plot(geneTree, xlab="", sub="",cex=0.3);
# Set the minimum module size
minModuleSize = 20;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);

table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf("./scWGCNA/Gene-module-tree.pdf", width = 8, height = 4)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")
dev.off()

#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower)

colnames(diss1) =rownames(diss1) =colnames(datExpr)[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
pdf("./scWGCNA/Gene-module-tree-rest.pdf", width = 8, height = 4)
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()
diag(diss1) = NA;


sizeGrWindow(7,7)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))

#Extract modules
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=colnames(datExpr)[which(dynamicColors==color)]
  write.csv(module, paste0("./scWGCNA/module_",color, ".csv"), row.names=FALSE, col.names=FALSE,quote=FALSE)
}

#Quantify module similarity by eigengene correlation.
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = orderMEs(MEList$eigengenes)

pdf("./scWGCNA/module-module-relationship.pdf", width = 8, height = 6)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
dev.off()

design=model.matrix(~0+ meta$Epithelial.annotation)
colnames(design)=levels(meta$Epithelial.annotation)
# Recalculate MEs with color labels
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

pdf("./scWGCNA/Module-trait-relationships.pdf", width = 7, height = 6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

save.image(file="./scWGCNA/WGCNA.RData")


modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

StemA = as.data.frame(design[,6])
names(StemA) = "StemA"
geneTraitSignificance = as.data.frame(cor(datExpr, StemA, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(StemA), sep="");
names(GSPvalue) = paste("p.GS.", names(StemA), sep="");

module = "red"
column = match(module, modNames);
moduleGenes = dynamicColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Stem_cell_A",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)#这个图作用不�?

# Select module
module = "yellow"
# Select module probes
probes = colnames(TOM)
inModule = (dynamicColors==module)
modProbes = probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[modProbes, modProbes]
dimnames(modTOM) = list(modProbes, modProbes)
 
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("./scWGCNA/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("./scWGCNA/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = dynamicColors[inModule]
)

# Select module
module = "green"

probes = colnames(TOM) 
inModule = (dynamicColors==module)
modProbes = probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[modProbes, modProbes]
dimnames(modTOM) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("./scWGCNA/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("./scWGCNA/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = dynamicColors[inModule]
)

save.image(file = "./scWGCNA/WGCNA.RData")

#WGCNA
#1.modules (DAVID)
#module_yellow
dat<-read.csv("./scWGCNA/DAVID/module_yellow.csv", header = T) %>% subset(FDR<0.05) %>% arrange(Category,Count)
rownames(dat)<-dat$Term

#GO
go<-subset(dat, Category %in% c("GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT"))
top<-go%>% group_by(Category) %>% top_n(n=5, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,12)
top$Term<- factor(top$Term, levels = top$Term)

library(RColorBrewer)
p1<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO Term", y="Count", title="yellow:GO Term Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))
ggsave("./scWGCNA/DAVID/module_yellow_go_enrichment.pdf", p1, width = 6.5, height = 3.5)

#kegg
kegg<-subset(dat, Category =="KEGG_PATHWAY")
top<-kegg%>% group_by(Category) %>% top_n(n=5, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,10)
top$Term<- factor(top$Term, levels = top$Term)
p2<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO Term", y="Count", title="yellow:KEGG Pathway Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))
ggsave("./scWGCNA/DAVID/module_yellow_kegg_enrichment.pdf", p2, width = 5, height = 2)

#module_green
dat<-read.csv("./scWGCNA/DAVID/module_green.csv", header = T) %>% subset(FDR<0.05) %>% arrange(Category,Count)
rownames(dat)<-dat$Term

#GO
go<-subset(dat, Category %in% c("GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT"))
top<-go%>% group_by(Category) %>% top_n(n=5, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,12)
top$Term<- factor(top$Term, levels = top$Term)

library(RColorBrewer)
p1<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO Term", y="Count", title="green:GO Term Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))
ggsave("./scWGCNA/DAVID/module_green_go_enrichment.pdf", p1, width = 5, height = 3.5)

#kegg
kegg<-subset(dat, Category =="KEGG_PATHWAY")
top<-kegg%>% group_by(Category) %>% top_n(n=5, wt=Count) %>% data.frame()
top$Term<-substring(top$Term ,10)
top$Term<- factor(top$Term, levels = top$Term)
p2<-ggplot2::ggplot(top, aes(Term, Count)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  geom_col(aes(fill= FDR))+
  scale_fill_gradient(low=brewer.pal(9,"Blues")[7], high=brewer.pal(9,"Blues")[4])+
  labs(x="GO Term", y="Count", title="green:KEGG Pathway Enrichment") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6.5))
ggsave("./scWGCNA/DAVID/module_green_kegg_enrichment.pdf", p2, width = 5, height = 1.6)
