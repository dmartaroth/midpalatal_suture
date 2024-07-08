# ## ######################################## ## #
#                 ANNOTATION SCRNA-SEQ MPS       #
# ## ######################################## ## #

# Updated on: June 8, 2023
# Code should be updated for reproducible file paths and saving

#Load RDS object from PopsicleR preprocessing
mps<-readRDS("data-output/RDS objects/mps.Rds")

#plot generation
library(Seurat)
library(tidyverse)
library(popsicleR)
library(dplyr)
library(patchwork)
library(ggplot2)
library(dittoSeq)
library(scCustomize)
library(Nebulosa)
library(effects)
library(cowplot)

pdf("figures/umap_top3dot_mps.pdf",width=20,height = 7)
p1<-dittoDimPlot(mps, "ident", do.label = TRUE,labels.repel = TRUE,labels.size = 2.5,labels.highlight = FALSE,opacity=0.9,do.ellipse = TRUE,
                 color.panel = c("darkolivegreen2", "orange", "purple", "lightcoral", "skyblue","maroon2","slateblue","gold","dodgerblue3","plum1","darkseagreen3","pink","violetred4","tomato1","orchid3","darkolivegreen3")) +labs(title = 'E15 midpalatal suture')

#Plot top 3 markers
cluster_markers <- FindAllMarkers(
  mps,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1.0,
  method='MAST'
)

top3<-cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

p2<-DotPlot(
  object = mps, features =   top3$gene,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15 midpalatal suture")+
  theme(plot.title = element_text(color="blue",size=12,face="bold"),plot.subtitle=element_text(color="red",face="italic"))

plot_grid(p1,p2)
dev.off()

pdf("figures/umap_mps.pdf",width=9,height = 7)
p1
dev.off()

#Find top genes for each cluster for annotation
#Export top 50 markers for each cluster as csv
cluster0.markers <- FindMarkers(mps, ident.1 = 0, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster0.markers),n=50),file="data-output/cluster0.markers.csv")
##mesenchyme

cluster1.markers <- FindMarkers(mps, ident.1 = 1, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster1.markers),n=50),file="data-output/cluster1.markers.csv")
##mesenchyme

cluster2.markers <- FindMarkers(mps, ident.1 = 2, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster2.markers),n=50),file="data-output/cluster2.markers.csv")
##mesenchyme

cluster3.markers <- FindMarkers(mps, ident.1 = 3, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster3.markers),n=50),file="data-output/cluster3.markers.csv")
##epithelium

cluster4.markers <- FindMarkers(mps, ident.1 = 4, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster4.markers),n=50),file="data-output/cluster4.markers.csv")
##neural

cluster5.markers <- FindMarkers(mps, ident.1 = 5, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster5.markers),n=50),file="data-output/cluster5.markers.csv")
##erythrocytes

cluster6.markers <- FindMarkers(mps, ident.1 = 6, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster6.markers),n=50),file="data-output/cluster6.markers.csv")
##"cerebral cortex" and "brain came up; unknown, but suspect neural

cluster7.markers <- FindMarkers(mps, ident.1 = 7, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster7.markers),n=50),file="data-output/cluster7.markers.csv")
##vascular

cluster8.markers <- FindMarkers(mps, ident.1 = 8, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster8.markers),n=50),file="data-output/cluster8.markers.csv")
##immune

cluster9.markers <- FindMarkers(mps, ident.1 = 9, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster9.markers),n=50),file="data-output/cluster9.markers.csv")
##neural

cluster10.markers <- FindMarkers(mps, ident.1 = 10, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster10.markers),n=50),file="data-output/cluster10.markers.csv")
##immune

cluster11.markers <- FindMarkers(mps, ident.1 = 11, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster11.markers),n=50),file="data-output/cluster11.markers.csv")
##neural

cluster12.markers <- FindMarkers(mps, ident.1 = 12, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster12.markers),n=50),file="data-output/cluster12.markers.csv")
##smooth muscle? synapses? predict muscle

cluster13.markers <- FindMarkers(mps, ident.1 = 13, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster13.markers),n=50),file="data-output/cluster13.markers.csv")
#neural?

cluster14.markers <- FindMarkers(mps, ident.1 = 14, min.pct = 0.25,only.pos = TRUE)
write.csv(head(rownames(cluster14.markers),n=50),file="data-output/cluster14.markers.csv")
##ciliated cells


#Predicted cluster identities
cluster<-c(0,1:14)
identity<-c('mesenchyme','mesenchyme','mesenchyme','epithelium','neural','eruthrocytes','neural','vascular','immune','neural','immune','neural','muscle','neural','cilia')
df<-data.frame(cluster,identity)
print(df)
pdf("figures/cluster_ident.pdf",width=3,height = 5)
p<-tableGrob(df)
grid.arrange(p)
dev.off()

#plot known suture markers
pdf("figures/DotPlot_markergenes.pdf",width=6,height = 4)
"marker_genes"<-c("Prrx1","Gli1","Axin2","Crabp1","Runx2","Sp7","Dmp1","Sost","Phex","Acp5","Nfatc1","Acan","Col2a1","Pecam1","Vwf")
DotPlot(
  object = mps, features =   marker_genes,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15 midpalatal suture")+
  theme(plot.title = element_text(color="blue",size=12,face="bold"),plot.subtitle=element_text(color="red",face="italic"))
dev.off()

pdf("figures/clustered_DotPlot_markergenes.pdf",width=10,height = 10)
Clustered_DotPlot(seurat_object = mps, features = marker_genes, k = 6)
Clustered_DotPlot(seurat_object = mps, features = top3$gene, k = 10)
dev.off()

#Generate plots of known markers
pdf("figures/FeaturePlot_osteocyte.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Dmp1","Sost"),pt.size = 0.5)
dev.off()

pdf("figures/FeaturePlot_osteoblast.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Runx2","Sp7"),pt.size = 0.5)
dev.off()

pdf("figures/FeaturePlot_osteoclast.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Acp5","Nfatc1"),pt.size = 0.5)
dev.off()

pdf("figures/FeaturePlot_chondrocyte.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Acan","Col2a1"),pt.size = 0.5)
dev.off()

pdf("figures/FeaturePlot_vasculature.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Pecam1","Vwf"),pt.size = 0.5)
dev.off()

pdf("figures/FeaturePlot_progenitors.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Prrx1","Gli1"),pt.size = 0.5)
dev.off()

pdf("figures/FeaturePlot_proprioceptors.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Runx3","Ntrk3"),pt.size = 0.5)
dev.off()

pdf("figures/FeaturePlot_neuronal_mechanoreceptors.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Ret","Ntrk2"),pt.size = 0.5)
dev.off()

pdf("figures/FeaturePlot_proneurogenic.pdf",width=9,height = 5)
FeaturePlot(mps,features=c("Sox2","Sox5"),pt.size = 0.5)
dev.off()


df<-data.frame(cluster,identity)
print(df)

#Assign new cluster identities based on cluster markers and UMAP relationships
levels(mps)
new.cluster.ids<-c("mes.1","mes.3","mes.2","epith","neu.2","eryth","neu.3","vasc","imm.1","neu.5","imm.2","neu.4","musc","neu.1","cilia")
names(new.cluster.ids)<-levels(mps)
mps$new.cluster.ids <- Idents(mps)
mps<-RenameIdents(mps,new.cluster.ids)
mps@meta.data$seurat_clusters <-  mps@meta.data$new.cluster.ids
DimPlot(mps)
Idents(mps)<-factor(x=Idents(mps),levels=sort(levels(mps)))
levels(mps)
DimPlot(mps)

pdf("figures/umap_mps_annotated_clusters.pdf",width=8,height = 8)
dittoDimPlot(mps, "ident", do.label = TRUE,labels.repel = TRUE,labels.size = 2.5,labels.highlight = FALSE,opacity=0.9,do.ellipse = TRUE,
             color.panel = c("darkolivegreen2", "orange", "purple", "lightcoral", "skyblue","maroon2","slateblue","gold","dodgerblue3","plum1","darkseagreen3","pink","violetred4","tomato1","orchid3","darkolivegreen3")) +labs(title = 'E15 midpalatal suture')
dev.off()

#save renamed Rds file
saveRDS(mps,file="data-output/mps_annotated.Rds")
