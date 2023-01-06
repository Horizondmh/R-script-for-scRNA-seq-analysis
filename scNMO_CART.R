## QC & AddClonotypeData----------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
gc()
library(Seurat)
library(tidyselect)
library(dplyr)
library(DoubletFinder)
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))
  tcr <- tcr[!duplicated(tcr$barcode), ]
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
  tcr <- tcr[, c(2,1,3)]
  rownames(tcr) <- tcr[,1]
  tcr[,1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}
dir1_1 = c('./data/matrix_/B_013_pre/')
dir1_2 = c('./data/tcr_bcr_matrix/B_013_pre_TCR/')
dir1_3 = c('./data/tcr_bcr_matrix/B_013_pre_BCR/')
sample_name1 <- c('B_013_pre')
scRNAlist <- list()
for(i in 1:N){
  counts <- Read10X(data.dir = dir1_1[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample_name1[i], min.cells = 3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]]$log10GenesPerUMI <- log10(scRNAlist[[i]]$nFeature_RNA)/log10(scRNAlist[[i]]$nCount_RNA)
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt < 10)
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nCount_RNA>= 1000  &
                             nCount_RNA< 40000 &
                             nFeature_RNA>= 200  &
                             log10GenesPerUMI > 0.7)
  ## preprocessing
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- ScaleData(scRNAlist[[i]], verbose = FALSE)
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], verbose = FALSE)
  scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], npcs = 30, verbose = FALSE)
  scRNAlist[[i]] <- RunUMAP(scRNAlist[[i]], reduction = "pca", dims = 1:30)
  scRNAlist[[i]] <- FindNeighbors(scRNAlist[[i]], reduction = "pca", dims = 1:30)
  scRNAlist[[i]] <- FindClusters(scRNAlist[[i]], resolution = 0.8)
}
for(i in 1:N){
  table(scRNAlist[[i]]$orig.ident)
  Doubletrate = ncol(scRNAlist[[i]])*8*1e-6
  sweep.res.list <- paramSweep_v3(scRNAlist[[i]], PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- scRNAlist[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sample14@meta.data$ClusteringResults
  nExp_poi <- round(Doubletrate*nrow(scRNAlist[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  scRNAlist[[i]] <- doubletFinder_v3(scRNAlist[[i]], PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  scRNAlist[[i]]$doubFind_res = scRNAlist[[i]]@meta.data %>% select(contains('DF.classifications'))
  scRNAlist[[i]]$doubFind_score = scRNAlist[[i]]@meta.data %>% select(contains('pANN'))
  scRNAlist[[i]] = subset(scRNAlist[[i]], doubFind_res == "Singlet")
  table(scRNAlist[[i]]$orig.ident)
  scRNAlist[[i]] <- add_clonotype(dir1_2[i], scRNAlist[[i]],"t")
  scRNAlist[[i]] <- add_clonotype(dir1_3[i], scRNAlist[[i]],"b")
  scRNAlist[[i]] <- subset(scRNAlist[[i]], cells = colnames(scRNAlist[[i]])[!(!is.na(scRNAlist[[i]]$t_clonotype_id) & !is.na(scRNAlist[[i]]$b_clonotype_id))])
}
scRNA_harmony <- merge(scRNAlist1, y = c(scRNAlist2))
table(scRNA_harmony$orig.ident)
scRNA_harmony <- NormalizeData(scRNA_harmony)
scRNA_harmony <- CellCycleScoring(object = scRNA_harmony,
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes,
                                  g1.features = cc.genes$g1.genes)
scRNA_harmony <- FindVariableFeatures(scRNA_harmony)
scRNA_harmony <- ScaleData(scRNA_harmony)
scRNA_harmony <- RunPCA(scRNA_harmony, verbose=FALSE)
library(harmony)
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
plot1 <- DimPlot(scRNA_harmony, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA_harmony, ndims=50, reduction="pca")
plotc <- plot1+plot2
plotc
pc.num = 1:30
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num)
#scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = pc.num)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = pc.num)
scRNA_harmony <- FindClusters(scRNA_harmony, reduction = "harmony", resolution = 0.8)
DimPlot(scRNA_harmony, reduction = "umap", label=T)
#DimPlot(scRNA_harmony, reduction = "tsne", label=T)
scRNA_harmony@meta.data$tcr = scRNA_harmony@meta.data$t_clonotype_id
scRNA_harmony@meta.data$tcr[!is.na(scRNA_harmony@meta.data$tcr)] <- "TCR"
scRNA_harmony@meta.data$bcr = scRNA_harmony@meta.data$b_clonotype_id
scRNA_harmony@meta.data$bcr[!is.na(scRNA_harmony@meta.data$bcr)] <- "BCR"
DimPlot(scRNA_harmony, reduction = "umap", label=T, group.by = "bcr")
saveRDS(scRNA_harmony, "scRNA_HC_PRE_3M.rds")


## Cell_clustering----------------------------------------------------------------------------------------------------------------------------------------------------
pdf(file="markerBubble_1.pdf",width=30,height=11)
cluster10Marker=c("CD3D","CD3E","CD4","CD8A","MKI67",
                  "CCR7","TCF7","LEF1","SELL","CCR6","NR4A1","NCR3","KLRB1","GPR183","IL7R","CD27",
                  "PRF1","GZMA","GZMB","GZMK","NKG7","FCGR3A","KLRG1",
                  "ICOS","PDCD1","LAG3","HAVCR2","CD200","CTLA4","ENTPD1","ITGAE","EOMES","IRF8",
                  "FOXP3","IL2RA","TRGV9","TRDV2",
                  "NCAM1","KLRF1","KLRC1","CD160")
DotPlot(object = scRNA_all, features = cluster10Marker)
dev.off()
pdf(file="markerBubble_2.pdf",width=15,height=11)
cluster10Marker=c("MS4A1","IGHD","NEIL1","CD27","CD38","ITGAX","TBX21","XBP1","PRDM1","IGHM","COCH","MKI67","IGHA1","IGHG1","MZB1")
DotPlot(object = scRNA_all, features = cluster10Marker,col.min = -1)
dev.off()
pdf(file="markerBubble_3.pdf",width=15,height=11)
cluster10Marker=c("LYZ","CD14","CD83","FCGR3A","C1QA","C1QB","C1QC","CSF1R","TREM2","APOE","TYROBP","CX3CR1","SLC2A5","P2RY13","P2RY12",
                  "CD1C","LILRA4","PF4","CD68","MARCO","FUT4","ITGAM","MME","CXCR2","SELL")
DotPlot(object = scRNA_all, features = cluster10Marker,col.min = -1)
dev.off()
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)
new.cluster.ids <- c("CD8+ T cells",
                     "CD4+ T cells",
                     "Monocytes",
                     "NK cells",
                     "CD4+ T cells",
                     "CD8+ T cells",
                     "B cells",
                     "CD8+ T cells",
                     "B cells",
                     "Monocytes",
                     "Monocytes",
                     "Monocytes",
                     "Monocytes",
                     "CD8+ T cells",
                     "NK cells",
                     "Proliferating cells",
                     "CD4+ T cells",
                     "Plasma cells",
                     "Monocytes",
                     "Monocytes",
                     "Monocytes",
                     "CD4+ T cells",
                     "B cells",
                     "Monocytes",
                     "Monocytes",
                     "Monocytes",
                     "Unknown",
                     "Unknown",
                     "Monocytes",
                     "CD8+ T cells")
scRNA_all@active.ident <- plyr::mapvalues(
  x = scRNA_all@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_all@meta.data$label = scRNA_all@active.ident
table(scRNA_all$orig.ident)
scRNA_all = subset(scRNA_all, label != "Unknown")
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T)
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T,group.by = "tcr")
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T,group.by = "bcr")
scRNA_all@active.ident = factor(scRNA_all@active.ident,
                                levels = c("CD4+ T cells","CD8+ T cells","NK cells","Proliferating cells",
                                           "B cells","Plasma cells","Monocytes"))
scRNA_all$label = scRNA_all@active.ident
scRNA_HC = subset(scRNA_all, SampleType == "HC")
scRNA_HC_PBMC = subset(scRNA_HC, SampleOrig == "blood")
scRNA_HC_CSF = subset(scRNA_HC, SampleOrig == "csf")
table(scRNA_HC_PBMC$label)
table(scRNA_HC_CSF$label)
# Fisher's exact test
data = read.csv("Fisher_exact_test.csv",header = T,row.names = 1)
fisher.test(data,alternative = "two.sided")$p.value
data = read.csv("BH_correction.csv",header = T,row.names = 1)
data$BH =p.adjust(data$Raw.p,method = "BH")
data$BH


## B_cell_reclustering------------------------------------------------------------------------------------------------------------------------------------------------
scRNA_B_P = subset(scRNA_all, label == "B cells"|label == "Plasma cells"|label == "Proliferating cells")
scRNA_B_P = subset(scRNA_B_P, bcr == "BCR")
pc.num = 1:30
scRNA_B_P <- RunUMAP(scRNA_B_P, reduction = "harmony", dims = pc.num)
scRNA_B_P <- FindNeighbors(scRNA_B_P, reduction = "harmony", dims = 1:30)
scRNA_B_P <- FindClusters(scRNA_B_P, reduction = "harmony", resolution = 1.5)
DimPlot(scRNA_B_P, label = T,reduction = "umap")
DimPlot(scRNA_B_P, label = T,reduction = "umap")
dir1_3 = c('./data/tcr_bcr_matrix_2/B_HC_2_BCR/')
bcr36 <- read.csv(paste(dir1_3[36],"filtered_contig_annotations.csv", sep=""))
bcr36 = bcr36[bcr36$chain == "IGH",]
bcr36 <- bcr36[!duplicated(bcr36$barcode), ]
bcr36 <- bcr36[,c("barcode", "c_gene")]
colnames(bcr36) = c("cell","c_gene")
bcr36$cell = paste(bcr36$cell,"_19_2", sep="")
bcr = rbind(bcr1,bcr2,bcr3,bcr4,bcr5,bcr6,bcr7,bcr8,bcr9,bcr10,bcr11,bcr12,
            bcr13,bcr14,bcr15,bcr16,bcr18,bcr19,bcr20,bcr21,bcr22,bcr23,bcr24,
            bcr25,bcr26,bcr27,bcr29,bcr31,bcr32,bcr34,bcr35,bcr36)
saveRDS(bcr,"bcr_isotype.rds")
bcr = readRDS("bcr_isotype.rds")
scRNA_B_P@meta.data$id  <- 1:nrow(scRNA_B_P@meta.data)
scRNA_B_P$cell = rownames(scRNA_B_P@meta.data)
scRNA_B_P@meta.data = merge(scRNA_B_P@meta.data, bcr, by = "cell",all.x = T)
scRNA_B_P@meta.data = scRNA_B_P@meta.data[order(scRNA_B_P@meta.data$id), ]
rownames(scRNA_B_P@meta.data) = scRNA_B_P$cell
sub_014 = subset(scRNA_B_P, orig.ident == "B_014_pre")
meta = sub_014@meta.data
colnames(meta)
meta = meta[,c(15,93,94)]
scRNA_B_P$c_gene[scRNA_B_P$c_gene == ""]<-"NA"
scRNA_B_P$Isotype = scRNA_B_P$c_gene
table(scRNA_B_P$Isotype)
scRNA_B_P$Isotype[scRNA_B_P$Isotype == "IGHE"]<-"NA"
scRNA_B_P$Isotype = factor(scRNA_B_P$Isotype, levels = c("IGHM","IGHD","IGHA1","IGHA2",
                                                         "IGHG1","IGHG2","IGHG3","IGHG4","NA"))
col6 = c("#98C9DD","#207CB5","#A6D38E","#37A849","#F69595","#EB2A2A","#FCBA71","#F58229","#C4C2C2")
library(ggplot2)
DimPlot(scRNA_B_P, reduction = "umap", cols = col6, group.by = "Isotype",pt.size = 0.2)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank())
pdf(file="markerBubble_1.pdf",width=30,height=11)
cluster10Marker=c("CD3D","CD3E","CD4","CD8A","MKI67",
                  "CCR7","TCF7","LEF1","SELL","CCR6","NR4A1","NCR3","KLRB1","GPR183","IL7R","CD27",
                  "PRF1","GZMA","GZMB","GZMK","NKG7","FCGR3A","KLRG1",
                  "ICOS","PDCD1","LAG3","HAVCR2","CD200","CTLA4","ENTPD1","ITGAE","EOMES","IRF8",
                  "FOXP3","IL2RA","TRGV9","TRDV2",
                  "NCAM1","KLRF1","KLRC1","CD160")
DotPlot(object = scRNA_B_P, features = cluster10Marker)
dev.off()
pdf(file="markerBubble_2.pdf",width=15,height=11)
cluster10Marker=c("MS4A1","IGHD","NEIL1","CD27","CD38","ITGAX","TBX21","XBP1","PRDM1","IGHM","COCH","MKI67","IGHA1","IGHG1","MZB1")
DotPlot(object = scRNA_B_P, features = cluster10Marker,col.min = -1)
dev.off()
pdf(file="markerBubble_3.pdf",width=15,height=11)
cluster10Marker=c("LYZ","CD14","CD83","FCGR3A","C1QA","C1QB","C1QC","CSF1R","TREM2","APOE","TYROBP","CX3CR1","SLC2A5","P2RY13","P2RY12",
                  "CD1C","LILRA4","PF4","CD68","MARCO","FUT4","ITGAM","MME","CXCR2","SELL")
DotPlot(object = scRNA_B_P, features = cluster10Marker,col.min = -1)
dev.off()
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
new.cluster.ids <- c("Naive B",
                     "Naive B",
                     "Immature B",
                     "Naive B",
                     "Naive B",
                     "Switched memory B",
                     "Non-switched memory B",
                     "Switched memory B",
                     "Naive B",
                     "Switched memory B",
                     "Naive B",
                     "Immature B",
                     "Non-switched memory B",
                     "Plasma cells/Plasmablasts",
                     "Immature B",
                     "Double negative B",
                     "Switched memory B",
                     "Double negative B",
                     "Non-switched memory B",
                     "Naive B",
                     "Non-switched memory B",
                     "Plasma cells/Plasmablasts",
                     "Double negative B",
                     "Plasma cells/Plasmablasts",
                     "Unknown",
                     "Unknown")
scRNA_B_P@active.ident <- plyr::mapvalues(
  x = scRNA_B_P@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_B_P@meta.data$label = scRNA_B_P@active.ident
scRNA_B_P = subset(scRNA_B_P, label != "Unknown")
# DimPlot(scRNA_B_P, reduction = "umap",label=T, repel = T)
DimPlot(scRNA_B_P, reduction = "umap",label=T, repel = T)
scRNA_B_P@active.ident = factor(scRNA_B_P@active.ident,
                                levels = c("Immature B","Naive B","Non-switched memory B",
                                           "Switched memory B","Plasma cells/Plasmablasts","Double negative B"))
scRNA_B_P$label = scRNA_B_P@active.ident
scRNA_B_P_HC_BL = subset(scRNA_B_P, SampleType != "3M")
col = c("#999899","#91C463","#5586AD","#9BC3DF","#FCB863","#F2D2D5")
library(ggplot2)
DimPlot(scRNA_B_P_HC_BL, reduction = "umap",repel = T, group.by = "label", cols = col,pt.size = 0.8)+
  NoLegend()+labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(size=15))+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
library(viridis)
library(ggplot2)
pdf(file="markerBubble_Bcells.pdf",width=10,height=4.5)
cluster10Marker=c("MS4A1","IGHD","CD19","IGHM","CD27","COCH","XBP1","PRDM1","IGHA1","IGHG1","MZB1","ITGAX","TBX21")
DotPlot(object = scRNA_B_P_HC_BL, features = cluster10Marker)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()
dev.off()
library(stringr)
scRNA_B_P_HC_BL$Isotype = str_replace_na(scRNA_B_P_HC_BL$Isotype)
scRNA_B_P_HC_BL$Isotype = factor(scRNA_B_P_HC_BL$Isotype,levels = c(
  "IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","NA"
))
col6 = c("#98C9DD","#207CB5","#A6D38E","#37A849","#F69595","#EB2A2A","#FCBA71","#F58229","#C4C2C2")
library(ggplot2)
DimPlot(scRNA_B_P_HC_BL, reduction = "umap", cols = col6, group.by = "Isotype",pt.size = 0.6)+
  NoLegend()+labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(size=15))+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
scRNA_B_P_HC = subset(scRNA_B_P, SampleType == "BL")
sub_immature = subset(scRNA_B_P_HC, label == "Immature B")
sub_naive = subset(scRNA_B_P_HC, label == "Naive B")
sub_nonswitched = subset(scRNA_B_P_HC, label == "Non-switched memory B")
sub_switched = subset(scRNA_B_P_HC, label == "Switched memory B")
sub_plasmablast = subset(scRNA_B_P_HC, label == "Plasma cells/Plasmablasts")
sub_dn = subset(scRNA_B_P_HC, label == "Double negative B")
table(sub_dn$Isotype)
clonotype = scRNA_B_P_HC_BL$b_clonotype_id
clonotype_df = data.frame(cell = names(clonotype),type = clonotype,orig.ident = scRNA_B_P_HC_BL$orig.ident)
clonotype_df = na.omit(clonotype_df)
for(i in 1:nrow(clonotype_df)){
  clonotype_df$type[i] = paste(clonotype_df$orig.ident[i],clonotype_df$type[i])
}
celltype = scRNA_B_P_HC_BL$label
celltype_df = data.frame(cell = names(celltype),celltype = celltype)
celltype_clonotype_df = merge(clonotype_df,celltype_df,all.y = T)
clononumber_df = data.frame(type = names(table(celltype_clonotype_df$type)),clononumber = table(celltype_clonotype_df$type))
clononumber_df = clononumber_df[,-2]
colnames(clononumber_df) = c("type","clononumber")
celltype_clonotype_clononumber_df = inner_join(celltype_clonotype_df,clononumber_df)
#celltype_clonotype_clononumber_df = merge(celltype_clonotype_df,clononumber_df,by = "type",all.x = T)
celltype_clonotype_clononumber_df$clonostate = ifelse(celltype_clonotype_clononumber_df$clononumber >= 5, ">=5","<5")
scRNA_B_P_HC_BL@meta.data$cell = rownames(scRNA_B_P_HC_BL@meta.data)
clononumber_df = celltype_clonotype_clononumber_df[,c(1,5)]
clononumber_df[,2] = as.numeric(clononumber_df[,2])
scRNA_B_P_HC_BL@meta.data = inner_join(scRNA_B_P_HC_BL@meta.data,clononumber_df)
#scRNA_B_P_HC_BL@meta.data = merge(scRNA_B_P_HC_BL@meta.data,clononumber_df,by = "cell",all.x = T)
clonostate_df = celltype_clonotype_clononumber_df[,c(1,6)]
scRNA_B_P_HC_BL@meta.data = inner_join(scRNA_B_P_HC_BL@meta.data,clonostate_df)
#scRNA_B_P_HC_BL@meta.data = merge(scRNA_B_P_HC_BL@meta.data,clonostate_df,by = "cell",all.x = T)
rownames(scRNA_B_P_HC_BL@meta.data) = scRNA_B_P_HC_BL@meta.data$cell
library(ggplot2)
DimPlot(scRNA_B_P_HC_BL, reduction = "umap", group.by = "clonostate",cols = c("#D1D0E4","#463179"),pt.size = 0.8)+
  NoLegend()+labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(size=15))+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
scRNA_B_P_BL = subset(scRNA_B_P_HC_BL, SampleType == "BL")
scRNA_P_BL = subset(scRNA_B_P_BL, label == "Plasma cells/Plasmablasts")
VlnPlot(scRNA_P_BL, "TNFRSF17", group.by = "clonostate")
VlnPlot(scRNA_P_BL, "TNFRSF17", group.by = "clononumber")
p = VlnPlot(scRNA_P_BL, "TNFRSF17", group.by = "clonostate")
data = p$data
scRNA_P_BL = NormalizeData(scRNA_P_BL)
all.genes = rownames(scRNA_P_BL)
scRNA_P_BL = ScaleData(scRNA_P_BL, features = all.genes)
TNFRSF17 = scRNA_P_BL@assays$RNA@scale.data[rownames(scRNA_P_BL@assays$RNA@scale.data) == "TNFRSF17",]
TNFRSF17_df = data.frame(cell = names(TNFRSF17),TNFRSF17_expression = TNFRSF17)
meta = scRNA_P_BL@meta.data
meta$cell = rownames(meta)
meta2 = inner_join(meta, TNFRSF17_df)
colnames(meta2)
meta3 = meta2[,c(97,98,99)]
write.csv(meta3, "ppb_BCMA_clononumber.csv")
library(immunarch)
dir1_3 = c('./data/tcr_bcr_matrix_2/B_015_pre_BCR/')
mydata <- read.csv(paste(dir1_3[24],"clonotypes.csv", sep=""))
colnames(mydata) = c("clonotype_id","Clones","proportion","CDR3.aa","CDR3.nt")
div_gini <- repDiversity(mydata, "gini")
div_gini
mydata2 = mydata[mydata$Clones!=1,]
div_gini2<- repDiversity(mydata2, "gini")
div_gini2
library(ggplot2)
library(ggExtra)
data = read.csv("gini.csv",row.names = 1)
piris <- ggplot(data, aes(all_gini, poly_gini)) +
  geom_point()+
  xlim(-0.2,1)+
  ylim(-0.2,1)
ggMarginal(piris, groupFill = F)
library(COSG)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)
library(COSG)
Idents(scRNA_B_P_HC_BL) = "label"
table(scRNA_B_P_HC_BL$label)
sample = subset(scRNA_B_P_HC_BL,downsample = 654)
marker_cosg <- cosg(
  sample,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=20)
marker_cosg$names
sample@active.ident = factor(sample@active.ident,
                             levels = c("Immature B","Naive B","Non-switched memory B",
                                        "Switched memory B","Plasma cells/Plasmablasts","Double negative B"))
sample$label = sample@active.ident
mat <- GetAssayData(sample, slot = "counts")
cluster_info <- sort(sample$label)
gene = c(marker_cosg$names[,1],marker_cosg$names[,2],marker_cosg$names[,3],
         marker_cosg$names[,4],marker_cosg$names[,5],marker_cosg$names[,6])
mat <- as.matrix(mat[gene, names(cluster_info)])
mat = log2(mat + 1)
mark_gene <- c("TCL1A","FCER2","IGHD","IGHM","TCL1A","IGHM","CCR7","FCER2","IGHD",
               "TNFRSF13B","COCH","AIM2","TNFRSF13B","XBP1","MZB1","JCHAIN","TNFRSF17","PRDM1","NXPH4")
gene_pos <- which(rownames(mat) %in% mark_gene)

row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = mark_gene))

SampleType_info = sample$SampleType
SampleType_df = data.frame(cell = names(SampleType_info),SampleType = SampleType_info)
SampleOrig_info = sample$SampleOrig
SampleOrig_df = data.frame(cell = names(SampleOrig_info),SampleOrig = SampleOrig_info)
cluster_df = data.frame(cell = names(cluster_info),cluster = cluster_info)
cycle_info = sample$Phase
cycle_df = data.frame(cell = names(cycle_info),cycle = cycle_info)
col_df = inner_join(cluster_df, SampleType_df)
col_df = inner_join(col_df,SampleOrig_df)
col_df = inner_join(col_df,cycle_df)
rownames(col_df) = col_df$cell
col_df = col_df[,-1]
table(col_df$cluster)
table(col_df$SampleType)
table(col_df$SampleOrig)
table(col_df$cycle)
#col = c("#999899","#91C463","#5586AD","#9BC3DF","#FCB863","#F2D2D5")
colours <- list(
  "cluster"=c("Immature B"="#999899",
              "Naive B"="#91C463",
              "Non-switched memory B"="#5586AD",
              "Switched memory B"="#9BC3DF",
              "Plasma cells/Plasmablasts"="#FCB863",
              "Double negative B"="#F2D2D5"),
  "SampleType"=c("BL"="#52679F","HC"="#B2C3E0"),
  "SampleOrig"=c("blood"="#CA7D6E","csf"="#99C5DF"),
  "cycle"=c("G1"="#6131A0","G2M"="#D76AFF","S"="#F2D1A8"))
# colours <- list(
#   "cluster"=c("blood Naive B cells"="#F1B9B8","csf Naive B cells"="#99D3E1",
#               "blood Non-switched memory B cells"="#CEB5AC","csf Non-switched memory B cells"="#D9D8D8",
#               "blood Switched memory B cells"="#D2F1A2","csf Switched memory B cells"="#EDE9AF",
#               "blood Plasma cells"="#CFD9F5","csf Plasma cells"="#D1CE60",
#               "blood Plasmablasts"="#72B95A","csf Plasmablasts"="#E4984C"),
#   "type"=c("B_012_pre"="#52679F","B_013_pre"="#52679F","B_014_pre"="#52679F","B_015_pre"="#52679F",
#            "B_016_pre"="#52679F","B_017_pre"="#52679F","C_012_pre"="#52679F","C_013_pre"="#52679F",
#            "C_014_pre"="#52679F","C_015_pre"="#52679F","C_016_pre"="#52679F","C_017_pre"="#52679F",
#            "B_HC_2"="#B2C3E0","B_HC_3"="#B2C3E0","C_HC_2"="#B2C3E0","C_HC_3"="#B2C3E0"),
#   "cycle"=c("G1"="#6131A0","G2M"="#E79D76","S"="#F2D1A8"))
col_anno = HeatmapAnnotation(df=col_df, 
                             which="col",
                             col=colours, 
                             annotation_width=unit(c(1, 4), "cm"), 
                             gap=unit(1, "mm"))
col_heatmap = colorRamp2(c(-0.5,0.2,2),c("#37538B","white","#7C262C"))
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        right_annotation = row_anno,
        top_annotation = col_anno,
        #column_split = cluster_info,
        col = col_heatmap,
        column_title = NULL,
        use_raster = F)
scRNA_B_P@active.ident = factor(scRNA_B_P@active.ident,
                                levels = c("Immature B cells","Naive B cells","Non-switched memory B cells","Switched memory B cells",
                                           "Plasmablasts","Plasma cells","Double negative B cells"))
scRNA_B_P$label = scRNA_B_P@active.ident
library(dplyr)
colnames(scRNA_B_P@meta.data)
sub_1 = subset(scRNA_B_P,orig.ident == "B_017_3m")
b = sub_1@meta.data[,c(15,89)]
b = na.omit(b)
duplicate_name = b %>% group_by(b_clonotype_id) %>% summarise(freq = n()) %>% filter(freq > 1) %>% select(b_clonotype_id)
duplicate_data = b[b$b_clonotype_id %in% duplicate_name$b_clonotype_id, ]
non_duplicate_data = b[!b$b_clonotype_id %in% duplicate_name$b_clonotype_id, ]
table(non_duplicate_data$label)
table(duplicate_data$label)
# Fisher's exact test
data = read.csv("Fisher_exact_test.csv",header = T,row.names = 1)
fisher.test(data,alternative = "two.sided")$p.value
data = read.csv("BH_correction.csv",header = T,row.names = 1)
data$BH =p.adjust(data$Raw.p,method = "BH")
data$BH
FeaturePlot(scRNA_B_P_HC_BL, features = "TNFRSF17",pt.size = 0.4)
FeaturePlot(scRNA_B_P_HC_BL, features = "IGHG1",pt.size = 0.4)
scRNA_B_P_HC_BL$label = scRNA_B_P_HC_BL@active.ident
#col = c("#DAB8DA","#9BC3DF","#789CBF","#B9B9B9","#7EB793","#F3754E","#FCB863")
col = c("#999899","#91C463","#5586AD","#9BC3DF","#FCB863","#F2D2D5")
VlnPlot(scRNA_B_P_HC_BL, features = "TNFRSF17", group.by = "label",pt.size = 0,cols = col)
p = VlnPlot(scRNA_B_P_HC_BL, features = "TNFRSF17", group.by = "label",pt.size = 0,cols = col)
data = p$data
scRNA_B_P_HC_BL_plasmablast = subset(scRNA_B_P_HC_BL, label == "Plasma cells/Plasmablasts")
VlnPlot(scRNA_B_P_HC_BL_plasmablast, features = "TNFRSF17", group.by = "SampleType",pt.size = 1)
scRNA_B_P_HC_BL_plasmablast = subset(scRNA_B_P_HC_BL, label == "Plasma cells/Plasmablasts")
col = c("#98C9DD","#207CB5","#A6D38E","#37A849","#F69595","#EB2A2A","#FCBA71","#F58229","#C4C2C2")
VlnPlot(scRNA_B_P_HC_BL_plasmablast, features = "TNFRSF17", group.by = "Isotype",pt.size = 1,cols = col)
p = VlnPlot(scRNA_B_P_HC_BL_plasmablast, features = "TNFRSF17", group.by = "Isotype",pt.size = 0,cols = col)
data = p$data
scRNA_B_P_BL = subset(scRNA_B_P, SampleType == "BL")
scRNA_B_P_BL_plasmablast = subset(scRNA_B_P_BL, label == "Plasma cells/Plasmablasts")
col = c("#98C9DD","#207CB5","#A6D38E","#37A849","#F69595","#EB2A2A","#FCBA71","#F58229","#C4C2C2")
VlnPlot(scRNA_B_P_BL_plasmablast, features = "TNFRSF17", group.by = "Isotype",pt.size = 1,cols = col)
p = VlnPlot(scRNA_B_P_BL_plasmablast, features = "TNFRSF17", group.by = "Isotype",pt.size = 0,cols = col)
data = p$data
library(tidyverse)
library(gghalves)
# data = gather(data)
# data = na.omit(data)
colnames(data) = c("TNFRSF17","Label")
variable <- c("Immature B","Naive B","Non-switched memory B","Switched memory B",
              "Plasma cells/Plasmablasts","Double negative B")
my_sort <-factor(variable,levels = variable)
ggplot(data)+
  geom_half_violin(aes(as.numeric(factor(Label,levels = my_sort))+0.1,
                       TNFRSF17,fill=factor(Label,levels = my_sort)),
                   side = 'r',cex=0.8)+ 
  geom_boxplot(aes(as.numeric(factor(Label,levels = my_sort))+0.1,
                   TNFRSF17,fill=factor(Label,levels = my_sort)),
               outlier.colour="black",width=0.1,cex=0.8)+
  geom_jitter(aes(as.numeric(factor(Label,levels = my_sort))-0.2,
                  TNFRSF17,color=factor(Label,levels = my_sort)), 
              width = 0.1,size=0.5) +
  scale_fill_manual(values = c("#999899","#91C463","#5586AD","#9BC3DF","#FCB863","#F2D2D5"))+
  scale_color_manual(values = c("#999899","#91C463","#5586AD","#9BC3DF","#FCB863","#F2D2D5"),guide='none')+
  scale_x_continuous(breaks = unique(as.numeric(factor(data$Label,levels = my_sort))), 
                     labels = unique(factor(data$Label,levels = my_sort)))+
  theme_test(base_size = 15)+
  labs(x=NULL,y='GSVA score')+
  scale_y_continuous(limits = c(0,3.5))+
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(face = 'italic'))
sub_ppb = subset(scRNA_B_P, label == "Plasma cells/Plasmablasts")
sub_ppb_012 = subset(sub_ppb, orig.ident == "B_014_pre")
meta = sub_ppb_012@meta.data
colnames(meta)
meta2 = meta[,c(1,2,15,95)]
B_012_pre = read.csv("consensus_annotations_B_014_pre.csv")
B_012_pre = B_012_pre[B_012_pre$chain == "IGH",]
B_012_3m = read.csv("consensus_annotations_B_014_3m.csv")
B_012_3m = B_012_3m[B_012_3m$chain == "IGH",]
merge = full_join(B_012_pre, B_012_3m,by="cdr3_nt_length")
merge$j_gene_iden = ifelse(merge$j_gene.x == merge$j_gene.y, "identical", "not-identical")
merge = subset(merge, j_gene_iden == 'identical')
for (i in 1:nrow(merge)) {
  a = merge$cdr3_nt.x[i]
  a1 = strsplit(a,split = "")
  a2 = unlist(a1)
  b = merge$cdr3_nt.y[i]
  b1 = strsplit(b,split = "")
  b2 = unlist(b1)
  merge$overlap[i] = sum(a2 == b2)/nchar(a)
}
write.csv(merge, "overlap_for_014.csv")
sub_012_pre = subset(scRNA_B_P,orig.ident == "C_017_pre")
sub_012_pre = subset(sub_012_pre, bcr == "BCR")
sub_012_pre$IGH = sub_012_pre$b_cdr3s_aa
sub_012_pre$cell = rownames(sub_012_pre@meta.data)
colnames(sub_012_pre@meta.data)
clonotype = sub_012_pre$b_clonotype_id
clonotype_df = data.frame(cell = names(clonotype),type = clonotype)
clononumber_df = data.frame(type = names(table(clonotype_df$type)),clononumber = table(clonotype_df$type))
clononumber_df = clononumber_df[,-2]
colnames(clononumber_df) = c("type","clononumber")
clononumber_df$clonofreq = clononumber_df$clononumber/nrow(clonotype_df)
clonotype_clononumber_df = merge(clonotype_df,clononumber_df,by = "type",all.x = T)
sub_012_pre@meta.data = inner_join(sub_012_pre@meta.data,clonotype_clononumber_df)
meta_012_pre = sub_012_pre@meta.data
colnames(meta_012_pre)
meta_012_pre = meta_012_pre[,c(1,15,16,89,94,97,99,100)]
for(i in 1:nrow(meta_012_pre)){
  meta_012_pre$IGH[i] = strsplit(meta_012_pre$b_cdr3s_aa[i], ";", fixed= T)[[1]][1]
}
head(meta_012_pre)
sub_012_3m = subset(scRNA_B_P,orig.ident == "C_017_3m")
sub_012_3m = subset(sub_012_3m, bcr == "BCR")
sub_012_3m$IGH = sub_012_3m$b_cdr3s_aa
sub_012_3m$cell = rownames(sub_012_3m@meta.data)
colnames(sub_012_3m@meta.data)
clonotype = sub_012_3m$b_clonotype_id
clonotype_df = data.frame(cell = names(clonotype),type = clonotype)
clononumber_df = data.frame(type = names(table(clonotype_df$type)),clononumber = table(clonotype_df$type))
clononumber_df = clononumber_df[,-2]
colnames(clononumber_df) = c("type","clononumber")
clononumber_df$clonofreq = clononumber_df$clononumber/nrow(clonotype_df)
clonotype_clononumber_df = merge(clonotype_df,clononumber_df,by = "type",all.x = T)
sub_012_3m@meta.data = inner_join(sub_012_3m@meta.data,clonotype_clononumber_df)
meta_012_3m = sub_012_3m@meta.data
colnames(meta_012_3m)
meta_012_3m = meta_012_3m[,c(1,15,16,89,94,97,99,100)]
for(i in 1:nrow(meta_012_3m)){
  meta_012_3m$IGH[i] = strsplit(meta_012_3m$b_cdr3s_aa[i], ";", fixed= T)[[1]][1]
}
head(meta_012_3m)
IGH_replication = intersect(meta_012_pre$IGH,meta_012_3m$IGH)
# meta_012_pre$REP = ifelse(meta_012_pre$IGH %in% IGH_replication, "rep", "no_rep")
# meta_012_pre_rep = meta_012_pre[meta_012_pre$REP == "rep",]
# meta_012_pre_norep = meta_012_pre[meta_012_pre$REP == "no_rep",]
# meta_012_3m$REP = ifelse(meta_012_3m$IGH %in% IGH_replication, "rep", "no_rep")
# meta_012_3m_rep = meta_012_3m[meta_012_3m$REP == "rep",]
# meta_012_3m_norep = meta_012_3m[meta_012_3m$REP == "no_rep",]
# write.csv(meta_012_3m_norep,"meta_014_3m_norep.csv")
write.csv(meta_012_pre,"meta_017_pre.csv")
write.csv(meta_012_3m,"meta_017_3m.csv")
# Fisher's exact test
data = read.csv("Fisher_exact_test.csv",header = T,row.names = 1)
fisher.test(data,alternative = "two.sided")$p.value
library(patchwork)
scRNA_B_P_HC_BL$SampleType = factor(scRNA_B_P_HC_BL$SampleType, levels = c("BL","HC"))
scRNA_P = subset(scRNA_B_P_HC_BL, label == "Plasma cells/Plasmablasts")
sub_IgGplasma = subset(scRNA_P,Isotype == "IGHG1"|Isotype == "IGHG2"|Isotype == "IGHG3"|Isotype == "IGHG4")
sub_IgAplasma = subset(scRNA_P,Isotype == "IGHA1"|Isotype == "IGHA2")
sub_IgMplasma = subset(scRNA_P,Isotype == "IGHM")
Marker1=read.csv("immunoglobulin_production.csv")
Marker2=read.csv("immunoglobulin_mediated_immune_response.csv")
Marker3=read.csv("B cell activation.csv")
Marker4=read.csv("B cell receptor signaling pathway.csv")
Marker = c(Marker1,Marker2,Marker3,Marker4)
# cols = c('#4676B4', '#FDCF58',"D83429")
library(ggplot2)
p1 = DotPlot(object = sub_IgGplasma, features = Marker,group.by = "SampleType",cols = c('#4676B4', '#FDCF58',"D83429"),scale.min = 0,scale.max = 40)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p2 = DotPlot(object = sub_IgAplasma, features = Marker,group.by = "SampleType",cols = c('#4676B4', '#FDCF58',"D83429"),scale.min = 0,scale.max = 40)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p3 = DotPlot(object = sub_IgMplasma, features = Marker,group.by = "SampleType",cols = c('#4676B4', '#FDCF58',"D83429"),scale.min = 0,scale.max = 40)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p1/p2/p3
table(scRNA_B_P$orig.ident)
scRNA_P = subset(scRNA_B_P, label == "Plasma cells/Plasmablasts")
table(scRNA_P$orig.ident)
SampleOrig = scRNA_B_P$SampleOrig
SampleOrig_df = data.frame(cell = names(SampleOrig),type = SampleOrig,label = scRNA_B_P$label)
for(i in 1:nrow(SampleOrig_df)){
  SampleOrig_df$SampleClass[i] = paste(SampleOrig_df$type[i],SampleOrig_df$label[i])
}
scRNA_B_P@meta.data$cell = rownames(scRNA_B_P@meta.data)
SampleOrig_df = SampleOrig_df[,c(1,4)]
scRNA_B_P@meta.data = inner_join(scRNA_B_P@meta.data,SampleOrig_df)
rownames(scRNA_B_P@meta.data) = scRNA_B_P@meta.data$cell
DimPlot(scRNA_B_P, group.by = "SampleClass")


## Monocle------------------------------------------------------------------------------------------------------------------------------------------------------------
library(monocle)
library(tidyverse)
library(dplyr)
library(patchwork)
Idents(scRNA_B_P) = scRNA_B_P$SampleClass
table(Idents(scRNA_B_P))
levels(scRNA_B_P)
data <- as(as.matrix(scRNA_B_P@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_B_P@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#build new cell data set
test <- newCellDataSet(data,phenoData =pd,featureData =fd,expressionFamily=negbinomial.size(),lowerDetectionLimit=0.5)
test <- estimateSizeFactors(test)
test <- estimateDispersions(test)
test <- detectGenes(test, min_expr = 0.1)
print(head(fData(test)))
expressed_genes = row.names(subset(fData(test),num_cells_expressed >= 10))
Sys.time()
diff <- differentialGeneTest(test[expressed_genes,],
                             fullModelFormulaStr = "~label")
Sys.time()
deg = subset(diff, qval < 0.01)
deg = deg[order(deg$qval, decreasing = F),]
#write.table(deg, "test.monocle.DEG.xls", col.names = T, row.names = F, sep = "\t", quote = F)
ordergene = rownames(deg)
test = setOrderingFilter(test, ordergene)
#ordergene = row.names(deg)[order(deg$qval)][1:2000]
test=reduceDimension(test,method = "DDRTree",max_components = 2) 
test=orderCells(test)
plot_cell_trajectory(test,color_by = "label",cell_size = 1.2)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  scale_color_manual(breaks = c("Immature B cells","Naive B cells","Non-switched memory B cells","Switched memory B cells",
                                "Plasma cells","Plasmablasts"), 
                     values=c("#F1B9B8","#99D3E1","#CEB5AC","#D9D8D8","#D2F1A2","#EDE9AF","#CFD9F5"))
plot_cell_trajectory(test,color_by = "State")+ theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
plot_cell_trajectory(test,color_by = "Pseudotime")+ theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ facet_wrap("~SampleType", nrow = 1)
library(viridis)
plot(plot_cell_trajectory(test, show_cell_names = F, color_by = "Pseudotime") + scale_color_viridis_c())+ theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
plot_cell_trajectory(test,color_by = "SampleClass",cell_size = 1)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  scale_color_manual(breaks = c("blood Immature B cells","csf Immature B cells",
                                "blood Naive B cells","csf Naive B cells",
                                "blood Non-switched memory B cells","csf Non-switched memory B cells",
                                "blood Switched memory B cells","csf Switched memory B cells",
                                "blood Plasma cells","csf Plasma cells",
                                "blood Plasmablasts","csf Plasmablasts"), 
                     values=c("#F5F5DC","#F0F8FF","#FFF0F5","#B0C4DE","#FFC0CB","#87CEFA","#FF69B4","#00BFFF","#DC143C","#0000FF","#8B008B","#000080"))
plot_cell_trajectory(test,color_by = "clonostate")+ theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
plot_cell_trajectory(test,color_by = "SampleType")+ theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
plot_cell_trajectory(test,color_by = "label")+facet_wrap(~label,nrow=1)
plot_cell_trajectory(test, color_by = "label",cell_size = 1) +
  facet_wrap(~State, nrow = 1)
plotdf = pData(test_NMO)
plotdf = subset(plotdf,SampleType == "NMO")
library(ggridges)
library(RColorBrewer)
library(scales)
ggplot(plotdf, aes(x=Pseudotime,y=CellClass,fill=CellClass))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+ 
  scale_fill_manual(values=c("#FFF0F5","#FFC0CB","#DC143C","#8B008B","#B0C4DE","#87CEFA","#0000FF","#000080"))
ggsave("ridgeplot_NMO.pdf",width = 23,height = 7,units = "cm")


## Pathway_enrichment_analysis----------------------------------------------------------------------------------------------------------------------------------------
scRNA_B_P_BL = subset(scRNA_B_P_HC_BL, SampleType == "BL")
scRNA_B_BL = subset(scRNA_B_P_BL, SampleType != "Plasma cells/Plasmablasts")
Idents(scRNA_B_BL) = "SampleOrig"
sce.markers <- FindAllMarkers(object = scRNA_B_BL,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
write.csv(sce.markers, "sce.markers_Bcells_pbmc_csf.csv")
library(GSVA)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(msigdbr)
human <- msigdbr(species = "Homo sapiens")
table(human$gs_subcat)
human_GO_bp = msigdbr(species = "Homo sapiens",
                      category = "C5",
                      subcategory = "GO:BP") %>%
  dplyr::select(gs_name,gene_symbol)
human_GO_bp_Set = human_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)
scRNA_B_P_BL = subset(scRNA_B_P, SampleType == "BL")
scRNA_PPB_BL = subset(scRNA_B_P_BL, label != "Plasma cells/Plasmablasts")
DefaultAssay(scRNA_PPB_BL) = "RNA"
scRNA_PPB_BL = NormalizeData(scRNA_PPB_BL)
Idents(scRNA_PPB_BL) = "orig.ident"
expr = AverageExpression(scRNA_PPB_BL, assays = "RNA", slot = "data")[[1]]
expr = expr[rowSums(expr)>0,]
expr = as.matrix(expr)
gsva.res = gsva(expr, human_GO_bp_Set, method = "ssgsea")
gsva_gobp = data.frame(Genesets = rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva_gobp, "B_gsva_GOBP.csv", row.names = F)
library(limma)
gsva_gobp = read.csv("B_gsva_GOBP.csv",row.names = 1)
group <- factor(c(rep("PBMC", 5), rep("CSF", 5)), levels = c('PBMC', 'CSF'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(gsva_gobp)
design
# PBMC VS CSF
compare <- makeContrasts(PBMC - CSF, levels=design)
fit <- lmFit(gsva_gobp, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1,n=Inf)
head(Diff)
write.csv(Diff,"diff_B_GOBP.csv")
library(stringr)
Diff$id = rownames(Diff)
Diff$id <- str_replace(Diff$id , "GOBP_","")
Diff$id <- gsub("_"," ",Diff$id)
Diff$id <- tolower(Diff$id)
write.csv(Diff,"diff_PPB.csv")
Diff = read.csv("diff_B_selected.csv",row.names = 1)
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
library(ggplot2)
library(ggthemes)
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, B_PBMC versus B_CSF') +
  guides(fill=F)+
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
high1 <- nrow(dat_plot)
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') +
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + 
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + 
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 
p


## cellchat-----------------------------------------------------------------------------------------------------------------------------------------------------------
library(ComplexHeatmap)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(NMF)
library(reticulate)
library(umapr)
library(patchwork)
Idents(scRNA_pbmc) = "label"
scRNA_pbmc_hc = subset(scRNA_pbmc, SampleType == "HC")
scRNA_pbmc_bl = subset(scRNA_pbmc, SampleType == "BL")
data.input = scRNA_pbmc_bl@assays$RNA@data
meta = scRNA_pbmc_bl@meta.data
unique(meta$label)
meta$label = droplevels(meta$label, exclude = setdiff(levels(meta$label),unique(meta$label)))
levels(meta$label)
cellchat <- createCellChat(object = scRNA_pbmc_bl, meta = meta, group.by = "label")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
netVisual_bubble(cellchat, sources.use = c(10), targets.use = c(1:9,11), remove.isolate = FALSE)#> Comparing communications on a single object
netVisual_bubble(cellchat, sources.use = c(1:9,11), targets.use = c(10), remove.isolate = FALSE)#
p = netVisual_bubble(cellchat, sources.use = c(1:9,11), targets.use = c(10), remove.isolate = FALSE)#> Comparing communications on a single object
data = p$data
data.dir = './cellchat_pbmc_ppb'
dir.create(data.dir)
setwd(data.dir)
cellchat.HC = readRDS("cellchat_HC.rds")
cellchat.BL = readRDS("cellchat_BL.rds")
object.list = list(HC = cellchat.HC, BL = cellchat.BL)
cellchat = mergeCellChat(object.list, add.names = names(object.list))
netVisual_bubble(cellchat, sources.use = c(1:8,11), targets.use = c(10), comparison = c(1,2), angle.x = 45)#> Comparing communications on a single object
p = netVisual_bubble(cellchat, sources.use = c(10), targets.use = c(1:8,11), comparison = c(1,2), angle.x = 45)#> Comparing communications on a single object
data = p$data
write.csv(data, "data_HC_BL_ppb_all_bubbleplot.csv")
scRNA_pbmc_ppb = subset(scRNA_pbmc, label == "Plasma cells")
VlnPlot(scRNA_pbmc_ppb, features = "CD74", group.by = "SampleType")
scRNA_pbmc_cd14 = subset(scRNA_pbmc, label == "CD4+ T cells")
VlnPlot(scRNA_pbmc_cd14, features = "CD4", group.by = "SampleType")
input = read.csv("data_netVisual_chord.csv", header = T)
input$source_short = substr(input$source, 1,4)
input$ligand = paste("L", input$ligand, sep = "_")
input$ligand = paste(input$source_short, input$ligand, sep = "_")
input$target_short = substr(input$target, 1,4)
input$receptor = paste("T", input$receptor, sep = "_")
input$receptor = paste(input$target_short, input$receptor, sep = "_")
input = input[,c(3,4,5)]
write.csv(input, "input.csv")
library(circlize)
library(reshape2)
input = read.csv("input.csv")
order = read.csv("order.csv")
order = order[,1]
chordDiagram(input,order = order,
             col = c(rep("#8CACC4",8),rep("#001E3F",2),
                     rep("#8CACC4",0),rep("#001E3F",2),
                     rep("#8CACC4",0),rep("#001E3F",1),
                     rep("#8CACC4",0),rep("#001E3F",2),
                     rep("#8CACC4",0),rep("#001E3F",1),
                     rep("#8CACC4",0),rep("#001E3F",1),
                     rep("#8CACC4",2),rep("#001E3F",1)),
             annotationTrack = c("name", "grid"))




## QC & AddClonotypeData----------------------------------------------------------------------------------------------------------------------------------------------
library(Seurat)
library(tidyselect)
library(dplyr)
library(DoubletFinder)
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))
  tcr <- tcr[!duplicated(tcr$barcode), ]
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
  tcr <- tcr[, c(2,1,3)]
  rownames(tcr) <- tcr[,1]
  tcr[,1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}
dir1_1 = c('./data/matrix/B_013_7d/')
dir1_2 = c('./data/tcr_bcr_matrix/B_013_7d_TCR/')
dir1_3 = c('./data/tcr_bcr_matrix/B_013_7d_BCR/')
sample_name1 <- c('B_013_7d')
scRNAlist <- list()
for(i in 1:N){
  data <- Read10X(data.dir = dir1_1[i])
  rownames(x = data[["Antibody Capture"]]) <- gsub(pattern = "TotalSeq_C0056", replacement = "BCMA",
                                                   x = rownames(x = data[["Antibody Capture"]]))
  rownames(x = data[["Antibody Capture"]]) <- gsub(pattern = "TotalSeq_C0072", replacement = "CD4",
                                                   x = rownames(x = data[["Antibody Capture"]]))
  rownames(x = data[["Antibody Capture"]]) <- gsub(pattern = "TotalSeq_C0080", replacement = "CD8A",
                                                   x = rownames(x = data[["Antibody Capture"]]))
  rownames(x = data[["Antibody Capture"]]) <- gsub(pattern = "TotalSeq_C0988", replacement = "BCMA_FITC",
                                                   x = rownames(x = data[["Antibody Capture"]]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data[["Gene Expression"]], project=sample_name1[i], min.cells = 3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]]$log10GenesPerUMI <- log10(scRNAlist[[i]]$nFeature_RNA)/log10(scRNAlist[[i]]$nCount_RNA)
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt < 10)
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nCount_RNA>= 1000  &
                             nCount_RNA< 40000 &
                             nFeature_RNA>= 200  &
                             log10GenesPerUMI > 0.7)
  scRNAlist[[i]][["ADT"]] <- CreateAssayObject(data[["Antibody Capture"]][, colnames(x = scRNAlist[[i]])])
  DefaultAssay(scRNAlist[[i]]) <- 'ADT'
  VariableFeatures(scRNAlist[[i]]) <- rownames(scRNAlist[[i]][["ADT"]])
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]], normalization.method = 'CLR', margin = 2) %>% 
    ScaleData() %>% RunPCA(reduction.name = 'apca')
  DefaultAssay(scRNAlist[[i]]) <- 'RNA'
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]]) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
}
for(i in 1:N){
  table(scRNAlist[[i]]$orig.ident)
  Doubletrate = ncol(scRNAlist[[i]])*8*1e-6
  sweep.res.list <- paramSweep_v3(scRNAlist[[i]], PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- scRNAlist[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sample14@meta.data$ClusteringResults
  nExp_poi <- round(Doubletrate*nrow(scRNAlist[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  scRNAlist[[i]] <- doubletFinder_v3(scRNAlist[[i]], PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  scRNAlist[[i]]$doubFind_res = scRNAlist[[i]]@meta.data %>% select(contains('DF.classifications'))
  scRNAlist[[i]]$doubFind_score = scRNAlist[[i]]@meta.data %>% select(contains('pANN'))
  scRNAlist[[i]] = subset(scRNAlist[[i]], doubFind_res == "Singlet")
  table(scRNAlist[[i]]$orig.ident)
}
for(i in 1:N){
  scRNAlist[[i]] <- add_clonotype(dir1_2[i], scRNAlist[[i]],"t")
  scRNAlist[[i]] <- add_clonotype(dir1_3[i], scRNAlist[[i]],"b")
  scRNAlist[[i]] <- subset(scRNAlist[[i]], cells = colnames(scRNAlist[[i]])[!(!is.na(scRNAlist[[i]]$t_clonotype_id) & !is.na(scRNAlist[[i]]$b_clonotype_id))])
}
table(!is.na(scRNAlist[[15]]$t_clonotype_id),!is.na(scRNAlist[[15]]$b_clonotype_id))

scRNA_harmony <- merge(scRNAlist1, y = c(scRNAlist2))
table(scRNA_harmony$orig.ident)
scRNA_harmony <- CellCycleScoring(object = scRNA_harmony,
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes,
                                  g1.features = cc.genes$g1.genes)
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
DefaultAssay(scRNA_harmony) <- 'ADT'
VariableFeatures(scRNA_harmony) <- rownames(scRNA_harmony[["ADT"]])
scRNA_harmony <- NormalizeData(scRNA_harmony, normalization.method = 'CLR', margin = 2)%>% ScaleData() %>% RunPCA(reduction.name = 'apca')
library(harmony)
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
scRNA_harmony <- FindMultiModalNeighbors(scRNA_harmony,
                                         reduction.list = list("pca", "apca"),
                                         dims.list = list(1:50, 1:4),
                                         modality.weight.name = "RNA.weight")
scRNA_harmony <- RunUMAP(scRNA_harmony, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
scRNA_harmony <- FindClusters(scRNA_harmony, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = 'pca', dims = 1:30, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = 'apca', dims = 1:4, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
p3 <- DimPlot(scRNA_harmony, reduction = 'rna.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(scRNA_harmony, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4
scRNA_harmony@meta.data$tcr = scRNA_harmony@meta.data$t_clonotype_id
scRNA_harmony@meta.data$tcr[!is.na(scRNA_harmony@meta.data$tcr)] <- "TCR"
scRNA_harmony@meta.data$bcr = scRNA_harmony@meta.data$b_clonotype_id
scRNA_harmony@meta.data$bcr[!is.na(scRNA_harmony@meta.data$bcr)] <- "BCR"
DimPlot(scRNA_harmony, reduction = "wnn.umap", group.by = "tcr")
scRNA_harmony = subset(scRNA_harmony, tcr == "TCR")
FeaturePlot(scRNA_harmony, features = "BCMA-FITC",reduction = "adt.umap",cols = c("lightgrey","#E5D963","#FCC40E","#FDA900","#FD7A00","#F44F1D","#FA1D0E","#FF1000"))
DefaultAssay(scRNA_harmony) <- 'RNA'
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
DefaultAssay(scRNA_harmony) <- 'ADT'
VariableFeatures(scRNA_harmony) <- rownames(scRNA_harmony[["ADT"]])
scRNA_harmony <- NormalizeData(scRNA_harmony, normalization.method = 'CLR', margin = 2)%>% ScaleData() %>% RunPCA(reduction.name = 'apca')
library(harmony)
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
scRNA_harmony <- FindMultiModalNeighbors(scRNA_harmony,
                                         reduction.list = list("pca", "apca"),
                                         dims.list = list(1:50, 1:4),
                                         modality.weight.name = "RNA.weight")
scRNA_harmony <- RunUMAP(scRNA_harmony, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
scRNA_harmony <- FindClusters(scRNA_harmony, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = 'pca', dims = 1:30, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = 'apca', dims = 1:4, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
p3 <- DimPlot(scRNA_harmony, reduction = 'rna.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(scRNA_harmony, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4
DimPlot(scRNA_harmony, reduction = "wnn.umap", label = T)
FeaturePlot(scRNA_harmony, features = "BCMA-FITC",reduction = "wnn.umap",cols = c("lightgrey","#E5D963","#FCC40E","#FDA900","#FD7A00","#F44F1D","#FA1D0E","#FF1000","#984EA3","#772B84","#FA1D0E","#FA1D0E"))
scRNA_harmony <- FindMultiModalNeighbors(scRNA_harmony,
                                         reduction.list = list("pca", "apca"),
                                         dims.list = list(1:50, 1:4),
                                         modality.weight.name = "RNA.weight")
scRNA_harmony <- RunUMAP(scRNA_harmony, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
scRNA_harmony <- FindClusters(scRNA_harmony, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
DimPlot(scRNA_harmony, reduction = "wnn.umap", label = T)


## Cell_clustering----------------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(scRNA_CAR) <- 'RNA'
pdf(file="markerBubble_1.pdf",width=30,height=11)
cluster10Marker=c("CD3D","CD3E","CD4","CD8A","MKI67",
                  "CCR7","TCF7","LEF1","SELL","CCR6","NR4A1","NCR3","KLRB1","GPR183","IL7R","CD27",
                  "PRF1","GZMA","GZMB","GZMK","NKG7","FCGR3A","KLRG1",
                  "ICOS","PDCD1","LAG3","HAVCR2","CD200","CTLA4","ENTPD1","ITGAE","EOMES","IRF8",
                  "FOXP3","IL2RA","TRGV9","TRDV2",
                  "NCAM1","KLRF1","KLRC1","CD160","CT103A-CAR")#T细胞marker #https://www.jianshu.com/p/0127c9b380c9
DotPlot(object = scRNA_CAR, features = cluster10Marker)
dev.off()
pdf(file="markerBubble_2.pdf",width=15,height=11)
cluster10Marker=c("MS4A1","IGHD","NEIL1","CD27","CD38","ITGAX","TBX21","XBP1","PRDM1","IGHM","COCH","MKI67","IGHA1","IGHG1","MZB1") #http://www.360doc.com/content/12/0121/07/76149697_988292774.shtml
DotPlot(object = scRNA_CAR, features = cluster10Marker,col.min = -1)
dev.off()
pdf(file="markerBubble_3.pdf",width=15,height=11)
cluster10Marker=c("LYZ","CD14","CD83","FCGR3A","C1QA","C1QB","C1QC","CSF1R","TREM2","APOE","TYROBP","CX3CR1","SLC2A5","P2RY13","P2RY12",
                  "CD1C","LILRA4","PF4","CD68","MARCO","FUT4","ITGAM","MME","CXCR2","SELL") #http://www.360doc.com/content/12/0121/07/76149697_988292774.shtml
DotPlot(object = scRNA_CAR, features = cluster10Marker,col.min = -1)
dev.off()
DefaultAssay(scRNA_CAR) <- 'ADT'
pdf(file="markerBubble_4.pdf",width=15,height=11)
cluster10Marker=c("CD4","CD8A","BCMA-FITC","BCMA") #http://www.360doc.com/content/12/0121/07/76149697_988292774.shtml
DotPlot(object = scRNA_CAR, features = cluster10Marker,col.min = -1)
dev.off()
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
new.cluster.ids <- c("CD4+ T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD4+ T",
                     "CD4/CD8 T",
                     "CD4+ T",
                     "CD4+ T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD4+ T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD8+ T",
                     "CD4+ T",
                     "Unknown",
                     "CD8+ T",
                     "CD4/CD8 T",
                     "CD8+ T",
                     "CD4+ T",
                     "Unknown",
                     "CD4/CD8 T",
                     "CD8+ T",
                     "CD4/CD8 T")
scRNA_CAR@active.ident <- plyr::mapvalues(
  x = scRNA_CAR@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_CAR@meta.data$label = scRNA_CAR@active.ident
scRNA_CAR = subset(scRNA_CAR, label != "Unknown")
DimPlot(scRNA_CAR,reduction = "wnn.umap",label = T, cols = c("#A6D38E","#FCBA71","#C0C0C0"))
scRNA_CAR$label = factor(scRNA_CAR$label,levels = c("CD4+ T","CD8+ T","CD4/CD8 T"))


## CD4_CART_reclustering----------------------------------------------------------------------------------------------------------------------------------------------
scRNA_CAR_CD4 = subset(scRNA_CAR, label == "CD4+ T" | label == "CD4/CD8 T")
DefaultAssay(scRNA_CAR_CD4) = "RNA"
scRNA_CAR_CD4 <- NormalizeData(scRNA_CAR_CD4)
scRNA_CAR_CD4 <- FindVariableFeatures(scRNA_CAR_CD4)
scRNA_CAR_CD4 <- ScaleData(scRNA_CAR_CD4)
scRNA_CAR_CD4 <- RunPCA(scRNA_CAR_CD4, verbose=FALSE)
library(harmony)
system.time({scRNA_CAR_CD4 <- RunHarmony(scRNA_CAR_CD4, group.by.vars = "orig.ident")})
plot1 <- DimPlot(scRNA_CAR_CD4, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA_CAR_CD4, ndims=50, reduction="pca")
plotc <- plot1+plot2
plotc
pc.num = 1:30
scRNA_CAR_CD4 <- RunUMAP(scRNA_CAR_CD4, reduction = "harmony", dims = pc.num)
scRNA_CAR_CD4 <- FindNeighbors(scRNA_CAR_CD4, reduction = "harmony", dims = pc.num)
scRNA_CAR_CD4 <- FindClusters(scRNA_CAR_CD4, reduction = "harmony", resolution = 0.5)
scRNA_CAR_CD4 <- FindClusters(scRNA_CAR_CD4, reduction = "harmony", resolution = 1)
DimPlot(scRNA_CAR_CD4, reduction = "umap", label=T)
saveRDS(scRNA_CAR_CD4,"scRNA_IP_PBMC_CITESEQ_CAR_CD4_0804.rds")
pdf(file="markerBubble_1.pdf",width=30,height=11)
cluster10Marker=c("CD3D","CD3E","CD4","CD8A","MKI67",
                  "CCR7","TCF7","LEF1","SELL","CCR6","NR4A1","NCR3","KLRB1","GPR183","IL7R","CD27",
                  "PRF1","GZMA","GZMB","GZMK","NKG7","FCGR3A","KLRG1",
                  "ICOS","PDCD1","LAG3","HAVCR2","CD200","CTLA4","ENTPD1","ITGAE","EOMES","IRF8",
                  "FOXP3","IL2RA","TRGV9","TRDV2",
                  "NCAM1","KLRF1","KLRC1","CD160","CT103A-CAR")
DotPlot(object = scRNA_CAR_CD4, features = cluster10Marker)
dev.off()
pdf(file="markerBubble_2.pdf",width=15,height=11)
cluster10Marker=c("MS4A1","IGHD","NEIL1","CD27","CD38","ITGAX","TBX21","XBP1","PRDM1","IGHM","COCH","MKI67","IGHA1","IGHG1","MZB1")
DotPlot(object = scRNA_CAR_CD4, features = cluster10Marker,col.min = -1)
dev.off()
pdf(file="markerBubble_3.pdf",width=15,height=11)
cluster10Marker=c("LYZ","CD14","CD83","FCGR3A","C1QA","C1QB","C1QC","CSF1R","TREM2","APOE","TYROBP","CX3CR1","SLC2A5","P2RY13","P2RY12",
                  "CD1C","LILRA4","PF4","CD68","MARCO","FUT4","ITGAM","MME","CXCR2","SELL")
DotPlot(object = scRNA_CAR_CD4, features = cluster10Marker,col.min = -1)
dev.off()
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
new.cluster.ids <- c("CD4+ Tem",
                     "Proliferating T cells",
                     "CD4+ Tcm",
                     "Proliferating T cells",
                     "CD4+ Tcm",
                     "CD4+ Tn",
                     "Proliferating T cells",
                     "CD4+ Tem",
                     "CD4+ Tem",
                     "Proliferating T cells",
                     "Unknown",
                     "Proliferating T cells",
                     "Treg",
                     "CD4+ Tem",
                     "Unknown",
                     "Proliferating T cells",
                     "CD4+ Tn",
                     "CD4+ Tem",
                     "Proliferating T cells")
scRNA_CAR_CD4@active.ident <- plyr::mapvalues(
  x = scRNA_CAR_CD4@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_CAR_CD4@meta.data$label = scRNA_CAR_CD4@active.ident
scRNA_CAR_CD4 = subset(scRNA_CAR_CD4, label != "Unknown")
DimPlot(scRNA_CAR_CD4, reduction = "umap", label=T)
scRNA_CAR_CD4@active.ident = factor(scRNA_CAR_CD4@active.ident, levels = c("CD4+ Tn","CD4+ Tcm","CD4+ Tem","Treg","Proliferating T cells"))
scRNA_CAR_CD4@meta.data$label = scRNA_CAR_CD4@active.ident
DimPlot(scRNA_CAR_CD4, reduction = "umap", label=T,group.by = "label",cols = c("#AA90AA","#AACDC5","#7EB793","#3383BC","#FCB863"))+NoLegend()
scRNA_CAR_CD4$SampleType = factor(scRNA_CAR_CD4$SampleType, levels = c("IP","7D","14D","21D","28D"))
col = c("#A4C8DB","#166BA0","#ABCE86","#4EA74C","#984EA3")
DimPlot(scRNA_CAR_CD4, label = T, reduction = "umap", group.by = "SampleType",cols = col)+NoLegend()
scRNA_CAR_CD4$SampleType2 = ifelse(scRNA_CAR_CD4$SampleType == "IP", "IP",
                                   ifelse(scRNA_CAR_CD4$SampleType == "7D"|scRNA_CAR_CD4$SampleType == "14D","Early","Late"))
scRNA_CAR_CD4$SampleType2 = factor(scRNA_CAR_CD4$SampleType2, levels = c("IP","Early","Late"))
sub_012 = subset(scRNA_CAR_CD4, orig.ident == "B_017_28d")
table(sub_012$label)
library(ggplot2)
FeaturePlot(scRNA_CAR_CD4, features = "CD27", reduction = "umap", cols = c("#BFBEBE","#7F69AF","#5165B1"))+NoLegend()
library(viridis)
library(ggplot2)
pdf(file="markerBubble_CD4CAR.pdf",width=9,height=2.5)
cluster10Marker=c("CCR7","TCF7","LEF1","IL7R","CD27","SELL","CD44","FOXP3","IL2RA","MKI67")
DotPlot(object = scRNA_CAR_CD4, features = cluster10Marker)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()
dev.off()
library(dplyr)
colnames(scRNA_CAR_CD4@meta.data)
sub_1 = subset(scRNA_CAR_CD4,orig.ident == "B_013_7d")
b = sub_1@meta.data[,c(12,88)]
b = na.omit(b)
duplicate_name = b %>% group_by(t_clonotype_id) %>% summarise(freq = n()) %>% filter(freq > 1) %>% select(t_clonotype_id)
duplicate_data = b[b$t_clonotype_id %in% duplicate_name$t_clonotype_id, ]
non_duplicate_data = b[!b$t_clonotype_id %in% duplicate_name$t_clonotype_id, ]
table(non_duplicate_data$label)
table(duplicate_data$label)
library(COSG)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)
Idents(scRNA_CAR_CD4) = "SampleType2"
table(scRNA_CAR_CD4$SampleType2)
scRNA_CAR_CD4 = NormalizeData(scRNA_CAR_CD4)
all.genes = rownames(scRNA_CAR_CD4)
scRNA_CAR_CD4 = ScaleData(scRNA_CAR_CD4, features = all.genes)
mat <- GetAssayData(scRNA_CAR_CD4, slot = "scale.data")
cluster_info <- sort(scRNA_CAR_CD4$SampleType2)
# gene = c(marker_cosg$names[,1][1:30],marker_cosg$names[,2],marker_cosg$names[,3])
gene = read.csv("CD4.csv")
gene = gene$Gene
mat <- as.matrix(mat[gene, names(cluster_info)])
scRNA_CAR_CD4$Patient = ifelse(scRNA_CAR_CD4$orig.ident == "CT103A_013"|scRNA_CAR_CD4$orig.ident == "B_013_7d"|scRNA_CAR_CD4$orig.ident == "B_013_14d"|scRNA_CAR_CD4$orig.ident == "B_013_21d"|scRNA_CAR_CD4$orig.ident == "B_013_28d","NMO_013",
                               ifelse(scRNA_CAR_CD4$orig.ident == "CT103A_014"|scRNA_CAR_CD4$orig.ident == "B_014_7d"|scRNA_CAR_CD4$orig.ident == "B_014_14d"|scRNA_CAR_CD4$orig.ident == "B_014_21d"|scRNA_CAR_CD4$orig.ident == "B_014_28d","NMO_014",
                                      ifelse(scRNA_CAR_CD4$orig.ident == "CT103A_015"|scRNA_CAR_CD4$orig.ident == "B_015_7d"|scRNA_CAR_CD4$orig.ident == "B_015_14d"|scRNA_CAR_CD4$orig.ident == "B_015_21d"|scRNA_CAR_CD4$orig.ident == "B_015_28d","NMO_015",
                                             ifelse(scRNA_CAR_CD4$orig.ident == "CT103A_016"|scRNA_CAR_CD4$orig.ident == "B_016_7d"|scRNA_CAR_CD4$orig.ident == "B_016_14d"|scRNA_CAR_CD4$orig.ident == "B_016_21d"|scRNA_CAR_CD4$orig.ident == "B_016_28d","NMO_016","NMO_017"))))
Patient_info = scRNA_CAR_CD4$Patient
Patient_df = data.frame(cell = names(Patient_info),Patient = Patient_info)
cluster_df = data.frame(cell = names(cluster_info),cluster = cluster_info)
col_df = inner_join(cluster_df, Patient_df)
rownames(col_df) = col_df$cell
col_df = col_df[,-1]
table(col_df$cluster)
table(col_df$Patient)
# colours <- list(
#   "cluster"=c("IP"="#A4C8DB",
#               "7D"="#166BA0",
#               "14D"="#ABCE86",
#               "21D"="#4EA74C",
#               "28D"="#984EA3"),
#   "Patient"=c("NMO_012"="#DAB8DA","NMO_013"="#9BC3DF","NMO_014"="#789CBF",
#               "NMO_015"="#B9B9B9","NMO_016"="#7EB793","NMO_017"="#F3754E"))
colours <- list(
  "cluster"=c("IP"="#A4C8DB",
              "Early"="#166BA0",
              "Late"="#ABCE86"),
  "Patient"=c("NMO_013"="#9BC3DF","NMO_014"="#789CBF",
              "NMO_015"="#B9B9B9","NMO_016"="#7EB793","NMO_017"="#F3754E"))
col_anno = HeatmapAnnotation(df=col_df, 
                             which="col", # 表示 column annotations
                             col=colours, 
                             annotation_width=unit(c(1, 4), "cm"), 
                             gap=unit(1, "mm"))
col_heatmap = colorRamp2(c(-2,0,2),c("#417D8B","white","#BB3E45"))
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        # right_annotation = row_anno,
        top_annotation = col_anno,
        #column_split = cluster_info,
        col = col_heatmap,
        column_title = NULL,
        use_raster = F)
clonotype = scRNA_CAR_CD4$t_clonotype_id
clonotype_df = data.frame(cell = names(clonotype),type = clonotype,orig.ident = scRNA_CAR_CD4$orig.ident)
clonotype_df = na.omit(clonotype_df)
for(i in 1:nrow(clonotype_df)){
  clonotype_df$type[i] = paste(clonotype_df$orig.ident[i],clonotype_df$type[i])
}
celltype = scRNA_CAR_CD4$label
celltype_df = data.frame(cell = names(celltype),celltype = celltype)
celltype_clonotype_df = merge(clonotype_df,celltype_df,all.y = T)
clononumber_df = data.frame(type = names(table(celltype_clonotype_df$type)),clononumber = table(celltype_clonotype_df$type))
clononumber_df = clononumber_df[,-2]
colnames(clononumber_df) = c("type","clononumber")
celltype_clonotype_clononumber_df = merge(celltype_clonotype_df,clononumber_df,by = "type",all.x = T)
celltype_clonotype_clononumber_df$clonostate = ifelse(celltype_clonotype_clononumber_df$clononumber >= 5, ">5","<5")
scRNA_CAR_CD4@meta.data$cell = rownames(scRNA_CAR_CD4@meta.data)
clononumber_df = celltype_clonotype_clononumber_df[,c(2,5)]
clononumber_df[,2] = as.numeric(clononumber_df[,2])
#scRNA_CAR_CD4@meta.data = merge(scRNA_CAR_CD4@meta.data,clononumber_df,by = "cell",all.x = T)
scRNA_CAR_CD4@meta.data = inner_join(scRNA_CAR_CD4@meta.data,clononumber_df)
clonostate_df = celltype_clonotype_clononumber_df[,c(2,6)]
#scRNA_CAR_CD4@meta.data = merge(scRNA_CAR_CD4@meta.data,clonostate_df,by = "cell",all.x = T)
scRNA_CAR_CD4@meta.data = inner_join(scRNA_CAR_CD4@meta.data,clonostate_df)
rownames(scRNA_CAR_CD4@meta.data) = scRNA_CAR_CD4@meta.data$cell
library(ggplot2)
DimPlot(scRNA_CAR_CD4, reduction = "umap", group.by = "clonostate",cols = c("#D1D0E4","#463179"),raster=FALSE)+
  NoLegend()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank())
# FeaturePlot(scRNA_CAR_CD4,features = "clononumber",max.cutoff = 500)+ 
#   theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
celltype_clonotype_clononumber_df$clonostate = ifelse(celltype_clonotype_clononumber_df$clononumber >= 10, ">=10",
                                                      ifelse(celltype_clonotype_clononumber_df$clononumber >= 5,"5-9",
                                                             ifelse(celltype_clonotype_clononumber_df$clononumber >= 2,"2-4","Unique")))
scRNA_CAR_CD4@meta.data$cell = rownames(scRNA_CAR_CD4@meta.data)
clononumber_df = celltype_clonotype_clononumber_df[,c(2,5)]
clononumber_df[,2] = as.numeric(clononumber_df[,2])
#scRNA_CAR_CD4@meta.data = merge(scRNA_CAR_CD4@meta.data,clononumber_df,by = "cell",all.x = T)
scRNA_CAR_CD4@meta.data = inner_join(scRNA_CAR_CD4@meta.data,clononumber_df)
clonostate_df = celltype_clonotype_clononumber_df[,c(2,6)]
#scRNA_CAR_CD4@meta.data = merge(scRNA_CAR_CD4@meta.data,clonostate_df,by = "cell",all.x = T)
scRNA_CAR_CD4@meta.data = inner_join(scRNA_CAR_CD4@meta.data,clonostate_df)
rownames(scRNA_CAR_CD4@meta.data) = scRNA_CAR_CD4@meta.data$cell
sub_10 = subset(scRNA_CAR_CD4, clonostate == "Unique")
table(sub_10$orig.ident)
DefaultAssay(scRNA_CAR_CD4) <- "RNA"
list = read.csv("genelist.csv")
cd_features <- list(list[,1])
Inscore <- AddModuleScore(scRNA_CAR_CD4,
                          features = cd_features,
                          ctrl = 100,
                          name = "CD_Features")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[93] <- 'Exhaustion_Score'
p = VlnPlot(Inscore,features = 'Exhaustion_Score', 
        pt.size = 0, adjust = 2,group.by = "SampleType",
        cols = c("#AB2A4F","#BC5572","#CD7F97","#DEABBA","#EFD5DD"))
p
data = p$data
table(scRNA_CAR_CD4$SampleType2)
table(scRNA_CAR_CD4$t_clonotype_id)
scRNA_CAR_CD4@meta.data$clonolabel = ifelse(scRNA_CAR_CD4$SampleType2 != "Late",'Other Timepoints',
                                        ifelse(scRNA_CAR_CD4$t_clonotype_id == "clonotype1" ,'Clonotype1',
                                            ifelse(scRNA_CAR_CD4$t_clonotype_id == "clonotype2" ,'Clonotype2',
                                                ifelse(scRNA_CAR_CD4$t_clonotype_id == "clonotype3" ,'Clonotype3',
                                                    ifelse(scRNA_CAR_CD4$t_clonotype_id == "clonotype4" ,'Clonotype4','Late') )
                                                      )))
scRNA_CAR_CD4$clonolabel = factor(scRNA_CAR_CD4$clonolabel,
                              levels = c("Late","Other Timepoints",
                                         "Clonotype1","Clonotype2",
                                         "Clonotype3","Clonotype4"))
table(scRNA_CAR_CD4$clonolabel)
col = c("#B2CFE5","#D2D1D2","#8A2A46","#D97351","#5384B6","#61549F")
library(ggplot2)
DimPlot(scRNA_CAR_CD4, reduction = "umap", group.by = "clonolabel",cols = col,pt.size = 1)+ NoLegend()


## Myeloid_reclustering-----------------------------------------------------------------------------------------------------------------------------------------------
scRNA_mono = subset(scRNA_all, label == "CD14+ monocytes"|label == "CD16+ monocytes"|
                               label == "Macrophages"|label == "pDCs"|label == "cDCs")
pc.num = 1:30
scRNA_mono <- RunUMAP(scRNA_mono, reduction = "harmony", dims = pc.num)
DimPlot(scRNA_mono, label =T)
scRNA_mono <- FindNeighbors(scRNA_mono, reduction = "harmony", dims = 1:30)
scRNA_mono <- FindClusters(scRNA_mono, reduction = "harmony", resolution = 1)
DimPlot(scRNA_mono, label = T,reduction = "umap")
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
new.cluster.ids <- c("CD14+ monocytes",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "CD16+ monocytes",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "cDCs",
                     "Macrophages",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "Unknown",
                     "CD14+ monocytes",
                     "Unknown",
                     "Macrophages",
                     "Macrophages",
                     "Unknown",
                     "Macrophages",
                     "pDCs",
                     "Macrophages",
                     "CD16+ monocytes",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown")
scRNA_mono@active.ident <- plyr::mapvalues(
  x = scRNA_mono@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_mono@meta.data$label = scRNA_mono@active.ident
scRNA_mono = subset(scRNA_mono, label != "Unknown")
DimPlot(scRNA_mono, reduction = "umap",label=T, repel = T)
scRNA_mono@active.ident = factor(scRNA_mono@active.ident, levels = c("CD14+ monocytes","CD16+ monocytes","Macrophages","cDCs","pDCs"))
scRNA_mono$label = scRNA_mono@active.ident
col = c("#C6CF50","#587F58","#9D0068","#3A6A89","#AA9BC5")
scRNA_mono_HC_BL = subset(scRNA_mono,SampleType!= "3M")
library(ggplot2)
DimPlot(scRNA_mono_HC_BL, reduction = "umap",label=T, repel = T, cols = col,pt.size = 0.8,group.by = "label")+
  NoLegend()+
  labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
col5 = c("#C986B6","#9BC99A")
DimPlot(scRNA_mono_HC_BL, reduction = "umap",label=T, repel = T, cols = col5,pt.size = 0.8,group.by = "SampleOrig")+
  NoLegend()+
  labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
sub_PBMC = subset(scRNA_mono_HC_BL, SampleOrig == "blood")
sub_CSF = subset(scRNA_mono_HC_BL, SampleOrig == "csf")
table(sub_PBMC$label)
table(sub_CSF$label)
Idents(scRNA_mono) = "seurat_clusters"
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,12,14,15,17,18,19,20)
new.cluster.ids <- c("CD14+ monocytes",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "CD16+ monocytes",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "cDCs",
                     "Macrophage-like",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "CD14+ monocytes",
                     "Microglia-like",
                     "Macrophage-like",
                     "Macrophage-like",
                     "pDCs",
                     "Microglia-like",
                     "CD16+ monocytes") 
scRNA_mono@active.ident <- plyr::mapvalues(
  x = scRNA_mono@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_mono@meta.data$label = scRNA_mono@active.ident
DimPlot(scRNA_mono, reduction = "umap",label=T, repel = T)
SampleOrig = scRNA_mono$SampleOrig
SampleOrig_df = data.frame(cell = names(SampleOrig),type = SampleOrig,label = scRNA_mono$label)
for(i in 1:nrow(SampleOrig_df)){
  SampleOrig_df$SampleClass[i] = paste(SampleOrig_df$type[i],SampleOrig_df$label[i])
}
scRNA_mono@meta.data$cell = rownames(scRNA_mono@meta.data)
SampleOrig_df = SampleOrig_df[,c(1,4)]
scRNA_mono@meta.data = inner_join(scRNA_mono@meta.data,SampleOrig_df)
rownames(scRNA_mono@meta.data) = scRNA_mono@meta.data$cell
DimPlot(scRNA_mono, group.by = "SampleClass")
table(scRNA_mono$SampleClass)
library(monocle)
library(tidyverse)
library(dplyr)
library(patchwork)
Idents(scRNA_mono) = scRNA_mono$SampleClass
table(Idents(scRNA_mono))
levels(scRNA_mono)
data <- as(as.matrix(scRNA_mono@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_mono@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#build new cell data set
test <- newCellDataSet(data,phenoData =pd,featureData =fd,expressionFamily=negbinomial.size(),lowerDetectionLimit=0.5)
test <- estimateSizeFactors(test)
test <- estimateDispersions(test)
test <- detectGenes(test, min_expr = 0.1)
print(head(fData(test)))
expressed_genes = row.names(subset(fData(test),num_cells_expressed >= 10))
Sys.time()
diff <- differentialGeneTest(test[expressed_genes,],
                             fullModelFormulaStr = "~label")
Sys.time()
deg = subset(diff, qval < 0.01)
deg = deg[order(deg$qval, decreasing = F),]
#write.table(deg, "test.monocle.DEG.xls", col.names = T, row.names = F, sep = "\t", quote = F)
ordergene = rownames(deg)
test = setOrderingFilter(test, ordergene)
#ordergene = row.names(deg)[order(deg$qval)][1:2000]
test=reduceDimension(test,method = "DDRTree",max_components = 2) 
test=orderCells(test)
library(COSG)
library(BuenColors)
library(circlize)
library(dplyr)
library(ComplexHeatmap)
scRNA_mono$CellClass = ifelse(scRNA_mono$SampleClass == "blood CD14+ monocytes"|scRNA_mono$SampleClass == "csf CD14+ monocytes","CD14+ monocytes",
                          ifelse(scRNA_mono$SampleClass == "blood CD16+ monocytes"|scRNA_mono$SampleClass == "csf CD16+ monocytes","CD16+ monocytes",
                              ifelse(scRNA_mono$SampleClass == "blood Macrophage-like","blood Macrophage-like",
                                  ifelse(scRNA_mono$SampleClass == "csf Macrophage-like","csf Macrophage-like","csf Microglia-like"))))
table(scRNA_mono$CellClass)
scRNA_mono$CellClass = factor(scRNA_mono$CellClass, levels = c("CD14+ monocytes","CD16+ monocytes","blood Macrophage-like","csf Macrophage-like","csf Microglia-like"))
scRNA_mono$SampleClass = factor(scRNA_mono$SampleClass, levels = c("blood CD14+ monocytes","blood CD16+ monocytes","blood Macrophage-like",
                                                                   "csf CD14+ monocytes","csf CD16+ monocytes","csf Macrophage-like","csf Microglia-like"))
Idents(scRNA_mono) = "SampleClass"
# sample = subset(scRNA_pre_3m_mono, downsample = 1000)
mat <- GetAssayData(scRNA_mono, slot = "counts")
cluster_info <- sort(scRNA_mono$SampleClass)
#gene = c(marker_cosg$names[,1],marker_cosg$names[,2],marker_cosg$names[,3],marker_cosg$names[,4],marker_cosg$names[,5],marker_cosg$names[,6])
gene = c("LYZ","S100A9","S100A8","FCN1","VCAN","S100A12",
         "CDKN1C","HES4","CKB","TCF7L2","FCGR3A",
         "CX3CR1","LYVE1","F13A1","MRC1",
         "P2RY12","TREM2","APOE","APOC1","GPNMB")
mat <- as.matrix(mat[gene, names(cluster_info)])
mat = log2(mat + 1)
SampleClass = scRNA_mono$SampleClass
SampleClass_df = data.frame(cell = names(SampleClass),SampleClass = SampleClass)
CellClass = scRNA_mono$CellClass
CellClass_df = data.frame(cell = names(CellClass),CellClass = CellClass)
col_df = inner_join(SampleClass_df, CellClass_df)
col_df = col_df[,-1]
col_df = col_df[order(col_df[,1]),]
# cycle_info = scRNA_mono$Phase
# cycle_df = data.frame(cell = names(cycle_info),cycle = cycle_info)
# col_df = join(cluster_df, type_df)
# col_df = join(col_df,cycle_df,by = "cell")
# rownames(col_df) = col_df$cell
# col_df = col_df[,-1]
# table(col_df$cluster)
# table(col_df$type)
# table(col_df$cycle)
colours <- list(
  "SampleClass"=c("blood CD14+ monocytes"="#FFC0CB","csf CD14+ monocytes"="#87CEFA",
                  "blood CD16+ monocytes"="#FF69B4","csf CD16+ monocytes"="#00BFFF",
                  "blood Macrophage-like"="#DC143C","csf Macrophage-like"="#0000FF",
                  "csf Microglia-like"="#000080"),
  "CellClass"=c("CD14+ monocytes"="#A4CCE1",
                "CD16+ monocytes"="#587F58",
                "blood Macrophage-like"="#C6CF50","csf Macrophage-like"="#1A547E",
                "csf Microglia-like"="#A389C1"))
#"type"=c("B_012_pre"="#52679F","B_013_pre"="#52679F","B_014_pre"="#52679F","B_015_pre"="#52679F",
#         "B_016_pre"="#52679F","B_017_pre"="#52679F","B_HC_2"="#B2C3E0","B_HC_3"="#B2C3E0","B_HC_4"="#B2C3E0"),
#"cycle"=c("G1"="#6131A0","G2M"="#E79D76","S"="#F2D1A8"))
col_anno = HeatmapAnnotation(df=col_df,
                             which="col",
                             col=colours,
                             annotation_width=unit(c(1, 4), "cm"),
                             gap=unit(1, "mm"))

# col_heatmap = colorRamp2(c(-1,1.5,3),c("#377EB8","white","#E41A1C"))
col_heatmap = colorRamp2(c(-1,0.2,3),c("#37538B","white","#7C262C"))
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        #right_annotation = row_anno,
        top_annotation = col_anno,
        #column_split = cluster_info,
        col = col_heatmap,
        column_title = NULL,
        use_raster = FALSE)
scRNA_mg = subset(scRNA_csf, label == "Macrophages")
scRNA_mg <- RunUMAP(scRNA_mg, reduction = "harmony", dims = 1:26)
scRNA_mg <- FindNeighbors(scRNA_mg, reduction = "harmony", dims = 1:26)
scRNA_mg <- FindClusters(scRNA_mg, reduction = "harmony", resolution = 1)
DimPlot(scRNA_mg, label = T,reduction = "umap")
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12)
new.cluster.ids <- c("Macrophage-like cells",
                     "Macrophage-like cells",
                     "Macrophage-like cells",
                     "Macrophage-like cells",
                     "Macrophage-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Macrophage-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Macrophage-like cells")
scRNA_mg@active.ident <- plyr::mapvalues(
  x = scRNA_mg@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_mg@meta.data$label = scRNA_mg@active.ident
DimPlot(scRNA_mg, reduction = "umap",label=T, repel = T)
scRNA_mg@active.ident = factor(scRNA_mg@active.ident, levels = c("Macrophage-like cells","Microglia-like cells"))
scRNA_mg$label = scRNA_mg@active.ident
col = c("#52195D","#FFE42A")
scRNA_mg_HC_BL = subset(scRNA_mg,SampleType != "3M")
library(ggplot2)
DimPlot(scRNA_mg_HC_BL, reduction = "umap",label=T, repel = T, cols = col,pt.size = 0.8,group.by = "label")+
  NoLegend()+
  labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_blank())+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
library(viridis)
library(ggplot2)
pdf(file="markerBubble.pdf",width=8,height=2.5)
cluster10Marker=c("LYVE1","F13A1",
                  "CX3CR1","TREM2","P2RY12","APOE","GPNMB","CTSD","APOC1","LGALS3")
DotPlot(object = scRNA_mg_HC_BL, features = cluster10Marker)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()
dev.off()
DefaultAssay(scRNA_mg_HC_BL) <- "RNA"
list = read.csv("genelist.csv")
cd_features <- list(list[,1])
Inscore <- AddModuleScore(scRNA_mg_HC_BL,
                          features = cd_features,
                          ctrl = 100,
                          name = "CD_Features")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[94] <- 'Inflammatory_Score'
VlnPlot(Inscore,features = 'Inflammatory_Score', 
        pt.size = 0, adjust = 2,group.by = "orig.ident")
library(ggplot2)
Inscore_3M = Inscore
Inscore_3M$Inflammatory_Score = ifelse(Inscore_3M$SampleType == "BL",Inscore_3M$Inflammatory_Score, -0.1)
mydata<- FetchData(Inscore_3M,vars = c("UMAP_1","UMAP_2","Inflammatory_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = Inflammatory_Score))+
  geom_point(size = 2)+
  xlim(-12,10)+
  ylim(-6,4)+
  scale_color_gradientn(values = seq(0,1,0.25),
                        limits=c(-0.1,0.4),
                        colours = c('#FDF5EE',"#F3CBB8","#E88F70","#D7523E","#A92E29","#611818"))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  NoLegend()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
library(GSVA)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(msigdbr)
human <- msigdbr(species = "Homo sapiens")
table(human$gs_cat)
table(human$gs_subcat)
human_GO_bp = msigdbr(species = "Homo sapiens",
                      category = "C2", 
                      subcategory = "CP:WIKIPATHWAYS") %>%
  dplyr::select(gs_name,gene_symbol)
human_GO_bp_Set = human_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)
table(scRNA_mg$orig.ident)
DefaultAssay(scRNA_mg) = "RNA"
scRNA_mg = NormalizeData(scRNA_mg)
Idents(scRNA_mg) = "orig.ident"
expr = AverageExpression(scRNA_mg, assays = "RNA", slot = "data")[[1]]
expr = expr[rowSums(expr)>0,]
expr = as.matrix(expr)
gsva.res = gsva(expr, human_GO_bp_Set, method = "ssgsea")
gsva_gobp = data.frame(Genesets = rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva_gobp, "MG_gsva_WIKIPATHWAYS.csv", row.names = F)


## CSOmap-------------------------------------------------------------------------------------------------------------------------------------------------------------
library(data.table)
table(scRNA_csf$orig.ident)
table(scRNA_csf$label)
scRNA_csf$label = gsub("B cells","BPPB",scRNA_csf$label)
scRNA_csf$label = gsub("Plasma cells","BPPB",scRNA_csf$label)
table(scRNA_csf$label)
sum(scRNA_csf_pre[["RNA"]]@data[,1])
TPM_csf_pre <- expm1(scRNA_csf_pre[["RNA"]]@data)
head(colSums(TPM_csf_pre))
TPM_csf_pre_df = as.data.frame(TPM_csf_pre)
fwrite(TPM_csf_pre_df, row.names = T, sep = "\t",file = "TPM.txt")

scRNA_csf_pre$cells = rownames(scRNA_csf_pre@meta.data)
colnames(scRNA_csf_pre@meta.data)
meta = scRNA_csf_pre@meta.data
colnames(meta)
meta = meta[,c(93,88)]
colnames(meta) = c("cells","labels")
fwrite(meta, row.names = F, sep = "\t",file = "label.txt")

#CSOmap <- function(DataSetName) {
library(plotly)
library(data.table)
# load data
TPMpath <- paste0("./data/", DataSetName, "/TPM.txt")
LRpath <- paste0("./data/", DataSetName, "/LR_pairs.txt")
Labelpath <- paste0("./data/", DataSetName, "/label.txt")
TPM <- read.table(TPMpath, header = TRUE, sep = "\t")
LR <- read.table(LRpath, header = FALSE, sep = "\t")
labelData <- read.table(Labelpath, header = TRUE, sep = "\t")
# create output path
dir.create(paste0("./results/"))
# set variables properly
genenames <- TPM$X
TPM <- TPM[,2:dim(TPM)[2]]
TPM[is.na(TPM)] = 0
cellnames <- colnames(TPM)
labels <- labelData$labels[match(cellnames, labelData$cells)]
labels[is.na(labels)] = "unlabeled"
standards <- unique(labels)
labelIx <- match(labels, standards)
cellCounts <- table(labelIx)
# find ligands and receptors TPM
ligandsIndex <- match(c(LR$V1, LR$V2), genenames)
receptorIndex <- match(c(LR$V2, LR$V1), genenames)
ligandsTPM <- as.matrix(TPM[ligandsIndex[!is.na(ligandsIndex)&!is.na(receptorIndex)],])
receptorTPM <- as.matrix(TPM[receptorIndex[!is.na(ligandsIndex)&!is.na(receptorIndex)],])
LRscores <- c(LR$V3, LR$V3)
LRscores <- LRscores[!is.na(ligandsIndex)&!is.na(receptorIndex)]
affinityMat <- t(ligandsTPM) %*% diag(LRscores) %*% receptorTPM
# discret affinityMat
for (i in 1:dim(affinityMat)) {
  affinityArray <- affinityMat[i,]
  affinityArray[i] = 0
  affinityArraySorted <- sort(affinityArray, decreasing = TRUE)
  affinityArray[affinityArray <= affinityArraySorted[50]] = 0
  affinityMat[i,] = affinityArray
}
# optimization
coords <- optimization(affinityMat)
coordsPath <- paste0("./results/", "/coordinates.txt")
row.names(coords) <- cellnames
colnames(coords) <- c('x', 'y', 'z')
# write coords
write.table(coords, coordsPath, quote = FALSE, sep = "\t")
# do statistical tests
dist <- as.matrix(dist(coords))
counts <- matrix(0, nrow = length(standards), ncol = length(standards))
colnames(counts) <- standards
rownames(counts) <- standards
k = 3 # default top k connections
topKs <- c()
diag(dist) <- Inf
for (i in 1:dim(dist)[1]) {
  distSorted <- sort(dist[i,])
  topKs <- c(topKs, distSorted[3])
}
topK <- median(topKs)
for (i in 1:dim(dist)[1]) {
  connects <- which(dist[i,]<=topK)
  for (j in connects) {
    counts[labelIx[i], labelIx[j]] = counts[labelIx[i], labelIx[j]] + 1
  }
}
diag(counts) <- diag(counts) / 2
# write counts
countsPath <- paste0("./results/", "/counts.txt")
write.table(counts, countsPath, quote = TRUE, sep = "\t")
# calculate stats
countsN <- (sum(counts) + sum(diag(counts)))/2
p_value <- copy(counts)
K = countsN
for (i in 1:dim(counts)[1]) {
  for (j in 1:dim(counts)[2]) {
    if (i == j) {
      M <- as.numeric(cellCounts[i]) * (as.numeric(cellCounts[i])-1) / 2
    } else {
      M <- as.numeric(cellCounts[i]) * (as.numeric(cellCounts[j]))
    }
    N <- sum(cellCounts) * (sum(cellCounts) - 1) / 2 - M
    p_value[i, j] <- phyper(counts[i, j], M, N, K, lower.tail = FALSE)
  }
}
q_value <- matrix(p.adjust(p_value), nrow = length(standards), ncol = length(standards))
colnames(q_value) <- standards
rownames(q_value) <- standards
# write p and q value
pvaluePath <- paste0("./results/", "/pvalue.txt")
qvaluePath <- paste0("./results/", "/qvalue.txt")
write.table(p_value, pvaluePath, quote = TRUE, sep = "\t")
write.table(q_value, qvaluePath, quote = TRUE, sep = "\t")
result = list()
result$coords = coords
result$counts = counts
result$pvalue = p_value
result$qvalue = q_value
# # plot 3d and heatmap
plot_ly(x=coords[,1], y=coords[,2], z=coords[,3], type="scatter3d", mode="markers", color=labels)
heatmap(p_value, scale = "none", Colv = NA, Rowv = NA)
heatmap(q_value, scale = "none", Colv = NA, Rowv = NA)
result
}

optimization <- function (affinityMat, initial_config = NULL, k = 3,  
                          max_iter = 1000, min_cost = 0, 
                          condition = "tight") 
{ # this function is inspired from tsne algorithm. We use similar gradient descent
  # method to optimize our target function specified in our paper.
  # condition can be loose or tight, we suggest using "loose" condition 
  # for dataset with over 10000 cells
  n = dim(affinityMat)[1]
  momentum = 0.5 # initial momentum
  final_momentum = 0.8 # final momentum
  mom_switch_iter = 250 # value to which momentum is changed
  epsilon = 1000 # initial learning rate
  min_gain = 0.01 # minimum gain for delta-bar-delta
  eps = 2^(-52)
  epoch = 100
  if (!is.null(initial_config) && is.matrix(initial_config)) {
    if (nrow(initial_config) != n | ncol(initial_config) != 
        k) {
      stop("initial_config argument does not match necessary configuration for X")
    }
    ydata = initial_config
  }
  else {
    ydata = matrix(rnorm(k * n), n)
  }
  P = 0.5 * (affinityMat + t(affinityMat))
  P[P < eps] <- eps
  P = P/sum(P)
  grads = matrix(0, nrow(ydata), ncol(ydata))
  incs = matrix(0, nrow(ydata), ncol(ydata))
  gains = matrix(1, nrow(ydata), ncol(ydata))
  for (iter in 1:max_iter) {
    
    sum_ydata = apply(ydata^2, 1, sum)
    d = sum_ydata + sweep(-2 * ydata %*% t(ydata), 2, -t(sum_ydata))
    num = 1/(1 + d)
    diag(num) = 0
    Q = num/sum(num)
    if (any(is.nan(num))) 
      message("NaN in grad. descent")
    Q[Q < eps] = eps
    P_Q = P - Q
    P_Q[P_Q > 0 & d <= 0.01] = -0.01;
    stiffnesses = 4 * (P - Q) * num
    for (i in 1:n) {
      grads[i, ] = apply(sweep(-ydata, 2, -ydata[i, ]) * 
                           stiffnesses[, i], 2, sum)
    }
    gains = ((gains + 0.2) * abs(sign(grads) != sign(incs)) + 
               gains * 0.8 * abs(sign(grads) == sign(incs)))
    gains[gains < min_gain] = min_gain
    incs = momentum * incs - epsilon * (gains * grads)
    ydata = ydata + incs
    ydata = sweep(ydata, 2, apply(ydata, 2, mean))
    if (iter == mom_switch_iter) 
      momentum = final_momentum
    if (iter%%epoch == 0) {
      cost = sum(apply(P * log((P + eps)/(Q + eps)), 1, 
                       sum))
      message("Iteration #", iter, " loss function cost is: ", 
              cost)
      if (cost < min_cost) 
        break
    }
    range = max(abs(ydata))
    if (condition == "tight") {
      if (range > 50 && iter%%10 == 0) {
        ydata = ydata * 50/range
      }
    } else {
      if (range > 50 && iter%%max_iter == 0) {
        ydata = ydata * 50/range
      }
    }
  }
  ydata
}
