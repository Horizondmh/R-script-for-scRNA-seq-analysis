rm(list = ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
gc()

library(Seurat)
library(tidyselect)
library(dplyr)

add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))
  
  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]
  
  # Only keep the barcode and clonotype columns.
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))
  
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
  
  # Reorder so barcodes are first column and set them as rownames.
  tcr <- tcr[, c(2,1,3)]
  rownames(tcr) <- tcr[,1]
  tcr[,1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}

dir1_1 = c('./data/matrix_/B_018_pre/')
dir1_2 = c('./data/tcr_bcr_matrix/B_018_pre_TCR/')
dir1_3 = c('./data/tcr_bcr_matrix/B_018_pre_BCR/')
sample_name1 <- c('B_018_pre')

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
  scRNAlist[[i]] <- FindClusters(scRNAlist[[i]], resolution = 0.5)
}


for(i in 1:N){
  table(scRNAlist[[i]]$orig.ident)
  Doubletrate = ncol(scRNAlist[[i]])*8*1e-6
  
  sweep.res.list <- paramSweep_v3(scRNAlist[[i]], PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- scRNAlist[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sample14@meta.data$ClusteringResults
  nExp_poi <- round(Doubletrate*nrow(scRNAlist[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  scRNAlist[[i]] <- doubletFinder_v3(scRNAlist[[i]], PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  ## save results
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
saveRDS(scRNA_harmony, "scRNA_HC_PRE.rds")



scRNA_all = readRDS("scRNA_HC_PRE.rds")
scRNA_all <- FindNeighbors(scRNA_all, reduction = "harmony", dims = 1:30)
scRNA_all <- FindClusters(scRNA_all, reduction = "harmony", resolution = 2.5)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42)
new.cluster.ids <- c("T_CD4_naive",
                     "T_CD4_mem",
                     "T_CD4_naive",
                     "T_CD8_mem",
                     "T_CD4_mem",
                     "T_CD8_tox",
                     "T_CD8_naive",
                     "Mono_CD14",
                     "NK",
                     "NK",
                     "Mono_CD14",
                     "B_naive",
                     "Mono_CD14",
                     "T_CD4_mem",
                     "B_mem",
                     "T_CD8_mem",
                     "Treg",
                     "Unknown",
                     "Mono_CD16",
                     "NK",
                     "Mono_CD14",
                     "cDCs",
                     "T_CD8_tox",
                     "T_CD8_tox",
                     "Mono_CD14",
                     "T_CD8_mem",
                     "T_CD8_tox",
                     "Mono_CD14",
                     "Cycling",
                     "T_CD8_naive",
                     "Unknown",
                     "Unknown",
                     "NKT",
                     "NKT",
                     "T_CD4_mem",
                     "Unknown",
                     "Plasma",
                     "T_CD4_mem",
                     "pDCs",
                     "Unknown",
                     "T_CD8_tox",
                     "Unknown",
                     "Unknown")
scRNA_all@active.ident <- plyr::mapvalues(
  x = scRNA_all@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_all@meta.data$label = scRNA_all@active.ident
table(scRNA_all$orig.ident)
scRNA_all = subset(scRNA_all, label != "Unknown")
scRNA_all <- RunUMAP(scRNA_all, reduction = "harmony", dims = 1:30, min.dist = 0.5, n.neighbors = 80L)
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T)
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T, group.by = "seurat_clusters")
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T,group.by = "tcr")
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T,group.by = "bcr")

sub_018 = subset(scRNA_all, orig.ident == "B_018_pre")
table(sub_018$label)
# Fisher's exact test
data = read.csv("Fisher_exact_test.csv",header = T,row.names = 1)
fisher.test(data,alternative = "two.sided")$p.value
data = read.csv("BH_correction.csv",header = T,row.names = 1)
data$BH =p.adjust(data$Raw.p,method = "BH")
data$BH

library(ggplot2)
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T,raster=FALSE)+
  NoLegend()+labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))

#group_by_tcr/bcr
scRNA_all@meta.data$tcr[is.na(scRNA_all@meta.data$tcr)] <- "NO TCR"
scRNA_all@meta.data$bcr[is.na(scRNA_all@meta.data$bcr)] <- "NO BCR"
scRNA_all@meta.data$tcrbcr = ifelse(scRNA_all$tcr == "TCR", "TCR",
                                    ifelse(scRNA_all$bcr == "BCR", "BCR","NA"))
scRNA_all$tcrbcr = factor(scRNA_all$tcrbcr, levels = c("BCR","TCR","NA"))
col4 = c("#99C5DF","#CA7D6E","#C4C2C2")
DimPlot(scRNA_all, reduction = "umap", cols = col4, group.by='tcrbcr',raster=FALSE)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank())

#group_by_Sample
col2 = c("#8AC7E3","#0D71A8","#E78C8B","#2E9845","#E61E2A",
         "#EF7C2C","#B6A2C6","#5D428F")
DimPlot(scRNA_all, reduction = "umap", cols = col2, group.by = "orig.ident",raster=FALSE)+
  NoLegend()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank())

scRNA_all@active.ident = factor(scRNA_all@active.ident, levels = c("T_CD4_naive","T_CD4_mem","Treg",
                        "Cycling","T_CD8_naive","T_CD8_mem","T_CD8_tox","NKT","NK",
                        "B_naive","B_mem","Plasma","pDCs","cDCs","Mono_CD16","Mono_CD14"))
scRNA_all$label = scRNA_all@active.ident
library(viridis)
library(ggplot2)
pdf(file="markerBubble.pdf",width=11,height=6)
cluster10Marker=c("CD3D","CD3E","CD4","CD8A",
                  "CCR7","TCF7","LEF1","GPR183","IL7R","CD27","FOXP3","IL2RA",
                  "MKI67","GZMK","PRF1","GZMA","GZMB","KLRC1","NCAM1","CD160",
                  "MS4A1","IGHD","IGHM","COCH","MZB1","XBP1","PRDM1",
                  "LILRA4","CD1C","FCGR3A","CD14","LYZ","S100A8")
DotPlot(object = scRNA_all, features = cluster10Marker)+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()
dev.off()

library(viridis)
library(ggpointdensity)
library(ggplot2)
library(Nebulosa)
scRNA_all_BL = subset(scRNA_all, SampleType == "BL")
scRNA_all_HC = subset(scRNA_all, SampleType == "HC")
data = cbind(Embeddings(object=scRNA_all_BL[['umap']]),FetchData(scRNA_all_BL,'SampleType'))
p = ggplot(data = data, 
           mapping = aes(x = UMAP_1, y = UMAP_2))+
  NoLegend()+
  geom_pointdensity()+
  scale_color_viridis()+ theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p

DefaultAssay(scRNA_all) <- "RNA"
list = read.csv("genelist.csv")
cd_features <- list(list[,1])
Inscore <- AddModuleScore(scRNA_all,
                          features = cd_features,
                          ctrl = 100,
                          name = "CD_Features")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[42] <- 'Inflammatory_Score'
VlnPlot(Inscore,features = 'Inflammatory_Score', 
        pt.size = 0, adjust = 2,group.by = "orig.ident")
library(ggplot2)
Inscore_3M = Inscore
Inscore_3M$Inflammatory_Score = ifelse(Inscore_3M$SampleType == "HC",Inscore_3M$Inflammatory_Score, -0.1)
mydata<- FetchData(Inscore_3M,vars = c("UMAP_1","UMAP_2","Inflammatory_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = Inflammatory_Score))+
  geom_point(size = 1)+
  xlim(-10,15)+
  ylim(-12,15)+
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

data = read.csv("Volcano_data.csv")
colnames(data) = c("label","qvalue","log2FoldChange")
library(ggrepel)
ggplot(data=data, aes(x=log2FoldChange, y=qvalue)) +
  geom_point(alpha = 0.8, size = 2.0) +
  #scale_colour_manual(name="",values=alpha(mycolor,0.7)) +
  labs(x="log2 fold change",  y="-log10 FDR q", title="HC vs NMO") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  #theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_hline(yintercept=-log10(0.05),linetype=2, color = 'black', size = 0.5) +
  geom_vline(xintercept=c(0),linetype=2, color = 'black', size = 0.5) +
  geom_text_repel(data=data,aes(x=log2FoldChange, y=qvalue,label=label),size=4) +
  xlim(-2,2)+
  ylim(0,200)+
  theme(plot.title = element_text(hjust = 0.5,size = 10, face = "bold"))

#cellchat
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
table(scRNA_all$SampleType)
table(scRNA_all$label)
Idents(scRNA_all) = "label"
scRNA_all_hc = subset(scRNA_all, SampleType == "HC")
scRNA_all_bl = subset(scRNA_all, SampleType == "BL")

data.input = scRNA_all_bl@assays$RNA@data
meta = scRNA_all_bl@meta.data
unique(meta$label)
meta$label = droplevels(meta$label, exclude = setdiff(levels(meta$label),unique(meta$label)))
levels(meta$label)
cellchat <- createCellChat(object = scRNA_all_bl, meta = meta, group.by = "label")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
# future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# netVisual_chord_gene(cellchat, sources.use = c(1:6), targets.use = c(1:6), slot.name = "netP", legend.pos.x = 10)#> Note: The second link end is drawn out of sector ' '.#> Note: The first link end is drawn out of sector 'MIF'.#> Note: The second link end is drawn out of sector ' '.#> Note: The first link end is drawn out of sector 'CXCL '.
netVisual_bubble(cellchat, sources.use = c(10), targets.use = c(1:9,11), remove.isolate = FALSE)#> Comparing communications on a single object
netVisual_bubble(cellchat, sources.use = c(1:9,11), targets.use = c(10), remove.isolate = FALSE)#
p = netVisual_bubble(cellchat, sources.use = c(1:9,11), targets.use = c(10), remove.isolate = FALSE)#> Comparing communications on a single object
data = p$data
#write.csv(data, "data_BL_bubbleplot.csv")
saveRDS(cellchat, "1118cellchat_BL.rds")

# data.dir = './cellchat_pbmc_ppb'
# dir.create(data.dir)
# setwd(data.dir)
cellchat.HC = readRDS("1118cellchat_HC.rds")
cellchat.BL = readRDS("1118cellchat_BL.rds")
object.list = list(HC = cellchat.HC, BL = cellchat.BL)
cellchat = mergeCellChat(object.list, add.names = names(object.list))
netVisual_bubble(cellchat, sources.use = c(10,11,12), targets.use = c(1:9,13:16), comparison = c(1,2), angle.x = 45)#> Comparing communications on a single object
p = netVisual_bubble(cellchat, sources.use = c(10,11,12), targets.use = c(1:9,13:16), comparison = c(1,2), angle.x = 45)#> Comparing communications on a single object
data = p$data
write.csv(data, "data_HC_BL_BPPB_other_bubbleplot.csv")

sub_B = subset(scRNA_all, label == "B_naive"|label == "B_mem"|label == "Plasma")
library(viridis)
library(ggplot2)
pdf(file="markerBubble.pdf",width=11,height=6)
cluster10Marker=c("ADGRE5","CLEC2B","COL9A3","HLA-DPA1","HLA-DRB1","HLA-E","ITGB2","MIF","NAMPT")
DotPlot(object = sub_B, features = cluster10Marker)+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()+coord_flip()
dev.off()

p = DotPlot(object = sub_B, features = cluster10Marker)+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()+coord_flip()
data1 = p$data
write.csv(data1, "B_subcluster.csv")

sub_other = subset(scRNA_all, label != "B_naive")
sub_other = subset(sub_other, label != "B_mem")
sub_other = subset(sub_other, label != "Plasma")
library(viridis)
library(ggplot2)
pdf(file="markerBubble.pdf",width=11,height=6)
cluster10Marker=c("ITGA5","ITGB1","CD74","CXCR4",
                  "CD44","ICAM1","CD4","KLRK1","KLRB1","CD55")
DotPlot(object = sub_other, features = cluster10Marker)+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()+coord_flip()
dev.off()

p = DotPlot(object = sub_other, features = cluster10Marker)+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()+coord_flip()
data2 = p$data
write.csv(data2, "Other_subcluster.csv")

sub_CD16 = subset(scRNA_all, label == "cDCs")
Idents(sub_CD16) = "SampleType"
sce.markers <- FindAllMarkers(object = sub_CD16,
                              only.pos = TRUE,
                              min.pct = 0,
                              logfc.threshold = 0)
write.csv(sce.markers,"sce.markers_cDCs.csv")

rt = read.csv(file = '1123B_subcluster.csv',header=T,sep=",")
head(rt)
rt$id = factor(rt$id, levels = c("B_naive","B_mem","Plasma"))
ggplot(rt,aes(x=id,y=features.plot)) + 
  geom_point(aes(size=pct.exp,color=avg_log2FC))+
  scale_size(rang=c(1,10))+
  scale_color_gradientn(values = seq(0,1,0.25),
                        limits=c(0,0.5),
                        colours = c('#D3D0D5',"#C5AFE0","#AA86EA","#3C1DFF","#502CFC","#2710FC"))+
  labs(
    color=expression(avg_log2FC),
    size="pct.exp",
    x="id"
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_bw()+
  theme(
    axis.text.x=element_text(angle=90,hjust=1),
    axis.text.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size=rel(1.5)),
    axis.title.y = element_blank()
  )
rt = read.csv(file = '1123Other_subcluster_selected.csv',header=T,sep=",")
head(rt)
rt$id = factor(rt$id, levels = c("T_CD4_naive","T_CD4_mem","Treg","Cycling",
                                 "T_CD8_naive","T_CD8_mem","T_CD8_tox","NKT","NK",
                                 "pDCs","cDCs","Mono_CD16","Mono_CD14"))
rt$features.plot = factor(rt$features.plot, 
                          levels = c("CD55","KLRB1","CD44","CD4-2","CD4-1","KLRK1","ICAM1",
                                 "CD74+CD44","CD74+CXCR4","ITGA5+ITGB1"))
ggplot(rt,aes(x=id,y=features.plot)) + 
  geom_point(aes(size=pct.exp,color=avg_log2FC))+
  scale_size(rang=c(1,10))+
  scale_color_gradientn(values = seq(0,1,0.25),
                        limits=c(0,1.52),
                        colours = c('#D3D0D5',"#C5AFE0","#AA86EA","#3C1DFF","#502CFC","#2710FC"))+
  labs(
    color=expression(avg_log2FC),
    size="pct.exp",
    x="id"
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_bw()+
  theme(
    axis.text.x=element_text(angle=90,hjust=1),
    axis.text.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size=rel(1.5)),
    axis.title.y = element_blank()
  )

sub_myeloid = subset(scRNA_all, label == "pDCs"|label == "cDCs"|label == "Mono_CD16"|label == "Mono_CD14")
table(sub_myeloid$orig.ident)
sub_myeloid_HC_BL_3M = subset(sub_myeloid, SampleType != "1M")
table(sub_myeloid_HC_BL_3M$orig.ident)
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
                      subcategory = "WIKIPATHWAYS") %>%
  dplyr::select(gs_name,gene_symbol)
human_GO_bp_Set = human_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)
DefaultAssay(sub_myeloid_HC_BL_3M) = "RNA"
sub_myeloid_HC_BL_3M = NormalizeData(sub_myeloid_HC_BL_3M)
Idents(sub_myeloid_HC_BL_3M) = "orig.ident"
expr = AverageExpression(sub_myeloid_HC_BL_3M, assays = "RNA", slot = "data")[[1]]
expr = expr[rowSums(expr)>0,]
expr = as.matrix(expr)
gsva.res = gsva(expr, human_GO_bp_Set, method = "ssgsea")
gsva_gobp = data.frame(Genesets = rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva_gobp, "Myeloid_gsva_WIKIPATHWAYS.csv", row.names = F)


library(Seurat)
library(tidyselect)
library(dplyr)
scRNA_B_P = subset(scRNA_all, label == "B cells"|label == "Plasma cells"|label == "Proliferating cells")
scRNA_B_P = subset(scRNA_B_P, bcr == "BCR")
pc.num = 1:30
scRNA_B_P <- RunUMAP(scRNA_B_P, reduction = "harmony", dims = pc.num)
scRNA_B_P <- FindNeighbors(scRNA_B_P, reduction = "harmony", dims = 1:30)
scRNA_B_P <- FindClusters(scRNA_B_P, reduction = "harmony", resolution = 1.5)
DimPlot(scRNA_B_P, label = T,reduction = "umap")
scRNA_B_P <- RunUMAP(scRNA_B_P, reduction = "harmony", dims = 1:30)
scRNA_B_P <- FindNeighbors(scRNA_B_P, reduction = "harmony", dims = 1:30)
scRNA_B_P <- FindClusters(scRNA_B_P, reduction = "harmony", resolution = 1.5)
# scRNA_B_P$tcr[is.na(scRNA_B_P$tcr)] = "NA"
# scRNA_B_P = subset(scRNA_B_P, tcr == "NA")
#DimPlot(scRNA_B_P, label = T,reduction = "tsne")
DimPlot(scRNA_B_P, label = T,reduction = "umap")
scRNA_B_P <- FindClusters(scRNA_B_P, reduction = "harmony", resolution = 3.5)
scRNA_B_P <- RunUMAP(scRNA_B_P, reduction = "harmony", dims = 1:25, min.dist = 2, n.neighbors = 200L)
DimPlot(scRNA_B_P, reduction = "umap",label=T, repel = T)

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
new.cluster.ids <- c("Naive B",
                     "Naive B",
                     "Immature B",
                     "Non-switched memory B",
                     "Naive B",
                     "Immature B",
                     "Naive B",
                     "Switched memory B",
                     "Non-switched memory B",
                     "Naive B",
                     "Switched memory B",
                     "Switched memory B",
                     "Switched memory B",
                     "Non-switched memory B",
                     "Immature B",
                     "Switched memory B",
                     "Naive B",
                     "Plasma cells/Plasmablasts",
                     "Naive B",
                     "Naive B",
                     "Naive B",
                     "Naive B",
                     "Plasma cells/Plasmablasts",
                     "Double negative B",
                     "Naive B",
                     "Double negative B",
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
col = c("#999899","#91C463","#5586AD","#9BC3DF","#FCB863","#F2D2D5")
library(ggplot2)
DimPlot(scRNA_B_P, reduction = "umap",repel = T, group.by = "label", cols = col,pt.size = 0.8)+
  NoLegend()+labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(size=15))+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))


sub_012_pre = subset(scRNA_B_P_HC_BL,orig.ident == "B_018_pre")
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
meta_012_pre = meta_012_pre[,c(1,15,16,37,43,46,48,49)]
for(i in 1:nrow(meta_012_pre)){
  meta_012_pre$IGH[i] = strsplit(meta_012_pre$b_cdr3s_aa[i], ";", fixed= T)[[1]][1]
}
head(meta_012_pre)

clonotype = scRNA_B_P$b_clonotype_id
clonotype_df = data.frame(cell = names(clonotype),type = clonotype,orig.ident = scRNA_B_P$orig.ident)
clonotype_df = na.omit(clonotype_df)
for(i in 1:nrow(clonotype_df)){
  clonotype_df$type[i] = paste(clonotype_df$orig.ident[i],clonotype_df$type[i])
}
celltype = scRNA_B_P$label
celltype_df = data.frame(cell = names(celltype),celltype = celltype)
celltype_clonotype_df = merge(clonotype_df,celltype_df,all.y = T)
clononumber_df = data.frame(type = names(table(celltype_clonotype_df$type)),clononumber = table(celltype_clonotype_df$type))
clononumber_df = clononumber_df[,-2]
colnames(clononumber_df) = c("type","clononumber")
celltype_clonotype_clononumber_df = inner_join(celltype_clonotype_df,clononumber_df)
#celltype_clonotype_clononumber_df = merge(celltype_clonotype_df,clononumber_df,by = "type",all.x = T)
celltype_clonotype_clononumber_df$clonostate = ifelse(celltype_clonotype_clononumber_df$clononumber >= 2, ">=2","<2")
scRNA_B_P@meta.data$cell = rownames(scRNA_B_P@meta.data)
clononumber_df = celltype_clonotype_clononumber_df[,c(1,5)]
clononumber_df[,2] = as.numeric(clononumber_df[,2])
scRNA_B_P@meta.data = inner_join(scRNA_B_P@meta.data,clononumber_df)
#scRNA_B_P@meta.data = merge(scRNA_B_P@meta.data,clononumber_df,by = "cell",all.x = T)
clonostate_df = celltype_clonotype_clononumber_df[,c(1,6)]
scRNA_B_P@meta.data = inner_join(scRNA_B_P@meta.data,clonostate_df)
#scRNA_B_P@meta.data = merge(scRNA_B_P@meta.data,clonostate_df,by = "cell",all.x = T)
rownames(scRNA_B_P@meta.data) = scRNA_B_P@meta.data$cell
library(ggplot2)
# scRNA_B_P_HC = subset(scRNA_B_P, SampleType == "HC")
# scRNA_B_P_HC = subset(scRNA_B_P, SampleType == "HC")
DimPlot(scRNA_B_P, reduction = "umap", group.by = "clonostate",cols = c("#D1D0E4","#463179"),pt.size = 0.8)+
  NoLegend()+labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(size=15))+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))

library(dplyr)
colnames(scRNA_B_P@meta.data)
sub_1 = subset(scRNA_B_P,orig.ident == "B_023_pre")
b = sub_1@meta.data[,c(15,37)]
b = na.omit(b)
duplicate_name = b %>% group_by(b_clonotype_id) %>% summarise(freq = n()) %>% filter(freq > 1) %>% select(b_clonotype_id)
duplicate_data = b[b$b_clonotype_id %in% duplicate_name$b_clonotype_id, ]
non_duplicate_data = b[!b$b_clonotype_id %in% duplicate_name$b_clonotype_id, ]
table(non_duplicate_data$label)
table(duplicate_data$label)
# # Fisher's exact test
# data = read.csv("Fisher_exact_test.csv",header = T,row.names = 1)
# fisher.test(data,alternative = "two.sided")$p.value
# data = read.csv("BH_correction.csv",header = T,row.names = 1)
# data$BH =p.adjust(data$Raw.p,method = "BH")
# data$BH

library(patchwork)
scRNA_B_P$SampleType = ifelse(scRNA_B_P$orig.ident == "B_HC_5"|scRNA_B_P$orig.ident == "B_HC_7"|scRNA_B_P$orig.ident == "B_HC_8", "HC",
                          ifelse(scRNA_B_P$orig.ident == "B_018_pre", "IMNM",
                              ifelse(scRNA_B_P$orig.ident == "B_021_pre"|scRNA_B_P$orig.ident == "B_022_pre", "MG","CIDP")))
table(scRNA_B_P$SampleType)
scRNA_B_P$SampleType = factor(scRNA_B_P$SampleType, levels = c("CIDP","MG","IMNM","HC"))
scRNA_P = subset(scRNA_B_P, label == "Plasma cells/Plasmablasts")
Marker=read.csv("genelist.csv")
DotPlot(object = scRNA_P, features = Marker,group.by = "SampleType",cols = c('#4676B4', '#FDCF58',"D83429"),scale.min = 0,scale.max = 40)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
Marker1=read.csv("immunoglobulin_production.csv")
Marker2=read.csv("immunoglobulin_mediated_immune_response.csv")
Marker3=read.csv("B cell activation.csv")
Marker4=read.csv("B cell receptor signaling pathway.csv")
Marker = c(Marker1,Marker2,Marker3,Marker4)
# cols = c('#4676B4', '#FDCF58',"D83429")
library(ggplot2)
DotPlot(object = scRNA_P, features = Marker,group.by = "SampleType",cols = c('#4676B4', '#FDCF58',"D83429"),scale.min = 0,scale.max = 40)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

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
                      category = "C2",
                      subcategory = "CP:WIKIPATHWAYS") %>%
  dplyr::select(gs_name,gene_symbol)
human_GO_bp_Set = human_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)
scRNA_B_P_BL_3M = subset(scRNA_B_P, SampleType == "BL"|SampleType == "3M")
DefaultAssay(scRNA_B_P_BL_3M) = "RNA"
scRNA_B_P_BL_3M = NormalizeData(scRNA_B_P_BL_3M)
Idents(scRNA_B_P_BL_3M) = "orig.ident"
expr = AverageExpression(scRNA_B_P_BL_3M, assays = "RNA", slot = "data")[[1]]
expr = expr[rowSums(expr)>0,]
expr = as.matrix(expr)
gsva.res = gsva(expr, human_GO_bp_Set, method = "ssgsea")
gsva_gobp = data.frame(Genesets = rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva_gobp, "B_gsva_WIKIPATHWAYS.csv", row.names = F)
library(limma)
gsva_gobp = read.csv("B_gsva_WIKIPATHWAYS.csv",row.names = 1)
group <- factor(c(rep("BL", 5), rep("X3M", 4)), levels = c('BL', 'X3M'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(gsva_gobp)
design
compare <- makeContrasts(X3M - BL, levels=design)
fit <- lmFit(gsva_gobp, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1,n=Inf)
head(Diff)
write.csv(Diff,"diff_B_WIKIPATHWAYS.csv")
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
  ylab('t value of GSVA score, B_BL versus B_3M') +
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
            hjust = 1,color = 'black')
p

library(ggplot2)
data=read.csv(file="diff_B_selected.csv",header=T,row.names=1,sep=",")
head(data)
data$Condition = ifelse(data$logFC>0,'UP','DOWN')
table(data$Condition)
data$label = rownames(data)
library(ggrepel) 
mycolor <- c("#5C75A5","#F5BE2B")
ggplot(data=data, aes(x=logFC, y=abs(t))) +
  geom_point(aes(color = Condition), alpha = 0.8, size = 8) +
  scale_colour_manual(name="",values=alpha(mycolor,0.7)) +
  labs(x="log2FC",  y="abs(t-value)", title="B_BL vs B_3M") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  #theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_hline(yintercept=c(2),linetype=2, color = 'black', size = 0.5) +
  geom_vline(xintercept=c(0),linetype=2, color = 'black', size = 0.5) +
  geom_text_repel(data=data,aes(x=logFC, y=abs(t),label=label),size=2) +
  xlim(-0.12,0.12)+
  ylim(1.8,5.2)+
  theme(plot.title = element_text(hjust = 0.5,size = 10, face = "bold"))

library(Seurat)##https://satijalab.org/seurat/articles/multimodal_vignette.html
library(tidyselect)
library(dplyr)
library(DoubletFinder)
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))
  
  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]
  
  # Only keep the barcode and clonotype columns.
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))
  
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
  
  # Reorder so barcodes are first column and set them as rownames.
  tcr <- tcr[, c(2,1,3)]
  rownames(tcr) <- tcr[,1]
  tcr[,1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}
dir1_1 = c('./data/matrix_11/CT103A_018/')
dir1_2 = c('./data/tcr_bcr_matrix_11/CT103A_018_TCR/')
dir1_3 = c('./data/tcr_bcr_matrix_11/CT103A_018_BCR/')
sample_name1 <- c('CT103A_018')
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
  ## preprocessing
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
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- scRNAlist[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sample14@meta.data$ClusteringResults
  nExp_poi <- round(Doubletrate*nrow(scRNAlist[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  scRNAlist[[i]] <- doubletFinder_v3(scRNAlist[[i]], PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  ## save results
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
table(!is.na(scRNAlist[[1]]$t_clonotype_id),!is.na(scRNAlist[[1]]$b_clonotype_id))

scRNA_harmony <- merge(scRNAlist1[[1]], y = c(scRNAlist1[[2]]))
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
DefaultAssay(scRNA_harmony) <- 'RNA'
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30, min.dist = 0.5, n.neighbors = 80L)
DimPlot(scRNA_harmony, reduction = "wnn.umap", label = T)
DimPlot(scRNA_harmony, reduction = "umap", label = T)

DefaultAssay(scRNA_harmony) <- 'ADT'
FeaturePlot(scRNA_harmony, features = "BCMA-FITC",reduction = "umap",cols = c("lightgrey","#FFD700","#FFA500","#FF8C00","#FF8C00","#FF7F50","#FF4500","#FF0000","#B22222","#8B0000"))
FeaturePlot(scRNA_harmony, features = "CD4",reduction = "umap",cols = c("lightgrey","lightgrey","#FFD700","#FCC40E","#FDA900","#FD7A00","#FD7A00","#FD7A00"))
FeaturePlot(scRNA_harmony, features = "CD8A",reduction = "umap",cols = c("lightgrey","lightgrey","#FFD700","#FCC40E","#FDA900","#FD7A00","#FD7A00","#FD7A00"))

DefaultAssay(scRNA_harmony) <- 'RNA'
FeaturePlot(scRNA_harmony, features = "CD8A",reduction = "umap",raster=F)
FeaturePlot(scRNA_harmony, features = "CD4",reduction = "umap",raster=F,cols = c("lightgrey","#A17CEB","#6A40F7","#0000FF","#00008B","#191970"))
FeaturePlot(scRNA_harmony, features = "CT103A-CAR",reduction = "umap",raster=F,cols = c("lightgrey","#6A40F7","#0000FF"))

DefaultAssay(scRNA_harmony) <- 'ADT'
sub_018_IP = subset(scRNA_harmony, orig.ident == "CT103A_018")
p = FeaturePlot(sub_018_IP, features = "BCMA-FITC",reduction = "adt.umap",cols = c("lightgrey","#FFD700","#FFA500","#FF8C00","#FF8C00","#FF7F50","#FF4500","#FF0000","#B22222","#8B0000","#4B0082"))
data2 = p$data
table(data2$BCMA.FITC)
data2$cell = rownames(data2)
head(sub_018_IP@meta.data)
sub_018_IP@meta.data$cell = rownames(sub_018_IP@meta.data)
sub_018_IP@meta.data = inner_join(sub_018_IP@meta.data, data2)
rownames(sub_018_IP@meta.data) = sub_018_IP@meta.data$cell
sub_018_IP$CAR = ifelse(sub_018_IP$BCMA.FITC <= 2, "no_car", "car")
table(sub_018_IP$CAR)
sub_018_IP_CAR = subset(sub_018_IP, CAR == "car")
sub_018_IP_noCAR = subset(sub_018_IP, CAR == "no_car")

scRNA_CAR = merge(sub_018_IP_CAR,
                  y = c(sub_019_IP_CAR,sub_021_IP_CAR,sub_022_IP_CAR,sub_023_IP_CAR,
                        scRNA_Hashtag_CAR,sub_023_1m_CAR))
DefaultAssay(scRNA_CAR) <- 'RNA'
table(scRNA_CAR$orig.ident)
scRNA_CAR <- NormalizeData(scRNA_CAR)
scRNA_CAR <- CellCycleScoring(object = scRNA_CAR,
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes,
                                  g1.features = cc.genes$g1.genes)
scRNA_CAR <- FindVariableFeatures(scRNA_CAR)
scRNA_CAR <- ScaleData(scRNA_CAR)
scRNA_CAR <- RunPCA(scRNA_CAR, verbose=FALSE)
library(harmony)
system.time({scRNA_CAR <- RunHarmony(scRNA_CAR, group.by.vars = "orig.ident")})
plot1 <- DimPlot(scRNA_CAR, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA_CAR, ndims=50, reduction="pca")
plotc <- plot1+plot2
plotc
pc.num = 1:30
scRNA_CAR <- RunUMAP(scRNA_CAR, reduction = "harmony", dims = pc.num)
#scRNA_CAR <- RunTSNE(scRNA_CAR, reduction = "harmony", dims = pc.num)
scRNA_CAR <- FindNeighbors(scRNA_CAR, reduction = "harmony", dims = pc.num)
scRNA_CAR <- FindClusters(scRNA_CAR, reduction = "harmony", resolution = 0.8)
DimPlot(scRNA_CAR, reduction = "umap", label=T)
#DimPlot(scRNA_CAR, reduction = "tsne", label=T)
scRNA_CAR@meta.data$tcr = scRNA_CAR@meta.data$t_clonotype_id
scRNA_CAR@meta.data$tcr[!is.na(scRNA_CAR@meta.data$tcr)] <- "TCR"
DimPlot(scRNA_CAR, reduction = "umap", label=T, group.by = "tcr")
scRNA_CAR <- FindClusters(scRNA_CAR, reduction = "harmony", resolution = 2)
DimPlot(scRNA_CAR, reduction = "umap", label=T)
table(scRNA_CAR$orig.ident)
pdf(file="markerBubble_1.pdf",width=30,height=11)
cluster10Marker=c("CD3D","CD3E","CD4","CD8A","MKI67",
                  "CCR7","TCF7","LEF1","SELL","CCR6","NR4A1","NCR3","KLRB1","GPR183","IL7R","CD27",
                  "PRF1","GZMA","GZMB","GZMK","NKG7","FCGR3A","KLRG1",
                  "ICOS","PDCD1","LAG3","HAVCR2","CD200","CTLA4","ENTPD1","ITGAE","EOMES","IRF8",
                  "FOXP3","IL2RA","TRGV9","TRDV2",
                  "NCAM1","KLRF1","KLRC1","CD160","CT103A-CAR")
DotPlot(object = scRNA_CAR, features = cluster10Marker)
dev.off()
pdf(file="markerBubble_2.pdf",width=15,height=11)
cluster10Marker=c("MS4A1","IGHD","NEIL1","CD27","CD38","ITGAX","TBX21","XBP1","PRDM1","IGHM","COCH","MKI67","IGHA1","IGHG1","MZB1") 
DotPlot(object = scRNA_CAR, features = cluster10Marker,col.min = -1)
dev.off()
pdf(file="markerBubble_3.pdf",width=15,height=11)
cluster10Marker=c("LYZ","CD14","CD83","FCGR3A","C1QA","C1QB","C1QC","CSF1R","TREM2","APOE","TYROBP","CX3CR1","SLC2A5","P2RY13","P2RY12",
                  "CD1C","LILRA4","PF4","CD68","MARCO","FUT4","ITGAM","MME","CXCR2","SELL") 
DotPlot(object = scRNA_CAR, features = cluster10Marker,col.min = -1)
dev.off()

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)
new.cluster.ids <- c("CD4+ T",
                     "CD8+ T",
                     "CD8+ T",
                     "CD4+ T",
                     "CD4+ T",
                     "CD4+ T",
                     "CD4/CD8 T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD4/CD8 T",
                     "CD4/CD8 T",
                     "CD8+ T",
                     "CD4/CD8 T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD4/CD8 T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD4/CD8 T",
                     "CD4/CD8 T",
                     "CD4+ T",
                     "CD8+ T",
                     "CD4+ T",
                     "CD4/CD8 T",
                     "CD4/CD8 T",
                     "CD8+ T",
                     "CD4/CD8 T",
                     "CD8+ T",
                     "CD4/CD8 T",
                     "CD4/CD8 T",
                     "CD4+ T",
                     "Unknown",
                     "CD4+ T",
                     "CD8+ T")
scRNA_CAR@active.ident <- plyr::mapvalues(
  x = scRNA_CAR@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_CAR@meta.data$label = scRNA_CAR@active.ident
scRNA_CAR = subset(scRNA_CAR, label != "Unknown")
DimPlot(scRNA_CAR,reduction = "umap",label = T, cols = c("#A6D38E","#FCBA71","#C0C0C0"))

library(Seurat)
library(tidyselect)
library(dplyr)
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
scRNA_CAR_CD4 <- FindClusters(scRNA_CAR_CD4, reduction = "harmony", resolution = 2)
scRNA_CAR_CD4 = subset(scRNA_CAR_CD4, seurat_clusters != "19")
scRNA_CAR_CD4 = subset(scRNA_CAR_CD4, seurat_clusters != "23")
scRNA_CAR_CD4 <- RunUMAP(scRNA_CAR_CD4, reduction = "harmony", dims = 1:30, min.dist = 1, n.neighbors = 150L)
scRNA_CAR_CD4 <- FindNeighbors(scRNA_CAR_CD4, reduction = "harmony", dims = pc.num)
scRNA_CAR_CD4 <- FindClusters(scRNA_CAR_CD4, reduction = "harmony", resolution = 2)
DimPlot(scRNA_CAR_CD4, reduction = "umap", label=T)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
new.cluster.ids <- c("CD4+ Tem",
                     "CD4+ Tcm",
                     "CD4+ Tcm",
                     "CD4+ Tcm",
                     "Proliferating T cells",
                     "CD4+ Tn",
                     "Proliferating T cells",
                     "CD4+ Tn",
                     "Proliferating T cells",
                     "Proliferating T cells",
                     "Proliferating T cells",
                     "Proliferating T cells",
                     "Proliferating T cells",
                     "Unknown",
                     "Proliferating T cells",
                     "Proliferating T cells",
                     "Treg",
                     "CD4+ Tem",
                     "Proliferating T cells",
                     "Proliferating T cells",
                     "Proliferating T cells",
                     "CD4+ Tcm",
                     "Proliferating T cells")
scRNA_CAR_CD4@active.ident <- plyr::mapvalues(
  x = scRNA_CAR_CD4@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)
scRNA_CAR_CD4@meta.data$label = scRNA_CAR_CD4@active.ident
scRNA_CAR_CD4 = subset(scRNA_CAR_CD4, label != "Unknown")
DimPlot(scRNA_CAR_CD4, reduction = "umap", label=T)
scRNA_CAR_CD4@active.ident = factor(scRNA_CAR_CD4@active.ident, levels = c("CD4+ Tn","CD4+ Tcm","CD4+ Tem","Treg","Proliferating T cells"))
scRNA_CAR_CD4$label = scRNA_CAR_CD4@active.ident
table(scRNA_CAR_CD4$orig.ident)
table(scRNA_CAR_CD4$SampleType)

DimPlot(scRNA_CAR_CD4, reduction = "umap", label=T,group.by = "label",cols = c("#AA90AA","#AACDC5","#7EB793","#3383BC","#FCB863"))+NoLegend()
scRNA_CAR_CD4$SampleType = factor(scRNA_CAR_CD4$SampleType, levels = c("IP","1M"))
col = c("#A4C8DB","#166BA0")
DimPlot(scRNA_CAR_CD4, label = T, reduction = "umap", group.by = "SampleType",cols = col)+NoLegend()

library(viridis)
library(ggplot2)
pdf(file="markerBubble_CD4CAR.pdf",width=8,height=2.5)
cluster10Marker=c("CCR7","TCF7","LEF1","IL7R","CD27","SELL","CD44","FOXP3","IL2RA","MKI67")
DotPlot(object = scRNA_CAR_CD4, features = cluster10Marker)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()
dev.off()

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
DimPlot(scRNA_CAR_CD4, label = T, reduction = "umap", group.by = "clonostate",cols = c("#D1D0E4","#463179"))+NoLegend()

table(scRNA_CAR_CD4$SampleType)
table(scRNA_CAR_CD4$t_clonotype_id)
scRNA_CAR_CD4@meta.data$clonolabel = ifelse(scRNA_CAR_CD4$SampleType != "1M",'Other Timepoints',
                                            ifelse(scRNA_CAR_CD4$t_clonotype_id == "clonotype1" ,'Clonotype1',
                                                   ifelse(scRNA_CAR_CD4$t_clonotype_id == "clonotype2" ,'Clonotype2',
                                                          ifelse(scRNA_CAR_CD4$t_clonotype_id == "clonotype3" ,'Clonotype3',
                                                                 ifelse(scRNA_CAR_CD4$t_clonotype_id == "clonotype4" ,'Clonotype4','1M') )
                                                   )))
scRNA_CAR_CD4$clonolabel = factor(scRNA_CAR_CD4$clonolabel,
                                  levels = c("1M","Other Timepoints",
                                             "Clonotype1","Clonotype2",
                                             "Clonotype3","Clonotype4"))
table(scRNA_CAR_CD4$clonolabel)
col = c("#B2CFE5","#D2D1D2","#8A2A46","#D97351","#5384B6","#61549F")
library(ggplot2)
DimPlot(scRNA_CAR_CD4, reduction = "umap", group.by = "clonolabel",cols = col,pt.size = 1)+ NoLegend()

DefaultAssay(scRNA_CAR_CD4) <- "RNA"
list = read.csv("genelist.csv")
cd_features <- list(list[,1])
Inscore <- AddModuleScore(scRNA_CAR_CD4,
                          features = cd_features,
                          ctrl = 100,
                          name = "CD_Features")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[83] <- 'Proliferation_Score'
VlnPlot(Inscore,features = 'Proliferation_Score',pt.size = 0, adjust = 2,group.by = "orig.ident")+NoLegend()
p = VlnPlot(Inscore,features = 'Proliferation_Score', 
            pt.size = 0, adjust = 2,group.by = "orig.ident")
p
data = p$data

library(COSG)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)
Idents(scRNA_CAR_CD4) = "SampleType"
table(scRNA_CAR_CD4$SampleType)
# sample = subset(scRNA_CAR_CD4,downsample = 5000)
marker_cosg <- cosg(
  scRNA_CAR_CD4,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=30)
marker_cosg$names
gene = c(marker_cosg$names[,1],marker_cosg$names[,2])
scRNA_CAR_CD4 = NormalizeData(scRNA_CAR_CD4)
all.genes = rownames(scRNA_CAR_CD4)
scRNA_CAR_CD4 = ScaleData(scRNA_CAR_CD4, features = all.genes)
mat <- GetAssayData(scRNA_CAR_CD4, slot = "scale.data")
cluster_info <- sort(scRNA_CAR_CD4$orig.ident)
mat <- as.matrix(mat[gene, names(cluster_info)])
SampleType_info = scRNA_CAR_CD4$SampleType
SampleType_df = data.frame(cell = names(SampleType_info),SampleType = SampleType_info)
cluster_df = data.frame(cell = names(cluster_info),cluster = cluster_info)
col_df = inner_join(cluster_df, SampleType_df)
rownames(col_df) = col_df$cell
col_df = col_df[,-1]
table(col_df$cluster)
table(col_df$SampleType)
colours <- list(
  "SampleType"=c("IP"="#A4C8DB",
                 "1M"="#166BA0"))
col_anno = HeatmapAnnotation(df=col_df, 
                             which="col", 
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
