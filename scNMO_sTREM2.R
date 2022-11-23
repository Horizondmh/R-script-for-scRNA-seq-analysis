rm(list = ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)
gc()
library(Seurat)
library(tidyselect)
library(dplyr)
library(DoubletFinder)
dir1_1 = c('./data/matrix/C_NMO_1/')
sample_name1 <- c('C_NMO_1')

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
}

scRNA_harmony <- merge(scRNAlist1,y = c(scRNAlist2))
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

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
new.cluster.ids <- c("CD4+ T cells",
                     "CD8+ T cells",
                     "CD4+ T cells",
                     "CD14+ monocytes",
                     "CD8+ T cells",
                     "B cells",
                     "NK cells",
                     "CD4+ T cells",
                     "Unknown",
                     "Macrophages",
                     "CD16+ monocytes",
                     "CD8+ T cells",
                     "cDCs",
                     "CD8+ T cells",
                     "Proliferating cells",
                     "Plasma cells",
                     "CD8+ T cells",
                     "pDCs",
                     "B cells",
                     "CD14+ monocytes",
                     "Macrophages",
                     "Unknown",
                     "CD8+ T cells",
                     "Unknown",
                     "CD8+ T cells",
                     "Unknown",
                     "Macrophages",
                     "CD14+ monocytes")

scRNA_all@active.ident <- plyr::mapvalues(
  x = scRNA_all@active.ident,
  from = current.cluster.ids,
  to = new.cluster.ids)

scRNA_all@meta.data$label = scRNA_all@active.ident
table(scRNA_all$orig.ident)
scRNA_all = subset(scRNA_all, label != "Unknown")
scRNA_all <- RunUMAP(scRNA_all, reduction = "harmony", dims = 1:30, min.dist = 0.3, n.neighbors = 30L)
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T)
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T)

table(scRNA_all$orig.ident)
scRNA_all@active.ident = factor(scRNA_all@active.ident, 
                                      levels = c("CD4+ T cells","CD8+ T cells","NK cells","Proliferating cells","B cells","Plasma cells",
                                                 "CD14+ monocytes","CD16+ monocytes","Macrophages","cDCs","pDCs"))
library(ggplot2)
DimPlot(scRNA_all, reduction = "umap",label=T, repel = T,raster=FALSE)+
  NoLegend()+labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
library(viridis)
library(ggplot2)
pdf(file="markerBubble_all.pdf",width=9,height=4)
cluster10Marker=c("CD3G","CD4","CD8A","KLRC1","MKI67","MS4A1","MZB1","S100A8","CD14","FCGR3A","LYVE1","MRC1",
                  "TREM2","APOE","CD1C","LILRA4")
DotPlot(object = scRNA_all, features = cluster10Marker)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()
dev.off()

library(viridis)
library(ggplot2)
pdf(file="markerBubble_PBMC.pdf",width=9,height=4)
cluster10Marker=c("CD3G","CD4","CD8A","KLRC1","MKI67","MS4A1","MZB1","S100A8","CD14","FCGR3A","LYVE1","MRC1",
                  "TREM2","APOE","CD1C","LILRA4")
DotPlot(object = scRNA_blood, features = cluster10Marker)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()
dev.off()

DefaultAssay(scRNA_csf) <- "RNA"
list = read.csv("genelist.csv")
cd_features <- list(list[,1])
Inscore <- AddModuleScore(scRNA_csf,
                          features = cd_features,
                          ctrl = 100,
                          name = "CD_Features")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[50] <- 'Inflammatory_Score'
VlnPlot(Inscore,features = 'Inflammatory_Score', 
        pt.size = 0, adjust = 2,group.by = "orig.ident")

meta = Inscore@meta.data[,c(48,50)]
write.csv(meta, "Inflammatory_Score_csf.csv")

scRNA_mono = subset(scRNA_all, label == "CD14+ monocytes"|label == "CD16+ monocytes"|
                               label == "Macrophages"|label == "pDCs"|label == "cDCs")
scRNA_mono <- RunUMAP(scRNA_mono, reduction = "harmony", dims = 1:20, min.dist = 0.3, n.neighbors = 30L)
DimPlot(scRNA_mono, label =T)
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
table(scRNA_mono$orig.ident)
scRNA_mono@active.ident = factor(scRNA_mono@active.ident, levels = c("CD14+ monocytes","CD16+ monocytes","Macrophages","cDCs","pDCs"))
scRNA_mono$label = scRNA_mono@active.ident
col = c("#C6CF50","#587F58","#9D0068","#3A6A89","#AA9BC5")
library(ggplot2)
DimPlot(scRNA_mono, reduction = "umap",label=T, repel = T, cols = col,pt.size = 0.8,group.by = "label")+
  NoLegend()+
  labs(x = "UMAP_1", y = "UMAP_2")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
table(scRNA_mono$orig.ident)
sub_PBMC = subset(scRNA_mono, SampleOrig == "blood")
sub_CSF = subset(scRNA_mono, SampleOrig == "csf")
table(sub_PBMC$label)
table(sub_CSF$label)

library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
library(GOplot)
library(org.Hs.eg.db)
library(clusterProfiler)

scRNA_mg = subset(scRNA_mono, label == "Macrophages")
scRNA_mg = subset(scRNA_mg, SampleOrig == "csf")
table(scRNA_mg$orig.ident)
scRNA_mg$SampleType = ifelse(scRNA_mg$orig.ident == "C_HC_5"|scRNA_mg$orig.ident == "C_HC_7"|
                             scRNA_mg$orig.ident == "C_HC_10", "HC","NMO")
scRNA_mg$SampleType = factor(scRNA_mg$SampleType, levels = c("HC","NMO"))

Idents(scRNA_mg) = "SampleType"
sce.markers <- FindAllMarkers(object = scRNA_mg,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

ID = bitr(sce.markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(ID) = c("gene","ENTREZID")
data = inner_join(ID,sce.markers)
geneList = data$ENTREZID

BP <- enrichGO(gene = geneList,  
               keyType = "ENTREZID",  
               OrgDb=org.Hs.eg.db,  
               ont = "BP",  
               pvalueCutoff = 0.01, 
               pAdjustMethod = "fdr", 
               minGSSize = 1,   
               maxGSSize = 500, 
               qvalueCutoff = 0.05,
               readable = TRUE) 
df <- BP@result
df_diff = df[df$qvalue < 0.05,]

go=data.frame(Category = rep("BP",587),ID = df_diff$ID,Term = df_diff$Description, Genes = gsub("/", ", ", df_diff$geneID), adj_pval = df_diff$p.adjust)

data$logFC = ifelse(data$cluster == "HC",-data$avg_log2FC,data$avg_log2FC)
data2 = data[,c(1,9)]
colnames(data2) = c("ID","logFC")
circ = circle_dat(go,data2)
write.csv(circ,"circ.csv")

circ_new = read.csv("circ_new.csv")
GOBubble(circ_new, labels = 1,table.legend =F)

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
scRNA_mg = subset(scRNA_all, label == "Macrophages")
scRNA_mg = subset(scRNA_mg, SampleOrig == "csf")
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

scRNA_mg = subset(scRNA_mono, label == "Macrophages")
scRNA_mg = subset(scRNA_mg, SampleOrig == "csf")
table(scRNA_mg$orig.ident)
scRNA_mg <- RunUMAP(scRNA_mg, reduction = "harmony", dims = 1:30)
scRNA_mg <- FindNeighbors(scRNA_mg, reduction = "harmony", dims = 1:30)
scRNA_mg <- FindClusters(scRNA_mg, reduction = "harmony", resolution = 1)
DimPlot(scRNA_mg, label = T,reduction = "umap")
table(scRNA_mg$orig.ident)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
new.cluster.ids <- c("Macrophage-like cells",
                     "Macrophage-like cells",
                     "Macrophage-like cells",
                     "Macrophage-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Macrophage-like cells",
                     "Macrophage-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Microglia-like cells",
                     "Macrophage-like cells",
                     "Microglia-like cells",
                     "Macrophage-like cells",
                     "Macrophage-like cells",
                     "Macrophage-like cells",
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
library(ggplot2)
DimPlot(scRNA_mg, reduction = "umap",label=T, repel = T, cols = col,pt.size = 0.8,group.by = "label")+
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
cluster10Marker=c("LYVE1","F13A1","MRC1","CD163","MARCO",
                  "CTSD","CSTB","LGALS3","CD9","ITGAX")
DotPlot(object = scRNA_mg, features = cluster10Marker)+ 
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_color_viridis()
dev.off()
sce.markers_mg = FindAllMarkers(object = scRNA_mg,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
