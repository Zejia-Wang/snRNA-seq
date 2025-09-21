library(multtest)
library(Seurat) 
library(dplyr)
library(mindr)
library(tidyverse)
seob_list <- list()
samples <- c('WT1', 'WT2','WT3','hm1','hm2','hm3')
for (sample in samples) {
  
  scrna_data <- Read10X(
    data.dir = str_c("data/", sample))
  
  seob <- CreateSeuratObject(
    counts = scrna_data,
    project = sample,
    min.cells = 3,  
    min.features = 200) 
  
  seob[['sample']] <- sample
  
  seob_list[[sample]] = seob  
}

dim(seob_list$WT1)
dim(seob_list$WT2)
dim(seob_list$WT3)
dim(seob_list$hm1)
dim(seob_list$hm2)
dim(seob_list$hm3)

seob <- merge(x = seob_list[[1]],
              y = seob_list[-1], 
              add.cell.ids = names(seob_list)) 
dim(seob)
seob@meta.data[1:4,] 

view(seob@meta.data) 
seob[["percent.mt"]] <-PercentageFeatureSet(seob,pattern = "^ATMG")
table(seob@meta.data$percent.mt)
seob[["percent.pt"]] <-PercentageFeatureSet(seob,pattern = "^ATCG")
table(seob@meta.data$percent.pt)

VlnPlot(seob,
        features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.pt"),
        group.by  = "sample",
        pt.size = 0.1,
        ncol = 4)

RidgePlot(object = seob, 
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.pt"),
          ncol = 1,
          group.by  = "sample")
滤
seob <- subset(seob,subset= nFeature_RNA >300 & nFeature_RNA <6000 & nCount_RNA <36000 & percent.pt <2 & percent.mt <5)
table(seob$orig.ident) 
dim(seob)

VlnPlot(seob,
        features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.pt"),
        group.by  = "sample",
        pt.size = 0, 
        ncol = 4)

RidgePlot(object = seob, 
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.pt"),
          ncol = 1,
          group.by  = "sample")

plot1 <- FeatureScatter(seob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(seob,feature1 = "nCount_RNA",feature2 = "percent.mt")
library(patchwork)
CombinePlots(plots = list(plot1,plot2))

seob_meta = seob@meta.data
class(seob_meta)
write.csv(seob_meta,"seob_meta.csv", row.names = T)


library(tidyverse)
library(patchwork)
library(Seurat)
table(seob$orig.ident) 
dim(seob)
min(seob$nCount_RNA)
sum(seob$nCount_RNA)

seob_list <- SplitObject(seob, split.by = "sample")
for(i in 1:length(seob_list)){
  seob <- seob_list[[i]]
  seob <- NormalizeData(seob, 
                        normalization.method = "LogNormalize")
  seob <- FindVariableFeatures(seob, 
                               selection.method = "vst", 
                               nfeatures = 3000) 
  seob_list[[i]] <- seob
}
rm(i, seob)

anchors <- FindIntegrationAnchors(
  object.list = seob_list, 
  normalization.method = "LogNormalize")

seob <- IntegrateData(anchorset = anchors)

DefaultAssay(seob) <- "integrated"

seob <- ScaleData(seob, features = rownames(seob))



seob <- RunPCA(seob)
ElbowPlot(seob, ndims = 100)
DimPlot(seob,
        reduction = "pca",
        group.by = "sample")

##t-SNE
seob <- RunTSNE(seob, 
                dims = 1:30) 

##UMAP
seob <- RunUMAP(seob, dims = 1:30)

p1 <- DimPlot(seob, 
              reduction = "pca", # pca, umap, tsne
              group.by = "sample",
              label = F)

p2 <- DimPlot(seob, 
              reduction = "tsne", # pca, umap, tsne
              group.by = "sample",
              label = F)

p3 <- DimPlot(seob, 
              reduction = "umap", # pca, umap, tsne
              group.by = "sample",
              label = F)

p1 + p2 + p3  + plot_layout(guides = "collect") & theme(legend.position = "top")


seob <- FindNeighbors(seob,
                      dims = 1:30)

seob <- FindClusters(seob,
                     resolution = 0.7, 
                     random.seed = 1) 

p1 <- DimPlot(seob, 
              reduction = "pca", # pca, umap, tsne
              group.by  = "seurat_clusters",
              label = T)

p2 <- DimPlot(seob, 
              reduction = "tsne", # pca, umap, tsne
              group.by  = "seurat_clusters",
              label = T)

p3 <- DimPlot(seob, 
              reduction = "umap", # pca, umap, tsne
              group.by  = "seurat_clusters",
              label = T)

p1 + p2 + p3



library(readr)
cell_marker <- read.table(file = "marker.txt", 
                          sep = "\t",
                          header = T) 
cell_marker

seob_list <- SplitObject(seob, split.by = "sample")

p3 <- DotPlot(seob_list[[1]], 
              assay = 'RNA',
              col.min = -1,
              col.max = 6,
              #dot.scale = 6,
              scale.min = 1,
              scale.max = 10,
              features = unique(cell_marker$Cell_Marker)) +
  coord_flip() +
  theme(axis.text = element_text(size = 8,
                                 angle = 0,
                                 hjust = 1))
pdf("marker.pdf",width = 8,height = 7)
p3
dev.off()
#####################
p1 <- DimPlot(seob, 
              reduction = "umap", # pca, umap, tsne
              group.by  = "seurat_clusters",
              split.by = "sample",
              label = T)

pdf("UMAP.pdf",width = 12,height = 7)
p1
dev.off()


p2 <- DotPlot(seob, 
              split.by = "sample",
              assay = 'RNA',
              features = unique(cell_marker$Cell_Marker)) +
  #coord_flip() + 
  theme(axis.text = element_text(size = 8,
                                 angle = 90,
                                 hjust = 1))




pdf("marker.pdf",width = 8,height = 25)
p2
dev.off()


#FeaturePlot
library(RColorBrewer)
p2_2 <- FeaturePlot(seob, 
            reduction = "umap", 
            #assay = "RNA",
            features = c("AT2G22850"), 
            #split.by = "sample",
            label = F,
            order = T,
            pt.size = 0.5) 

pdf("geneID.pdf",width = 15,height = 7)
p2_2
dev.off()

cluster2type <- c()

seob[['cell_type']] = unname(cluster2type[seob@meta.data$seurat_clusters])

p3 <- DimPlot(seob, 
        reduction = "umap", 
        group.by = "cell_type",
        label = TRUE, 
        split.by = "sample",
        pt.size = 1) 

pdf("ano_UMAP.pdf",width = 14,height = 7)
p3
dev.off()

seob_meta = seob@meta.data
class(seob_meta)
write.csv(seob_meta,"seob_meta.csv", row.names = T)




