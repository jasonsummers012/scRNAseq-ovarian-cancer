library(Seurat)
library(SeuratData)
library(ggplot2)
library(tidyverse)
library(rPanglaoDB)

#Create Seurat object
ovarian_cancer_data = Read10X_h5("data/17k_Ovarian_Cancer_scFFPE_count_filtered_feature_bc_matrix.h5")
ovarian_cancer = CreateSeuratObject(
  counts = ovarian_cancer_data,
  min.cells = 3,
  min.features = 200
)

#Filter out poor quality cells
ovarian_cancer[["percent.mt"]] = PercentageFeatureSet(
  ovarian_cancer,
  pattern = "^MT-"
)

ovarian_cancer = subset(
  ovarian_cancer,
  subset = nFeature_RNA > 200 &
  nFeature_RNA < 2500 &
  percent.mt < 5
)

#Normalize the data
ovarian_cancer = NormalizeData(
  ovarian_cancer,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

#Find variable features
ovarian_cancer = FindVariableFeatures(
  ovarian_cancer,
  selection.method = "vst",
  nfeatures = 2000
)

top10 = head(VariableFeatures(ovarian_cancer), 10)
plot1 = VariableFeaturePlot(ovarian_cancer)
plot2 = LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)

ggsave(
  "results/variable_features.png",
  plot2,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

#Scale the data
all_genes = rownames(ovarian_cancer)
ovarian_cancer = ScaleData(
  ovarian_cancer,
  features = all_genes
)

#Perform PCA
ovarian_cancer = RunPCA(
  ovarian_cancer,
  features = VariableFeatures(object = ovarian_cancer)
)

heatmap = DimHeatmap(
  ovarian_cancer,
  dims = 1,
  cells = 500,
  balanced = TRUE,
  fast = FALSE
)

ggsave(
  "results/pca_heatmap.png",
  heatmap,
  width = 20,
  height = 16,
  dpi = 600,
)

elbow_plot = ElbowPlot(ovarian_cancer)

ggsave(
  "results/elbow_plot.png",
  elbow_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

#Perform clustering
ovarian_cancer = FindNeighbors(
  ovarian_cancer,
  dims = 1:12
)

ovarian_cancer = FindClusters(
  ovarian_cancer,
  resolution = 0.5
)

#Perform non-linear dimensional reduction
ovarian_cancer = RunUMAP(
  ovarian_cancer,
  dims = 1:10
)

dim_plot = DimPlot(
  ovarian_cancer,
  reduction = "umap",
)

ggsave(
  "results/dim_plot.png",
  dim_plot,
  width = 10,
  height = 8,
  dpi = 300
)

#Find markers for every clusters
ovarian_markers = FindAllMarkers(
  ovarian_cancer,
  only.pos = TRUE,
  logfc.threshold = 1,
  min.pct = 0.25
)

ovarian_markers = ovarian_markers[ovarian_markers$p_val_adj < 0.05, ]

write.csv(ovarian_markers, "results/cluster_markers.csv")

#Save and load Seurat object
saveRDS(
  ovarian_cancer,
  file = "results/ovarian_cancer_final.rds"
)

ovarian_cancer = readRDS("results/ovarian_cancer_final.rds")

#Identify fibroblasts
fibroblast_markers = c("COL1A1", "COL3A1", "FN1", "ACTA2", "PDGFRA")

fibroblast_samples = getMarkers(include = fibroblast_markers)
human_fibroblast = fibroblast_samples[fibroblast_samples$Specie == "Homo sapiens", ]

fibroblast_clusters = ovarian_markers %>%
  filter(gene %in% fibroblast_markers) %>%
  group_by(cluster) %>%
  summarize(
    n_fibroblast_markers = n(),
    avg_log2FC = mean(avg_log2FC),
    .groups = "drop"
  ) %>%
  arrange(desc(n_fibroblast_markers), desc(avg_log2FC))

vln_plot = VlnPlot(
  ovarian_cancer,
  features = fibroblast_markers,
  group.by = "seurat_clusters",
  pt.size = 0.1
)

ggsave(
  "results/vlnplot_fibroblast_markers_cluster0.png",
  vln_plot,
  width = 10,
  height = 8,
  dpi = 300
)