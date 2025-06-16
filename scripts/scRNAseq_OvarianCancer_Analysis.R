library(Seurat)
library(SeuratData)
library(ggplot2)

#Create Seurat object
ovarian_cancer_data = Read10X_h5("data/17k_Ovarian_Cancer_scFFPE_count_filtered_feature_bc_matrix.h5")
ovarian_cancer = CreateSeuratObject(
  counts = ovarian_cancer_data,
  min.cells = 3,
  min.features = 200
)

#Filter out cells with >5% mitochondrial counts
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
ovarian_cancer = ScaleData(ovarian_cancer, features = all_genes)