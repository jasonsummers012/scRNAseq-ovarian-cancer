# scRNAseq-Ovarian-Cancer

## Project Overview

This repository contains a reproducible analysis pipeline for single-cell RNA sequencing (scRNA-seq) data from human ovarian cancer tissue, using the 10x Genomics “17k Human Ovarian Cancer (FFPE) Single Cell Gene Expression” dataset. The analysis is performed in R with the Seurat package and aims to explore cellular heterogeneity, identify major cell populations, and investigate gene expression patterns within the tumor microenvironment.

---

## Dataset

- **Source:** 10x Genomics
- **Sample:** Human ovarian papillary serous carcinoma (stage III-B), FFPE tissue
- **Cells:** ~17,000 single cells
- **Sequencing Platform:** 10x Genomics Chromium, Illumina NovaSeq 6000
- **Data Format:** Filtered feature/cell matrix (HDF5, `.h5`)
- **License:** CC BY 4.0

---

## Project Structure

scRNAseq_ovarian_cancer/
│
├── data/ # Raw and processed data files (not tracked in Git)
├── scripts/ # R scripts and analysis notebooks
│ └── OvarianCancer_scRNAseq_Analysis.R
├── results/ # Figures, tables, and output
├── docs/ # Documentation and notes
├── .gitignore
├── scRNAseq_ovarian_cancer.Rproj
└── README.md

---

## Analysis Workflow

1. **Data Import and Quality Control**
2. **Normalization and Feature Selection**
3. **Dimensionality Reduction and Clustering**
4. **Cell Type Annotation**
5. **Differential Expression and Visualization**
6. **Reporting and Interpretation**

---

## Requirements

- R (≥ 4.0)
- RStudio (recommended)
- Packages: Seurat, tidyverse, remotes, SeuratData

Install with:
install.packages("Seurat")
install.packages("tidyverse")
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")

---

## How to Run

1. **Clone the repository**
git clone https://github.com/jasonsummers012/scRNAseq-ovarian-cancer.git

text
2. **Download the dataset**
- Place the `.h5` file from 10x Genomics in the `data/` folder.
3. **Open the R project in RStudio**
- Open `scRNAseq_ovarian_cancer.Rproj`
4. **Run the main analysis script**
- Execute `scripts/OvarianCancer_scRNAseq_Analysis.R` step by step.

---

## Results

- Identification of major cell types in ovarian cancer tissue
- Visualization of cell clusters and gene expression patterns
- Insights into tumor microenvironment heterogeneity

---

## License

This project and its analysis scripts are released under the MIT License. The dataset is provided under the CC BY 4.0 license by 10x Genomics.

---

## References

- [10x Genomics: 17k Human Ovarian Cancer (FFPE) Single Cell Gene Expression](https://www.10xgenomics.com/datasets/17k-human-ovarian-cancer-scFFPE)
- [Seurat: Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
