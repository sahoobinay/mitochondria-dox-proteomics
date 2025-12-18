# 🧬 Proteomics Analysis: Mitochondrial Biogenesis Modulates Doxorubicin-Induced Cardiotoxicity

This repository contains the full R-based analysis pipeline and publication-quality figures for our study investigating how mitochondrial biogenesis influences doxorubicin (DOX)-induced cardiotoxicity in neonatal rat ventricular myocytes (NRVMs).

---

## 📖 Project Summary

Doxorubicin is a potent chemotherapeutic agent known for its dose-dependent cardiotoxic effects. Neonatal rat ventricular myocytes (NRVMs) were cultured for either 3 or 6 days to represent low and high mitochondrial biogenesis states, respectively, and treated with 0, 2, or 4 µM DOX. Using global proteomic profiling, we analyzed how mitochondrial content modulates proteome-wide DOX responses.

---

## 📁 Repository Structure

```bash
proteomics-analysis/
├── README.md                   # Project overview and usage instructions
├── LICENSE                     # License information (MIT)
├── data/
│   └── data.xlsx               # Processed input dataset (optional, or replace with a link)
├── scripts/
│   ├── Proteomics.R            # Complete preprocessing, PCA, and heatmap analysis
│   ├── Mitochondrial_Biogenesis_Proteins.R  # Subset and heatmap for mito-biogenesis proteins
│   ├── cluster.R               # Custom hierarchical clustering code
│   ├── top_20_fc.R            # Fold-change analysis between D6_0uM and D3_0uM
│   └── proteomics_dataset_specific.R        # Dataset-specific logic and annotations
├── figures/
│   ├── 1_PCA.pdf               # PCA of proteomics data
│   ├── 2_Heatmap.pdf           # Top 500 most variable proteins heatmap
│   ├── 3_D3_vs_D6_baseline.pdf # FC plot: D6_0uM vs D3_0uM
│   ├── 4_Baseline_FC.pdf       # Top 20 FC plot
│   └── 5_GO_baseline.pdf       # GO enrichment: mitochondrial response genes
```

---

## 🚀 How to Run the Analysis

### 📦 Install Dependencies

In R:
```r
install.packages(c("tidyverse", "pheatmap", "ggrepel", "readxl", "limma", "org.Rn.eg.db", "clusterProfiler", "AnnotationDbi", "ComplexHeatmap"))
```

### 🔁 Run the Pipeline

Run the scripts in this order:

1. `scripts/Proteomics.R` – log-transform, imputation, PCA, and global heatmap
2. `scripts/Mitochondrial_Biogenesis_Proteins.R` – identify and visualize mito-biogenesis proteins
3. `scripts/top_20_fc.R` – log2 fold-change analysis between D3 and D6
4. `scripts/cluster.R` – clustering helper functions
5. `scripts/proteomics_dataset_specific.R` – metadata annotations



---

## 📊 Key Outputs

- **PCA**: Distinguishes Day 3 vs Day 6 proteomes, shows dose-dependent shifts
- **Heatmap**: Top 500 proteins stratified by timepoint and DOX dose
- **Baseline Comparison**: D6_0uM vs D3_0uM log2 FC plot and volcano-style highlight
- **Top 20**: Mitochondrial-related proteins driving differential response
- **GO Enrichment**: Angiogenesis, mitochondrial dynamics, and cellular stress response

---

## 📜 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

## 🙋 Contact & Citation

For questions, collaborations, or reuse, contact:

**Dr. Binay K. Sahoo**  
Postdoctoral Fellow, Stanford School of Medicine  
[sbinay@stanford.edu]

Please cite this work if used in academic or translational projects. Citation coming soon (preprint under preparation).
