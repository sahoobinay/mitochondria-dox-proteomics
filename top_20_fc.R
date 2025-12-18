# --- Load libraries ---
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# --- Load data ---
df <- read.csv("top20_differentially_expressed_genes_D6.csv", row.names = 1)

# --- Focus on D6 samples only ---
d6 <- df %>% select(starts_with("D6_"))

# --- Compute log2 fold changes relative to 0 µM ---
d6$logFC_2uM <- d6$D6_2uM - d6$D6_0uM
d6$logFC_4uM <- d6$D6_4uM - d6$D6_0uM

# --- Identify top 20 most variable genes ---
d6$meanFC <- rowMeans(d6[, c("logFC_2uM", "logFC_4uM")])
d6$absFC <- abs(d6$meanFC)
top20_genes <- d6 %>%
  arrange(desc(absFC)) %>%
  slice_head(n = 20) %>%
  rownames()

# --- Extract and scale expression values ---
expr_top <- df[top20_genes, c("D6_0uM", "D6_2uM", "D6_4uM")]
expr_scaled <- t(scale(t(expr_top)))

# --- Create clean gene names on rows ---
rownames(expr_scaled) <- top20_genes  # ensures gene labels are used

# --- Create heatmap with labeled rows ---
pheatmap(
  df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Top 20 Differentially Expressed Proteins",
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(255),
  fontsize_row = 10,         # make gene names visible
  fontsize_col = 10,
  border_color = NA,
  angle_col = 45,            # tilt column labels
  treeheight_row = 15,
  treeheight_col = 15
)
write.csv(expr_top, "top20_differentially_expressed_genes_D6.csv",row.names = TRUE)
