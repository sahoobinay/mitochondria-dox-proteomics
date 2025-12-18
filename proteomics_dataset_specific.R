# 1) Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gggenes)
library(pheatmap)
library(limma)
library(tidyr)
# 2) Read the Excel into a tibble
raw_df <- read_excel("data.xlsx", sheet = 1)

# 2.1) Check that “Protein_Name” is spelled exactly like this
print(names(raw_df))
# e.g. [1] "Protein_Name" "3 day 0uM" …

# 2.2) Peek at the first few Protein_Name values
head(raw_df$Protein_Name, 10)
# You should see "D3Z9P5", "P30152", etc.

# 3) Force raw_df into a base data.frame so that rownames are supported
raw_df <- as.data.frame(raw_df)

# 3.1) Now set the rownames to the Protein_Name column:
rownames(raw_df) <- raw_df$Protein_Name

# 3.2) Double‐check that this worked:
head(rownames(raw_df), 10)
# You should see "D3Z9P5", "P30152", "Q6AYT4", … exactly as raw_df$Protein_Name

# 4) Drop the Protein_Name column so that only the six numeric columns remain
expr_df <- raw_df %>% dplyr::select(-Protein_Name)

# 4.1) Confirm you have exactly the six intensity columns
print(names(expr_df))
# [1] "3 day 0uM" "3 day 2uM" "3 day 4uM" "6 day 0uM" "6 day 2uM" "6 day 4uM"

# 5) Now convert expr_df (a data.frame with rownames) into a numeric matrix
#    Because expr_df is a data.frame, as.matrix(...) will preserve its row names.
expr_mat <- expr_df
expr_mat[] <- lapply(expr_mat, as.numeric)  # coerce each column to numeric
expr_mat <- as.matrix(expr_mat)
rownames(expr_mat) <- rownames(expr_df)

# 5.1) Verify that rownames(expr_mat) are the true Protein IDs
head(rownames(expr_mat), 10)
# You should see "D3Z9P5", "P30152", "Q6AYT4", "A0A0G2K5Z4", … (not 1, 2, 3, …)

# From here you can continue with filtering, imputation, log2, median‐centering, etc.

# 3) Filter proteins present in ≥ half (≥3) of the six columns
keep_idx <- rowSums(!is.na(expr_mat)) >= 3
expr_mat_filt <- expr_mat[keep_idx, ]
cat("Kept", nrow(expr_mat_filt), "proteins out of", nrow(expr_mat), "\n")

# 4) Impute remaining NAs with a small constant (1)
expr_mat_filt[is.na(expr_mat_filt)] <- 1

# 5) Log2-transform (+1 offset) and median-center each column
expr_log2 <- log2(expr_mat_filt + 1)
expr_norm <- sweep(expr_log2,
                   2,
                   apply(expr_log2, 2, median),
                   FUN = "-")

# 5.1) Quick boxplot to verify centering

boxplot(
    as.data.frame(expr_norm),
    las  = 2,
    main = "Median-centered log2 Intensities (per sample)",
    ylab = "Intensity"
  )

# 6) PCA on the six samples (columns)# 6) PCA on boxplot()the six samples (columns)
pca_res <- prcomp(t(expr_norm), center = TRUE, scale. = FALSE)
pc_scores <- as.data.frame(pca_res$x)
pc_scores$Sample <- rownames(pc_scores)  # e.g. "3 day 0uM", etc.


ggplot(pc_scores, aes(x = PC1, y = PC2, label = Sample,)) +
  geom_point(size = 3) +
  geom_text_repel(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA of 6 Samples",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2,2], 1), "%)")
  )


dist_mat <- as.matrix(dist(t(expr_norm)))
pheatmap(dist_mat,
         clustering_method = "complete",
         display_numbers    = FALSE,
         main               = "Distance")


group <- factor(colnames(expr_norm),
                levels = c("D3_0uM","D3_2uM","D3_4uM","D6_0uM","D6_2uM","D6_4uM"))

design <- model.matrix(~ 0 + group)  # group: factor(c("D3_0","D3_2",...))
colnames(design) <- levels(group)
contrast.D3_2_vs_D3_0 <- makeContrasts(
  D3_2_vs_D3_0 = D3_2uM - D3_0uM,
  levels       = design
)
fit <- lmFit(expr_norm, design)
contrast.D3_2_vs_D3_0 <- makeContrasts(D3_2uM - D3_0uM, levels = design)
fit2 <- contrasts.fit(fit, contrast.D3_2_vs_D3_0)
logFC_D3_2_vs_D3_0 <- fit2$coefficients[, "D3_2uM - D3_0uM"]

df_fc <- data.frame(
  Protein = rownames(expr_norm),
  logFC   = logFC_D3_2_vs_D3_0,
  row.names = NULL,
  stringsAsFactors = FALSE
)
head(df_fc)
# Volcano plot with ggplot2  
ggplot(df_fc, aes(x = logFC, y = -log10(abs(logFC)))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "D3_2 vs D3_0",
    x     = "log2FC D3_2uM vs D3_0uM",
    y     = "-log10 |log2FC|"
  )


# 7) K-means clustering on proteins across those six columns
k <- 4
set.seed(123)
km_res <- kmeans(expr_norm, centers = k, nstart = 25)
cluster_labels <- km_res$cluster
table(cluster_labels)

# 7.1) Prepare a long-format data.frame for trajectory plotting
traj_df <- expr_norm %>%
  as.data.frame() %>%
  mutate(Protein = rownames(expr_norm),
         Cluster = factor(cluster_labels)) %>%
  pivot_longer(
    cols = -c(Protein, Cluster),
    names_to  = "Condition",
    values_to = "Expression"
  )
traj_df$Condition <- factor(
  traj_df$Condition,
  levels = colnames(expr_norm),
  ordered = TRUE
)

ggplot(traj_df, aes(x = Condition, y = Expression, group = Protein, color = Cluster)) +
  geom_line(alpha = 0.2) +
  stat_summary(
    fun    = mean,
    geom   = "line",
    aes(group = Cluster),
    linewidth   = 1
  ) +
  facet_wrap(~ Cluster, ncol = 1) +
  theme_minimal() +
  labs(
    title = paste0("K-means (k=", k, ") Protein Clusters"),
    x     = "Condition",
    y     = "Median-centered log2(Intensity + 1)"
  ) +
  theme(legend.position = "none")

# 8) Heatmap of the top 100 most variable proteins (by variance across six columns)
protein_var <- apply(expr_norm, 1, var)
N <- 100
if (nrow(expr_norm) < N) N <- nrow(expr_norm)
topN_idx <- order(protein_var, decreasing = TRUE)[1:N]
topN_mat <- expr_norm[topN_idx, ]

# 8.1) Optional: annotate rows by their k-means cluster
annotation_row <- data.frame(
  Cluster = factor(cluster_labels[topN_idx])
)
rownames(annotation_row) <- rownames(topN_mat)

pheatmap(
  expr_mat,
  scale                    = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_row           = annotation_row,
  show_rownames            = FALSE,
  show_colnames            = TRUE,
  fontsize_row             = 6,
  fontsize_col             = 10,
  main                     = paste0("Heatmap: Top ", N, " Variable Proteins")
)

# 9) Compute descriptive log2-fold changes between 6 day vs 3 day at each concentration
#     – 0uM:   LFC_6_vs_3_at_0uM   = ("6 day 0uM" − "3 day 0uM")
#     – 2uM:   LFC_6_vs_3_at_2uM   = ("6 day 2uM" − "3 day 2uM")
#     – 4uM:   LFC_6_vs_3_at_4uM   = ("6 day 4uM" − "3 day 4uM")

lfc_df <- data.frame(
  Protein             = rownames(expr_norm),
  LFC_6_vs_3_0uM      = expr_norm[, "D6_0uM"] - expr_norm[, "D3_0uM"],
  LFC_6_vs_3_2uM      = expr_norm[, "D6_2uM"] - expr_norm[, "D3_2uM"],
  LFC_6_vs_3_4uM      = expr_norm[, "D6_4uM"] - expr_norm[, "D3_4uM"]
)

# 9.1) View the top 10 proteins by absolute LFC at 0uM
lfc_df %>%
  arrange(desc(abs(LFC_6_vs_3_0uM))) %>%
  slice(1:10)

# 9.2) Likewise, top 10 at 2uM or 4uM if desired:
lfc_df %>%
  arrange(desc(abs(LFC_6_vs_3_2uM))) %>%
  slice(1:10)

lfc_df %>%
  arrange(desc(abs(LFC_6_vs_3_4uM))) %>%
  slice(1:10)

# 10) (Optional) Plot a small volcano-style scatter for each concentration pair
#      Here we only have one sample per group—so these "fold changes" are purely descriptive:
lfc_df %>%
  ggplot(aes(x = LFC_6_vs_3_0uM, y = LFC_6_vs_3_2uM)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Scatter of LFC (6 day vs 3 day): 0 µM vs 2 µM",
    x     = "LFC @ 0 µM (6d vs 3d)",
    y     = "LFC @ 2 µM (6d vs 3d)"
  )
# 9) Finally, extract all proteins in “cluster 1”:
proteins_in_cluster1 <- names(hc_clusters)[hc_clusters == 3]
length(proteins_in_cluster1)   # how many proteins in cluster 1?
proteins_in_cluster1            # character vector of those protein IDs
# 1) Install (if needed) and load a lightweight Excel‐writer package
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
library(writexl)

# 2) Turn the character vector into a one‐column data.frame
cluster1_df <- data.frame(Protein_Name = proteins_in_cluster1)

# 3) Write it to an .xlsx file in your working directory
write_xlsx(cluster1_df, path = "cluster3_proteins.xlsx")

# 4) (Optional) Confirm the file exists
list.files(pattern = "cluster1_proteins.xlsx")
# Heatmap of upregulated proteins (no p-values, just showing expression)
# ---- Fold Change Calculation: D6-0uM vs D3-0uM ----
# Identify columns for each group (adjust if your sample naming differs)
d3_0_cols <- grep("D3.*0uM", colnames(expr_log2), value=TRUE)
d6_0_cols <- grep("D6.*0uM", colnames(expr_log2), value=TRUE)

# Calculate mean log2 expression for each group
d3_0_mean <- rowMeans(expr_log2[, d3_0_cols, drop=FALSE])
d6_0_mean <- rowMeans(expr_log2[, d6_0_cols, drop=FALSE])

# Compute log2 fold change (D6-0uM minus D3-0uM)
fc_D6vD3_0uM <- d6_0_mean - d3_0_mean

# Get proteins upregulated in D6-0uM (e.g., log2FC > 1)
up_increased <- names(fc_D6vD3_0uM)[fc_D6vD3_0uM > 1]
# Print or save as needed
print(up_increased)



top_fc <- sort(fc_D6vD3_0uM, decreasing = TRUE)[1:500]
barplot(top_fc, las = 2, cex.names = 0.7, main = "Top 20 FC: D6-0uM vs D3-0uM",
        ylab = "Log2 Fold Change")

# If using UniProt IDs for your proteins:
library(clusterProfiler)
library(org.Rn.eg.db)
entrez_ids <- bitr(up_increased, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
ego <- enrichGO(entrez_ids$ENTREZID, OrgDb = org.Rn.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
dotplot(ego, showCategory = 15) +
  ggtitle('GO Biological Procell Enrichment')

expr_up <- expr_log2[up_increased, , drop=FALSE]

# Draw heatmap
pheatmap(expr_up,
         
         scale = "row",
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of Upregulated Proteins (Day3_vs_Day6")

# End of script