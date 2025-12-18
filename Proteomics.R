# Load required packages
library(tidyverse)
library(limma)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(org.Rn.eg.db)  # For rat annotations
library(AnnotationDbi)
library(ggfortify)
library(ggrepel)

# 1. Load and Preprocess Data ------------------------------------------------
# Replace path with your actual file path
data_raw <- readxl::read_excel("data.xlsx", sheet = "Sheet1")

# Set Protein_ID as row names
data <- data_raw %>%
  column_to_rownames("Protein_Name") %>%
  as.matrix()

# Log2 transform
data_log <- log2(data + 1)

# Remove proteins with >50% missing values
keep <- rowMeans(!is.na(data_log)) >= 0.5
data_log <- data_log[keep, ]

# Impute remaining missing values (left-censored)
impute_left_censored <- function(x) {
  q <- quantile(x, probs = 0.01, na.rm = TRUE)
  x[is.na(x)] <- rnorm(sum(is.na(x)), mean = q, sd = 0.3)
  return(x)
}
data_imputed <- t(apply(data_log, 1, impute_left_censored))
write.csv(data_imputed, "data_imputed.csv", row.names = TRUE)
# 2. Design Matrix -----------------------------------------------------------
# Extract sample groups
groups <- c("D3_0uM", "D3_2uM", "D3_4uM", "D6_0uM", "D6_2uM", "D6_4uM")
design <- model.matrix(~ 0 + factor(groups))
colnames(design) <- groups

# 3. PCA ---------------------------------------------------------------------
pca <- prcomp(t(data_imputed), scale. = TRUE)
autoplot(pca, data = data.frame(group = groups), colour = 'group') +
  ggtitle("PCA of Samples") + theme_minimal()

# 4. Heatmap -----------------------------------------------------------------
# Get column names from the imputed matrix
sample_names <- colnames(data_imputed)
top_var <- head(order(apply(data_imputed, 1, var), decreasing = TRUE), 100)
# Derive annotation from sample names
annotation_col <- data.frame(
  Time = ifelse(grepl("^D3", sample_names), "Day3", "Day6"),
  Dose = sub(".*_(\\d+)uM", "\\1", sample_names)
)
rownames(annotation_col) <- sample_names

# Create the heatmap
pheatmap(data_imputed[top_var, ],
         annotation_col = annotation_col,
         show_rownames = TRUE,
         main = "Top 100 Most Variable Proteins")
# 5. Differential Expression (Limma) -----------------------------------------
fit <- lmFit(data_imputed, design)

# Set up contrasts
contrast.matrix <- makeContrasts(
  D3_2vs0 = D3_2uM - D3_0uM,
  D3_4vs0 = D3_4uM - D3_0uM,
  D6_2vs0 = D6_2uM - D6_0uM,
  D6_4vs0 = D6_4uM - D6_0uM,
  D6vsD3_baseline = D6_0uM - D3_0uM,
  D6vsD3_DOX4 = D6_4uM - D3_4uM,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, coef = "D6_4vs0", number = Inf, sort.by = "P")
# 6. Volcano plot ---------------------------------------------------
# Step 1: Compute fold-change
fc_D6_4vs0 <- data_imputed[, "D6_4uM"] - data_imputed[, "D6_0uM"]

# Step 2: Prepare data frame
fc_table <- data.frame(
  Protein = rownames(data_imputed),
  log2FC = fc_D6_4vs0
)
fc_table <- fc_table %>%
  mutate(significant = abs(log2FC) > 1)

# Step 3: Volcano plot with labels
ggplot(fc_table, aes(x = log2FC, y = 0, color = significant)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(data = subset(fc_table, significant),
                  aes(label = Protein),
                  max.overlaps = 25,
                  size = 3,
                  box.padding = 0.4,
                  point.padding = 0.3) +
  theme_minimal() +
  labs(title = "Log2 Fold Change: D6_4uM vs D6_0uM",
       x = "log2 Fold Change",
       y = "")
# 7. Enrichment Analysis (Optional) ------------------------------------------
# Convert protein IDs (assumes UniProt-style IDs)
uniprot_ids <- rownames(results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ])

# Use biomaRt or org.Rn.eg.db to map to Entrez
entrez_ids <- mapIds(org.Rn.eg.db,
                     keys = uniprot_ids,
                     column = "ENTREZID",
                     keytype = "UNIPROT",
                     multiVals = "first")

# Remove NAs
entrez_ids <- na.omit(entrez_ids)

# GO enrichment
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Rn.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment (Biological Processes)")

# 8. Save Differential Expression Table --------------------------------------
write.csv(results, "DE_proteins_D6_4vs0.csv", row.names = TRUE)