# 1) Load packages
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
library(readxl)
library(dplyr)

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
expr_df <- raw_df %>% select(-Protein_Name)

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

# 6) Example: filter, impute, log2, center
keep_idx <- rowSums(!is.na(expr_mat)) >= 3
expr_mat_filt <- expr_mat[keep_idx, ]
expr_mat_filt[is.na(expr_mat_filt)] <- 1
expr_log2 <- log2(expr_mat_filt + 1)
expr_norm <- sweep(expr_log2,
                   2,
                   apply(expr_log2, 2, median),
                   FUN = "-")

# 6.1) Check again that rownames survived
head(rownames(expr_norm), 10)
# Still the same protein IDs.

# 7) Select top 100 by variance
protein_var <- apply(expr_norm, 1, var)
top100_idx <- order(protein_var, decreasing = TRUE)[1:100]
top100_mat <- expr_norm[top100_idx, ]

# 7.1) Make sure rownames(top100_mat) are your IDs
head(rownames(top100_mat), 10)

# 8) Recompute the hierarchical clustering exactly as pheatmap does
row_dist <- dist(top100_mat, method = "euclidean")
row_hc   <- hclust(row_dist, method = "complete")

# 8.1) Cut the tree into k clusters (e.g. k = 4)
k <- 4
hc_clusters <- cutree(row_hc, k = k)

# 8.2) Now inspect the first few entries of hc_clusters
head(hc_clusters, 10)
# You should see something like:
# D3Z9P5   P30152   Q6AYT4   A0A0G2K5Z4   D3ZA55   P10361   A0A0G2K338  F1M8F6  P09606  A0A0G2JZT6
#       2        3        1           2        4        1           2        3        1            4

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
