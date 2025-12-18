library(org.Rn.eg.db)

# Your protein IDs (row names of expression matrix)
protein_ids <- rownames(data_imputed)

# Map UniProt → Entrez
uniprot_to_entrez <- mapIds(org.Rn.eg.db,
                            keys = protein_ids,
                            column = "ENTREZID",
                            keytype = "UNIPROT",
                            multiVals = "first")
# Filter out NA values
valid_entrez_ids <- uniprot_to_entrez[!is.na(uniprot_to_entrez)]
# Use GO terms related to mitochondrial biogenesis
go_terms <- c("GO:0007005", "GO:0016042")  # mitochondrion organization, mitochondrial biogenesis

# Get associated Entrez IDs
go_genes <- AnnotationDbi::select(org.Rn.eg.db,
                                  keys = go_terms,
                                  columns = "ENTREZID",
                                  keytype = "GO")
# Filter valid Entrez IDs to those associated with mitochondrial biogenesis
go_genes <- go_genes[go_genes$ENTREZID %in% valid_entrez_ids, ]
mito_entrez_ids <- unique(na.omit(go_genes$ENTREZID))
# Match to your protein IDs
mito_matched <- names(uniprot_to_entrez)[uniprot_to_entrez %in% mito_entrez_ids]

# Subset expression data
data_mito <- data_imputed[rownames(data_imputed) %in% mito_matched, ]
# Column annotation
annotation_col <- data.frame(
  Time = ifelse(grepl("^D3", colnames(data_mito)), "Day3", "Day6"),
  Dose = sub(".*_(\\d+)uM", "\\1", colnames(data_mito))
)
rownames(annotation_col) <- colnames(data_mito)
library(pheatmap)

pheatmap(data_mito,
         scale = "row",
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize = 5,
         main = "Mitochondrial Biogenesis Proteins")

write.csv(data_mito, "mitochondrial_biogenesis_proteins.csv", row.names = TRUE)
