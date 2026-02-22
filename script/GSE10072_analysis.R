# ==========================================================
# GSE10072 Transcriptomic Analysis
# Adenocarcinoma vs Normal Lung Tissue
# ==========================================================

# -----------------------------
# 1. Load Libraries
# -----------------------------
library(GEOquery)
library(limma)
library(hgu133a.db)
library(AnnotationDbi)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(umap)
library(clusterProfiler)
library(org.Hs.eg.db)

# -----------------------------
# 2. Download Dataset
# -----------------------------
options(download.file.method = "libcurl")
gset <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
ex <- exprs(gset)

# -----------------------------
# 3. Define Groups
# -----------------------------
group_info <- pData(gset)[["source_name_ch1"]]
groups <- make.names(group_info)
gset$group <- factor(groups)

# -----------------------------
# 4. Design Matrix & Limma
# -----------------------------
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

contrast_matrix <- makeContrasts(
  Adenocarcinoma.of.the.Lung - Normal.Lung.Tissue,
  levels = design
)

fit <- lmFit(ex, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

# -----------------------------
# 5. Annotate Probe → Gene
# -----------------------------
probe_ids <- rownames(topTableResults)

gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

# -----------------------------
# 6. Count Significant Genes
# -----------------------------
sig <- topTableResults[
  topTableResults$adj.P.Val < 0.01 &
    abs(topTableResults$logFC) > 1,
]

n_total <- nrow(sig)
n_up <- sum(sig$logFC > 1)
n_down <- sum(sig$logFC < -1)

cat("Total significant:", n_total, "\n")
cat("Upregulated:", n_up, "\n")
cat("Downregulated:", n_down, "\n")

# -----------------------------
# 7. Boxplot
# -----------------------------
group_colors <- ifelse(
  gset$group == "Adenocarcinoma.of.the.Lung",
  "black",
  "red"
)

png("boxplot_GSE10072.png", width = 1000, height = 800)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Expression Distribution",
  ylab = "Expression Value (log2)"
)
legend(
  "topright",
  legend = c("Cancer", "Normal"),
  fill = c("black", "red"),
  bty = "n"
)
dev.off()

# -----------------------------
# 8. Density Plot
# -----------------------------
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

expr_long$Group <- ifelse(
  expr_long$Group == "Adenocarcinoma.of.the.Lung",
  "Cancer",
  "Normal"
)

density_plot <- ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Expression Distribution",
    x = "Expression Value (log2)",
    y = "Density"
  )

ggsave("density_plot_GSE10072.png", density_plot)

# -----------------------------
# 9. UMAP
# -----------------------------
umap_input <- t(ex)
umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

umap_df$Group <- ifelse(
  umap_df$Group == "Adenocarcinoma.of.the.Lung",
  "Cancer",
  "Normal"
)

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot",
    x = "UMAP 1",
    y = "UMAP 2"
  )

ggsave("umap_GSE10072.png", umap_plot)

# -----------------------------
# 10. Volcano Plot
# -----------------------------
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val
)

volcano_data$status <- "Not Significant"

volcano_data$status[
  volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01
] <- "Upregulated (Cancer)"

volcano_data$status[
  volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01
] <- "Downregulated (Cancer)"

volcano_plot <- ggplot(volcano_data,
                       aes(x = logFC,
                           y = -log10(adj.P.Val),
                           color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c(
    "Downregulated (Cancer)" = "blue",
    "Not Significant" = "grey",
    "Upregulated (Cancer)" = "red"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot DEG Kanker Paru",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

ggsave("volcano_plot_GSE10072.png", volcano_plot)

# -----------------------------
# 11. Heatmap (Top 50)
# -----------------------------
topTableResults <- topTableResults[order(topTableResults$adj.P.Val), ]
top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,
  top50$SYMBOL
)

rownames(mat_heatmap) <- gene_label

annotation_col <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)

annotation_col$Group <- ifelse(
  annotation_col$Group == "Adenocarcinoma.of.the.Lung",
  "Cancer",
  "Normal"
)

png("heatmap_GSE10072.png", width = 1000, height = 800)
pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)
dev.off()

# -----------------------------
# 12. Enrichment Analysis
# -----------------------------
sig_genes <- na.omit(sig$SYMBOL)

gene_entrez <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ggsave("GO_GSE10072.png",
       dotplot(ego, showCategory = 10))

ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

ggsave("KEGG_GSE10072.png",
       dotplot(ekegg, showCategory = 10))

# -----------------------------
# 13. Save Results
# -----------------------------
write.csv(topTableResults,
          "GSE10072_DEG_results.csv",
          row.names = FALSE)

write.csv(sig,
          "DEG_GSE10072_significant_only.csv",
          row.names = FALSE)

write.csv(as.data.frame(ego),
          "GO_Enrichment_GSE10072.csv",
          row.names = FALSE)

write.csv(as.data.frame(ekegg),
          "KEGG_Enrichment_GSE10072.csv",
          row.names = FALSE)

cat("Analysis completed successfully.\n")