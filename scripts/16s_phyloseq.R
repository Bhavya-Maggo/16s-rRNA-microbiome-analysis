# ─────────────────────────────────────────────
# Load libraries
# ─────────────────────────────────────────────
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(DESeq2)
library(microbiome)
library(ggtree)
library(ape)
library(RColorBrewer)
library(biomformat)
library(pheatmap)


# ─────────────────────────────────────────────
# File paths
# ─────────────────────────────────────────────
biom_file <- "E:/metagenome/exports/feature-table.biom"
tax_file  <- "E:/metagenome/exports/taxonomy.tsv"
tree_file <- "E:/metagenome/exports/tree.nwk"
meta_file <- "E:/metagenome/metadata.tsv"
out_dir   <- "E:/metagenome/R_results"

# ─────────────────────────────────────────────
# Load BIOM table
# ─────────────────────────────────────────────
otu <- read_biom(biom_file)
otu_mat <- as(biom_data(otu), "matrix")

otu_tab <- otu_table(otu_mat, taxa_are_rows = TRUE)

# ─────────────────────────────────────────────
# Taxonomy parsing
# ─────────────────────────────────────────────
tax_df <- read.table(tax_file, header = TRUE, sep = "\t", row.names = 1)

tax_split <- strsplit(as.character(tax_df$Taxon), ";\\s*")

tax_mat <- t(sapply(tax_split, function(x) {
  x <- gsub("[a-z]__", "", trimws(x)) 
  length(x) <- 7
  x
}))

colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(tax_mat) <- rownames(tax_df)

tax_tab <- tax_table(as.matrix(tax_mat))

# ─────────────────────────────────────────────
# Metadata
# ─────────────────────────────────────────────
meta <- read.table(meta_file, header = TRUE, sep = "\t", row.names = 1, comment.char = "")
meta_tab <- sample_data(meta)

# Grouping column must match the column used in the bash
group_var <- "time" 

# ─────────────────────────────────────────────
# Tree
# ─────────────────────────────────────────────
tree <- read.tree(tree_file)

# ─────────────────────────────────────────────
# Build phyloseq object
# ─────────────────────────────────────────────
ps <- phyloseq(otu_tab, tax_tab, meta_tab, phy_tree(tree))
cat("Initial:", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")

# Extract tree from phyloseq
tree <- phy_tree(ps)

# Get ASV IDs from tree
asv_ids <- tree$tip.label

# Get taxonomy table
tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)

# Create labels: ASV number + Genus  Family  Phylum
asv_labels <- data.frame(
  ASV = asv_ids,
  Label = paste0("ASV", seq_along(asv_ids))
) %>%
  left_join(tax, by = "ASV") %>%
  mutate(Label_tax = ifelse(!is.na(Genus) & Genus != "",
                            paste0(Label, " - ", Genus),
                            ifelse(!is.na(Family) & Family != "",
                                   paste0(Label, " - ", Family),
                                   paste0(Label, " - ", Phylum))))

# Assign these combined labels to the tree tips
tree$tip.label <- asv_labels$Label_tax[match(tree$tip.label, asv_labels$ASV)]

p_tree <- ggtree(tree) +
  geom_tiplab(aes(label=label), size = 4, hjust = 0) +  # left-align labels
  theme_tree2() +
  ggtitle("Phylogenetic Tree (ASV + Taxonomy)")

ggsave(file.path(out_dir, "p_tree.png"), plot = p_tree,
       width = 12, height = 12)

# ─────────────────────────────────────────────
# Filtering
# ─────────────────────────────────────────────
ps <- subset_taxa(ps,
                  !is.na(Phylum) &
                    Kingdom %in% c("Bacteria", "Archaea") &
                    !is.na(Family) &
                    !is.na(Class)
)

cat("After filtering:", ntaxa(ps), "ASVs\n")

# ─────────────────────────────────────────────
# Relative abundance
# ─────────────────────────────────────────────
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))

# ─────────────────────────────────────────────
# Phylum barplot
# ─────────────────────────────────────────────
ps.phylum <- tax_glom(ps.rel, "Phylum")

# Identify rare taxa
abund <- taxa_sums(ps.phylum)
rare_taxa <- names(abund[abund < 0.01 * sum(abund)])

p_phyla <- plot_bar(ps.phylum, fill = "Phylum") +
  facet_wrap(as.formula(paste("~", group_var))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(out_dir, "p_phyla.png"), 
       plot = p_phyla, width = 10, height = 6)

# ─────────────────────────────────────────────
# Rarefaction
# ─────────────────────────────────────────────
set.seed(711)
rare_depth <- 3000   # confirm from alpha-rarefaction.qzv

ps.rare <- rarefy_even_depth(ps, sample.size = rare_depth,
                             rngseed = 711, trimOTUs = TRUE)

# ─────────────────────────────────────────────
# DESeq2
# ─────────────────────────────────────────────
ps_deseq <- ps
sample_data(ps_deseq)$group <- factor(sample_data(ps_deseq)[[group_var]])

dds <- phyloseq_to_deseq2(ps_deseq, ~ group)
dds <- DESeq(dds, fitType = "parametric")
res <- results(dds, alpha = 0.05)

res_df <- as.data.frame(res)
res_df$ASV <- rownames(res_df)

tax_r <- as.data.frame(tax_table(ps_deseq))
tax_r$ASV <- rownames(tax_r)

res_df <- merge(res_df, tax_r, by = "ASV")

sig_df <- res_df[!is.na(res_df$padj) &
                   res_df$padj < 0.05 &
                   abs(res_df$log2FoldChange) > 1, ]

write.csv(sig_df, file.path(out_dir, "DESeq2_sig.csv"), row.names = FALSE)
cat("Significant ASVs:", nrow(sig_df), "\n")

# ─────────────────────────────────────────────
# Volcano plot
# ─────────────────────────────────────────────

# Create labels: Genus > Family > ASV
res_df2 <- merge(res_df, tax, by = "ASV", all.x = TRUE)

res_df2 <- res_df2 %>%
  mutate(label = coalesce(as.character(Genus.y), as.character(Family.y), ASV))

res_df2$status <- "ns"
res_df2$status[res_df2$padj < 0.05 & res_df2$log2FoldChange >  1] <- "up"
res_df2$status[res_df2$padj < 0.05 & res_df2$log2FoldChange < -1] <- "down"

# Top 5 up-regulated
top_up <- res_df2 %>%
  filter(status == "up") %>%
  arrange(-log2FoldChange) %>%
  slice_head(n = 5)

# Top 5 down-regulated
top_down <- res_df2 %>%
  filter(status == "down") %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 5)

# Combine top labels
top_labels <- bind_rows(top_up, top_down)

p_vol <- ggplot(res_df2, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(data = top_labels, 
                  aes(label = label),
                  size = 3,
                  box.padding = 0.3,
                  max.overlaps = 10) +   # prevents too many overlapping labels
  scale_color_manual(values = c(up = "firebrick", down = "steelblue", ns = "grey60")) +
  theme_bw() +
  labs(title = "DESeq2 Volcano Plot", x = "log2 Fold Change", y = "-log10(padj)")

ggsave(file.path(out_dir, "p_vol.png"), plot = p_vol,
       width = 10, height = 6)

# ═══════════════════════════════════════════════════════════════
# HEATMAP
# ═══════════════════════════════════════════════════════════════

if (nrow(sig_df) > 0) {
  
  # Extract counts for significant ASVs only
  sig_asvs <- sig_df$ASV
  ps.sig <- prune_taxa(sig_asvs, ps)
  ps.sig.rel <- transform_sample_counts(ps.sig, function(x) x / sum(x))
  
  mat_sig <- as.matrix(otu_table(ps.sig.rel))
  if (!taxa_are_rows(ps.sig.rel)) mat_sig <- t(mat_sig)
  
  # Label rows with Genus (fallback to Family → Order if missing)
  tax_sig <- as.data.frame(tax_table(ps.sig.rel))
  rownames(mat_sig) <- ifelse(
    !is.na(tax_sig$Genus) & tax_sig$Genus != "",
    paste0(tax_sig$Genus, " (", rownames(mat_sig), ")"),
    ifelse(
      !is.na(tax_sig$Family) & tax_sig$Family != "",
      paste0(tax_sig$Family, " (", rownames(mat_sig), ")"),
      rownames(mat_sig)
    )
  )
  
  meta_df2 <- as(sample_data(ps.sig.rel), "data.frame")
  anno_col3 <- data.frame(
    Group = meta_df2[[group_var]],
    row.names = rownames(meta_df2)
  )
  
  mat_sig_log <- log10(mat_sig + 1e-6)
  
  p_heatmap_deseq <- pheatmap(
    mat_sig_log,
    annotation_col   = anno_col3,
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    scale            = "row",
    color            = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
    border_color     = NA,
    fontsize_row     = 7,
    fontsize_col     = 8,
    main             = paste0("Significant ASVs (DESeq2 padj<0.05, |LFC|>1), n=", nrow(sig_df)),
    filename         = file.path(out_dir, "heatmap.png"),
    width            = 12,
    height           = max(6, nrow(sig_df) * 0.3)   # auto-scale height by ASV count
  )
  
} else {
  cat("No significant ASVs found — skipping DESeq2 heatmap.\n")
}
