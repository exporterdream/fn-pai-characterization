#!/usr/bin/env Rscript
# DESeq2 differential expression across growth phases (untreated, n=9)
# Pseudo-counts = mean_coverage * gene_length / 75, rounded to integers.
suppressMessages(library(DESeq2))

cov <- read.delim("coverage_mean.tsv", row.names = 1, check.names = FALSE)
ann <- read.delim("gene_annotation.tsv", stringsAsFactors = FALSE)
meta <- read.delim("sample_metadata.tsv", stringsAsFactors = FALSE)
rpkm <- read.delim("rpkm.tsv", row.names = 1, check.names = FALSE)

untreated <- meta$gsm[meta$treatment == "untreated"]
tex <- meta$gsm[meta$treatment == "tex_plus"]

# Filter: RPKM >= 0.5 in >= 3 of 9 untreated samples
keep <- rowSums(rpkm[, untreated] >= 0.5) >= 3
cat("Genes passing filter:", sum(keep), "/", nrow(rpkm), "\n")

run_deseq <- function(samples, label) {
  sub_cov <- cov[keep, samples]
  len <- ann$length[match(rownames(sub_cov), ann$locus_tag)]
  pseudo <- round(sweep(as.matrix(sub_cov), 1, len / 75, `*`))
  storage.mode(pseudo) <- "integer"

  cd <- meta[match(samples, meta$gsm), ]
  cd$phase <- factor(cd$phase, levels = c("E", "M", "S"))
  rownames(cd) <- cd$gsm

  dds <- DESeqDataSetFromMatrix(pseudo, cd, design = ~ phase)
  dds <- DESeq(dds, quiet = TRUE)

  contrasts <- list(M_vs_E = c("phase", "M", "E"),
                    S_vs_E = c("phase", "S", "E"),
                    S_vs_M = c("phase", "S", "M"))
  out <- list()
  for (nm in names(contrasts)) {
    res <- as.data.frame(results(dds, contrast = contrasts[[nm]]))
    res$locus_tag <- rownames(res)
    res <- merge(res, ann[, c("locus_tag", "product", "gene", "PAI")], by = "locus_tag")
    res <- res[order(res$padj), ]
    write.table(res, paste0("de_", label, "_", nm, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
    cat(label, nm, ": ", n_sig, "DE genes (FDR<0.05)\n")
    out[[nm]] <- res
  }
  out
}

de_unt <- run_deseq(untreated, "untreated")
de_tex <- run_deseq(tex, "texplus")

# PAI-focused summary (untreated)
pai <- ann[!is.na(ann$PAI) & ann$PAI != "", ]
pai_de <- data.frame(locus_tag = pai$locus_tag, PAI = pai$PAI, product = pai$product)
for (nm in c("M_vs_E", "S_vs_E", "S_vs_M")) {
  r <- de_unt[[nm]]
  m <- match(pai_de$locus_tag, r$locus_tag)
  pai_de[[paste0(nm, "_log2FC")]] <- r$log2FoldChange[m]
  pai_de[[paste0(nm, "_padj")]] <- r$padj[m]
}
write.table(pai_de, "de_PAI_untreated.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("PAI DE table written.\n")

#!/usr/bin/env Rscript
# Genome-wide Spearman co-expression + PAI enrichment tests
# PRIMARY:     paired within-vs-between Wilcoxon (are PAI genes more connected
#              to each other than to the rest of the genome?)
# SENSITIVITY: permutation test (within-PAI |rho| vs size-matched random sets)
# Both tests run on each individual island AND on the combined set of all PAI genes.
# No bnlearn. No Fisher exact. No hierarchical clustering modules.

rpkm <- read.delim("rpkm.tsv", row.names = 1, check.names = FALSE)
ann <- read.delim("gene_annotation.tsv", stringsAsFactors = FALSE)
meta <- read.delim("sample_metadata.tsv", stringsAsFactors = FALSE)

untreated <- meta$gsm[meta$treatment == "untreated"]
keep <- rowSums(rpkm[, untreated] >= 0.5) >= 3
genes <- rownames(rpkm)[keep]
cat("Filtered genes:", length(genes), "\n")

logx <- log2(as.matrix(rpkm[genes, untreated]) + 1)

# Spearman correlation matrix (genes x genes) on n=9 untreated samples
cmat <- cor(t(logx), method = "spearman")
cat("Correlation matrix:", nrow(cmat), "x", ncol(cmat), "\n")
saveRDS(cmat, "spearman_matrix_untreated.rds")

pai <- ann[!is.na(ann$PAI) & ann$PAI != "", ]
pai_genes <- intersect(pai$locus_tag, genes)
non_pai_genes <- setdiff(genes, pai_genes)
cat("PAI genes passing filter:", length(pai_genes), "/ 75\n")

upper_tri <- function(m) m[upper.tri(m)]

# ============================================================
# PRIMARY: Paired within-vs-between Wilcoxon signed-rank test
# ============================================================
# For each PAI gene: mean |rho| with same-group genes vs mean |rho| with
# all non-PAI genes.  Tests whether genes in a group are more tightly
# connected to each other than to the genomic background.
# Run on each individual island AND on the combined set of all PAI genes.

paired_test <- function(group_genes, group_name) {
  within_vals <- sapply(group_genes, function(g) {
    mean(abs(cmat[g, setdiff(group_genes, g)]))
  })
  between_vals <- sapply(group_genes, function(g) {
    mean(abs(cmat[g, non_pai_genes]))
  })
  diff <- within_vals - between_vals
  wt <- wilcox.test(within_vals, between_vals, paired = TRUE,
                    alternative = "greater")
  data.frame(PAI = group_name, n_genes = length(group_genes),
             mean_within_absrho = round(mean(within_vals), 4),
             mean_between_absrho = round(mean(between_vals), 4),
             mean_diff = round(mean(diff), 4),
             n_genes_within_gt_between = sum(diff > 0),
             wilcoxon_p = signif(wt$p.value, 4))
}

# Individual islands
pai1_genes <- intersect(pai$locus_tag[pai$PAI == "PAI1"], genes)
pai2_genes <- intersect(pai$locus_tag[pai$PAI == "PAI2"], genes)
pai3_genes <- intersect(pai$locus_tag[pai$PAI == "PAI3"], genes)

paired_res <- rbind(
  paired_test(pai1_genes, "PAI1"),
  paired_test(pai2_genes, "PAI2"),
  paired_test(pai3_genes, "PAI3"),
  paired_test(pai_genes, "ALL_COMBINED")   # union of all PAI genes
)
print(paired_res)
write.table(paired_res, "paired_wilcoxon.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("Paired Wilcoxon test written (PRIMARY test).\n")

# ============================================================
# SENSITIVITY: Permutation test (within-group |rho| vs random gene sets)
# ============================================================
# Compares observed mean within-group |rho| to the null distribution from
# 10,000 size-matched random gene sets drawn from all expressed genes.
# Run on each individual island AND on the combined set.

set.seed(42)
perm_test <- function(group_genes, group_name, n_perm = 10000) {
  k <- length(group_genes)
  obs <- mean(abs(upper_tri(cmat[group_genes, group_genes, drop = FALSE])))
  perm <- replicate(n_perm, {
    rg <- sample(genes, k)
    mean(abs(upper_tri(cmat[rg, rg, drop = FALSE])))
  })
  p <- (sum(perm >= obs) + 1) / (n_perm + 1)
  data.frame(PAI = group_name, n_genes = k, obs_mean_absrho = round(obs, 4),
             perm_mean = round(mean(perm), 4), perm_p = signif(p, 4))
}

perm_res <- rbind(
  perm_test(pai1_genes, "PAI1"),
  perm_test(pai2_genes, "PAI2"),
  perm_test(pai3_genes, "PAI3"),
  perm_test(pai_genes, "ALL_COMBINED")
)
print(perm_res)
write.table(perm_res, "permutation_test.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("Permutation test written (SENSITIVITY test).\n")

cat("Done.\n")
