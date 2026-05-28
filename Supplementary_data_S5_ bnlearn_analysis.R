# ============================================================================
# CO-EXPRESSION ANALYSIS OF ALL THREE PATHOGENICITY ISLANDS (PAI1, PAI2, PAI3)
# F. nucleatum ATCC 25586 — Supplementary Data S5
#
# Methods:
#   PAI1 (10 genes): Pearson pairwise correlation with BH-FDR correction
#   -r values range from -1 to +1. Positive r indicates co-expression.
#   -BH-FDR correction controls the false discovery rate across 45 tests.\n")
    # -Pearson r > ~0.47 reaches p < 0.05 (uncorrected)=significant

#   PAI2 (32 genes): bnlearn Bayesian network (Hill-Climbing, BIC-Gaussian, 500 bootstrap)
#   PAI3 (33 genes): bnlearn Bayesian network (Hill-Climbing, BIC-Gaussian, 500 bootstrap)
# 0.5 = edge appears in majority of bootstrap samples. Proportion where edge has this specific direction.
#> 0.5 = direction is consistent in majority of samples.
#An edge A -> B suggests that expression of gene A is predictive of expression of gene B, given all other genes in the network. [I don't really conside direction in main analysis that match]

#bnlearn is more like a support, not really main logic chain. See below to see why it is only support!!!
#CORRELATION ≠ CAUSATION: Edges represent statistical dependencies, not necessarily direct physical interactions or causal relationships.
# SAMPLE SIZE: With only 18 samples, results are hypothesis-generating.
#Bootstrap analysis helps but cannot fully compensate for small n.
# Edge directions should be interpreted cautiously. I don't really consider direction, this is why.
# I want to show Edges supported by BOTH bnlearn AND STRING (convergent evidence from two independent methods)
#
# Data Source: GSE161360 (GEO) - F. nucleatum ATCC 25586 RNA-seq
#   18 samples: 3 growth phases (E/M/S) × 2 treatments (untreated/TEX+) × 3 replicates
#   WIG coverage files → strand-specific gene-level coverage → log2(|coverage| + 1)
#
# Significance thresholds:
#   PAI1 Pearson: BH-FDR < 0.05
#   PAI2/PAI3 bnlearn: bootstrap strength > 0.5 AND direction > 0.5
#
# Seed: 42 (all random operations)
# ============================================================================

# STEP 0: SET WORKING DIRECTORY TO THE SCRIPT'S OWN FOLDER
script_dir <- tryCatch({
    dirname(sys.frame(1)$ofile)          # works via source()
}, error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("--file=", "", file_arg)))
    } else {
        getwd()                          # fallback
    }
})
setwd(script_dir)
cat("Working directory set to:", getwd(), "\n")

# ----------------------------------------------------------------------------
# STEP 1: INSTALL AND LOAD REQUIRED PACKAGES
# ----------------------------------------------------------------------------

# Install if not already installed (uncomment as needed):
# install.packages("bnlearn", repos = "https://cloud.r-project.org/")

# library(bnlearn)
# cat("bnlearn version:", as.character(packageVersion("bnlearn")), "\n")

# ----------------------------------------------------------------------------
# STEP 2: LOAD EXPRESSION DATA FOR ALL THREE PAIs
# ----------------------------------------------------------------------------
# Expression matrices were pre-computed from GSE161360 WIG files.
# Any code from a random AI can help you open zip file, or you can open yourself
# extraction pipeline can be done by python, what I do is generated excel file:

#   Rows: 18 samples (named by growth phase, treatment, replicate)
#   Columns: PAI genes (FN#### format)
#   Values: log2(|strand-specific coverage| + 1)



pai1_expr <- read.csv("pai1_expression_log2.csv", row.names = 1)
pai2_expr <- read.csv("pai2_expression_log2.csv", row.names = 1)
pai3_expr <- read.csv("pai3_expression_log2.csv", row.names = 1)

cat("PAI1:", nrow(pai1_expr), "samples x", ncol(pai1_expr), "genes\n")
cat("PAI2:", nrow(pai2_expr), "samples x", ncol(pai2_expr), "genes\n")
cat("PAI3:", nrow(pai3_expr), "samples x", ncol(pai3_expr), "genes\n")

cat("\n")
cat("============================================================================\n")
cat("PART A: PAI1 — PEARSON PAIRWISE CORRELATION\n")
cat("============================================================================\n")

# Compute all pairwise Pearson correlations
genes_pai1 <- colnames(pai1_expr)
n_genes <- length(genes_pai1)
n_pairs <- n_genes * (n_genes - 1) / 2

cat("Computing", n_pairs, "pairwise Pearson correlations for", n_genes, "PAI1 genes...\n")

# Initialize results data frame
pai1_results <- data.frame(
    gene1     = character(n_pairs),
    gene2     = character(n_pairs),
    r         = numeric(n_pairs),
    p_value   = numeric(n_pairs),
    stringsAsFactors = FALSE
)

idx <- 1
for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
        g1 <- genes_pai1[i]
        g2 <- genes_pai1[j]
        ct <- cor.test(pai1_expr[[g1]], pai1_expr[[g2]], method = "pearson")
        pai1_results[idx, ] <- list(g1, g2, ct$estimate, ct$p.value)
        idx <- idx + 1
    }
}

# BH-FDR correction across all 45 tests
pai1_results$p_adj_BH <- p.adjust(pai1_results$p_value, method = "BH")
pai1_results$significant_FDR05 <- pai1_results$p_adj_BH < 0.05

# Sort by r (descending)
pai1_results <- pai1_results[order(-pai1_results$r), ]

cat("\nPAI1 Pearson correlation results:\n")
cat("  Total pairs tested:", n_pairs, "\n")
cat("  Significant at FDR < 0.05:", sum(pai1_results$significant_FDR05), "\n")

cat("\nTop 10 pairs by Pearson r:\n")
print(head(pai1_results[, c("gene1","gene2","r","p_value","p_adj_BH","significant_FDR05")], 10))

cat("\nAll significant pairs (FDR < 0.05):\n")
sig_pai1 <- pai1_results[pai1_results$significant_FDR05, ]
print(sig_pai1[, c("gene1","gene2","r","p_value","p_adj_BH")])

cat("\nFN1885 (hemolysin) pairs:\n")
hem_pairs <- pai1_results[pai1_results$gene1 == "FN1885" | pai1_results$gene2 == "FN1885", ]
print(hem_pairs[, c("gene1","gene2","r","p_value","p_adj_BH","significant_FDR05")])

# Save PAI1 results
# install.packages("openxlsx")
library(openxlsx)
write.xlsx(pai1_results, "pai1_pearson_correlations.xlsx", rowNames = FALSE)
cat("\nPAI1 results saved: pai1_pearson_correlations.xlsx\n")


cat("\n")
cat("============================================================================\n")
cat("PART B: PAI2 — bnlearn BAYESIAN NETWORK\n")
cat("============================================================================\n")

pai2_df <- as.data.frame(pai2_expr)

cat("PAI2 data: ", nrow(pai2_df), "samples x", ncol(pai2_df), "genes\n")
cat("\nExpression summary:\n")
print(summary(as.vector(as.matrix(pai2_df))))

# --- Structure learning ---
cat("\nLearning PAI2 network structure (Hill-Climbing, BIC-Gaussian)...\n")
set.seed(42)
pai2_bn <- hc(pai2_df, score = "bic-g")
cat("PAI2 network: ", length(nodes(pai2_bn)), "nodes,", narcs(pai2_bn), "edges\n")

# --- Bootstrap analysis ---
cat("\nRunning PAI2 bootstrap analysis (500 iterations)...\n")
cat("This may take several minutes...\n")
set.seed(42)
pai2_boot <- boot.strength(
    pai2_df,
    R = 500,
    algorithm = "hc",
    algorithm.args = list(score = "bic-g")
)
cat("PAI2 bootstrap complete.\n")

# Filter significant edges
pai2_sig <- pai2_boot[pai2_boot$strength > 0.5 & pai2_boot$direction > 0.5, ]
pai2_sig <- pai2_sig[order(-pai2_sig$strength), ]

cat("\nPAI2 significant edges (strength > 0.5, direction > 0.5):\n")
cat("  Total significant edges:", nrow(pai2_sig), "\n")
cat("\nTop 10 edges by bootstrap strength:\n")
print(head(pai2_sig, 10))

# Save PAI2 results
write.csv(pai2_boot, "pai2_bnlearn_all_edges.csv", row.names = FALSE)
write.xlsx(pai2_sig,  "pai2_bnlearn_significant_edges.xlsx", row.names = FALSE)
cat("\nPAI2 results saved\n")

# same pipeline, name change, it work

cat("\n")
cat("============================================================================\n")
cat("PART C: PAI3 — bnlearn BAYESIAN NETWORK\n")
cat("============================================================================\n")

pai3_df <- as.data.frame(pai3_expr)

cat("PAI3 data:", nrow(pai3_df), "samples x", ncol(pai3_df), "genes\n")
cat("\nExpression summary:\n")
print(summary(as.vector(as.matrix(pai3_df))))

# --- Structure learning ---
cat("\nLearning PAI3 network structure (Hill-Climbing, BIC-Gaussian)...\n")
set.seed(42)
pai3_bn <- hc(pai3_df, score = "bic-g")
cat("PAI3 network:", length(nodes(pai3_bn)), "nodes,", narcs(pai3_bn), "edges\n")

# --- Bootstrap analysis ---
cat("\nRunning PAI3 bootstrap analysis (500 iterations)...\n")
cat("This may take several minutes...\n")
set.seed(42)
pai3_boot <- boot.strength(
    pai3_df,
    R = 500,
    algorithm = "hc",
    algorithm.args = list(score = "bic-g")
)
cat("PAI3 bootstrap complete.\n")

# Filter significant edges
pai3_sig <- pai3_boot[pai3_boot$strength > 0.5 & pai3_boot$direction > 0.5, ]
pai3_sig <- pai3_sig[order(-pai3_sig$strength), ]

cat("\nPAI3 significant edges (strength > 0.5, direction > 0.5):\n")
cat("  Total significant edges:", nrow(pai3_sig), "\n")
cat("\nTop 10 edges by bootstrap strength:\n")
print(head(pai3_sig, 10))

# Save PAI3 results
write.csv(pai3_boot, "pai3_bnlearn_all_edges.csv", row.names = FALSE)
write.xlsx(pai3_sig,  "pai3_bnlearn_significant_edges.xlsx", row.names = FALSE)
cat("\nPAI3 results saved\n")

# ============================================================================
# STEP 3: SUMMARY
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("SUMMARY OF CO-EXPRESSION ANALYSIS\n")
cat("============================================================================\n")
cat("\n")
cat(sprintf("PAI1 (Pearson correlation, n=18 samples, 10 genes):\n"))
cat(sprintf("  Pairs tested:          %d\n", n_pairs))
cat(sprintf("  Significant (FDR<0.05): %d\n", sum(pai1_results$significant_FDR05)))
cat(sprintf("  Strongest pair:        %s–%s (r = %.3f, FDR = %.2e)\n",
    pai1_results$gene1[1], pai1_results$gene2[1],
    pai1_results$r[1], pai1_results$p_adj_BH[1]))
cat("\n")
cat(sprintf("PAI2 (bnlearn, n=18 samples, 32 genes):\n"))
cat(sprintf("  Significant edges:     %d\n", nrow(pai2_sig)))
cat(sprintf("  Strongest edge:        %s→%s (strength = %.3f)\n",
    pai2_sig$from[1], pai2_sig$to[1], pai2_sig$strength[1]))
cat("\n")
cat(sprintf("PAI3 (bnlearn, n=18 samples, 33 genes):\n"))
cat(sprintf("  Significant edges:     %d\n", nrow(pai3_sig)))
cat(sprintf("  Strongest edge:        %s→%s (strength = %.3f)\n",
    pai3_sig$from[1], pai3_sig$to[1], pai3_sig$strength[1]))

