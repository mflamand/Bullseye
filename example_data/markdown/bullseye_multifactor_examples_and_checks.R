# Run from the directory holding both R scripts:
# Rscript bullseye_multifactor_examples_and_checks.R
source("bullseye_functions.R")
set.seed(2026)
meta <- expand.grid(Replicate = seq_len(3), Batch = c("B1", "B2"),
                     Treatment = c("Vehicle", "Drug"), Genotype = c("WT", "KO"),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
meta$Batch <- factor(meta$Batch)
meta$Treatment <- factor(meta$Treatment, levels = c("Vehicle", "Drug"))
meta$Genotype <- factor(meta$Genotype, levels = c("WT", "KO"))
rownames(meta) <- paste0("S", seq_len(nrow(meta)))
ns <- 30L
X <- model.matrix(~Batch + Genotype * Treatment, meta)
beta <- setNames(rep(0, ncol(X)), colnames(X))
beta["(Intercept)"] <- -2.5
beta["GenotypeKO"] <- 0.5
beta["TreatmentDrug"] <- 0.3
beta["GenotypeKO:TreatmentDrug"] <- 0.4
p <- plogis(drop(X %*% beta))
cv <- matrix(sample(150:300, ns*nrow(meta), replace = TRUE), ns)
mt <- matrix(rbinom(length(cv), size = as.vector(cv),
                     prob = rep(p, each = ns)), ns)
rownames(cv) <- rownames(mt) <- paste0("site", seq_len(ns))
colnames(cv) <- colnames(mt) <- rownames(meta)
cv <- as.data.frame(cv); mt <- as.data.frame(mt)

# 1. Additive adjustment for batch; test only genotype.
a <- bullseye(cv, mt, meta, design = ~Batch + Genotype,
               group_var = "Genotype", test_coefficients = "GenotypeKO",
               link = "bin", diagnostics = TRUE, verbose = FALSE)

# 2. Full interaction model: test individual coefficients.
interaction <- bullseye(cv, mt, meta, design = ~Batch + Genotype * Treatment,
                         group_var = "Genotype",
                         test_coefficients = c("GenotypeKO", "GenotypeKO:TreatmentDrug"),
                         link = "quasibinomial", verbose = FALSE)

# 3. Explicit linear contrasts: genotype effects within each treatment,
#    and their difference (the interaction).
tests <- list(
  KO_vs_WT_Vehicle = c(GenotypeKO = 1),
  KO_vs_WT_Drug = c(GenotypeKO = 1, "GenotypeKO:TreatmentDrug" = 1),
  Interaction = c("GenotypeKO:TreatmentDrug" = 1)
)
args <- list(data_cov = cv, data_mut = mt, colData = meta,
              design = ~Batch + Genotype * Treatment, group_var = "Genotype",
              contrasts = tests, link = "bin", verbose = FALSE)
r1 <- do.call(bullseye, c(args, list(numCores = 1, batchSize = 7)))
r2 <- do.call(bullseye, c(args, list(numCores = 2, batchSize = 3, queueFactor = 3)))
stopifnot(isTRUE(all.equal(r1, r2, tolerance = 1e-9)))
for (label in names(tests)) {
  stopifnot(isTRUE(all.equal(r1[[paste0(label, "_FDR")]],
                             p.adjust(r1[[paste0(label, "_pvalue")]], "BH"))))
}
# Compare first site's custom contrast to ordinary glm using identical
# pseudocount-adjusted counts. This validates the matrix covariance calculation.
d <- meta
pc <- exp(mean(log(as.numeric(unlist(cv[1, ], use.names = FALSE)) + 0.5))) * 0.005
d$cov <- as.numeric(unlist(cv[1, ], use.names = FALSE)) + pc
d$mut <- as.numeric(unlist(mt[1, ], use.names = FALSE)) + pc
mod <- suppressWarnings(glm(mut/cov ~ Batch + Genotype * Treatment,
                             weights = cov, family = binomial(), data = d))
v <- setNames(rep(0, ncol(X)), colnames(X))
v[c("GenotypeKO", "GenotypeKO:TreatmentDrug")] <- 1
expected <- sum(v * coef(mod))
expected_se <- sqrt(drop(t(v) %*% vcov(mod) %*% v))
stopifnot(isTRUE(all.equal(r1$KO_vs_WT_Drug_coefficient[1], unname(expected), tolerance=1e-8)),
          isTRUE(all.equal(r1$KO_vs_WT_Drug_SE[1], unname(expected_se), tolerance=1e-8)))

# 4. Missing coverage that destroys site-specific design rank.
cv_bad <- cv; mt_bad <- mt
remove <- which(meta$Batch == "B2")
cv_bad[1, remove] <- 0; mt_bad[1, remove] <- 0
bad <- bullseye(cv_bad, mt_bad, meta, design = ~Batch + Genotype,
                 group_var = "Genotype", link = "bin", verbose = FALSE)
stopifnot(bad$fit_status[1] == "rank_deficient_site_design")

# 5. Optional beta-binomial parity check; may take time.
if (requireNamespace("aod", quietly = TRUE)) {
  args$data_cov <- cv[1:6, , drop = FALSE]
  args$data_mut <- mt[1:6, , drop = FALSE]
  args$link <- "betabin"
  bb1 <- do.call(bullseye, c(args, list(numCores = 1, batchSize = 2)))
  bb2 <- do.call(bullseye, c(args, list(numCores = 2, batchSize = 1)))
  stopifnot(isTRUE(all.equal(bb1, bb2, tolerance = 1e-8)))
  print(attr(bb1, "run_summary"))
  if (!any(is.finite(bb1$KO_vs_WT_Drug_pvalue)))
    warning("No valid beta-binomial tests in smoke check: inspect diagnostics before use.")
} else message("Beta-binomial check skipped: aod is not installed.")
message("All executed assertions passed.")

# Inspect exact coefficient/contrast definitions and fitting outcomes:
print(colnames(attr(r1, "design_matrix")))
print(attr(r1, "contrast_matrix"))
print(attr(r1, "run_summary"))
