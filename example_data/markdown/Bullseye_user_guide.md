# Bullseye user guide

## Scope and validation status

This guide describes the latest multifactor edition generated in this conversation,
including underscore-formatted output names and the corrected apply dispatch.
It does not describe the original single-factor script.

Bullseye tests site-level differences in editing using mutation and coverage counts.
It supports fixed-effect additive models, numeric covariates, interactions, and
user-specified one-dimensional linear contrasts. It does not support mixed-effects
syntax, offsets, omnibus multi-degree-of-freedom tests, or automatic marginal means.

## 1. Installation and loading

R's base/recommended packages handle all models except beta-binomial fitting.
Install aod once if you plan to use link = "betabin":

```r
install.packages("aod")
source("bullseye_functions.R")
```

After replacing the script, source it again to replace BOTH Bullseye() and its
helper functions. Loading only the main function can leave incompatible helpers
from an older revision in the session.

## 2. Input objects

### data_cov: coverage counts

A data.frame with one row per site and one column per sample. Row names are unique
site identifiers. An optional numeric column `n` holds a positive site multiplier.

### data_mut: mutation/edit counts

A data.frame with the same site row names in exactly the same order. Sample
columns must contain counts satisfying 0 <= mutation <= coverage. If an `n` column
is supplied, it must match the coverage table's multiplier values.

### colData: sample metadata

A data.frame with unique sample IDs as row names. These names must match count-table
sample columns. Every variable used in the design must exist in this table.
Metadata may be in a different sample order: the script selects count columns in
colData order. Count-table samples absent from colData are not analyzed.

Example layout:

```
coverage:              mutations:             metadata:
       WT1 WT2 KO1 KO2        WT1 WT2 KO1 KO2          Genotype Batch
site_A 100 120 110 130  site_A  5   8  20  26   WT1    WT       B1
site_B  80  90 100 110  site_B  4   3  10  11   WT2    WT       B2
                                             KO1    KO       B1
                                             KO2    KO       B2
```

This tiny layout illustrates structure only; it is not a recommendation for
replication or a demonstration of adequate power.

Use factors for categorical predictors and numeric vectors for continuous
predictors. Character metadata columns are converted to factors. Missing design
metadata are rejected, not silently omitted.

```r
colData$Genotype <- factor(colData$Genotype, levels = c("WT", "KO"))
colData$Batch <- factor(colData$Batch)
# If Age is a continuous numeric variable, leave it numeric.
```

A numeric group_var is rejected: convert a numerically coded grouping variable
to a factor explicitly. Set reference levels before calling Bullseye.

## 3. Quick start: genotype effect adjusted for batch

```r
source("bullseye_functions.R")
colData$Genotype <- factor(colData$Genotype, levels = c("WT", "KO"))
colData$Batch <- factor(colData$Batch)

res <- bullseye(
  data_cov = data_cov,
  data_mut = data_mut,
  colData = colData,
  design = ~ Batch + Genotype,
  group_var = "Genotype",
  test_coefficients = "GenotypeKO",
  link = "betabin",
  numCores = 4,
  batchSize = 50,
  queueFactor = 8,
  diagnostics = TRUE,
  prefilter = FALSE,
  all_values = FALSE
)

head(res)
table(res$fit_status)
attr(res, "run_summary")
```

The model-based genotype effect includes batch adjustment. The WT/KO pooled
editing summaries and KO_log2FoldChange do not.

## 4. Complete function interface

```r
bullseye(
  data_cov, data_mut, colData = NULL,
  filter.col = NULL, filter.with = NULL,
  min = 0.05, filter.sites = NULL, dataset = NULL,
  design = NULL, min.cov = 10, link = "betabin", numCores = 1,
  batchSize = 50L, queueFactor = 8L, diagnostics = FALSE,
  prefilter = FALSE, min_replicates = 2L,
  betabin_control = list(maxit = 2000), betabin_phi_ini = NULL,
  all_values = FALSE, verbose = TRUE, group_var = NULL,
  test_coefficients = NULL, contrasts = NULL
)
```

| Argument | Meaning |
|---|---|
| data_cov, data_mut | Site-by-sample count data.frames described above. |
| colData | Sample metadata; row names identify samples. |
| design | One-sided fixed-effects model formula. |
| group_var | Categorical metadata column for descriptive summaries, output filters and coverage prefilter. Required explicitly for multifactor designs. |
| test_coefficients | Exact model.matrix column names to test. |
| contrasts | Named list of named numeric coefficient weights. Use instead of test_coefficients. |
| link | "betabin", "bin", or "quasibinomial"; this argument selects a family, not an arbitrary link function. All use a logit mean model here. |
| numCores | Positive integer worker count; 1 runs serially. Workers may be capped by the number of batches. |
| batchSize | Maximum sites per worker task; default 50. |
| queueFactor | Maximum queued batches per worker in each bounded queue; default 8. |
| diagnostics | Include per-site timing/convergence/warning columns. |
| prefilter | Enable coverage-only eligibility before fitting; default FALSE. |
| min_replicates | When prefilter is TRUE, require at least this many covered samples in EVERY group_var level; default 2. |
| min.cov | Per-sample coverage cutoff, multiplied by site n. Coverage below cutoff is masked. |
| min | Editing threshold for called samples and filter.sites, on the fraction-times-n scale. |
| filter.sites | NULL, "all" (at least one summary group reaches min), or "called" (every summary group reaches min). Applied to pooled editing AFTER FDR. |
| filter.col, filter.with | Keep metadata rows where colData[[filter.col]] == filter.with. |
| dataset | Optional table with chr, start, end, strand; restrict to matching chr_start_end_strand site IDs. |
| betabin_control | Named optimizer-control list passed to aod::betabin; default maxit 2000. |
| betabin_phi_ini | NULL or one starting beta-binomial dispersion strictly between 0 and 1. |
| all_values | Include detailed count/coverage/called-sample summaries; default FALSE. |
| verbose | Print counts and fit-status summary. |

## 5. Model specification and test selection

### Single factor

```r
design = ~ Genotype
group_var = "Genotype"
test_coefficients = "GenotypeKO"
```

For a single categorical variable, group_var can be omitted. Explicit selection
is still helpful for clarity.

### Additive adjustment

```r
design = ~ Batch + Genotype
design = ~ Batch + Age + Genotype
```

Batch and Age affect the fitted genotype test but not descriptive pooled summaries.
The script rejects globally rank-deficient designs, for example perfect
confounding between batch and genotype.

### Inspect exact coefficient names

```r
colnames(model.matrix(~ Batch + Genotype, data = colData))
```

Select using R's original names, e.g. "GenotypeKO", NOT the formatted output name
"Genotype_KO". The actual names depend on factor levels and contrast coding.
If neither test_coefficients nor contrasts is supplied, all non-intercept design
columns are tested, including nuisance terms such as batch and age.

### Interaction model

```r
colData$Treatment <- factor(colData$Treatment, levels = c("Vehicle", "Drug"))
form <- ~ Batch + Genotype * Treatment
colnames(model.matrix(form, data = colData))
```

Under ordinary treatment coding with WT and Vehicle references:
- GenotypeKO: KO vs WT in Vehicle.
- TreatmentDrug: Drug vs Vehicle in WT.
- GenotypeKO:TreatmentDrug: difference of genotype log-odds effects between Drug
  and Vehicle; it is not the KO-vs-WT effect within Drug by itself.

### Named linear contrasts

```r
tests <- list(
  KO_vs_WT_Vehicle = c(GenotypeKO = 1),
  KO_vs_WT_Drug = c(GenotypeKO = 1, "GenotypeKO:TreatmentDrug" = 1),
  Interaction = c("GenotypeKO:TreatmentDrug" = 1)
)

res <- Bullseye(
  data_cov, data_mut, colData,
  design = ~ Batch + Genotype * Treatment,
  group_var = "Genotype",
  contrasts = tests,
  link = "betabin",
  numCores = 4,
  diagnostics = TRUE
)
```

A contrast is L'beta, with SE = sqrt(L' V L). Omitted coefficient weights are zero.
Use contrasts OR test_coefficients, not both. Custom labels are preserved verbatim.
Each contrast is one test, not an omnibus test of multiple independent effects.

### No-intercept designs

A formula such as ~0 + Genotype is accepted. Its individual coefficients describe
level-specific log odds, not group differences. Use an explicit difference contrast
when comparing groups, e.g. c(GenotypeKO=1, GenotypeWT=-1).

## 6. Choosing a family: implementation differences

- bin: binomial GLM with dispersion fixed at 1.
- quasibinomial: GLM with estimated dispersion.
- betabin: aod beta-binomial fit with one common-over-samples, site-specific
  dispersion parameter. This is not a mixed-effects model with a subject random
  intercept. It is typically more computationally demanding.

The original pseudocount is retained: 0.005 times the geometric mean of
(coverage + 0.5), added equally to mutation and coverage in covered samples.
Descriptive editing summaries do not use this pseudocount. Binomial fits can
record expected noninteger-success warnings because of it.

## 7. Filtering, missingness and eligibility

The order of operations is:
1. Subset samples and requested sites; validate design and counts.
2. Mask coverage below min.cov*n and matching mutation entries to zero.
3. Optionally apply coverage prefiltering.
4. Drop zero-coverage samples separately for each site and fit the full design.
5. Correct valid p-values by BH separately for each requested test.
6. Apply output-only filter.sites.

Missing coverage is masked. Missing mutation counts where coverage passes the
threshold are rejected. Full design metadata must be complete.

With prefilter=TRUE, a site must have min_replicates covered samples in EVERY
summary group. Excluded rows have NA tests and status prefiltered_coverage; they
remain in the returned table unless output filtering removes them. This does not
guarantee sufficient information across every covariate stratum.

At every site, the covered-sample design must be full rank and have positive
residual df. Otherwise ALL its tests return NA. The implementation deliberately
does not recover individual estimable contrasts from a rank-deficient site.

filter.sites="all" means at least one group reaches min; "called" means every
group reaches min. These rules use pooled editing, not model-adjusted estimates.
Changing output filtering does not recalculate FDR. Enabling prefiltering can
change the set of valid tests and therefore FDR.

## 8. Result columns

### Always-present/descriptive columns

| Column | Meaning |
|---|---|
| site | Input site row name. |
| n | Site multiplier; defaults to 1, not replicate count. |
| WT, KO, etc. | Pooled mutation/coverage fraction times n, after masking. |
| KO_log2FoldChange | log2(KO pooled editing / reference pooled editing). Unadjusted. |
| baseEdit | Mean per-sample editing across all groups among samples reaching min; 0 if none. Not reference editing. |
| fit_status | Outcome of fitting and validation. |

Pooled editing is coverage-weighted, whereas baseEdit is an unweighted mean of
called sample fractions. Values are ordinary fractions when n=1; n!=1 changes the
scale and can yield values above 1. A zero denominator can yield NA or an infinite
empirical fold change. The site multiplier cancels in the fold-change ratio.

### Per coefficient or named contrast

| Suffix | Meaning |
|---|---|
| _coefficient | Model-based natural-log odds coefficient or specified linear contrast. |
| _SE | Its standard error on the same scale. |
| _pvalue | One-df chi-square Wald p-value for coefficient/contrast equal to zero. |
| _FDR | BH-adjusted p-value across valid sites for this requested test. |

Factor output prefixes insert an underscore, e.g. Genotype_KO. Interaction
components are treated similarly: Genotype_KO:Treatment_Drug. Internal names are
unchanged. Exact factor-level tokens are reformatted; unusual/non-treatment
contrast encodings can retain R's original names. Ambiguous/colliding labels are
not silently overwritten. Custom contrast labels are unchanged.

For additive genotype contrasts, exp(coefficient) gives the adjusted odds ratio.
Do not confuse it with the editing-fraction ratio. Interaction coefficients are
log ratios of odds ratios. Pooled log2FoldChange remains unadjusted even when the
coefficient adjusts for batch, age, or other covariates.

FDR is corrected separately for each requested coefficient/contrast, NOT jointly
across all contrasts. It includes only valid returned p-values, before output
filtering; excluded/failed tests are NA. Frequent failures require investigation.
Quasibinomial tests also use chi-square Wald p-values, not finite-df t/F tests.

### all_values=TRUE

For every group g:
- g_sum_cov: retained coverage sum.
- g_mut: retained mutation sum.
- g_min_cov: NUMBER OF SAMPLES with positive retained coverage, not minimum coverage.
- g_called_sites: NUMBER OF SAMPLES reaching the editing threshold, not distinct sites.
- g_average_called: mean editing in those called samples, or NA when none.

### diagnostics=TRUE

- fit_code: convergence code (0 means reported convergence).
- fit_evaluations: aod objective evaluation count from its iterations slot; NA for GLMs.
- fit_phi: beta-binomial dispersion; NA for other families in this implementation.
- fit_singular_hessian: aod singular-Hessian flag.
- fit_elapsed_sec: site-level elapsed fitting/preparation time, excluding queue wait.
- fit_warning: signalled warnings/errors.
- fit_optimizer_output: captured printed aod output, including some internally handled warnings.
- fit_n_samples: covered samples retained at the site.
- fit_residual_df: retained samples minus model-matrix columns.

### fit_status values

ok; prefiltered_coverage; insufficient_residual_df;
rank_deficient_site_design; optimizer_not_converged; singular_hessian;
invalid_dispersion; invalid_covariance; partially_testable;
no_testable_contrast; fit_failed.

Status ok means the implemented numerical checks passed, not that the model is
biologically correct or guaranteed to reach a global optimum. Invalid fits do not
produce substituted p-values.

## 9. Metadata attributes

```r
attr(res, "design_matrix")   # global sample-by-coefficient matrix
attr(res, "contrast_matrix") # coefficient-by-requested-test weights
attr(res, "test_name_map")   # internal_name -> output_prefix
attr(res, "run_summary")     # site counts, statuses, valid tests, reference, scope
```

The reference summary group is the first level of group_var. Model references
are set by the actual model-matrix contrast coding; they need not coincide with
simple first-level interpretations under custom coding.

## 10. Parallel performance and diagnostics

```r
res <- Bullseye(data_cov, data_mut, colData,
                design=~Batch+Genotype, group_var="Genotype",
                test_coefficients="GenotypeKO", numCores=4,
                batchSize=50, queueFactor=8, diagnostics=TRUE)

slow <- res[order(res$fit_elapsed_sec, decreasing=TRUE), ]
head(slow[, c("site", "fit_status", "fit_elapsed_sec", "fit_evaluations", "fit_phi")], 20)
table(res$fit_status)
```

PSOCK workers provide cross-platform execution. Small batches are dynamically
scheduled within bounded queues; queues still synchronize at their boundaries.
Full count/editing matrices remain in parent RAM. More workers are not guaranteed
to improve runtime. Benchmark representative data and inspect difficult fits.

Hessian computation is required for the current covariance-based tests and remains
enabled. Lowering maxit can increase failures; it is not a wall-clock timeout.
There are no automatic retries, model fallback or hard fit time limits.

## 11. Inspecting and exporting results

```r
head(res)
attr(res, "test_name_map")

significant <- res[
  !is.na(res$Genotype_KO_FDR) & res$Genotype_KO_FDR < 0.05,
  , drop=FALSE
]

# Descriptive percentage-point difference for ordinary n=1 sites:
# For general n, divide by n to recover the pooled mutation-fraction difference.
res$KO_delta_mutation_pp <- 100 * (res$KO - res$WT) / res$n

write.csv(res, "bullseye_results.csv", row.names=FALSE)
saveRDS(res, "bullseye_results.rds")
```

CSV preserves columns but not the custom design/contrast/run-summary attributes.
RDS preserves the R object and those attributes.

## 12. First-run validation and troubleshooting

Start with a small subset, preserving identical row ordering in both count tables:

```r
ii <- seq_len(min(100L, nrow(data_cov)))
args <- list(data_cov=data_cov[ii,,drop=FALSE],
             data_mut=data_mut[ii,,drop=FALSE], colData=colData,
             design=~Batch+Genotype, group_var="Genotype",
             test_coefficients="GenotypeKO", link="bin", verbose=FALSE)
a <- do.call(Bullseye, c(args, list(numCores=1, batchSize=7)))
b <- do.call(Bullseye, c(args, list(numCores=2, batchSize=5)))
all.equal(a, b, tolerance=1e-8)
```

Leave diagnostics FALSE for this comparison, since timing differs between runs.
Validate beta-binomial fits separately and compare selected sites to direct
model fits with identical masking and pseudocounts.

| Issue | Action |
|---|---|
| 'jobs' of mode 'function' not found | An older dispatcher/helper is loaded. Source the complete latest script; .be_batch must use design_matrix, not X, as its forwarded argument. |
| Unknown test coefficient | Inspect colnames(model.matrix(design, colData)); use internal names without the new display underscore. |
| Global rank deficiency | Check confounding, redundant terms, empty levels and number of samples. |
| Site rank deficiency | Coverage masking removed samples needed for full-model estimation. Inspect site-level coverage and design. |
| Many beta-binomial failures | Enable diagnostics and inspect code, warning text, phi and covariance failures before trusting FDR. |
| Unexpected WT/KO proportions | Check group_var, factor levels, coverage masking and the n multiplier. |
| Output lacks detailed counts | Set all_values=TRUE. |
| Runtime remains high | Inspect per-site time/evaluations; numerical optimization can dominate despite parallel scheduling. |

Do not treat successful execution alone as validation of a new statistical pipeline.
