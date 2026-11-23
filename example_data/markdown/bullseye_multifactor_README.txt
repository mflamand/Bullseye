BULLSEYE MULTIFACTOR EDITION

FILES
bullseye_functions.R: source this file.
bullseye_multifactor_examples_and_checks.R: executable synthetic-data examples and
checks; run from the directory containing both R scripts.

SUPPORTED MODELS
Fixed-effect one-sided formulas: ~Genotype, ~Batch+Genotype,
~Batch+Genotype+Age, ~Batch+Genotype*Treatment, ~0+Genotype.
Numeric covariates, interactions, standard model.matrix contrast coding and named
linear contrasts are supported. No random effects, offsets, or omnibus multi-df
factor tests. Every reported test is a 1-df chi-square Wald test, including for
quasibinomial fits (no finite-df t/F correction).
Categorical variables must be factors (character columns are converted); numeric
variables are continuous. Set factor reference levels explicitly before fitting.

GROUPING VS MODEL
For multifactor models supply group_var, e.g. group_var="Genotype".
It determines only descriptive pooled editing summaries, log2FoldChange,
output-only filter.sites and optional coverage prefiltering. It is not inferred
from formula order. Pooled editing/log2FoldChange are UNADJUSTED. Fitted
coefficients, SE, p-values and FDR reflect the complete design.
With interactions, a main coefficient is conditional on other factor reference
levels / numeric zero; it is not an overall averaged effect. An interaction
coefficient is a difference of log-odds effects, not a simple group difference.

SELECTING TESTS
Default: all model-matrix coefficients except (Intercept), including nuisance
terms and interactions. Prefer test_coefficients=c("GenotypeKO") for targeted
analysis. See colnames(model.matrix(design, colData)) for exact names; formula
factor ordering and contrast coding affect these names.
Alternatively supply contrasts=list(my_test=c(GenotypeKO=1, ...)). Each vector
contains weights on exact model-matrix columns. Missing weights are zero.
Only one of test_coefficients and contrasts may be supplied. Contrasts return
L'beta with SE=sqrt(L' V L). No automatic marginal means are computed.

OUTPUT NAMING CHANGE
Model result prefixes now use full coefficient names or user contrast labels:
GenotypeKO_coefficient, GenotypeKO_SE, GenotypeKO_pvalue, GenotypeKO_FDR.
These replace the old KO_coefficient-style names for statistical results.
Descriptive group columns WT/KO and KO_log2FoldChange retain their group labels.
SE columns are included by default. No adjusted log2FoldChange is fabricated.
Single-factor calls can omit group_var if its variable is categorical.

RETAINED CONTROLS
link = betabin, bin, quasibinomial; no betareg alias.
numCores=1, batchSize=50, queueFactor=8; cross-platform PSOCK dynamic scheduling.
Only bounded queues of count slices are sent; the small global sample design
matrix is reused in tasks. Full count/editing matrices remain in parent RAM.
all_values=FALSE hides per-group sum_cov, mut, min_cov, called_sites,
average_called. diagnostics=TRUE adds time, code, evaluations, dispersion,
Hessian status, warnings/output, covered sample count and residual df.
betabin_control=list(maxit=2000), betabin_phi_ini=NULL configure aod fitting.
Strict convergence and positive-definite fixed covariance are required.
The original pseudocount (added equally to coverage and mutation) is retained.
No model fallback, automatic retries or hard wall-clock timeouts.

PREFILTER AND RANK POLICY
prefilter=FALSE by default. If TRUE, require min_replicates=2 covered samples in
EVERY group_var level. This is deliberately more conservative than the previous
version's reference-plus-one-comparison eligibility rule. It is coverage-only,
not an editing filter. It is necessary but not sufficient for full design rank.
Globally confounded/rank-deficient designs are rejected before fitting.
At each site, zero/masked coverage samples are dropped; the complete global model
matrix must remain full rank with positive residual df. Otherwise all site tests
are NA. This intentionally conservative policy does not rescue individually
estimable contrasts from a rank-deficient site design. More covariates can
therefore increase NA counts. No missing metadata are silently dropped.

FDR AND SUMMARY
BH correction separately across valid fitted sites for EACH requested coefficient
or contrast, before output-only filter.sites. No correction across the collection
of contrasts. Filtered/failed fits are NA; high failure rates need investigation.
Enabling prefilter changes the testing family. Attributes design_matrix,
contrast_matrix and run_summary document the fitted design and testing scope.

VALIDATION
Generated and checked for string/delimiter balance in Python only. R is not
installed in the creation environment: no R parse or numerical validation has
been executed. Run the supplied assertions and inspect diagnostics on real data
before production analysis. This is a substantial revision, not a tested release.
