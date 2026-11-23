# Bullseye: multifactor fixed-effect editing analysis.
# Supported families: betabin (requires aod), bin, quasibinomial.
# Display names insert underscores between factor variables and levels.
# Internal model coefficient names remain unchanged for selecting tests.

.be_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x) || x > .Machine$integer.max)
    stop(name, " must be a positive integer.")
  as.integer(x)
}
.be_slot <- function(x, name, default = NA) {
  if (name %in% methods::slotNames(x)) methods::slot(x, name) else default
}
.be_tests <- function(X, test_coefficients = NULL, contrasts = NULL) {
  cn <- colnames(X)
  if (!is.null(contrasts) && !is.null(test_coefficients))
    stop("Specify test_coefficients OR contrasts, not both.")
  if (is.null(contrasts)) {
    if (is.null(test_coefficients)) test_coefficients <- setdiff(cn, "(Intercept)")
    if (!is.character(test_coefficients) || !length(test_coefficients) ||
        anyNA(test_coefficients) || anyDuplicated(test_coefficients) ||
        any(!test_coefficients %in% cn))
      stop("test_coefficients must be unique design-matrix column names.")
    L <- diag(length(cn))[, match(test_coefficients, cn), drop = FALSE]
    dimnames(L) <- list(cn, test_coefficients)
  } else {
    if (!is.list(contrasts) || !length(contrasts) || is.null(names(contrasts)) ||
        anyNA(names(contrasts)) || any(!nzchar(names(contrasts))) ||
        anyDuplicated(names(contrasts)))
      stop("contrasts must be a uniquely named list of named numeric vectors.")
    L <- matrix(0, length(cn), length(contrasts), dimnames = list(cn, names(contrasts)))
    for (j in seq_along(contrasts)) {
      v <- contrasts[[j]]
      if (!is.numeric(v) || !length(v) || any(!is.finite(v)) || is.null(names(v)) ||
          anyNA(names(v)) || anyDuplicated(names(v)) || any(!names(v) %in% cn) || all(v == 0))
        stop("Contrast weights must be finite and named with known design columns.")
      L[names(v), j] <- v
    }
  }
  L
}

# Rename exact factor/level tokens, not arbitrary substrings. Custom contrast
# labels are retained verbatim. Non-treatment encodings stay as R names them.
.be_display_labels <- function(labels, colData, custom_contrasts = FALSE) {
  if (custom_contrasts) return(labels)
  raw <- nice <- character()
  for (variable in names(colData)) {
    if (!is.factor(colData[[variable]])) next
    lv <- levels(colData[[variable]])
    raw <- c(raw, paste0(variable, lv))
    nice <- c(nice, paste0(variable, "_", lv))
  }
  vapply(labels, function(label) {
    parts <- strsplit(label, ":", fixed = TRUE)[[1L]]
    for (i in seq_along(parts)) {
      hits <- which(raw == parts[i])
      if (length(hits) == 1L) parts[i] <- nice[hits]
    }
    paste(parts, collapse = ":")
  }, character(1), USE.NAMES = FALSE)
}

.be_fit <- function(cov, mut, X, L, link, control, phi_ini) {
  k <- ncol(L)
  out <- list(coefficient = rep(NA_real_, k), se = rep(NA_real_, k),
              pvalue = rep(NA_real_, k), status = "fit_failed")
  dg <- list(fit_code = NA_integer_, fit_evaluations = NA_integer_,
             fit_phi = NA_real_, fit_singular_hessian = NA, fit_elapsed_sec = NA_real_,
             fit_warning = "", fit_optimizer_output = "", fit_n_samples = NA_integer_,
             fit_residual_df = NA_integer_)
  start <- proc.time()[["elapsed"]]; warns <- character()
  out <- tryCatch(withCallingHandlers({
    keep <- is.finite(cov) & cov > 0
    x <- X[keep, , drop = FALSE]; cv <- cov[keep]; mt <- mut[keep]
    dg$fit_n_samples <- nrow(x); dg$fit_residual_df <- nrow(x)-ncol(x)
    if (nrow(x) <= ncol(x)) {
      out$status <- "insufficient_residual_df"
    } else if (qr(x)$rank < ncol(x)) {
      out$status <- "rank_deficient_site_design"
    } else {
      pc <- exp(mean(log(cv+0.5)))*0.005
      cv <- cv+pc; mt <- mt+pc
      ok <- TRUE
      if (link == "betabin") {
        internal <- sprintf("BE_x%04d", seq_len(ncol(x)))
        d <- as.data.frame(x); names(d) <- internal
        d$BE_cov <- cv; d$BE_mut <- mt
        f <- stats::reformulate(internal, response = "cbind(BE_mut, BE_cov-BE_mut)",
                                intercept = FALSE)
        txt <- utils::capture.output({
          mod <- aod::betabin(f, random = ~1, data = d, phi.ini = phi_ini,
                              warnings = TRUE, hessian = TRUE, control = control)
        })
        dg$fit_optimizer_output <- paste(txt, collapse = " | ")
        dg$fit_code <- as.integer(.be_slot(mod, "code", NA_integer_))[1L]
        dg$fit_evaluations <- as.integer(.be_slot(mod, "iterations", NA_integer_))[1L]
        dg$fit_phi <- as.numeric(.be_slot(mod, "random.param", NA_real_))[1L]
        dg$fit_singular_hessian <- as.logical(.be_slot(mod, "singular.hessian", NA))[1L]
        if (is.na(dg$fit_code) || dg$fit_code != 0L) {
          out$status <- "optimizer_not_converged"; ok <- FALSE
        } else if (isTRUE(dg$fit_singular_hessian)) {
          out$status <- "singular_hessian"; ok <- FALSE
        } else if (!is.finite(dg$fit_phi) || dg$fit_phi < 0 || dg$fit_phi >= 1) {
          out$status <- "invalid_dispersion"; ok <- FALSE
        }
        if (ok) {
          b0 <- aod::coef(mod); V0 <- aod::vcov(mod)
          if (!all(internal %in% names(b0)) || !all(internal %in% rownames(V0)) ||
              !all(internal %in% colnames(V0)))
            stop("Unexpected aod coefficient names; cannot safely match design columns.")
          b <- unname(b0[internal]); V <- unname(V0[internal, internal, drop = FALSE])
        }
      } else {
        fam <- if (link == "bin") stats::binomial() else stats::quasibinomial()
        mod <- stats::glm.fit(x = x, y = mt/cv, weights = cv, family = fam)
        dg$fit_code <- if (isTRUE(mod$converged)) 0L else 1L
        if (!isTRUE(mod$converged)) {
          out$status <- "optimizer_not_converged"; ok <- FALSE
        } else {
          class(mod) <- c("glm", "lm")
          b <- unname(stats::coef(mod)); V <- unname(stats::vcov(mod))
        }
      }
      if (ok) {
        if (any(!is.finite(b)) || any(!is.finite(V)) ||
            inherits(try(chol((V+t(V))/2), silent = TRUE), "try-error")) {
          out$status <- "invalid_covariance"
        } else {
          est <- drop(crossprod(L, b)); var <- colSums(L*(V %*% L))
          valid <- is.finite(est) & is.finite(var) & var > 0
          out$coefficient[valid] <- est[valid]; out$se[valid] <- sqrt(var[valid])
          out$pvalue[valid] <- stats::pchisq(est[valid]^2/var[valid], 1, lower.tail = FALSE)
          out$status <- if (all(valid)) "ok" else if (any(valid)) {
            "partially_testable"
          } else "no_testable_contrast"
        }
      }
    }
    out
  }, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {
    warns <<- c(warns, paste("ERROR:", conditionMessage(e)))
    out$status <- "fit_failed"; out
  })
  dg$fit_elapsed_sec <- proc.time()[["elapsed"]]-start
  dg$fit_warning <- paste(unique(warns), collapse = " | ")
  out$diagnostics <- as.data.frame(dg, stringsAsFactors = FALSE)
  out
}

# Do not rename design_matrix to X: that conflicts with the apply APIs.
.be_batch <- function(job, design_matrix, L, link, control, phi_ini) {
  n <- length(job$rows); k <- ncol(L)
  b <- se <- p <- matrix(NA_real_, n, k)
  status <- character(n); dg <- vector("list", n)
  for (i in seq_len(n)) {
    z <- .be_fit(cov = job$cov[i, ], mut = job$mut[i, ], X = design_matrix,
                  L = L, link = link, control = control, phi_ini = phi_ini)
    b[i, ] <- z$coefficient; se[i, ] <- z$se; p[i, ] <- z$pvalue
    status[i] <- z$status; dg[[i]] <- z$diagnostics
  }
  list(rows = job$rows, b = b, se = se, p = p, status = status,
       diagnostics = do.call(rbind, dg))
}

bullseye <- function(data_cov, data_mut, colData = NULL,
                     filter.col = NULL, filter.with = NULL,
                     min = 0.05, filter.sites = NULL, dataset = NULL,
                     design = NULL, min.cov = 10, link = "betabin", numCores = 1,
                     batchSize = 50L, queueFactor = 8L, diagnostics = FALSE,
                     prefilter = FALSE, min_replicates = 2L,
                     betabin_control = list(maxit = 2000), betabin_phi_ini = NULL,
                     all_values = FALSE, verbose = TRUE, group_var = NULL,
                     test_coefficients = NULL, contrasts = NULL) {
  link <- match.arg(link, c("betabin", "bin", "quasibinomial"))
  numCores <- .be_int(numCores, "numCores"); batchSize <- .be_int(batchSize, "batchSize")
  queueFactor <- .be_int(queueFactor, "queueFactor")
  min_replicates <- .be_int(min_replicates, "min_replicates")
  for (nm in c("diagnostics", "prefilter", "all_values", "verbose")) {
    z <- get(nm)
    if (!is.logical(z) || length(z) != 1L || is.na(z)) stop(nm, " must be TRUE/FALSE.")
  }
  if (!is.numeric(min) || length(min) != 1L || !is.finite(min) || min < 0 || min > 1 ||
      !is.numeric(min.cov) || length(min.cov) != 1L || !is.finite(min.cov) || min.cov < 0)
    stop("Invalid min or min.cov.")
  if (!is.null(filter.sites) && (length(filter.sites) != 1L || is.na(filter.sites) ||
      !filter.sites %in% c("all", "called"))) stop("Invalid filter.sites.")
  if (!is.data.frame(colData) || !nrow(colData) || anyDuplicated(rownames(colData)))
    stop("Provide colData with unique sample row names.")
  if (!is.null(filter.col)) {
    if (length(filter.col) != 1L || !filter.col %in% names(colData) ||
        length(filter.with) != 1L || is.na(filter.with)) stop("Invalid sample filter.")
    take <- !is.na(colData[[filter.col]]) & colData[[filter.col]] == filter.with
    colData <- colData[take, , drop = FALSE]
  }
  if (!nrow(colData)) stop("No samples remain.")
  if (!inherits(design, "formula") || length(design) != 2L)
    stop("design must be a one-sided fixed-effects formula.")
  if (!length(all.vars(design)) || any(!all.vars(design) %in% names(colData)))
    stop("All design variables must exist in colData.")
  dt <- paste(deparse(design), collapse = "")
  if (grepl("|", dt, fixed = TRUE) || grepl("offset(", dt, fixed = TRUE))
    stop("Random-effects syntax and offsets are unsupported.")
  if (is.null(group_var)) {
    vv <- all.vars(design)
    if (length(vv) != 1L) stop("Specify group_var explicitly for a multifactor design.")
    group_var <- vv
  }
  if (!is.character(group_var) || length(group_var) != 1L || is.na(group_var) ||
      !group_var %in% names(colData)) stop("Invalid group_var.")
  if (anyNA(colData[[group_var]])) stop("group_var contains missing values.")
  if (is.numeric(colData[[group_var]]) && !is.factor(colData[[group_var]]))
    stop("Convert group_var to a factor explicitly.")
  colData[[group_var]] <- droplevels(as.factor(colData[[group_var]]))
  for (nm in names(colData)) {
    if (is.character(colData[[nm]])) colData[[nm]] <- factor(colData[[nm]])
    if (is.factor(colData[[nm]])) colData[[nm]] <- droplevels(colData[[nm]])
  }
  lev <- levels(colData[[group_var]])
  if (length(lev) < 2L) stop("group_var needs at least two observed levels.")
  reserved <- c("site", "n", "baseEdit", "fit_status", "fit_code", "fit_evaluations",
                 "fit_phi", "fit_singular_hessian", "fit_elapsed_sec", "fit_warning",
                 "fit_optimizer_output", "fit_n_samples", "fit_residual_df")
  if (any(lev %in% reserved)) stop("Group levels collide with reserved output names.")
  ref <- lev[1L]
  mf <- stats::model.frame(design, data = colData, na.action = stats::na.fail)
  if (!is.null(stats::model.offset(mf))) stop("Offsets are unsupported.")
  X <- stats::model.matrix(design, data = mf)
  if (nrow(X) != nrow(colData) || !ncol(X) || any(!is.finite(X))) stop("Invalid design matrix.")
  if (qr(X)$rank < ncol(X)) stop("Global design is rank deficient; check confounding.")
  if (nrow(X) <= ncol(X)) stop("Design requires positive residual degrees of freedom.")
  L <- .be_tests(X, test_coefficients, contrasts)
  display <- .be_display_labels(colnames(L), colData, custom_contrasts = !is.null(contrasts))
  if (anyDuplicated(display)) stop("Formatted test labels collide; use distinct custom contrast labels.")
  test_map <- data.frame(internal_name = colnames(L), output_prefix = display,
                          stringsAsFactors = FALSE)
  if (!is.list(betabin_control) || (length(betabin_control) &&
      (is.null(names(betabin_control)) || anyNA(names(betabin_control)) ||
       any(!nzchar(names(betabin_control)))))) stop("betabin_control must be a named list.")
  if (!is.null(betabin_phi_ini) && (!is.numeric(betabin_phi_ini) ||
      length(betabin_phi_ini) != 1L || !is.finite(betabin_phi_ini) ||
      betabin_phi_ini <= 0 || betabin_phi_ini >= 1))
    stop("betabin_phi_ini must be NULL or strictly between 0 and 1.")
  if (link == "betabin" && !requireNamespace("aod", quietly = TRUE))
    stop("Install aod for link='betabin'.")
  if (!is.data.frame(data_cov) || !is.data.frame(data_mut) ||
      !identical(rownames(data_cov), rownames(data_mut)) || anyDuplicated(rownames(data_cov)))
    stop("Counts must be data.frames with identical unique site row names.")
  samples <- rownames(colData)
  if (!all(samples %in% names(data_cov)) || !all(samples %in% names(data_mut)))
    stop("Some colData samples are absent from count tables.")
  if (!is.null(dataset)) {
    if (!all(c("chr", "start", "end", "strand") %in% names(dataset)))
      stop("dataset requires chr, start, end, strand.")
    ids <- paste(dataset$chr, dataset$start, dataset$end, dataset$strand, sep = "_")
    take <- rownames(data_cov) %in% ids
    data_cov <- data_cov[take, , drop = FALSE]; data_mut <- data_mut[take, , drop = FALSE]
  }
  ns <- nrow(data_cov)
  if (!ns) stop("No sites remain.")
  sites <- rownames(data_cov)
  n <- if ("n" %in% names(data_cov)) data_cov$n else rep(1, ns)
  if (!is.numeric(n) || any(!is.finite(n)) || any(n <= 0)) stop("n must be positive finite.")
  if ("n" %in% names(data_mut) && !isTRUE(all.equal(n, data_mut$n)))
    stop("Count tables have inconsistent n columns.")
  cv <- as.matrix(data_cov[, samples, drop = FALSE]); mt <- as.matrix(data_mut[, samples, drop = FALSE])
  if (!is.numeric(cv) || !is.numeric(mt) || any(is.infinite(cv)) || any(is.infinite(mt)) ||
      any(cv < 0, na.rm = TRUE) || any(mt < 0, na.rm = TRUE) || any(mt > cv, na.rm = TRUE))
    stop("Counts must be finite nonnegative numbers or NA, with mut <= cov.")
  low <- is.na(cv) | sweep(cv, 1L, min.cov*n, "<")
  if (any(is.na(mt) & !low)) stop("Mutation counts are NA where coverage passes threshold.")
  cv[low] <- 0; mt[low] <- 0
  edit <- mt/cv*n; edit[!is.finite(edit)] <- NA_real_
  called <- edit; called[!is.na(called) & called < min] <- NA_real_
  groups <- split(seq_along(samples), colData[[group_var]])
  coverage <- vapply(groups, function(ii) rowSums(cv[, ii, drop = FALSE] > 0), numeric(ns))
  coverage <- matrix(coverage, nrow = ns, dimnames = list(sites, lev))
  eligible <- rep(TRUE, ns)
  if (prefilter) eligible <- apply(coverage >= min_replicates, 1L, all)
  rows <- which(eligible); nt <- ncol(L)
  B <- SE <- P <- matrix(NA_real_, ns, nt, dimnames = list(sites, colnames(L)))
  status <- rep("prefiltered_coverage", ns)
  dg <- data.frame(fit_code = rep(NA_integer_, ns), fit_evaluations = rep(NA_integer_, ns),
                    fit_phi = rep(NA_real_, ns), fit_singular_hessian = rep(NA, ns),
                    fit_elapsed_sec = rep(NA_real_, ns), fit_warning = rep("", ns),
                    fit_optimizer_output = rep("", ns), fit_n_samples = rep(NA_integer_, ns),
                    fit_residual_df = rep(NA_integer_, ns), stringsAsFactors = FALSE)
  nw <- if (length(rows)) as.integer(min(numCores, ceiling(length(rows)/batchSize))) else 0L
  cl <- NULL
  if (nw > 1L) {
    cl <- parallel::makePSOCKcluster(nw)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c(".be_slot", ".be_fit", ".be_batch"),
                             envir = environment(bullseye))
  }
  if (verbose) message("Sites: ", ns, "; eligible: ", length(rows),
                        "; design columns: ", ncol(X), "; tests: ", nt)
  pos <- 1L
  while (pos <= length(rows)) {
    nj <- as.integer(min(as.double(queueFactor)*max(1L, nw), ceiling((length(rows)-pos+1)/batchSize)))
    jobs <- vector("list", nj)
    for (j in seq_len(nj)) {
      end <- min(length(rows), pos+batchSize-1); rr <- rows[seq.int(pos, end)]
      jobs[[j]] <- list(rows = rr, cov = cv[rr, , drop = FALSE], mut = mt[rr, , drop = FALSE])
      pos <- end+1
    }
    if (is.null(cl)) {
      fits <- lapply(X = jobs, FUN = .be_batch, design_matrix = X, L = L,
                      link = link, control = betabin_control, phi_ini = betabin_phi_ini)
    } else {
      fits <- parallel::parLapplyLB(cl = cl, X = jobs, fun = .be_batch,
                                    design_matrix = X, L = L, link = link,
                                    control = betabin_control, phi_ini = betabin_phi_ini,
                                    chunk.size = 1L)
    }
    for (z in fits) {
      B[z$rows, ] <- z$b; SE[z$rows, ] <- z$se; P[z$rows, ] <- z$p
      status[z$rows] <- z$status; dg[z$rows, ] <- z$diagnostics[, names(dg), drop = FALSE]
    }
    rm(jobs, fits, z)
  }
  out <- data.frame(site = sites, n = n, stringsAsFactors = FALSE)
  for (g in lev) {
    ii <- groups[[g]]
    total_cv <- rowSums(cv[, ii, drop = FALSE]); total_mt <- rowSums(mt[, ii, drop = FALSE])
    pooled <- total_mt/total_cv*n; pooled[!is.finite(pooled)] <- NA_real_
    if (g %in% names(out)) stop("Group names collide with output columns.")
    out[[g]] <- pooled
    if (all_values) {
      cols <- paste0(g, c("_sum_cov", "_mut", "_min_cov", "_called_sites", "_average_called"))
      if (any(cols %in% c(names(out), lev))) stop("Group detail names collide with output columns.")
      out[[cols[1L]]] <- total_cv; out[[cols[2L]]] <- total_mt; out[[cols[3L]]] <- coverage[, g]
      out[[cols[4L]]] <- rowSums(edit[, ii, drop = FALSE] >= min, na.rm = TRUE)
      a <- rowMeans(called[, ii, drop = FALSE], na.rm = TRUE); a[is.nan(a)] <- NA_real_
      out[[cols[5L]]] <- a
    }
  }
  for (g in lev[-1L]) {
    nm <- paste0(g, "_log2FoldChange")
    if (nm %in% names(out)) stop("Fold-change output name collision.")
    out[[nm]] <- log2(out[[g]]/out[[ref]])
  }
  for (j in seq_len(nt)) {
    nn <- paste0(display[j], c("_coefficient", "_SE", "_pvalue", "_FDR"))
    if (any(nn %in% c(names(out), reserved))) stop("Test names collide with output columns.")
    out[[nn[1L]]] <- B[, j]; out[[nn[2L]]] <- SE[, j]; out[[nn[3L]]] <- P[, j]
    q <- rep(NA_real_, ns); valid <- is.finite(P[, j])
    q[valid] <- stats::p.adjust(P[valid, j], method = "BH"); out[[nn[4L]]] <- q
  }
  nc <- rowSums(!is.na(called))
  out$baseEdit <- ifelse(nc > 0, rowSums(called, na.rm = TRUE)/nc, 0)
  out$fit_status <- status
  if (diagnostics) out <- cbind(out, dg)
  if (anyDuplicated(names(out))) stop("Duplicate output columns.")
  if (!is.null(filter.sites)) {
    ge <- as.matrix(out[, lev, drop = FALSE]); pass <- !is.na(ge) & ge >= min
    take <- if (filter.sites == "called") rowSums(pass) == length(lev) else rowSums(pass) > 0
    out <- out[take, , drop = FALSE]
  }
  rownames(out) <- NULL
  attr(out, "design_matrix") <- X; attr(out, "contrast_matrix") <- L
  attr(out, "test_name_map") <- test_map
  attr(out, "run_summary") <- list(input_sites = ns, prefiltered_sites = sum(!eligible),
    submitted_sites = length(rows), output_sites = nrow(out), fit_status_counts = table(status),
    valid_tests_per_comparison = colSums(is.finite(P)), group_var = group_var,
    reference = ref, design = deparse(design),
    FDR_scope = "BH separately per requested test, valid fits before output filtering",
    prefilter = prefilter, min_replicates = if (prefilter) min_replicates else NA_integer_)
  if (verbose) message("Fit status: ", paste(names(table(status)), table(status), sep = "=", collapse = ", "))
  out
}
