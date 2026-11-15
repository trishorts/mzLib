# Reference values for SETTEST1 S4b's share models (QuantProject DEF-DIFF-SHARE-MODEL, DATA-DEFINITIONS v3.15, ruling S1):
#   a single:  log2 value ~ peptide + (design) + marker:(design) + (1 | sample)
#   a set:     log2 value ~ peptide + member:(design) + member:marker:(design) + (1 | sample) + (1 | sample:member)
# REML. The share effect is the contrast on the marker: terms (a set: their equal-weight average over the members), its
# se from the fit and its df Satterthwaite's (lmerTest::contest).
#
# Run ONCE, by hand, to (re)generate the fixtures beside this script:
#   E:/Programs/R-4.6.1/bin/Rscript.exe make_settest_share_model_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# The fixed part is an explicit matrix, as make_settest_model_fixtures.R builds it: per member, an intercept and the
# design's other columns, then marker x (the design's columns but the intercept), then one indicator for each of the
# member's peptides but the first in ordinal order. marker:intercept is left out: marker is constant within a peptide, so it
# is the sum of the marker peptides' own terms.
#
# Robust weighting (GR-8, GR-22): msqrob2 1.20.0's reweighting run to convergence, the loop of make_settest_model_fixtures.R
# unchanged (weights in the lmer call, each round refitted from scratch, cap 100).

suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(MASS) })
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", sprintf("%.17g", v))
put <- function(df, name) write.table(df, file.path(here, name), sep = "\t", quote = FALSE, row.names = FALSE)
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhobeg = 0.2, rhoend = 1e-12, maxfun = 1e5),
                    check.conv.singular = "ignore", check.conv.grad = "ignore", calc.derivs = FALSE)
maxitRob <- 100; tolRob <- 1e-6

robustFit <- function(f, d) {
  model <- lmerTest::lmer(f, data = d, REML = TRUE, control = ctrl)
  sseOld <- model@devcomp$cmp["pwrss"]
  it <- 0; converged <- FALSE
  while (it < maxitRob) {
    it <- it + 1
    res <- resid(model)
    d$w <- MASS::psi.huber(res / mad(res, 0))
    model <- lmerTest::lmer(f, data = d, weights = w, REML = TRUE, control = ctrl)
    sse <- model@devcomp$cmp["pwrss"]
    if (abs(sseOld - sse) / sseOld <= tolRob) { converged <- TRUE; break }
    sseOld <- sse
  }
  list(model = model, weights = d$w, iterations = it, converged = converged)
}

seed <- 20261115
set.seed(seed)

scenario <- function(name) {
  set <- startsWith(name, "set")
  context <- endsWith(name, "context")
  n <- if (context) 12 else if (endsWith(name, "missing")) 10 else 8
  samples <- data.frame(sample = sprintf("s%02d", 1:n), group = rep(c("young", "old"), length.out = n),
                        treatment = if (context) rep(c("a", "b"), each = n / 2) else "none")
  members <- if (!set) list(c(2, 3)) else if (context) list(c(2, 3), c(1, 2), c(2, 2), c(1, 3)) else list(c(2, 3), c(1, 2), c(2, 2))
  missing <- if (endsWith(name, "missing")) 0.12 else if (context) 0.05 else 0
  list(name = name, set = set, context = context, samples = samples, members = members, missing = missing)
}

simulate <- function(sc) {
  n <- nrow(sc$samples); old <- as.numeric(sc$samples$group == "old"); trt <- as.numeric(sc$samples$treatment == "b")
  shared <- rnorm(n, 0, 0.35)                                   # (1 | sample)
  rows <- NULL
  for (m in seq_along(sc$members)) {
    mu <- rnorm(1, 22, 1.2); change <- rnorm(1, 0.5, 0.4); share <- rnorm(1, -0.6, 0.3)
    te <- if (sc$context) rnorm(1, 0, 0.5) else 0; mte <- if (sc$context) rnorm(1, 0, 0.3) else 0
    own <- if (sc$set) rnorm(n, 0, 0.2) else rep(0, n)          # (1 | sample:member)
    counts <- sc$members[[m]]
    for (kind in c("marker", "reference")) {
      for (k in seq_len(counts[if (kind == "marker") 1 else 2])) {
        pep <- rnorm(1, 0, 0.8)
        mk <- kind == "marker"
        y <- mu + pep + change * old + te * trt + (if (mk) share * old + mte * trt else 0) + shared + own + rnorm(n, 0, 0.2)
        if (sc$missing > 0) y[runif(n) < sc$missing] <- NA
        rows <- rbind(rows, data.frame(member = sprintf("G%d", m), peptide = sprintf("G%d%s%d", m, if (mk) "M" else "R", k),
                                       kind = kind, sample = sc$samples$sample, y = y))
      }
    }
  }
  if (endsWith(sc$name, "missing")) rows$y[rows$kind == "marker" & rows$sample == "s03"] <- NA   # no marker value in one sample
  rows[!is.na(rows$y), ]
}

run <- function(sc) {
  d <- simulate(sc)
  d <- merge(d, sc$samples, by = "sample")
  d <- d[order(d$member, d$peptide, d$sample, method = "radix"), ]
  members <- sort(unique(d$member), method = "radix")
  X <- NULL; names <- NULL
  for (m in members) {
    is <- as.numeric(d$member == m); mk <- is * (d$kind == "marker")
    X <- cbind(X, is, is * (d$group == "old")); names <- c(names, paste0(m, ":intercept"), paste0(m, ":old"))
    if (sc$context) { X <- cbind(X, is * (d$treatment == "b")); names <- c(names, paste0(m, ":treatmentb")) }
    X <- cbind(X, mk * (d$group == "old")); names <- c(names, paste0(m, ":marker:old"))
    if (sc$context) { X <- cbind(X, mk * (d$treatment == "b")); names <- c(names, paste0(m, ":marker:treatmentb")) }
  }
  for (m in members) {
    peps <- sort(unique(d$peptide[d$member == m]), method = "radix")
    for (p in peps[-1]) { X <- cbind(X, as.numeric(d$peptide == p)); names <- c(names, p) }
  }
  colnames(X) <- names
  d$sm <- paste(d$sample, d$member, sep = ":")
  f <- if (sc$set) y ~ 0 + X + (1 | sample) + (1 | sm) else y ~ 0 + X + (1 | sample)
  L <- numeric(ncol(X)); L[grep(":marker:old$", names)] <- 1 / length(members)
  out <- function(fit, tag, extra = list()) {
    ct <- contest(fit, L, joint = FALSE)
    vc <- as.data.frame(VarCorr(fit))
    vsm <- if (sc$set) vc$vcov[vc$grp == "sm"] else NA
    put(do.call(data.frame, c(list(estimate = fmt(ct$Estimate), se = fmt(ct$`Std. Error`), df = fmt(ct$df),
                                   sigma2 = fmt(sigma(fit)^2), var_sample = fmt(vc$vcov[vc$grp == "sample"]),
                                   var_sample_member = fmt(vsm), reml = fmt(REMLcrit(fit))), extra,
                              list(members = length(members), peptides = length(unique(d$peptide)), observations = nrow(d)))),
        sprintf("settest_share_model_%s_%s.tsv", sc$name, tag))
    cat(sc$name, tag, ": estimate ", ct$Estimate, " se ", ct$`Std. Error`, " df ", ct$df, "\n", sep = "")
  }
  put(sc$samples, sprintf("settest_share_model_%s_samples.tsv", sc$name))
  put(data.frame(member = d$member, peptide = d$peptide, kind = d$kind, sample = d$sample, y = fmt(d$y)),
      sprintf("settest_share_model_%s_values.tsv", sc$name))
  out(lmerTest::lmer(f, data = d, REML = TRUE, control = ctrl), "result")
  rb <- robustFit(f, d)
  stopifnot(rb$converged)
  put(data.frame(member = d$member, peptide = d$peptide, sample = d$sample, weight = fmt(rb$weights)),
      sprintf("settest_share_model_%s_robust_weights.tsv", sc$name))
  out(rb$model, "robust_result", list(iterations = rb$iterations, down_weighted = sum(rb$weights < 1)))
}

for (s in c("single_balanced", "single_missing", "single_context", "set_balanced", "set_context")) run(scenario(s))

writeLines(c(
  sprintf("generated_utc %s UTC", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")),
  sprintf("R %s", R.version.string),
  sprintf("lme4 %s; lmerTest %s; MASS %s", packageVersion("lme4"), packageVersion("lmerTest"), packageVersion("MASS")),
  sprintf("seed %d", seed),
  "script make_settest_share_model_fixtures.R",
  "lmer REML on an explicit fixed matrix (per member: intercept, old[, treatment b], marker:old[, marker:treatment b]; peptides within member), single (1|sample), set (1|sample) + (1|sample:member), bobyqa rhoend 1e-12; the share effect and its Satterthwaite df via lmerTest::contest",
  "robust (*_robust_*) = msqrob2 1.20.0 .robust_fitting rules (psi.huber k = 1.345 on resid / mad(res, 0); pwrss tol 1e-6), cap 100, each round refitted from scratch with weights in the call, as make_settest_model_fixtures.R"),
  file.path(here, "PROVENANCE_settest_share_model.txt"))
