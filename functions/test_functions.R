# Checks that the stand-alone functions in REE_functions.R give exactly the app's results.
# Run from the repository root:  Rscript functions/test_functions.R
suppressMessages({library(ggplot2); library(MASS)})
app <- new.env()
owd <- setwd("app")
for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = app)
setwd(owd)
source("functions/REE_functions.R")

n_pass <- 0L
check <- function(desc, cond) {
  if (!isTRUE(cond)) stop("FAILED: ", desc, call. = FALSE)
  n_pass <<- n_pass + 1L
}
same <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  isTRUE(all(is.na(a) == is.na(b))) && isTRUE(all(abs(a[!is.na(a)] - b[!is.na(b)]) <= tol * pmax(1, abs(b[!is.na(b)]))))
}

# test data: the app's example (molar), reference data (mass and molar), and gappy/odd patterns
ex <- read.csv("app/data/example_user_data.csv", stringsAsFactors = FALSE)
ref <- local({ owd <- setwd("app"); on.exit(setwd(owd)); app$load_reference() })
ref <- ref[rowSums(!is.na(ref[, app$REE])) > 0, ]
set.seed(1)
refs <- ref[sample(nrow(ref), min(600, nrow(ref))), ]
gappy <- ex
gappy$Pr[c(1, 5)] <- NA; gappy$Tb[2] <- NA; gappy$Tm[c(3, 6)] <- NA; gappy$Lu[4] <- NA; gappy$La[7] <- NA
gappy$Ho[8] <- NA; gappy$Er[9] <- NA; gappy$Gd[10] <- gappy$Gd[10] * 5; gappy$Ce[11] <- 0
datasets <- list(example = list(d = ex, b = "molar"), gappy = list(d = gappy, b = "molar"),
                 reference = list(d = refs, b = refs$basis))

# ---- REE_ratios vs the app's prepare_user_data(), every option combination ----------------------
opts <- expand.grid(norm = c(names(app$NORMALISERS), "custom"), lree = c("LaPrNd", "PrNdSm"),
                    gap = c("interpolate", "groupmean", "none"), gd = c("measured", "ignore"),
                    gd_method = c("BD1996", "BDgeo", "SmNd", "TbDy", "BLambdaR"),
                    ce_form = c("geoPrNd", "arithLaPr", "geoLaPr"), stringsAsFactors = FALSE)
set.seed(2)
opts <- opts[sample(nrow(opts), 400), ]          # 400 random combinations (of 3240)
opts$gd_baseline <- sample(c(1, 1, 1.1, 1.4), nrow(opts), replace = TRUE)
pairs <- c(HREE_LREE = "x", MREE_MREEstar = "y", Ce_Cestar = "Ce_Ceanom", La_Lastar = "La_Laanom",
           Eu_Eustar = "Eu_Euanom", Gd_Gdstar = "Gd_Gdstar", Gd_anth = "Gd_anth", REE_rough = "REE_rough",
           n_estimated = "n_estimated", n_REE_measured = "n_REE_measured")
for (dn in names(datasets)) {
  d <- datasets[[dn]]$d; b <- datasets[[dn]]$b
  for (i in seq_len(nrow(opts))) {
    o <- opts[i, ]
    nc <- if (o$norm == "custom") app$NORMALISERS$NASC_McLennan1989 * 1.1 else NULL
    cf <- app$default_cfg(norm = o$norm, norm_custom = nc, lree = o$lree, gap = o$gap, gd = o$gd,
                          gd_method = o$gd_method, gd_baseline = o$gd_baseline, ce_form = o$ce_form)
    dd <- d; dd$basis <- NULL
    a <- tryCatch(app$prepare_user_data(dd, "concentrations", b, cf)$data, error = function(e) e)
    s <- tryCatch(REE_ratios(dd, basis = b, norm = o$norm, norm_custom = nc, lree = o$lree, gap = o$gap,
                             gd = o$gd, gd_method = o$gd_method, gd_baseline = o$gd_baseline,
                             ce_form = o$ce_form, append = FALSE), error = function(e) e)
    if (inherits(a, "error") || inherits(s, "error")) {
      check(sprintf("%s / combination %d: both fail or both succeed", dn, i), inherits(a, "error") && inherits(s, "error"))
      next
    }
    s <- s[a$row, ]                                # the app drops empty rows; compare the same rows
    for (k in names(pairs)) check(sprintf("%s / combination %d / %s", dn, i, k), same(s[[k]], a[[pairs[[k]]]]))
  }
}
cat(sprintf("REE_ratios: %d option combinations x %d data sets agree with the app\n", nrow(opts), length(datasets)))

# ---- screening ------------------------------------------------------------------------------------
s <- REE_ratios(gappy, rough_max = 0.05, min_ree = 14, append = FALSE)
check("screening flags too few REE", all(s$screen[s$n_REE_measured < 14 & is.finite(s$HREE_LREE)] == "excluded (too few REE)"))
check("screening flags roughness", all(s$screen[s$n_REE_measured == 14 & is.finite(s$REE_rough) & s$REE_rough > 0.05] == "excluded (roughness)"))
s2 <- REE_ratios(gappy, req_ree = "Pr", append = FALSE)
check("required REE", all(s2$screen[is.na(gappy$Pr)] == "excluded (too few REE)"))

# ---- REE_lambda vs the app ------------------------------------------------------------------------
for (dn in names(datasets)) {
  d <- datasets[[dn]]$d; b <- datasets[[dn]]$b
  for (g in c("measured", "ignore")) {
    cf <- app$default_cfg(gd = g); dd <- d; dd$basis <- NULL
    a <- app$prepare_user_data(dd, "concentrations", b, cf, lambda = TRUE)$data
    s <- REE_lambda(dd, basis = b, gd = g, append = FALSE)[a$row, ]
    for (k in c("lambda0", "lambda1", "lambda2", "lambda3", "lambda4", "tau1", "tau2", "tau3", "tau4", "lambda_rough"))
      check(sprintf("lambda %s / %s / %s", dn, g, k), same(s[[k]], a[[k]], 1e-10))
  }
}
cat("REE_lambda: agrees with the app (Gd measured and ignored, 3 data sets)\n")

# ---- REE_spectrum vs the app ----------------------------------------------------------------------
for (nm in c(names(app$NORMALISERS))) for (n2 in list(NULL, "Nd", "Yb")) {
  cf <- app$default_cfg(norm = nm); dd <- ex; dd$basis <- NULL
  a <- app$ree_spectra(app$prepare_user_data(dd, "concentrations", "molar", cf)$conc, "molar", cf, n2)
  s <- REE_spectrum(dd, basis = "molar", norm = nm, norm2 = n2, colour = "group", plot = FALSE)
  for (el in app$REE) {
    check(sprintf("spectrum %s %s", nm, el), same(s$normalised[[paste0("norm_", el)]], a$n[, el]))
    if (!is.null(n2)) check(sprintf("spectrum norm2 %s %s", nm, el), same(s$normalised[[paste0("norm2_", el)]], a$n2[, el]))
  }
}
p <- REE_spectrum(ex, colour = "water_depth_m", x = "radius", plot = FALSE)$plot
check("spectrum figure builds (numeric colour, radius axis)", inherits(ggplot2::ggplot_build(p), "ggplot_built"))
p <- REE_spectrum(gappy, colour = "group", connect = FALSE, plot = FALSE)$plot
check("spectrum figure builds (groups, broken lines)", inherits(ggplot2::ggplot_build(p), "ggplot_built"))
cat("REE_spectrum: normalised values agree with the app; figures build\n")
cat(sprintf("All %d checks passed.\n", n_pass))
