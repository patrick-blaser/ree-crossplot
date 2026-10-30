# =============================================================================
# REE_functions.R - stand-alone R functions behind the REE cross-plot app
#
#   REE_ratios()    HREE/LREE, MREE/MREE*, Ce/Ce*, La/La*, Eu/Eu*, Gd/Gd*, anthropogenic Gd
#                   and REE roughness
#   REE_lambda()    lambda (orthogonal polynomial) shape parameters after Anenburg and
#                   Williams (2022), base-R port of BLambdaR, with lambda roughness
#   REE_spectrum()  normalised REE patterns (spectra), table and ggplot2 figure
#   REE_read()      read CSV / TSV / PANGAEA .tab files (metadata header skipped, citation kept)
#
# The functions give the same results as the app (https://patrick-blaser.github.io/ree-crossplot/)
# with the same options and defaults. They only need base R; REE_spectrum() also needs ggplot2.
# Usage:   source("REE_functions.R")
#
# Author: Patrick Blaser. Licence: MIT (see LICENSE in the repository).
# If you use these functions, please cite:
#   Blaser, P., Gutjahr, M., Poeppelmeier, F., Frank, M., Kaboth-Bahr, S., Lippold, J., 2020.
#   Labrador Sea bottom water provenance and REE exchange during the past 35,000 years.
#   Earth Planet. Sci. Lett. 542, 116299. https://doi.org/10.1016/j.epsl.2020.116299
# and, where used, the method references given in the function documentation below.
# =============================================================================

# ---- constants ---------------------------------------------------------------------------------

.REE <- c("La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu")

.ATOMIC_MASS <- c(La = 138.905, Ce = 140.116, Pr = 140.908, Nd = 144.242, Sm = 150.36,
                  Eu = 151.964, Gd = 157.25, Tb = 158.925, Dy = 162.500, Ho = 164.930,
                  Er = 167.259, Tm = 168.934, Yb = 173.045, Lu = 174.967)
.Y_MASS <- 88.906

# normalising references, ppm (mass basis)
.NORMALISERS <- lapply(list(
  PAAS_Pourmand2012 = c(44.56, 88.25, 10.15, 37.32, 6.884, 1.215, 6.043, 0.8914, 5.325, 1.053,
                        3.075, 0.451, 3.012, 0.4386),
  PAAS_McLennan1989 = c(38.2, 79.6, 8.83, 33.9, 5.55, 1.08, 4.66, 0.774, 4.68, 0.991,
                        2.85, 0.405, 2.82, 0.433),
  NASC_McLennan1989 = c(32, 73, 7.9, 33, 5.7, 1.24, 5.2, 0.85, 5.8, 1.04, 3.4, 0.5, 3.1, 0.48),
  EUS_Bau2018       = c(44.3, 88.5, 10.6, 39.5, 7.30, 1.48, 6.34, 0.944, 5.86, 1.17,
                        3.43, 0.492, 3.26, 0.485),
  CI_Pourmand2012   = c(0.2469, 0.6321, 0.0959, 0.4854, 0.1556, 0.0599, 0.2093, 0.0378, 0.2577,
                        0.0554, 0.1667, 0.0261, 0.1694, 0.0256),
  CI_ONeill2016     = c(0.2472, 0.6308, 0.0950, 0.4793, 0.1542, 0.0592, 0.2059, 0.0375, 0.2540,
                        0.0554, 0.1645, 0.0258, 0.1684, 0.0251)), function(v) setNames(v, .REE))
.NORM_LABELS <- c(PAAS_Pourmand2012 = "PAAS - Pourmand et al. (2012)",
                  PAAS_McLennan1989 = "PAAS - McLennan (1989)",
                  EUS_Bau2018       = "EUS (European Shale) - Bau et al. (2018)",
                  NASC_McLennan1989 = "NASC - McLennan (1989)",
                  CI_Pourmand2012   = "CI chondrite - Pourmand et al. (2012)",
                  CI_ONeill2016     = "CI chondrite - O'Neill (2016)")
# Y of the references (spectra only; NA = no value)
.Y_NORM <- c(PAAS_Pourmand2012 = 27.31, PAAS_McLennan1989 = 27, EUS_Bau2018 = 31.9,
             NASC_McLennan1989 = NA_real_, CI_Pourmand2012 = NA_real_, CI_ONeill2016 = NA_real_)

.LREE_SETS <- list(LaPrNd = c("La", "Pr", "Nd"), PrNdSm = c("Pr", "Nd", "Sm"))
.MREE_SET <- c("Gd", "Tb", "Dy")
.HREE_SET <- c("Tm", "Yb", "Lu")
.REE_Z <- setNames(c(57:60, 62:71), .REE)          # atomic numbers (Pm, 61, omitted)

# effective ionic radii (pm) of the 3+ ions, 6-coordinated (Shannon, 1976), for the spectra
.RADII_PM <- c(Y = 90.0, La = 103.2, Ce = 101.0, Pr = 99.0, Nd = 98.3, Sm = 95.8, Eu = 94.7,
               Gd = 93.5, Tb = 92.3, Dy = 91.2, Ho = 90.1, Er = 89.0, Tm = 88.0, Yb = 86.8, Lu = 86.1)

# ---- input handling ----------------------------------------------------------------------------

#' Map column names to REE symbols ("La", "La_ppm", "Nd (ppb)", "Yb [mg/kg]" ...); ratios such as
#' "Ce/Ce*" or "La/Yb" are not taken for the element
.REE_COL_PATTERN <- "^[A-Za-z]{2}(?![A-Za-z]|\\s*/)"
.ree_detect_columns <- function(nms) {
  hit <- grepl(.REE_COL_PATTERN, nms, perl = TRUE)
  sym <- rep(NA_character_, length(nms))
  sym[hit] <- regmatches(nms, regexpr(.REE_COL_PATTERN, nms, perl = TRUE))
  idx <- match(tolower(sym), tolower(.REE))
  out <- setNames(nms[!is.na(idx)], .REE[idx[!is.na(idx)]])
  out[!duplicated(names(out))]
}

.ree_detect_y <- function(nms) {
  hit <- grep("^Y(?![A-Za-z]|\\s*/)", nms, perl = TRUE, value = TRUE)
  if (length(hit)) hit[1] else NULL
}

#' Unit basis from the units in the REE column names ("La [mg/kg]", "Nd_pmol/kg" ...)
.ree_units_basis <- function(cols) {
  u <- vapply(cols, function(nm) {
    b <- regmatches(nm, regexpr("\\[[^]]+\\]|\\([^)]+\\)", nm))
    if (length(b)) gsub("^[[(]|[])]$", "", b) else sub("^[A-Za-z]{2}[ _.]*", "", nm)
  }, character(1))
  cls <- ifelse(grepl("mol", u, ignore.case = TRUE), "molar",
         ifelse(grepl("ppm|ppb|ppt|g/g|g/kg|g/t|g g|mg|\u00b5g|ug|ng|pg", u, ignore.case = TRUE), "mass", NA_character_))
  if (all(is.na(cls))) return(NA_character_)
  if (any(is.na(cls)) || length(unique(cls)) > 1) return("mixed")
  cls[1]
}

.ree_number <- function(v) {
  if (is.numeric(v)) return(as.numeric(v))
  suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(v)), fixed = TRUE)))
}

#' Build the concentration table (Y + La..Lu) and the measured-value matrix from user data.
#' `cols` optionally maps element symbols to column names, e.g. c(La = "La139", Ce = "Ce140").
.ree_prepare <- function(data, cols = NULL) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  det <- .ree_detect_columns(names(data))
  if (!is.null(cols)) {
    bad <- setdiff(unname(cols), names(data))
    if (length(bad)) stop("Column(s) not found in `data`: ", paste(bad, collapse = ", "))
    det[names(cols)] <- cols
  }
  ycol <- if (!is.null(cols) && "Y" %in% names(cols)) cols[["Y"]] else .ree_detect_y(names(data))
  conc <- as.data.frame(matrix(NA_real_, nrow(data), length(.REE), dimnames = list(NULL, .REE)))
  for (el in intersect(names(det), .REE)) conc[[el]] <- .ree_number(data[[det[[el]]]])
  conc$Y <- if (!is.null(ycol)) .ree_number(data[[ycol]]) else NA_real_
  cm <- as.matrix(conc[, c("Y", .REE)])
  meas <- !is.na(cm) & cm > 0
  list(conc = conc, meas = meas, mapping = det[intersect(.REE, names(det))], ycol = ycol)
}

#' Concentration basis per row: a scalar, one value per row, or the name of a column of `data`;
#' NULL uses a column named "basis" (any case) if there is one, otherwise the units in the REE
#' column names ("La [mg/kg]" -> mass, "La [pmol/kg]" -> molar), otherwise "mass".
.ree_basis <- function(data, basis) {
  if (is.null(basis)) {
    bc <- names(data)[tolower(names(data)) == "basis"]
    if (length(bc)) basis <- data[[bc[1]]] else {
      ub <- .ree_units_basis(.ree_detect_columns(names(data)))
      if (identical(ub, "mixed")) warning("The REE columns have different units; `basis` = \"mass\" assumed. Check the units or give `basis`.")
      basis <- if (ub %in% c("mass", "molar")) ub else "mass"
    }
  } else if (length(basis) == 1 && basis %in% names(data)) basis <- data[[basis]]
  b <- rep_len(tolower(trimws(as.character(basis))), nrow(data))
  if (any(!b %in% c("mass", "molar")))
    stop("`basis` must be 'mass' or 'molar' (got: ", paste(unique(b[!b %in% c("mass", "molar")]), collapse = ", "), ").")
  b
}

.ree_norm_vector <- function(norm, norm_custom) {
  if (identical(norm, "custom")) {
    v <- norm_custom
    if (is.null(v) || is.null(names(v)) || any(!is.finite(v[.REE])) || any(v[.REE] <= 0))
      stop("`norm_custom` must be a named vector with positive values for all 14 REE (La-Lu), in ppm.")
    return(v[.REE])
  }
  if (!norm %in% names(.NORMALISERS))
    stop("Unknown `norm`: ", norm, ". Choose one of ", paste(c(names(.NORMALISERS), "custom"), collapse = ", "), ".")
  .NORMALISERS[[norm]]
}

.ree_norm_label <- function(norm) if (identical(norm, "custom")) "own reference" else unname(.NORM_LABELS[norm])

# ---- core calculations (identical to the app) ----------------------------------------------------

.ree_normalise <- function(conc, basis, nv) {
  m <- as.matrix(conc[, .REE, drop = FALSE]); storage.mode(m) <- "double"
  if (nrow(m) == 0) return(m)
  ref <- rbind(mass = nv[.REE], molar = nv[.REE] / .ATOMIC_MASS[.REE])
  n <- m / ref[basis, , drop = FALSE]
  n[!is.finite(n) | n <= 0] <- NA_real_
  dimnames(n) <- dimnames(m)
  n
}

.ree_zline <- function(x, el, a, b) {
  w <- (.REE_Z[[el]] - .REE_Z[[a]]) / (.REE_Z[[b]] - .REE_Z[[a]])
  unname(exp(log(x[, a]) + w * (log(x[, b]) - log(x[, a]))))
}

.ree_fill <- function(n) {
  rules <- list(
    Pr = function(x) { v <- .ree_zline(x, "Pr", "La", "Nd"); i <- is.na(v)
                       if (any(i)) v[i] <- .ree_zline(x[i, , drop = FALSE], "Pr", "Nd", "Sm"); v },
    Tb = function(x) .ree_zline(x, "Tb", "Gd", "Dy"),
    Ho = function(x) .ree_zline(x, "Ho", "Dy", "Er"),
    Tm = function(x) { v <- .ree_zline(x, "Tm", "Er", "Yb"); i <- is.na(v)
                       if (any(i)) v[i] <- .ree_zline(x[i, , drop = FALSE], "Tm", "Ho", "Yb"); v },
    Er = function(x) .ree_zline(x, "Er", "Ho", "Tm"),
    Lu = function(x) .ree_zline(x, "Lu", "Tm", "Yb"))
  filled <- matrix(FALSE, nrow(n), ncol(n), dimnames = dimnames(n))
  for (pass in 1:2) for (el in names(rules)) {
    miss <- is.na(n[, el]); if (!any(miss)) next
    v <- rules[[el]](n); ok <- miss & is.finite(v)
    n[ok, el] <- v[ok]; filled[ok, el] <- TRUE
  }
  list(n = n, filled = filled)
}

.ree_gd_star <- function(n, method) {
  switch(method,
    BD1996 = unname(0.33 * n[, "Sm"] + 0.67 * n[, "Tb"]),
    BDgeo  = unname(n[, "Sm"]^(1/3) * n[, "Tb"]^(2/3)),
    SmNd   = unname(n[, "Sm"]^2 / n[, "Nd"]),
    TbDy   = unname(n[, "Tb"]^2 / n[, "Dy"]),
    stop("Unknown `gd_method`: ", method))
}

.ree_roughness <- function(n) {
  d <- function(a, b, c) (n[, a] * 2 / (n[, b] + n[, c]) - 1)^2
  unname(sqrt(d("Nd", "Pr", "Sm") + d("Tb", "Gd", "Dy") + d("Dy", "Tb", "Ho") + d("Ho", "Dy", "Er") +
              d("Er", "Ho", "Tm") + d("Tm", "Er", "Yb") + d("Yb", "Tm", "Lu")))
}

.ree_check_columns <- function(present, lree, gap, gd, gd_method) {
  lre <- .LREE_SETS[[lree]]
  gd_inputs <- switch(gd_method, BD1996 = c("Sm", "Tb"), BDgeo = c("Sm", "Tb"), SmNd = c("Sm", "Nd"),
                      TbDy = c("Tb", "Dy"), BLambdaR = .REE)
  ignore <- identical(gd, "ignore")
  need <- list()
  if (identical(gap, "interpolate")) {
    need <- c(need, as.list(if (lree == "LaPrNd") c("La", "Nd") else c("Nd", "Sm")), "Yb")
    need <- c(need, if (!ignore) list("Gd") else as.list(gd_inputs))
    need <- c(need, list("Dy"), list(c("Er", "Ho", "Tm")))
  } else if (identical(gap, "groupmean")) {
    need <- list(lre, if (ignore) c("Tb", "Dy") else .MREE_SET, .HREE_SET)
    if (ignore) need <- c(need, as.list(gd_inputs))
  } else {
    need <- c(as.list(lre), as.list(.HREE_SET), as.list(setdiff(.MREE_SET, if (ignore) "Gd")))
    if (ignore) need <- c(need, as.list(gd_inputs))
  }
  miss <- vapply(unique(need), function(alt) if (any(alt %in% present)) NA_character_ else
    if (length(alt) == 1) alt else paste0("one of ", paste(alt, collapse = "/")), character(1))
  unique(miss[!is.na(miss)])
}

# ---- lambda (orthogonal polynomial) fit: base-R port of BLambdaR ---------------------------------

.LAMBDA_RADII <- c(La = 1.160, Ce = 1.143, Pr = 1.126, Nd = 1.109, Sm = 1.079, Eu = 1.066,
                   Gd = 1.053, Tb = 1.040, Dy = 1.027, Ho = 1.015, Er = 1.004, Tm = 0.994,
                   Yb = 0.985, Lu = 0.977)                 # O'Neill (2016), Angstrom

.LAMBDA <- local({
  radii <- .LAMBDA_RADII; N <- length(radii)
  a <- sum(radii); b <- sum(radii^2); c <- sum(radii^3); d <- sum(radii^4)
  e <- sum(radii^5); f <- sum(radii^6); g <- sum(radii^7)
  bet <- mean(radii)
  v <- solve(matrix(c(-a, N, -b, a), 2, byrow = TRUE), c(-b, -c))
  gam <- Re(polyroot(c(v[2], -v[1], 1)))
  v <- solve(matrix(c(-b, a, -N, -c, b, -a, -d, c, -b), 3, byrow = TRUE), c(-c, -d, -e))
  del <- Re(polyroot(c(-v[3], v[2], -v[1], 1)))
  v <- solve(matrix(c(-c, b, -a, N, -d, c, -b, a, -e, d, -c, b, -f, e, -d, c), 4, byrow = TRUE),
             c(-d, -e, -f, -g))
  eps <- Re(polyroot(c(v[4], -v[3], v[2], -v[1], 1)))
  f1 <- function(x) x - bet
  f2 <- function(x) (x - gam[1]) * (x - gam[2])
  f3 <- function(x) (x - del[1]) * (x - del[2]) * (x - del[3])
  f4 <- function(x) (x - eps[1]) * (x - eps[2]) * (x - eps[3]) * (x - eps[4])
  z <- (57:71)[-5]
  fit <- lm(z ~ poly(radii, 7))
  list(f = list(f1, f2, f3, f4), fmat = cbind(f1(radii), f2(radii), f3(radii), f4(radii)),
       trn = function(r) predict(fit, newdata = data.frame(radii = r)))
})

.ree_tetrad <- function(x) (1 - x^2) / 2 + sqrt(((1 - x^2)^2) / 4)
.ree_tf <- function(x, t) .ree_tetrad((x - (55.25 + 3.5 * t)) / 1.75)

.ree_blambda <- function(conc, basis, anomalies = c(La = FALSE, Ce = TRUE, Eu = TRUE, Gd = FALSE),
                         fit_factors = c(l4 = TRUE, t1 = TRUE, t2 = TRUE, t3 = TRUE, t4 = TRUE),
                         min_points = 10) {
  L <- .LAMBDA; radii <- .LAMBDA_RADII; norm <- .NORMALISERS$CI_ONeill2016
  an <- c(La = FALSE, Ce = TRUE, Eu = TRUE, Gd = FALSE); an[names(anomalies)] <- anomalies
  ff <- c(l4 = TRUE, t1 = FALSE, t2 = FALSE, t3 = FALSE, t4 = FALSE); ff[names(fit_factors)] <- fit_factors
  m <- as.matrix(conc[, .REE, drop = FALSE]); storage.mode(m) <- "double"
  m[!is.finite(m) | m <= 0] <- NA
  mol <- basis == "molar"
  if (any(mol)) m[mol, ] <- sweep(m[mol, , drop = FALSE], 2, .ATOMIC_MASS[.REE], "*")
  nm <- c("lambda0", "lambda1", "lambda2", "lambda3", "lambda4", "tau1", "tau2", "tau3", "tau4")
  npar <- 1 + length(nm) + length(.REE)
  excl <- c(if (an[["Gd"]]) 7, if (an[["Eu"]]) 6, if (an[["Ce"]]) 2, if (an[["La"]]) 1)
  fit_one <- function(ppm) {
    k <- log(ppm / norm)
    if (!(sum(!is.na(k[3:4])) > 0 && sum(!is.na(k[7:9])) > 0 && sum(!is.na(k[12:14])) > 0))
      return(rep(NA_real_, npar))
    rad.fit <- radii
    if (length(excl)) { k <- k[-excl]; rad.fit <- rad.fit[-excl] }
    if (sum(is.finite(k)) < min_points) return(rep(NA_real_, npar))
    X <- cbind(L$f[[1]](rad.fit), L$f[[2]](rad.fit), L$f[[3]](rad.fit))
    if (ff[["l4"]]) X <- cbind(X, L$f[[4]](rad.fit))
    zz <- L$trn(rad.fit)
    for (t in 1:4) if (ff[[paste0("t", t)]]) X <- cbind(X, .ree_tf(zz, t))
    ok <- is.finite(k)
    A <- cbind(1, X)[ok, , drop = FALSE]
    fit <- lm.fit(A, k[ok])
    cf <- fit$coefficients; cf[is.na(cf)] <- NA_real_
    tn <- paste0("t", 1:4)
    pos <- 1 + 3 + (if (ff[["l4"]]) 1 else 0)
    res <- c(cf[1:4], if (ff[["l4"]]) cf[5] else NA)
    j <- pos
    for (t in 1:4) if (ff[[tn[t]]]) { j <- j + 1; res <- c(res, cf[j]) } else res <- c(res, NA_real_)
    res <- unname(res)
    r0 <- ifelse(is.na(res), 0, res)
    zr <- L$trn(radii)
    pat <- exp(r0[1] + as.vector(L$fmat %*% r0[2:5]) +
                 r0[6] * .ree_tf(zr, 1) + r0[7] * .ree_tf(zr, 2) + r0[8] * .ree_tf(zr, 3) + r0[9] * .ree_tf(zr, 4))
    rss <- sum(fit$residuals^2); nn <- sum(ok); pp <- sum(!is.na(cf))
    tss <- sum((k[ok] - mean(k[ok]))^2)
    adj <- if (nn - pp > 0 && tss > 0) 1 - (rss / (nn - pp)) / (tss / (nn - 1)) else NA_real_
    c(adj, res, ppm / (pat * norm))
  }
  out <- t(apply(m, 1, function(r) tryCatch(fit_one(r), error = function(e) rep(NA_real_, npar))))
  if (nrow(m) == 0) out <- matrix(NA_real_, 0, npar)
  out <- as.data.frame(out)
  names(out) <- c("R2", nm, paste0(.REE, "_MFR"))
  mf <- out[, paste0(c("Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"), "_MFR"), drop = FALSE]
  out$lambda_rough <- sqrt(rowSums((mf - 1)^2)) / 7
  out
}

# ---- screening (as in the app's "Screening" section) ------------------------------------------

.ree_screen <- function(x, y, rough, meas, rough_max, min_ree, req_ree) {
  st <- rep("ok", length(x))
  st[!is.finite(x) | !is.finite(y)] <- "not computable (missing values)"
  if ((min_ree %||% 0) > 0 || length(req_ree)) {
    bad <- rowSums(meas[, .REE, drop = FALSE]) < (min_ree %||% 0)
    if (length(req_ree)) {
      if (any(!req_ree %in% .REE)) stop("`req_ree` must contain REE symbols (La..Lu).")
      bad <- bad | !apply(meas[, req_ree, drop = FALSE], 1, all)
    }
    st[bad] <- "excluded (too few REE)"
  }
  if (!is.null(rough_max)) st[st == "ok" & is.finite(rough) & rough > rough_max] <- "excluded (roughness)"
  st
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# =================================================================================================
# User functions
# =================================================================================================

#' Read a data file: CSV (comma or semicolon with decimal commas), TSV/TXT, or a PANGAEA .tab file
#'
#' For PANGAEA files the metadata block ("/* DATA DESCRIPTION: ... */") is skipped and the dataset
#' citation, DOI and licence are returned in attr(x, "source"). Column names are kept as they are
#' (e.g. "La [mg/kg]"), which the other functions recognise, including the units.
#' @param path file path
#' @return data.frame
REE_read <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines)) lines[1] <- sub("^\\ufeff", "", lines[1])
  first <- which(nzchar(trimws(lines)))[1]
  if (is.na(first)) stop("The file is empty.")
  src <- NULL
  if (startsWith(trimws(lines[first]), "/*")) {
    end <- which(trimws(lines) == "*/" | grepl("\\*/\\s*$", lines))
    end <- end[end >= first][1]
    if (is.na(end)) stop("The file starts with a '/*' header that is never closed ('*/').")
    h <- lines[first:end]
    key <- sub("^([^\t:]+):\t.*$", "\\1", h); key[!grepl("^[^\t:]+:\t", h)] <- NA_character_
    field <- function(k) {
      i <- which(key == k)[1]; if (is.na(i)) return(NA_character_)
      j <- i
      while (j < length(h) && is.na(key[j + 1]) && startsWith(h[j + 1], "\t")) j <- j + 1
      trimws(gsub("\\s+", " ", paste(sub("^[^\t]*\t", "", h[i:j]), collapse = " ")))
    }
    cit <- sub("[ ,;]+$", "", sub("\\s+In: .*$", "", field("Citation")))
    doi <- regmatches(cit, regexpr("10\\.[0-9]{4,}/[^ ,;]+", cit))
    src <- list(type = "PANGAEA", citation = cit, doi = if (length(doi)) sub("[.,]$", "", doi) else NA_character_,
                licence = field("License"))
    lines <- lines[-seq_len(end)]
    first <- which(nzchar(trimws(lines)))[1]
    if (is.na(first)) stop("The file contains no data table after its header.")
  }
  lines <- lines[first:length(lines)]
  hdr <- lines[1]
  seps <- c(`,` = lengths(regmatches(hdr, gregexpr(",", hdr))), `;` = lengths(regmatches(hdr, gregexpr(";", hdr))),
            `\t` = lengths(regmatches(hdr, gregexpr("\t", hdr))))
  sep <- if (!is.null(src)) "\t" else names(which.max(seps))
  con <- textConnection(lines, encoding = "UTF-8"); on.exit(close(con))
  df <- utils::read.table(con, header = TRUE, sep = sep, dec = if (sep == ";") "," else ".",
                          quote = if (is.null(src)) "\"" else "", comment.char = "", stringsAsFactors = FALSE,
                          check.names = FALSE, na.strings = c("NA", "", "n.d.", "nd", "NaN"), strip.white = TRUE)
  if (!is.null(src)) attr(df, "source") <- src
  df
}

#' REE ratios, anomalies, anthropogenic Gd and REE roughness
#'
#' Same calculation and options as the REE cross-plot app (sidebar sections "Your data",
#' "Normalisation & ratios", "Screening" and "Anthropogenic Gd").
#'
#' @param data data.frame with REE concentrations. Columns are recognised when their name starts
#'   with the element symbol ("La", "La_ppm", "Nd (pmol/kg)" ...); see `cols`.
#' @param basis "mass" (e.g. ppm, ppb, ng/g) or "molar" (e.g. pmol/kg); a single value, one value
#'   per row, or the name of a column. NULL (default) uses a column named "basis" if present,
#'   otherwise the units in the REE column names (e.g. "La [mg/kg]", as in PANGAEA files),
#'   otherwise "mass".
#' @param norm normalising reference: "PAAS_Pourmand2012" (default), "PAAS_McLennan1989",
#'   "NASC_McLennan1989", "EUS_Bau2018", "CI_Pourmand2012", "CI_ONeill2016" or "custom".
#' @param norm_custom with norm = "custom": named vector of the 14 REE in ppm (mass basis).
#' @param lree LREE definition: "LaPrNd" (La + Pr + Nd; Martin et al. 2010, Blaser et al. 2020)
#'   or "PrNdSm" (Pr + Nd + Sm).
#' @param gap missing Pr, Tb, Ho, Er, Tm, Lu: "interpolate" (log-linear line through the
#'   neighbours as a function of atomic number; Barrat and Bayon 2024), "groupmean" (mean of the
#'   measured members of each REE group instead of sums) or "none".
#' @param ce_form Ce/Ce*: "geoPrNd" (Ce Nd / Pr^2; Barrat and Bayon 2024, default), "arithLaPr"
#'   (2 Ce / (La + Pr)) or "geoLaPr" (Ce / sqrt(La Pr)).
#' @param gd "measured" (use measured Gd) or "ignore" (replace Gd by gd_baseline * Gd*, for
#'   samples with anthropogenic Gd).
#' @param gd_method Gd*: "BD1996" (0.33 Sm + 0.67 Tb; Bau and Dulski 1996), "BDgeo"
#'   (Sm^(1/3) Tb^(2/3)), "SmNd" (Sm^2/Nd), "TbDy" (Tb^2/Dy) or "BLambdaR" (polynomial fit
#'   without Gd, Ce and Eu; experimental).
#' @param gd_baseline natural Gd/Gd* (1 = none); anthropogenic Gd = Gd * max(0, 1 - baseline / (Gd/Gd*)).
#' @param rough_max optional REE roughness limit (the app uses 0.3); rougher samples are flagged
#'   in `screen`. Samples without a roughness value are not flagged.
#' @param min_ree optional minimum number of measured REE (of 14); fewer -> flagged in `screen`.
#' @param req_ree optional REE that must be measured, e.g. c("Ce", "Pr").
#' @param cols optional named vector mapping element symbols to column names, e.g.
#'   c(La = "La139", Ce = "Ce140"), for columns that are not recognised automatically.
#' @param append return the input columns together with the results (default) or the results only.
#' @param settings add columns stating the calculation settings (as in the app's results table).
#' @return data.frame, one row per input row: row, HREE_LREE, MREE_MREEstar, Ce_Cestar, La_Lastar,
#'   Eu_Eustar, Gd_Gdstar, Gd_anth (input units), REE_rough, n_estimated, n_REE_measured, basis,
#'   screen ("ok", "excluded (too few REE)", "excluded (roughness)", "not computable (missing values)").
#'   All anomalies are computed from the normalised pattern; Gd_Gdstar always uses measured Gd.
#' @examples
#'   res <- REE_ratios(mydata, basis = "molar", norm = "EUS_Bau2018", rough_max = 0.3)
#'   plot(res$HREE_LREE, res$MREE_MREEstar, log = "xy")
REE_ratios <- function(data, basis = NULL, norm = "PAAS_Pourmand2012", norm_custom = NULL,
                       lree = c("LaPrNd", "PrNdSm"), gap = c("interpolate", "groupmean", "none"),
                       ce_form = c("geoPrNd", "arithLaPr", "geoLaPr"), gd = c("measured", "ignore"),
                       gd_method = c("BD1996", "BDgeo", "SmNd", "TbDy", "BLambdaR"), gd_baseline = 1,
                       rough_max = NULL, min_ree = 0, req_ree = NULL, cols = NULL,
                       append = TRUE, settings = FALSE) {
  lree <- match.arg(lree); gap <- match.arg(gap); ce_form <- match.arg(ce_form)
  gd <- match.arg(gd); gd_method <- match.arg(gd_method)
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  pr <- .ree_prepare(data, cols)
  miss <- .ree_check_columns(names(pr$mapping), lree, gap, gd, gd_method)
  if (length(miss)) stop("Missing REE column(s): ", paste(miss, collapse = ", "),
                         ". Column names must start with the element symbol, or use `cols`.")
  b <- .ree_basis(data, basis)
  nv <- .ree_norm_vector(norm, norm_custom)
  conc <- pr$conc
  n <- .ree_normalise(conc, b, nv)
  nr <- nrow(n)
  base <- if (is.null(gd_baseline) || !is.finite(gd_baseline)) 1 else gd_baseline

  # Gd anomaly of the measured pattern and anthropogenic Gd
  gs <- if (identical(gd_method, "BLambdaR")) {
    lg <- .ree_blambda(conc, b, anomalies = c(Gd = TRUE),
                       fit_factors = c(l4 = TRUE, t1 = FALSE, t2 = FALSE, t3 = FALSE, t4 = FALSE),
                       min_points = 8)$Gd_MFR
    unname(n[, "Gd"] / lg)
  } else .ree_gd_star(n, gd_method)
  anom <- unname(n[, "Gd"] / gs)
  gd_anth <- ifelse(is.finite(anom), conc$Gd * pmax(0, 1 - base / anom), NA_real_)
  if (identical(gd, "ignore")) n[, "Gd"] <- gs * base

  rough <- .ree_roughness(n)
  filled <- matrix(FALSE, nr, ncol(n), dimnames = dimnames(n))
  if (identical(gap, "interpolate")) { f <- .ree_fill(n); n <- f$n; filled <- f$filled }

  grp <- if (identical(gap, "groupmean")) {
    function(els) { v <- rowMeans(n[, els, drop = FALSE], na.rm = TRUE); v[is.nan(v)] <- NA_real_; v }
  } else function(els) rowSums(n[, els, drop = FALSE])
  L <- grp(.LREE_SETS[[lree]]); H <- grp(.HREE_SET); M <- grp(.MREE_SET)
  n_est <- if (identical(gap, "groupmean")) {
    unname(rowSums(is.na(n[, unique(c(.LREE_SETS[[lree]], .HREE_SET, .MREE_SET)), drop = FALSE])))
  } else unname(rowSums(filled[, c("Pr", "Tb", "Tm", "Lu"), drop = FALSE]))

  out <- data.frame(
    row = seq_len(nr),
    HREE_LREE = unname(H / L),
    MREE_MREEstar = unname(M / ((H + L) / 2)),
    Ce_Cestar = unname(switch(ce_form,
                              arithLaPr = 2 * n[, "Ce"] / (n[, "La"] + n[, "Pr"]),
                              geoLaPr   = n[, "Ce"] / sqrt(n[, "La"] * n[, "Pr"]),
                              geoPrNd   = n[, "Ce"] * n[, "Nd"] / n[, "Pr"]^2)),
    La_Lastar = unname(n[, "La"] * n[, "Nd"]^2 / n[, "Pr"]^3),
    Eu_Eustar = unname(n[, "Eu"] / sqrt(n[, "Sm"] * n[, "Gd"])),
    Gd_Gdstar = anom,
    Gd_anth = unname(gd_anth),
    REE_rough = rough,
    n_estimated = as.integer(n_est),
    n_REE_measured = as.integer(rowSums(pr$meas[, .REE, drop = FALSE])),
    basis = b, stringsAsFactors = FALSE)
  out$screen <- .ree_screen(out$HREE_LREE, out$MREE_MREEstar, out$REE_rough, pr$meas,
                            rough_max, min_ree, req_ree)
  if (settings) {
    out$normalisation <- .ree_norm_label(norm)
    out$LREE_definition <- c(LaPrNd = "La-Pr-Nd", PrNdSm = "Pr-Nd-Sm")[[lree]]
    out$gap_method <- c(interpolate = "interpolate", groupmean = "group mean", none = "none")[[gap]]
    out$Gd_handling <- if (identical(gd, "ignore")) "ignored (replaced by estimate)" else "measured"
    out$Gd_star_method <- gd_method
    out$Gd_baseline <- base
    out$Ce_anomaly <- c(geoPrNd = "Ce Nd / Pr^2", arithLaPr = "2 Ce / (La + Pr)", geoLaPr = "Ce / sqrt(La Pr)")[[ce_form]]
  }
  attr(out, "column_mapping") <- pr$mapping
  if (append) {
    keep <- setdiff(names(data), names(out))
    res <- cbind(data[, keep, drop = FALSE], out, stringsAsFactors = FALSE)
    attr(res, "column_mapping") <- pr$mapping
    return(res)
  }
  out
}

#' Lambda shape parameters (orthogonal polynomials of the ionic radius) and lambda roughness
#'
#' Base-R port of BLambdaR (Anenburg and Williams, 2022; https://lambdar.rses.anu.edu.au/blambdar/),
#' as in the app's lambda space: patterns normalised to CI chondrite (O'Neill, 2016), log-
#' transformed, fitted with lambda0..lambda4 and the tetrad terms tau1..tau4. Ce and Eu (and Gd
#' with gd = "ignore") are left out of the fit. The app plots -lambda1 (x) against -lambda2 (y).
#'
#' @param data,basis,cols as in REE_ratios(). The normalisation does not apply here.
#' @param gd "measured" (Gd included in the fit) or "ignore" (Gd left out), as in the app.
#' @param anomalies optional override of the elements left out of the fit, e.g.
#'   c(La = TRUE, Ce = TRUE, Eu = TRUE, Gd = FALSE).
#' @param tetrads fit the four tetrad terms (TRUE, as in the app) or polynomials only.
#' @param min_points minimum number of usable elements after leaving out the anomalous ones
#'   (default 10: the full fit has 9 parameters).
#' @param rough_max optional lambda roughness limit (the app uses 0.03); see REE_ratios().
#' @param min_ree minimum number of measured REE for `screen` (the app uses 12 in lambda space).
#' @param req_ree optional REE that must be measured.
#' @param mfr also return the measured/fitted ratios <El>_MFR of every element.
#' @param append return the input columns together with the results (default).
#' @return data.frame: row, R2 (adjusted), lambda0..lambda4, tau1..tau4, lambda_rough
#'   (scatter of Tb-Lu around the fit), n_REE_measured, basis, screen, and optionally <El>_MFR.
REE_lambda <- function(data, basis = NULL, gd = c("measured", "ignore"), anomalies = NULL,
                       tetrads = TRUE, min_points = 10, rough_max = NULL, min_ree = 12,
                       req_ree = NULL, mfr = FALSE, cols = NULL, append = TRUE) {
  gd <- match.arg(gd)
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  pr <- .ree_prepare(data, cols)
  if (!length(pr$mapping)) stop("No REE columns found in `data` (names must start with the element symbol, or use `cols`).")
  b <- .ree_basis(data, basis)
  an <- c(Gd = identical(gd, "ignore"))
  if (!is.null(anomalies)) an <- c(an[setdiff(names(an), names(anomalies))], anomalies)
  ff <- if (tetrads) c(l4 = TRUE, t1 = TRUE, t2 = TRUE, t3 = TRUE, t4 = TRUE) else
    c(l4 = TRUE, t1 = FALSE, t2 = FALSE, t3 = FALSE, t4 = FALSE)
  bl <- .ree_blambda(pr$conc, b, anomalies = an, fit_factors = ff, min_points = min_points)
  out <- data.frame(row = seq_len(nrow(bl)),
                    bl[, c("R2", "lambda0", "lambda1", "lambda2", "lambda3", "lambda4",
                           "tau1", "tau2", "tau3", "tau4", "lambda_rough")],
                    n_REE_measured = as.integer(rowSums(pr$meas[, .REE, drop = FALSE])),
                    basis = b, stringsAsFactors = FALSE)
  out$screen <- .ree_screen(-out$lambda1, -out$lambda2, out$lambda_rough, pr$meas,
                            rough_max, min_ree, req_ree)
  if (mfr) out <- cbind(out, bl[, paste0(.REE, "_MFR")])
  rownames(out) <- NULL
  attr(out, "column_mapping") <- pr$mapping
  if (append) {
    keep <- setdiff(names(data), names(out))
    res <- cbind(data[, keep, drop = FALSE], out, stringsAsFactors = FALSE)
    attr(res, "column_mapping") <- pr$mapping
    return(res)
  }
  out
}

#' Normalised REE patterns (spectra): table and figure
#'
#' As the app's Spectra tab: concentrations normalised to the chosen reference, optionally
#' divided by one element (secondary normalisation), plotted against the element (atomic-number
#' order, gap at Pm) or the ionic radius (Shannon, 1976, 3+, 6-coordinated; Y as a separate point).
#'
#' @param data,basis,norm,norm_custom,cols as in REE_ratios(); a column starting with "Y" is used
#'   for yttrium (only references with a Y value: PAAS and EUS).
#' @param norm2 optional element for the secondary normalisation, e.g. "Nd".
#' @param keep optional logical vector selecting the rows to plot (e.g.
#'   REE_ratios(...)$screen == "ok"); the returned table always has every row.
#' @param colour optional column used for the colours: categorical (one colour per value) or
#'   numeric (continuous scale, e.g. age or depth). NULL = one colour.
#' @param id optional column used as sample label in the long table.
#' @param x "element" or "radius".
#' @param log logarithmic y axis.
#' @param show_y show Y (if available).
#' @param connect connect points across missing elements (TRUE) or break the line there (FALSE).
#' @param points draw points on top of the lines.
#' @param legend,line_size,alpha,base_size appearance.
#' @param plot print the figure (TRUE) or only return it.
#' @return invisibly, a list: `plot` (ggplot object), `normalised` (one row per input row, columns
#'   norm_<El> and, with norm2, norm2_<El>), `long` (plotting table).
REE_spectrum <- function(data, basis = NULL, norm = "PAAS_Pourmand2012", norm_custom = NULL,
                         norm2 = NULL, keep = NULL, colour = NULL, id = NULL,
                         x = c("element", "radius"), log = TRUE, show_y = TRUE, connect = TRUE,
                         points = TRUE, legend = TRUE, line_size = 0.5, alpha = 0.8, base_size = 14,
                         cols = NULL, plot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("REE_spectrum() needs the ggplot2 package.")
  x <- match.arg(x)
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  pr <- .ree_prepare(data, cols)
  if (!length(pr$mapping)) stop("No REE columns found in `data`.")
  b <- .ree_basis(data, basis)
  nv <- .ree_norm_vector(norm, norm_custom)
  yn <- if (identical(norm, "custom")) {
    if (!is.null(norm_custom) && "Y" %in% names(norm_custom) && is.finite(norm_custom[["Y"]]) &&
        norm_custom[["Y"]] > 0) norm_custom[["Y"]] else NA_real_
  } else unname(.Y_NORM[norm])
  els <- c("Y", .REE)

  # normalisation (molar data taken as pmol/kg; the absolute scale is kept, as in the app)
  m <- as.matrix(pr$conc[, els, drop = FALSE]); storage.mode(m) <- "double"
  nvy <- c(Y = yn, nv); mass <- c(Y = .Y_MASS, .ATOMIC_MASS)
  ref <- rbind(mass = nvy[els], molar = nvy[els] * 1e9 / mass[els])
  n <- m / ref[ifelse(b == "molar", "molar", "mass"), , drop = FALSE]
  n[!is.finite(n) | n <= 0] <- NA_real_
  colnames(n) <- els
  n2 <- NULL
  if (!is.null(norm2)) {
    if (!norm2 %in% els) stop("`norm2` must be one of ", paste(els, collapse = ", "), ".")
    n2 <- n / n[, norm2]; n2[!is.finite(n2) | n2 <= 0] <- NA_real_
  }
  vals <- if (is.null(n2)) n else n2

  # table with every row
  outc <- if (any(is.finite(n[, "Y"]))) els else .REE
  tab <- data.frame(row = seq_len(nrow(n)), basis = b, normalisation = .ree_norm_label(norm),
                    secondary_normalisation = if (is.null(norm2)) "none" else norm2, stringsAsFactors = FALSE)
  nn <- n[, outc, drop = FALSE]; colnames(nn) <- paste0("norm_", outc); tab <- cbind(tab, nn)
  if (!is.null(n2)) { n2m <- n2[, outc, drop = FALSE]; colnames(n2m) <- paste0("norm2_", outc); tab <- cbind(tab, n2m) }

  # long table for plotting
  keep <- if (is.null(keep)) rep(TRUE, nrow(n)) else as.logical(rep_len(keep, nrow(n)))
  keep[is.na(keep)] <- FALSE
  sample <- if (!is.null(id) && id %in% names(data)) as.character(data[[id]]) else sprintf("row%03d", seq_len(nrow(n)))
  cv <- if (!is.null(colour)) {
    if (!colour %in% names(data)) stop("`colour` column not found: ", colour)
    v <- data[[colour]]; vn <- .ree_number(v)
    if (is.numeric(v) || (mean(is.finite(vn)) >= 0.5 && length(unique(vn[is.finite(vn)])) > 2)) vn else
      { v <- as.character(v); v[is.na(v) | v == ""] <- "(no group)"; v }
  } else rep("my data", nrow(n))
  numeric_col <- is.numeric(cv)
  idx <- which(keep)
  pos <- c(Y = 0, La = 1, Ce = 2, Pr = 3, Nd = 4, Sm = 6, Eu = 7, Gd = 8, Tb = 9, Dy = 10,
           Ho = 11, Er = 12, Tm = 13, Yb = 14, Lu = 15)
  long <- data.frame(sid = rep(idx, times = length(els)), sample = rep(sample[idx], times = length(els)),
                     .c = rep(cv[idx], times = length(els)), element = rep(els, each = length(idx)),
                     value = as.vector(vals[idx, els, drop = FALSE]), stringsAsFactors = FALSE)
  long$pos <- unname(pos[long$element]); long$radius <- unname(.RADII_PM[long$element])

  # figure (same construction as the app)
  d <- long
  if (!show_y || !is.finite(yn)) d <- d[d$element != "Y", ]
  d$xv <- if (x == "radius") d$radius else d$pos
  d <- d[is.finite(d$xv), ]
  d$ok <- is.finite(d$value) & d$value > 0
  if (!numeric_col) {
    lv <- unique(as.character(cv))
    key <- vapply(lv, function(s) {
      parts <- regmatches(s, gregexpr("[0-9]+|[^0-9]+", s))[[1]]
      paste(ifelse(grepl("^[0-9]+$", parts),
                   formatC(suppressWarnings(as.numeric(parts)), width = 15, flag = "0", format = "f", digits = 0),
                   tolower(parts)), collapse = "")
    }, character(1))
    d$.c <- factor(d$.c, levels = lv[order(key, lv)])
  }
  ld <- d[d$element != "Y" | x == "element", ]
  ld <- if (connect) ld[ld$ok, ] else { ld$value[!ld$ok] <- NA_real_; ld }
  ld <- ld[order(ld$sid, ld$xv), ]
  pd <- d[d$ok & (points | d$element == "Y"), ]
  user_cols <- c("#E69F00", "#D55E00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#999999")
  fill_cols <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
  ref_short <- sub(" - .*$", "", .ree_norm_label(norm))
  ylab <- if (is.null(norm2)) sprintf("REE / %s", ref_short) else
    sprintf("(REE / %s) / (%s / %s)", ref_short, norm2, ref_short)

  p <- ggplot2::ggplot(d[d$ok, ], ggplot2::aes(xv, value))
  if (nrow(ld)) p <- p + ggplot2::geom_path(data = ld, ggplot2::aes(group = sid, colour = .c),
                                           linewidth = line_size, alpha = alpha, na.rm = TRUE)
  if (nrow(pd)) p <- p + ggplot2::geom_point(data = pd, ggplot2::aes(colour = .c),
                                            size = line_size * 3 + 0.6, alpha = alpha)
  if (numeric_col) {
    p <- p + ggplot2::scale_colour_gradientn(colours = fill_cols, na.value = "grey70", name = colour %||% "value")
  } else {
    ng <- nlevels(d$.c)
    cl <- if (ng <= length(user_cols)) user_cols[seq_len(ng)] else grDevices::hcl.colors(ng, "Dark 3")
    p <- p + ggplot2::scale_colour_manual(values = setNames(cl, levels(d$.c)), name = NULL, drop = FALSE)
  }
  p <- p + ggplot2::labs(y = ylab, x = if (x == "radius") "ionic radius (pm)" else NULL) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                   axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black"))
  if (x == "radius") p <- p + ggplot2::scale_x_reverse() else {
    br <- c(pos[c("Y", "La", "Ce", "Pr", "Nd")], Pm = 5, pos[c("Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu")])
    if (!show_y || !is.finite(yn)) br <- br[names(br) != "Y"]
    p <- p + ggplot2::scale_x_continuous(breaks = unname(br), labels = names(br))
  }
  if (log) p <- p + ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides = "l")
  if (!legend || (!numeric_col && nlevels(d$.c) <= 1 && is.null(colour))) p <- p + ggplot2::theme(legend.position = "none")
  if (plot) print(p)
  long$.c <- NULL
  if (!is.null(colour)) long[[colour]] <- rep(cv[idx], times = length(els))
  invisible(list(plot = p, normalised = tab, long = long))
}
