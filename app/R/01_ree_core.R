# ---------------------------------------------------------------------------
# REE normalisation and ratio calculations
#
# Default definitions (Martin et al., 2010; Blaser et al., 2020):
#   HREE/LREE     = (Tm + Yb + Lu) / (La + Pr + Nd)
#   MREE/MREE*    = (Gd + Tb + Dy) / mean(HREE, LREE)      [HREE, LREE = sums above]
#   Ce/Ce*        = 2 * Ce / (La + Pr)
# with every element first normalised to a shale (default PAAS, Pourmand et al., 2012).
# Options (see `default_cfg()`): normalisation standard, LREE definition (La-Pr-Nd or
# Pr-Nd-Sm), gap handling, and how a (possibly anthropogenic) Gd is treated.
# ---------------------------------------------------------------------------

REE <- c("La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu")

# Standard atomic weights (g/mol)
ATOMIC_MASS <- c(La = 138.905, Ce = 140.116, Pr = 140.908, Nd = 144.242, Sm = 150.36,
                 Eu = 151.964, Gd = 157.25, Tb = 158.925, Dy = 162.500, Ho = 164.930,
                 Er = 167.259, Tm = 168.934, Yb = 173.045, Lu = 174.967)

# ---- normalisation standards (ppm = ug/g, mass basis) ----------------------
mk <- function(...) setNames(c(...), REE)
NORMALISERS <- list(
  PAAS_Pourmand2012 = mk(44.56, 88.25, 10.15, 37.32, 6.884, 1.215, 6.043, 0.8914, 5.325, 1.053,
                         3.075, 0.451, 3.012, 0.4386),
  PAAS_McLennan1989 = mk(38.2, 79.6, 8.83, 33.9, 5.55, 1.08, 4.66, 0.774, 4.68, 0.991,
                         2.85, 0.405, 2.82, 0.433),
  NASC_McLennan1989 = mk(32, 73, 7.9, 33, 5.7, 1.24, 5.2, 0.85, 5.8, 1.04, 3.4, 0.5, 3.1, 0.48),
  EUS_Bau2018       = mk(44.3, 88.5, 10.6, 39.5, 7.30, 1.48, 6.34, 0.944, 5.86, 1.17,
                         3.43, 0.492, 3.26, 0.485),
  CI_Pourmand2012   = mk(0.2469, 0.6321, 0.0959, 0.4854, 0.1556, 0.0599, 0.2093, 0.0378, 0.2577,
                         0.0554, 0.1667, 0.0261, 0.1694, 0.0256),
  CI_ONeill2016     = mk(0.2472, 0.6308, 0.0950, 0.4793, 0.1542, 0.0592, 0.2059, 0.0375, 0.2540,
                         0.0554, 0.1645, 0.0258, 0.1684, 0.0251)
)
NORMALISER_LABELS <- c(
  PAAS_Pourmand2012 = "PAAS - Pourmand et al. (2012)",
  PAAS_McLennan1989 = "PAAS - McLennan (1989)",
  EUS_Bau2018       = "EUS (European Shale) - Bau et al. (2018)",
  NASC_McLennan1989 = "NASC - McLennan (1989)",
  CI_Pourmand2012   = "CI chondrite - Pourmand et al. (2012)",
  CI_ONeill2016     = "CI chondrite - O'Neill (2016)",
  custom            = "Own reference (upload one row)"
)
DEFAULT_NORM <- "PAAS_Pourmand2012"
rm(mk)

# Kept for the end-member/regression tests: PAAS as used for the published figure
PAAS_PPM <- NORMALISERS$PAAS_McLennan1989

# Scale used when a molar end-member (seawater) is expressed in mass units for mixing with a
# mass-based end-member (nodule). Arbitrary weighting that places the bend of the authigenic
# array exactly as in the published figure (pmol/kg * g/mol * 1e-3).
SEAWATER_TO_MASS_SCALE <- 1e-3

# Yttrium (only used for the spectra tab; not part of the ratios). Y of the shale standards as
# used in P. Blaser's REE_spectrum() function; NA = no value available.   # VERIFY against sources
Y_MASS <- 88.906
Y_NORM <- c(PAAS_Pourmand2012 = 27.31, PAAS_McLennan1989 = 27, EUS_Bau2018 = 31.9,
            NASC_McLennan1989 = NA_real_, CI_Pourmand2012 = NA_real_, CI_ONeill2016 = NA_real_)
SPEC_ELEMENTS <- c("Y", REE)

#' Y of the chosen normalising reference (NA if not available)
norm_y <- function(cfg = default_cfg()) {
  if (identical(cfg$norm, "custom")) {
    v <- cfg$norm_custom
    return(if (!is.null(v) && "Y" %in% names(v) && is.finite(v[["Y"]]) && v[["Y"]] > 0) v[["Y"]] else NA_real_)
  }
  unname(Y_NORM[cfg$norm])
}

#' Like ree_normalise(), but returns Y + La..Lu (Y is NA when the table has no Y column or the
#' reference has no Y value) and keeps the absolute scale (molar data: pmol/kg).
ree_normalise_y <- function(conc, basis = "mass", cfg = default_cfg()) {
  cc <- as.data.frame(conc, stringsAsFactors = FALSE)
  for (el in setdiff(SPEC_ELEMENTS, names(cc))) cc[[el]] <- NA_real_
  m <- as.matrix(cc[, SPEC_ELEMENTS, drop = FALSE]); storage.mode(m) <- "double"
  basis <- rep_len(tolower(trimws(as.character(basis))), nrow(m))
  if (nrow(m) == 0) return(m)
  nv <- c(Y = norm_y(cfg), norm_ppm(cfg)); mass <- c(Y = Y_MASS, ATOMIC_MASS[REE])
  # molar data are taken as pmol/kg: ppm (mg/kg) -> pmol/kg is x 1e9 / atomic mass (absolute scale
  # matters here, unlike for the ratios, where it cancels)
  ref <- rbind(mass = nv[SPEC_ELEMENTS], molar = nv[SPEC_ELEMENTS] * 1e9 / mass[SPEC_ELEMENTS])
  n <- m / ref[ifelse(basis == "molar", "molar", "mass"), , drop = FALSE]
  n[!is.finite(n) | n <= 0] <- NA_real_
  dimnames(n) <- dimnames(m)
  n
}

LREE_SETS <- list(LaPrNd = c("La", "Pr", "Nd"), PrNdSm = c("Pr", "Nd", "Sm"))
MREE_SET <- c("Gd", "Tb", "Dy")
HREE_SET <- c("Tm", "Yb", "Lu")

ratio_cols <- c("HREE_LREE", "MREE_MREEstar", "Ce_Ceanom", "Eu_Euanom")

#' All user-selectable calculation options in one list.
#' @field norm id in NORMALISERS or "custom" (then `norm_custom` = named ppm vector)
#' @field lree "LaPrNd" (Martin et al. 2010) or "PrNdSm"
#' @field gap "interpolate" (log-linear neighbours), "groupmean" (mean of the measured members
#'   of each REE group) or "none"
#' @field gd "measured" or "ignore" (replace Gd by baseline * Gd*)
#' @field gd_method "BD1996" (Bau and Dulski), "SmNd" (Sm^2/Nd), "TbDy" (mean of Tb and Dy),
#'   "BLambdaR" (polynomial fit without Gd, Ce, Eu)
#' @field gd_baseline natural Gd/Gd* that is divided out (1 = none)
default_cfg <- function(...) {
  modifyList(list(norm = DEFAULT_NORM, norm_custom = NULL, lree = "LaPrNd", gap = "interpolate",
                  gd = "measured", gd_method = "BD1996", gd_baseline = 1), list(...))
}

#' Normalising vector (ppm, mass basis) for a configuration
norm_ppm <- function(cfg = default_cfg()) {
  if (identical(cfg$norm, "custom")) {
    v <- cfg$norm_custom
    if (is.null(v) || any(!is.finite(v[REE])) || any(v[REE] <= 0)) {
      stop("Own reference: all 14 REE (La-Lu) with positive values are required.")
    }
    return(v[REE])
  }
  NORMALISERS[[cfg$norm]]
}

norm_label <- function(cfg) if (identical(cfg$norm, "custom")) "own reference" else
  unname(NORMALISER_LABELS[cfg$norm])

#' Normalise a table of REE concentrations.
#' @param conc data.frame/matrix with columns La..Lu
#' @param basis "mass" or "molar"; scalar or one value per row
#' @param norm named vector of normalising ppm values (mass basis); molar basis divides by the
#'   atomic masses, which is identical to converting the data to mass first.
ree_normalise <- function(conc, basis = "mass", norm = NORMALISERS[[DEFAULT_NORM]]) {
  m <- as.matrix(conc[, REE, drop = FALSE])
  storage.mode(m) <- "double"
  basis <- rep_len(tolower(trimws(as.character(basis))), nrow(m))
  if (any(!basis %in% c("mass", "molar"))) {
    stop("`basis` must be 'mass' or 'molar' (got: ",
         paste(unique(basis[!basis %in% c("mass", "molar")]), collapse = ", "), ").")
  }
  if (nrow(m) == 0) return(m)
  ref <- rbind(mass = norm[REE], molar = norm[REE] / ATOMIC_MASS[REE])
  n <- m / ref[basis, , drop = FALSE]
  n[!is.finite(n) | n <= 0] <- NA_real_
  dimnames(n) <- dimnames(m)
  n
}

#' Fill gaps in normalised REE patterns by log-linear interpolation of neighbours.
#' Pr, Tb, Ho, Tm are geometric means of their neighbours (Pr falls back on extrapolation from
#' Nd-Sm when La is missing, Tm on Ho-Yb when Er is missing, Er on Ho-Tm); Lu is extrapolated
#' from the Tm-Yb slope. Returns the filled matrix and a logical matrix flagging estimates.
ree_fill <- function(n) {
  rules <- list(
    Pr = function(x) { v <- sqrt(x[, "La"] * x[, "Nd"]); i <- is.na(v)
                       if (any(i)) v[i] <- x[i, "Nd"]^2 / x[i, "Sm"]; v },
    Tb = function(x) sqrt(x[, "Gd"] * x[, "Dy"]),
    Ho = function(x) sqrt(x[, "Dy"] * x[, "Er"]),
    Tm = function(x) { v <- sqrt(x[, "Er"] * x[, "Yb"]); i <- is.na(v)
                       if (any(i)) v[i] <- sqrt(x[i, "Ho"] * x[i, "Yb"]); v },
    Er = function(x) sqrt(x[, "Ho"] * x[, "Tm"]),
    Lu = function(x) x[, "Yb"]^2 / x[, "Tm"]
  )
  filled <- matrix(FALSE, nrow(n), ncol(n), dimnames = dimnames(n))
  for (pass in 1:2) {
    for (el in names(rules)) {
      miss <- is.na(n[, el])
      if (!any(miss)) next
      v <- rules[[el]](n)
      ok <- miss & is.finite(v)
      n[ok, el] <- v[ok]
      filled[ok, el] <- TRUE
    }
  }
  list(n = n, filled = filled)
}

#' Geogenic ("natural") Gd, Gd*, from PAAS-normalised neighbours
#' BD1996: Bau and Dulski (1996); SmNd: Gd* = Sm^2/Nd (Kulaksiz and Bau, 2007; Barrat et al. 2026)
gd_star <- function(n, method = "BD1996") {
  switch(method,
    BD1996 = unname(0.33 * n[, "Sm"] + 0.67 * n[, "Tb"]),
    SmNd   = unname(n[, "Sm"]^2 / n[, "Nd"]),
    TbDy   = unname((n[, "Tb"] + n[, "Dy"]) / 2),
    stop("Unknown Gd method: ", method))
}

#' REE roughness: scatter of the Nd-Lu pattern around locally smooth neighbours (normalised
#' values; measured elements only, so interpolated gaps cannot make a pattern look smooth).
#' NA unless Pr, Nd, Sm, Gd, Tb, Dy, Ho, Er, Tm, Yb and Lu are all available.
ree_roughness <- function(n) {
  d <- function(a, b, c) (n[, a] * 2 / (n[, b] + n[, c]) - 1)^2
  unname(sqrt(d("Nd", "Pr", "Sm") + d("Tb", "Gd", "Dy") + d("Dy", "Tb", "Ho") + d("Ho", "Dy", "Er") +
              d("Er", "Ho", "Tm") + d("Tm", "Er", "Yb") + d("Yb", "Tm", "Lu")))
}

#' Ratios from a matrix of normalised REE
#' @param filled optional logical matrix of estimated values (for `n_estimated`)
ratios_from_normalised <- function(n, filled = NULL, cfg = default_cfg()) {
  lre <- LREE_SETS[[cfg$lree]]
  grp <- if (identical(cfg$gap, "groupmean")) {
    function(els) { v <- rowMeans(n[, els, drop = FALSE], na.rm = TRUE); v[is.nan(v)] <- NA_real_; v }
  } else function(els) rowSums(n[, els, drop = FALSE])
  L <- grp(lre); H <- grp(HREE_SET); M <- grp(MREE_SET)
  out <- data.frame(
    HREE_LREE     = unname(H / L),
    MREE_MREEstar = unname(M / ((H + L) / 2)),
    Ce_Ceanom     = unname(2 * n[, "Ce"] / (n[, "La"] + n[, "Pr"])),
    Eu_Euanom     = unname(n[, "Eu"] / sqrt(n[, "Sm"] * n[, "Gd"])),
    row.names = NULL
  )
  out$n_estimated <- if (is.null(filled)) 0L else
    unname(rowSums(filled[, c("Pr", "Tb", "Tm", "Lu"), drop = FALSE]))
  out
}

#' Compute HREE/LREE, MREE/MREE*, Ce/Ce*, Eu/Eu*, Gd/Gd*, Gd_anth and REE roughness from raw
#' concentrations.
#' @param conc table with columns La..Lu (units: see `basis`)
#' @param basis "mass"/"molar", scalar or per row
#' @param cfg see `default_cfg()`
#' @param lambda_gd_anom optional Gd/Gd* from the polynomial (BLambdaR) fit; only used for
#'   gd_method = "BLambdaR"
#' @return data.frame; `Gd_anth` is in the unit of the input concentrations
ree_ratios <- function(conc, basis = "mass", cfg = default_cfg(), lambda_gd_anom = NULL) {
  n <- ree_normalise(conc, basis, norm_ppm(cfg))
  nr <- nrow(n)
  if (nr == 0) return(data.frame(HREE_LREE = numeric(), MREE_MREEstar = numeric(),
                                 Ce_Ceanom = numeric(), Eu_Euanom = numeric(), n_estimated = integer(),
                                 Gd_Gdstar = numeric(), Gd_anth = numeric(), REE_rough = numeric()))
  # Gd anomaly of the measured pattern
  gs <- if (identical(cfg$gd_method, "BLambdaR")) {
    if (is.null(lambda_gd_anom)) rep(NA_real_, nr) else unname(n[, "Gd"] / lambda_gd_anom)
  } else gd_star(n, cfg$gd_method)
  anom <- unname(n[, "Gd"] / gs)
  gd_conc <- as.numeric(as.matrix(conc[, "Gd", drop = FALSE]))
  base <- if (is.null(cfg$gd_baseline) || !is.finite(cfg$gd_baseline)) 1 else cfg$gd_baseline
  gd_anth <- ifelse(is.finite(anom), gd_conc * pmax(0, 1 - base / anom), NA_real_)
  if (identical(cfg$gd, "ignore")) n[, "Gd"] <- gs * base       # NA stays NA (missing)
  rough <- ree_roughness(n)
  filled <- matrix(FALSE, nr, ncol(n), dimnames = dimnames(n))
  if (identical(cfg$gap, "interpolate")) {
    f <- ree_fill(n); n <- f$n; filled <- f$filled
  }
  out <- ratios_from_normalised(n, filled, cfg)
  if (identical(cfg$gap, "groupmean")) {
    used <- unique(c(LREE_SETS[[cfg$lree]], HREE_SET, MREE_SET))
    out$n_estimated <- unname(rowSums(is.na(n[, used, drop = FALSE])))
  }
  out$Gd_Gdstar <- anom
  out$Gd_anth <- unname(gd_anth)
  out$REE_rough <- rough
  out
}

#' Which REE are needed to obtain both plot axes?
#' @return character vector describing what is missing (empty if everything needed is present)
check_ree_columns <- function(present, cfg = default_cfg()) {
  need <- list()   # each entry: alternatives (at least one must be present)
  lre <- LREE_SETS[[cfg$lree]]
  gd_inputs <- switch(cfg$gd_method, BD1996 = c("Sm", "Tb"), SmNd = c("Sm", "Nd"), TbDy = c("Tb", "Dy"),
                      BLambdaR = REE)
  ignore <- identical(cfg$gd, "ignore")
  if (identical(cfg$gap, "interpolate")) {
    need <- c(need, as.list(if (cfg$lree == "LaPrNd") c("La", "Nd") else c("Nd", "Sm")), "Yb")
    if (!ignore) need <- c(need, list("Gd")) else need <- c(need, as.list(gd_inputs))
    need <- c(need, list("Dy"), list(c("Er", "Ho", "Tm")))
  } else if (identical(cfg$gap, "groupmean")) {
    need <- list(lre, if (ignore) c("Tb", "Dy") else MREE_SET, HREE_SET)
    if (ignore) need <- c(need, as.list(gd_inputs))
  } else {
    need <- c(as.list(lre), as.list(HREE_SET), as.list(setdiff(MREE_SET, if (ignore) "Gd")))
    if (ignore) need <- c(need, as.list(gd_inputs))
  }
  miss <- vapply(unique(need), function(alt) if (any(alt %in% present)) NA_character_ else
    if (length(alt) == 1) alt else paste0("one of ", paste(alt, collapse = "/")), character(1))
  unique(miss[!is.na(miss)])
}

# Convert a molar end-member (pmol/kg) to the mass scale used for mixing
to_mass_scale <- function(conc, basis) {
  conc <- unlist(conc[REE])
  if (tolower(basis) == "molar") conc * ATOMIC_MASS[REE] * SEAWATER_TO_MASS_SCALE else conc
}
