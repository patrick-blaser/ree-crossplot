# ---------------------------------------------------------------------------
# Authigenic and detrital arrays, built from end-member compositions
# (data/endmembers.csv) rather than shipped as pre-computed curves. Mixing is linear in
# concentration, which is also linear in normalised space, so every normalisation standard works.
# ---------------------------------------------------------------------------

load_endmembers <- function(path = "data/endmembers.csv") {
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("NA", ""))
}

AUTH_R <- c(0, 10^seq(-4, 0, by = 0.04))
DETR_N <- 101

# normalised, gap-filled pattern of one end-member (1 x 14 matrix); gaps are always interpolated
endmember_normalised <- function(em_row, cfg = default_cfg()) {
  conc <- to_mass_scale(em_row, em_row$basis)
  n <- ree_normalise(as.data.frame(t(conc)), "mass", norm_ppm(cfg))
  ree_fill(n)$n
}

# gap-filled concentration pattern (mass scale) of one end-member; independent of the standard
endmember_conc <- function(em_row) {
  nrm <- NORMALISERS$PAAS_McLennan1989
  n <- ree_fill(ree_normalise(as.data.frame(t(to_mass_scale(em_row, em_row$basis))), "mass", nrm))$n
  n * matrix(nrm[REE], nrow(n), length(REE), byrow = TRUE)
}

mix_normalised <- function(a, b, r) outer(1 - r, a[1, ]) + outer(r, b[1, ])

.ratio_cfg <- function(cfg) { cfg$gap <- "interpolate"; cfg }

#' Authigenic array: mixing between the most HREE-rich seawater (r = 0) and the most
#' MREE-rich diagenetic Fe-Mn nodule (r = 1), following Du et al. (2016). Mixing is linear in
#' concentration with log-spaced r (as in the published figure).
authigenic_array <- function(em, cfg = default_cfg(), r = AUTH_R) {
  sw  <- em[em$id == "seawater_HREE", ]
  nod <- em[em$id == "nodule_CCFZ", ]
  n <- mix_normalised(endmember_normalised(sw, cfg), endmember_normalised(nod, cfg), r)
  colnames(n) <- REE
  cbind(data.frame(r = r), ratios_from_normalised(n, NULL, .ratio_cfg(cfg)))
}

.DETR_IDS <- c(GAC = "crust_GAC", MORB = "MORB_Gale2013", CPS = "crust_CPS")
.DETR_PAIRS <- list(c("GAC", "MORB"), c("MORB", "CPS"), c("GAC", "CPS"))

#' Detrital array: mixing curves GAC-MORB, MORB-CPS and GAC-CPS. Returns the three curves, the
#' closed polygon bounding them and the three end-members.
detrital_array <- function(em, cfg = default_cfg(), n = DETR_N) {
  r <- seq(0, 1, length.out = n)
  nrm <- lapply(.DETR_IDS, function(id) endmember_normalised(em[em$id == id, ], cfg))
  rc <- .ratio_cfg(cfg)
  curve <- function(p) {
    m <- mix_normalised(nrm[[p[1]]], nrm[[p[2]]], r); colnames(m) <- REE
    cbind(data.frame(curve = paste(p, collapse = "-"), r = r), ratios_from_normalised(m, NULL, rc))
  }
  curves <- do.call(rbind, lapply(.DETR_PAIRS, curve))
  poly <- rbind(curves[curves$curve == "GAC-MORB", ],
                curves[curves$curve == "MORB-CPS", ],
                curves[curves$curve == "GAC-CPS", ][n:1, ])
  ends <- do.call(rbind, lapply(names(nrm), function(k) {
    cbind(data.frame(end = k), ratios_from_normalised(nrm[[k]], NULL, rc))
  }))
  list(curves = curves, polygon = poly, endpoints = ends)
}

#' Lambda-space coordinates (-lambda1, -lambda2) of both arrays (BLambdaR fit of the mixed
#' concentration patterns). Computed once by data-raw/build_data.R and stored in
#' data/arrays_lambda.csv; independent of the chosen shale standard.
arrays_lambda_table <- function(em) {
  ids <- c("seawater_HREE", "nodule_CCFZ", unname(.DETR_IDS))
  conc <- setNames(lapply(ids, function(id) endmember_conc(em[em$id == id, ])), ids)
  mixc <- function(a, b, r) { m <- outer(1 - r, conc[[a]][1, ]) + outer(r, conc[[b]][1, ]); colnames(m) <- REE; m }
  one <- function(label, m, r) {
    b <- blambda(as.data.frame(m), "mass")
    data.frame(layer = label, r = r, lx = -b$lambda1, ly = -b$lambda2)
  }
  r <- seq(0, 1, length.out = DETR_N)
  rbind(one("authigenic", mixc("seawater_HREE", "nodule_CCFZ", AUTH_R), AUTH_R),
        one("GAC-MORB", mixc("crust_GAC", "MORB_Gale2013", r), r),
        one("MORB-CPS", mixc("MORB_Gale2013", "crust_CPS", r), r),
        one("GAC-CPS", mixc("crust_GAC", "crust_CPS", r), r))
}

#' Add lambda-space columns (lx, ly) to the arrays and keep ratio-space copies as (x, y)
attach_lambda <- function(arrays, lam) {
  a <- arrays$authigenic; l <- lam[lam$layer == "authigenic", ]
  a$lx <- l$lx[match(round(a$r, 8), round(l$r, 8))]; a$ly <- l$ly[match(round(a$r, 8), round(l$r, 8))]
  arrays$authigenic <- a
  cv <- arrays$detrital$curves
  key <- paste(lam$layer, round(lam$r, 8)); k2 <- paste(cv$curve, round(cv$r, 8))
  cv$lx <- lam$lx[match(k2, key)]; cv$ly <- lam$ly[match(k2, key)]
  arrays$detrital$curves <- cv
  p <- arrays$detrital$polygon; k3 <- paste(p$curve, round(p$r, 8))
  p$lx <- lam$lx[match(k3, key)]; p$ly <- lam$ly[match(k3, key)]
  arrays$detrital$polygon <- p
  arrays
}
