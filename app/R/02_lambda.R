# ---------------------------------------------------------------------------
# Lambda (orthogonal polynomial) description of REE patterns
#
# Base-R port of the BLambdaR calculation (Anenburg and Williams, 2022; BLambdaR,
# https://lambdar.rses.anu.edu.au), as used by P. Blaser. The pattern is normalised to chondrite
# (O'Neill, 2016), log-transformed and fitted by orthogonal polynomials of the ionic radius
# (lambda0..lambda4), optionally plus four tetrad terms (tau1..tau4). Elements with anomalies
# (Ce, Eu, optionally La, Gd) are left out of the fit; Y is not used.
# ---------------------------------------------------------------------------

LAMBDA_RADII <- c(La = 1.160, Ce = 1.143, Pr = 1.126, Nd = 1.109, Sm = 1.079, Eu = 1.066,
                  Gd = 1.053, Tb = 1.040, Dy = 1.027, Ho = 1.015, Er = 1.004, Tm = 0.994,
                  Yb = 0.985, Lu = 0.977)               # O'Neill (2016)
LAMBDA_NORM <- NORMALISERS$CI_ONeill2016


.LAMBDA <- local({
  radii <- LAMBDA_RADII
  N <- length(radii)
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
  z <- (57:71)[-5]                                   # atomic numbers without Pm
  fit <- lm(z ~ poly(radii, 7))
  list(f = list(f1, f2, f3, f4), fmat = cbind(f1(radii), f2(radii), f3(radii), f4(radii)),
       trn = function(r) predict(fit, newdata = data.frame(radii = r)))
})

.tetrad <- function(x) (1 - x^2) / 2 + sqrt(((1 - x^2)^2) / 4)
.tf <- function(x, t) .tetrad((x - (55.25 + 3.5 * t)) / 1.75)

#' Fit lambda/tau parameters to REE patterns.
#' @param conc table with columns La..Lu (any positive concentration unit; per-row `basis`)
#' @param anomalies elements excluded from the fit
#' @param fit_factors which of lambda4 and tau1..tau4 are fitted (lambda0..3 always are)
#' @return data.frame: R2, lambda0..lambda4, tau1..tau4, <El>_MFR (measured/fitted), lambda_rough
blambda <- function(conc, basis = "mass", anomalies = c(La = FALSE, Ce = TRUE, Eu = TRUE, Gd = FALSE),
                    fit_factors = c(l4 = TRUE, t1 = TRUE, t2 = TRUE, t3 = TRUE, t4 = TRUE),
                    norm = LAMBDA_NORM, min_points = 10) {
  # min_points: patterns with fewer usable elements (after leaving out the anomalous ones) give NA;
  # the default fit has 9 parameters, so 10 points leave one degree of freedom
  L <- .LAMBDA; radii <- LAMBDA_RADII
  an <- c(La = FALSE, Ce = TRUE, Eu = TRUE, Gd = FALSE); an[names(anomalies)] <- anomalies
  ff <- c(l4 = TRUE, t1 = FALSE, t2 = FALSE, t3 = FALSE, t4 = FALSE); ff[names(fit_factors)] <- fit_factors
  m <- as.matrix(conc[, REE, drop = FALSE]); storage.mode(m) <- "double"
  m[!is.finite(m) | m <= 0] <- NA
  basis <- rep_len(tolower(as.character(basis)), nrow(m))
  mol <- basis == "molar"
  if (any(mol)) m[mol, ] <- sweep(m[mol, , drop = FALSE], 2, ATOMIC_MASS[REE], "*")
  nm <- c("lambda0", "lambda1", "lambda2", "lambda3", "lambda4", "tau1", "tau2", "tau3", "tau4")
  npar <- 1 + length(nm) + length(REE)
  excl <- c(if (an[["Gd"]]) 7, if (an[["Eu"]]) 6, if (an[["Ce"]]) 2, if (an[["La"]]) 1)

  fit_one <- function(ppm) {
    k <- log(ppm / norm)
    if (!(sum(!is.na(k[3:4])) > 0 && sum(!is.na(k[7:9])) > 0 && sum(!is.na(k[12:14])) > 0))
      return(rep(NA_real_, npar))
    rad.fit <- radii
    if (length(excl)) { k <- k[-excl]; rad.fit <- rad.fit[-excl] }
    if (sum(is.finite(k)) < min_points) return(rep(NA_real_, npar))
    X <-cbind(L$f[[1]](rad.fit), L$f[[2]](rad.fit), L$f[[3]](rad.fit))
    if (ff[["l4"]]) X <- cbind(X, L$f[[4]](rad.fit))
    zz <- L$trn(rad.fit)
    for (t in 1:4) if (ff[[paste0("t", t)]]) X <- cbind(X, .tf(zz, t))
    ok <- is.finite(k)
    A <- cbind(1, X)[ok, , drop = FALSE]
    fit <- lm.fit(A, k[ok])
    cf <- fit$coefficients; cf[is.na(cf)] <- NA_real_
    tn <- paste0("t", 1:4)
    pos <- 1 + 3 + (if (ff[["l4"]]) 1 else 0)          # index of the last polynomial coefficient
    res <- c(cf[1:4], if (ff[["l4"]]) cf[5] else NA)
    j <- pos
    for (t in 1:4) if (ff[[tn[t]]]) { j <- j + 1; res <- c(res, cf[j]) } else res <- c(res, NA_real_)
    res <- unname(res)
    r0 <- ifelse(is.na(res), 0, res)
    zr <- L$trn(radii)
    pat <- exp(r0[1] + as.vector(L$fmat %*% r0[2:5]) +
                 r0[6] * .tf(zr, 1) + r0[7] * .tf(zr, 2) + r0[8] * .tf(zr, 3) + r0[9] * .tf(zr, 4))
    rss <- sum(fit$residuals^2); nn <- sum(ok); pp <- sum(!is.na(cf))
    tss <- sum((k[ok] - mean(k[ok]))^2)
    adj <- if (nn - pp > 0 && tss > 0) 1 - (rss / (nn - pp)) / (tss / (nn - 1)) else NA_real_
    c(adj, res, ppm / (pat * norm))
  }
  out <- t(apply(m, 1, function(r) tryCatch(fit_one(r), error = function(e) rep(NA_real_, npar))))
  if (nrow(m) == 0) out <- matrix(NA_real_, 0, npar)
  out <- as.data.frame(out)
  names(out) <- c("R2", nm, paste0(REE, "_MFR"))
  mf <- out[, paste0(c("Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"), "_MFR"), drop = FALSE]
  out$lambda_rough <- sqrt(rowSums((mf - 1)^2)) / 7
  out
}

#' Gd/Gd* from a polynomial fit that leaves out Gd, Ce and Eu (tetrads off for stability)
lambda_gd_anomaly <- function(conc, basis = "mass") {
  blambda(conc, basis, anomalies = c(Gd = TRUE), fit_factors = c(l4 = TRUE, t1 = FALSE, t2 = FALSE, t3 = FALSE, t4 = FALSE),
          min_points = 8)$Gd_MFR
}
