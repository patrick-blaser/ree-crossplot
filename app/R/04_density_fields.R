# ---------------------------------------------------------------------------
# Highest-density-region (HDR) fields from 2-D kernel density estimates
# ---------------------------------------------------------------------------

# bilinear interpolation of a kde2d grid at points (used for large samples, where the exact
# pairwise kernel sum would be slow and memory hungry)
.interp_grid <- function(k, x, y) {
  gx <- k$x; gy <- k$y
  ix <- pmin(pmax(findInterval(x, gx), 1), length(gx) - 1)
  iy <- pmin(pmax(findInterval(y, gy), 1), length(gy) - 1)
  tx <- (x - gx[ix]) / (gx[ix + 1] - gx[ix]); ty <- (y - gy[iy]) / (gy[iy + 1] - gy[iy])
  k$z[cbind(ix, iy)] * (1 - tx) * (1 - ty) + k$z[cbind(ix + 1, iy)] * tx * (1 - ty) +
    k$z[cbind(ix, iy + 1)] * (1 - tx) * ty + k$z[cbind(ix + 1, iy + 1)] * tx * ty
}

#' Polygon(s) enclosing `coverage` of the points, from a Gaussian 2-D KDE on linear axes.
#' @param adjust multiplier on the default (normal-reference) bandwidth
#' @return data.frame(x, y, piece) or NULL if fewer than `min_n` finite points
hdr_polygons <- function(x, y, coverage = 0.9, adjust = 1, n = 200, min_n = 8, exact_max = 1200) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < min_n) return(NULL)
  h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y)) * adjust
  if (any(!is.finite(h)) || any(h <= 0)) return(NULL)
  hh <- h / 4  # kde2d divides h by 4 internally; keep point densities consistent with it
  lims <- c(range(x) + c(-1, 1) * 4 * hh[1], range(y) + c(-1, 1) * 4 * hh[2])
  big <- length(x) > exact_max
  k <- MASS::kde2d(x, y, h = h, n = if (big) 150 else n, lims = lims)
  dens <- if (big) .interp_grid(k, x, y) else {
    ax <- stats::dnorm(outer(x, x, "-") / hh[1])   # density at the data points (exact, same kernel)
    ay <- stats::dnorm(outer(y, y, "-") / hh[2])
    rowSums(ax * ay) / (length(x) * hh[1] * hh[2])
  }
  thr <- stats::quantile(dens, probs = 1 - coverage, names = FALSE)
  cl <- grDevices::contourLines(k$x, k$y, k$z, levels = thr)
  if (length(cl) == 0) return(NULL)
  do.call(rbind, lapply(seq_along(cl), function(i) {
    data.frame(x = cl[[i]]$x, y = cl[[i]]$y, piece = i)
  }))
}
