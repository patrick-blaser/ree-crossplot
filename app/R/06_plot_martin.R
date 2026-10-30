# ---------------------------------------------------------------------------
# The MREE/MREE* vs HREE/LREE cross-plot (and its lambda-space variant)
# ---------------------------------------------------------------------------

FIELD_COLOURS <- c(forams = "#7B2D8E", `pore water` = "#4C9A2A", seawater = "#2166AC")
FIELD_NAMES   <- c(forams = "foraminifera", `pore water` = "pore water", seawater = "seawater")
# Okabe-Ito colour-blind-safe palette for user groups
USER_COLOURS <- c("#E69F00", "#D55E00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#999999")
USER_SHAPES  <- c(21, 23, 24, 22, 25)
FILL_COLOURS <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")   # viridis anchors

DEFAULT_XLIM <- c(0.3, 6.7)
DEFAULT_YLIM <- c(0.55, 2.08)

load_reference <- function(path = "data/reference_data.csv") {
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("NA", ""))
}

#' Add ratio-space coordinates (x, y) and the other derived quantities to the reference data
#' for a given configuration (normalisation, LREE definition, gap handling, Gd handling).
reference_ratios <- function(ref, cfg = default_cfg()) {
  r <- ree_ratios(ref[, REE], ref$basis, cfg)
  ref$x <- r$HREE_LREE
  ref$y <- r$MREE_MREEstar
  ref$n_estimated <- r$n_estimated
  ref$REE_rough <- r$REE_rough
  ref
}

.space_cols <- function(space) if (identical(space, "lambda")) c("lx", "ly") else c("x", "y")
.array_cols <- function(space) if (identical(space, "lambda")) c("lx", "ly") else c("HREE_LREE", "MREE_MREEstar")

# plot window -------------------------------------------------------------------------------
# Robust extent of user data: the 0.5-99.5 % range once there are at least 200 points, so a few
# extreme outliers do not squash the rest (they stay off-scale; the status line counts them).
.robust_range <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) >= 200) stats::quantile(v, c(0.005, 0.995), names = FALSE) else range(v)
}

#' Axis limits (in data units, before any log transform).
#' @param window "auto" (window holds the arrays, the field outlines and the bulk of the user
#'   data), "user" (zoom to the user data) or "arrays" (range of the two arrays only)
#' @param ext extent of the drawn field outlines, list(x = range, y = range), or NULL
#' @param log log10 axes (padding is applied in log space; ratio space only)
plot_limits <- function(ref, arrays, user = NULL, space = "ratio", window = "auto", ext = NULL, log = FALSE) {
  ac <- .array_cols(space); uc <- .space_cols(space)
  lam <- identical(space, "lambda"); lg <- isTRUE(log) && !lam
  tr <- if (lg) log10 else identity; inv <- if (lg) function(v) 10^v else identity
  arr <- rbind(setNames(arrays$authigenic[, ac], c("x", "y")),
               setNames(arrays$detrital$curves[, ac], c("x", "y")))
  arr <- arr[is.finite(arr$x) & is.finite(arr$y), ]
  pad <- function(r, f, top = f) r + c(-diff(r) * f, diff(r) * top)
  ax <- tr(range(arr$x)); ay <- tr(range(arr$y))
  if (identical(window, "arrays")) return(list(xlim = inv(pad(ax, 0.06)), ylim = inv(pad(ay, 0.06))))

  # base window: arrays + field outlines (or reference points when no outlines are drawn)
  bx <- ax; by <- ay
  if (!is.null(ext)) { bx <- range(bx, tr(ext$x)); by <- range(by, tr(ext$y)) }
  else if (!is.null(ref) && nrow(ref)) {
    rx <- ref[[uc[1]]]; ry <- ref[[uc[2]]]; ok <- is.finite(rx) & is.finite(ry)
    if (any(ok)) { bx <- range(bx, tr(.robust_range(rx[ok]))); by <- range(by, tr(.robust_range(ry[ok]))) }
  }
  bx <- pad(bx, 0.04); by <- pad(by, 0.04, top = 0.12)          # headroom for the legend
  if (!lam && !lg) {                                             # keep the layout of the published figure
    bx <- range(bx, DEFAULT_XLIM); by <- range(by, DEFAULT_YLIM)
  }
  xlim <- bx; ylim <- by
  if (!is.null(user)) {
    ux <- user[[uc[1]]]; uy <- user[[uc[2]]]
    ok <- is.finite(ux) & is.finite(uy) & (!lg | (ux > 0 & uy > 0))
    if (any(ok)) {
      ur <- .robust_range(tr(ux[ok])); vr <- .robust_range(tr(uy[ok]))
      if (identical(window, "user")) {
        px <- max(diff(ur) * 0.1, 0.02 * diff(bx)); py <- max(diff(vr) * 0.1, 0.02 * diff(by))
        xlim <- ur + c(-px, px); ylim <- vr + c(-py, py)
      } else {
        sx <- 0.01 * diff(bx); sy <- 0.01 * diff(by)
        xlim <- c(min(xlim[1], ur[1] - sx), max(xlim[2], ur[2] + sx))
        ylim <- c(min(ylim[1], vr[1] - sy), max(ylim[2], vr[2] + sy))
      }
    }
  }
  list(xlim = inv(xlim), ylim = inv(ylim))
}

.log_breaks <- function(lim, cand) cand[cand >= lim[1] & cand <= lim[2]]

#' Build the figure.
#' @param ref reference compilation with x, y (and lx, ly), already filtered
#' @param arrays list(authigenic, detrital) from authigenic_array()/detrital_array()
#' @param user NULL or data.frame(sample, group, x, y[, lx, ly, .fill])
#' @param opts list: space ("ratio"/"lambda"), log, window, fields, ref_points, authigenic,
#'   detrital, labels, coverage, adjust, point_size, alpha, user_fields, base_size,
#'   group_levels (fixes colours/shapes when groups are hidden), fill_label (when `.fill` is used)
plot_martin <- function(ref, arrays, user = NULL, opts = list()) {
  o <- modifyList(list(space = "ratio", log = FALSE, window = "auto", fields = TRUE, ref_points = TRUE,
                       authigenic = TRUE, detrital = TRUE, labels = TRUE, coverage = 0.95, adjust = 1.0,
                       point_size = 3.2, alpha = 0.9, user_fields = FALSE, zoom_user = FALSE,
                       base_size = 14, group_levels = NULL, fill_label = NULL, legend = TRUE), opts)
  if (isTRUE(o$zoom_user) && identical(o$window, "auto")) o$window <- "user"   # old option name
  lam <- identical(o$space, "lambda")
  lg <- isTRUE(o$log) && !lam
  tr <- if (lg) log10 else identity
  uc <- .space_cols(o$space); ac <- .array_cols(o$space)
  ext <- list(x = numeric(), y = numeric())     # extent of the drawn field outlines
  p <- ggplot2::ggplot()

  # --- reference fields and points -----------------------------------------
  ref_x <- ref[[uc[1]]]; ref_y <- ref[[uc[2]]]
  field_centre <- list()
  for (tp in names(FIELD_COLOURS)) {
    sel <- ref$type == tp & is.finite(ref_x) & is.finite(ref_y)
    if (!any(sel)) next
    s <- data.frame(x = ref_x[sel], y = ref_y[sel])
    col <- FIELD_COLOURS[[tp]]
    field_centre[[tp]] <- c(stats::median(s$x), stats::median(s$y))
    if (o$fields) {
      poly <- hdr_polygons(s$x, s$y, o$coverage, o$adjust)
      if (!is.null(poly)) {
        ext$x <- range(ext$x, poly$x); ext$y <- range(ext$y, poly$y)
        poly$x <- tr(poly$x); poly$y <- tr(poly$y)
        p <- p +
          ggplot2::geom_polygon(data = poly, ggplot2::aes(x, y, group = piece), fill = col,
                                alpha = 0.22, colour = NA) +
          ggplot2::geom_path(data = poly, ggplot2::aes(x, y, group = piece), colour = col,
                             linewidth = 0.6)
      }
    }
    if (o$ref_points) {
      shp <- c(forams = 2, `pore water` = 0, seawater = 16)[[tp]]
      s$x <- tr(s$x); s$y <- tr(s$y)
      p <- p + ggplot2::geom_point(data = s, ggplot2::aes(x, y), colour = col, shape = shp,
                                   size = 1.1, alpha = if (nrow(s) > 1000) 0.35 else 0.55, stroke = 0.4)
    }
  }

  # --- arrays ---------------------------------------------------------------
  have_det <- !is.null(arrays$detrital) && all(ac %in% names(arrays$detrital$curves))
  have_aut <- !is.null(arrays$authigenic) && all(ac %in% names(arrays$authigenic))
  if (o$detrital && have_det) {
    dp <- data.frame(x = tr(arrays$detrital$polygon[[ac[1]]]), y = tr(arrays$detrital$polygon[[ac[2]]]))
    cv <- data.frame(x = tr(arrays$detrital$curves[[ac[1]]]), y = tr(arrays$detrital$curves[[ac[2]]]),
                     curve = arrays$detrital$curves$curve)
    p <- p +
      ggplot2::geom_polygon(data = dp, ggplot2::aes(x, y), fill = "grey55", alpha = 0.30, colour = NA) +
      ggplot2::geom_path(data = cv, ggplot2::aes(x, y, group = curve), colour = "black", linewidth = 0.6)
    if (o$labels) {
      first <- function(nm, which) { d <- cv[cv$curve == nm, ]; d[if (which == "first") 1 else nrow(d), ] }
      lab <- rbind(cbind(first("GAC-MORB", "first")[, 1:2], label = "GAC"),
                   cbind(first("GAC-MORB", "last")[, 1:2], label = "MORB"),
                   cbind(first("MORB-CPS", "last")[, 1:2], label = "CPS"))
      if (!lam) lab$x[lab$label == "MORB"] <- lab$x[lab$label == "MORB"] + if (lg) 0.045 else 0.25
      p <- p +
        ggplot2::geom_label(data = lab, ggplot2::aes(x, y, label = label), fill = "grey80",
                            size = o$base_size * 0.30, label.padding = grid::unit(0.15, "lines"),
                            linewidth = 0.3)
      if (!lam) {
        p <- p + ggplot2::annotate("text", x = tr(3.3), y = tr(1.555), label = "detrital array",
                                   angle = if (lg) 6 else 9, size = o$base_size * 0.40, hjust = 0.5)
      } else {
        p <- p + ggplot2::annotate("text", x = mean(range(dp$x)), y = max(dp$y) + 0.03 * diff(range(dp$y)),
                                   label = "detrital array", size = o$base_size * 0.40, vjust = 0)
      }
    }
  }
  if (o$labels && o$fields) {
    if (lam) {
      fl <- do.call(rbind, lapply(names(field_centre), function(tp)
        data.frame(x = field_centre[[tp]][1], y = field_centre[[tp]][2], label = FIELD_NAMES[[tp]], tp = tp)))
      if (!is.null(fl)) p <- p + ggplot2::annotate("text", x = fl$x, y = fl$y, label = fl$label,
                                                   colour = unname(FIELD_COLOURS[fl$tp]),
                                                   size = o$base_size * 0.38, fontface = "bold")
    } else {
      p <- p + ggplot2::annotate("text", x = tr(c(3.6, 2.9, 5.7)), y = tr(c(1.30, 1.74, 0.86)),
                                 label = c("foraminifera", "pore\nwater", "seawater"),
                                 colour = unname(FIELD_COLOURS[c("forams", "pore water", "seawater")]),
                                 size = o$base_size * 0.38, lineheight = 0.85)
    }
  }

  # --- user data ------------------------------------------------------------
  if (!is.null(user) && nrow(user) > 0) {
    u <- user[is.finite(user[[uc[1]]]) & is.finite(user[[uc[2]]]), , drop = FALSE]
    u$X <- tr(u[[uc[1]]]); u$Y <- tr(u[[uc[2]]])
    lv <- if (is.null(o$group_levels)) unique(user$group) else o$group_levels
    u$group <- factor(u$group, levels = lv)
    ng <- length(lv)
    cols <- setNames(rep_len(USER_COLOURS, ng), lv)
    shps <- setNames(rep_len(USER_SHAPES, ng), lv)
    if (o$user_fields) {
      for (g in levels(droplevels(u$group))) {
        sg <- u[u$group == g, ]
        poly <- hdr_polygons(sg[[uc[1]]], sg[[uc[2]]], o$coverage, o$adjust)
        if (!is.null(poly)) {
          poly$x <- tr(poly$x); poly$y <- tr(poly$y)
          p <- p +
            ggplot2::geom_polygon(data = poly, ggplot2::aes(x, y, group = piece), fill = cols[[g]],
                                  alpha = 0.25, colour = NA) +
            ggplot2::geom_path(data = poly, ggplot2::aes(x, y, group = piece), colour = "grey20",
                               linewidth = 0.6)
        }
      }
    }
    if (".fill" %in% names(u) && any(is.finite(u$.fill))) {
      p <- p +
        ggplot2::geom_point(data = u, ggplot2::aes(X, Y, fill = .fill, shape = group),
                            colour = "black", size = o$point_size, alpha = o$alpha, stroke = 0.7) +
        ggplot2::scale_fill_gradientn(colours = FILL_COLOURS, na.value = "grey80",
                                      name = if (is.null(o$fill_label)) "value" else o$fill_label) +
        ggplot2::scale_shape_manual(values = shps, name = NULL, drop = TRUE) +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(fill = "grey75")),
                        fill = ggplot2::guide_colourbar(order = 2))
    } else {
      p <- p +
        ggplot2::geom_point(data = u, ggplot2::aes(X, Y, fill = group, shape = group),
                            colour = "black", size = o$point_size, alpha = o$alpha, stroke = 0.7) +
        ggplot2::scale_fill_manual(values = cols, name = NULL, drop = TRUE) +
        ggplot2::scale_shape_manual(values = shps, name = NULL, drop = TRUE)
    }
  }

  # the authigenic array is drawn last, on top of everything else
  if (o$authigenic && have_aut) {
    au <- data.frame(x = tr(arrays$authigenic[[ac[1]]]), y = tr(arrays$authigenic[[ac[2]]]))
    p <- p + ggplot2::geom_path(data = au, ggplot2::aes(x, y), colour = "black", linewidth = 1.3,
                                linetype = "dashed")
  }
  lims <- plot_limits(ref, arrays, user, o$space, o$window,
                      ext = if (o$fields && length(ext$x)) ext else NULL, log = lg)
  xl <- tr(lims$xlim); yl <- tr(lims$ylim)
  p <- p +
    ggplot2::coord_cartesian(xlim = xl, ylim = yl, expand = FALSE) +
    ggplot2::labs(x = if (lam) expression(-lambda[1]) else "HREE/LREE",
                  y = if (lam) expression(-lambda[2]) else "MREE/MREE*") +
    ggplot2::theme_classic(base_size = o$base_size) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                   axis.line = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(colour = "black"),
                   legend.position = "inside", legend.position.inside = c(0.995, 0.995),
                   legend.justification = c(1, 1),
                   legend.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.4),
                   legend.key.size = grid::unit(1.1, "lines"))
  if (!isTRUE(o$legend)) p <- p + ggplot2::theme(legend.position = "none")
  if (lg) {
    bx <- .log_breaks(lims$xlim, c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20))
    by <- .log_breaks(lims$ylim, c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4))
    p <- p + ggplot2::scale_x_continuous(breaks = log10(bx), labels = bx) +
      ggplot2::scale_y_continuous(breaks = log10(by), labels = by)
  } else {
    p <- p + ggplot2::scale_x_continuous(breaks = function(l) pretty(l, 7))
  }
  attr(p, "limits") <- lims            # (data units) lets the app count samples that lie off-scale
  p
}
