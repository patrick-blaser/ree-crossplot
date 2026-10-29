# ---------------------------------------------------------------------------
# REE patterns (spectra) of the user data: normalised concentrations against the element or the
# ionic radius. Idea and options after P. Blaser's REE_spectrum() function.
# ---------------------------------------------------------------------------

# position on the element axis: Y at 0, La..Lu in atomic-number order; Pm (5) stays empty
SPEC_POS <- c(Y = 0, La = 1, Ce = 2, Pr = 3, Nd = 4, Sm = 6, Eu = 7, Gd = 8, Tb = 9, Dy = 10,
              Ho = 11, Er = 12, Tm = 13, Yb = 14, Lu = 15)
SPEC_BREAKS <- c(SPEC_POS[c("Y", "La", "Ce", "Pr", "Nd")], Pm = 5, SPEC_POS[c("Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu")])

#' Effective ionic radii (pm) of the 3+ ions, 6-coordinated (Shannon 1976), from the table used by P. Blaser
load_radii <- function(path = "data/ionic_radii_3plus.csv") {
  r <- read.csv(path, stringsAsFactors = FALSE)
  setNames(r$radius_pm, r$element)
}

#' Normalised patterns of the user data.
#' @param conc table with Y (optional) and La..Lu
#' @param norm2 NULL or the element every pattern is divided by after the normalisation (secondary
#'   normalisation: removes absolute concentration differences)
#' @return list(n = normalised matrix, n2 = secondary-normalised matrix or NULL)
ree_spectra <- function(conc, basis = "mass", cfg = default_cfg(), norm2 = NULL) {
  n <- ree_normalise_y(conc, basis, cfg)
  n2 <- NULL
  if (!is.null(norm2) && norm2 %in% colnames(n)) {
    n2 <- n / n[, norm2]
    n2[!is.finite(n2) | n2 <= 0] <- NA_real_
  }
  list(n = n, n2 = n2)
}

#' Long table (one row per sample and element) for plotting
#' @param m matrix with the elements as columns
#' @param meta data.frame with one row per row of `m` (row, sample, group, colour variable ...)
spectra_long <- function(m, meta, radii) {
  nr <- nrow(m); els <- colnames(m)
  out <- cbind(meta[rep(seq_len(nr), times = length(els)), , drop = FALSE],
               element = rep(els, each = nr), value = as.vector(m))
  out$pos <- unname(SPEC_POS[out$element])
  out$radius <- unname(radii[out$element])
  rownames(out) <- NULL
  out
}

#' @param long output of spectra_long() (`value`, `element`, `pos`, `radius`, `sid`, `group`, optional `.col`)
#' @param o list: x ("element"/"radius"), log, show_y, connect, points, legend, line_size, alpha,
#'   ylab, colour_label (set when `.col` is numeric), group_levels, base_size
plot_spectra <- function(long, o = list()) {
  o <- modifyList(list(x = "element", log = TRUE, show_y = TRUE, connect = TRUE, points = TRUE, legend = TRUE,
                       line_size = 0.5, alpha = 0.8, ylab = "sample / reference", colour_label = NULL,
                       group_levels = NULL, base_size = 14), o)
  d <- long
  if (!o$show_y) d <- d[d$element != "Y", ]
  d$.c <- NA
  d$xv <- if (identical(o$x, "radius")) d$radius else d$pos
  d <- d[is.finite(d$xv), ]
  d$ok <- is.finite(d$value) & d$value > 0
  numeric_col <- ".col" %in% names(d) && is.numeric(d$.col)
  if (!numeric_col) {
    lv <- if (is.null(o$group_levels)) natural_sort(d$group) else o$group_levels
    d$group <- factor(d$group, levels = lv)
  }
  d$.c <- if (numeric_col) d$.col else d$group
  # lines: Y is part of the line on the element axis, a separate point on the radius axis
  ld <- d[d$element != "Y" | identical(o$x, "element"), ]
  # connect = FALSE keeps the missing values so that the line breaks at a missing element
  ld <- if (o$connect) ld[ld$ok, ] else { ld$value[!ld$ok] <- NA_real_; ld }
  ld <- ld[order(ld$sid, ld$xv), ]
  pd <- d[d$ok & (o$points | d$element == "Y"), ]
  p <- ggplot2::ggplot(d[d$ok, ], ggplot2::aes(xv, value))
  if (nrow(ld)) p <- p + ggplot2::geom_path(data = ld, ggplot2::aes(group = sid, colour = .c),
                                           linewidth = o$line_size, alpha = o$alpha, na.rm = TRUE)
  if (nrow(pd)) p <- p + ggplot2::geom_point(data = pd, ggplot2::aes(colour = .c),
                                            size = o$line_size * 3 + 0.6, alpha = o$alpha)
  if (numeric_col) {
    p <- p + ggplot2::scale_colour_gradientn(colours = FILL_COLOURS, na.value = "grey70",
                                             name = if (is.null(o$colour_label)) "value" else o$colour_label)
  } else {
    ng <- nlevels(d$group)
    cols <- if (ng <= length(USER_COLOURS)) USER_COLOURS[seq_len(ng)] else grDevices::hcl.colors(ng, "Dark 3")
    p <- p + ggplot2::scale_colour_manual(values = setNames(cols, levels(d$group)), name = NULL, drop = FALSE)
  }
  p <- p + ggplot2::labs(y = o$ylab, x = if (identical(o$x, "radius")) "ionic radius (pm)" else NULL) +
    ggplot2::theme_classic(base_size = o$base_size) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                   axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black"))
  if (identical(o$x, "radius")) p <- p + ggplot2::scale_x_reverse()
  else {
    br <- SPEC_BREAKS[if (o$show_y) names(SPEC_BREAKS) else setdiff(names(SPEC_BREAKS), "Y")]
    p <- p + ggplot2::scale_x_continuous(breaks = unname(br), labels = names(br))
  }
  if (isTRUE(o$log)) p <- p + ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides = "l")
  if (!isTRUE(o$legend)) p <- p + ggplot2::theme(legend.position = "none")
  p
}
