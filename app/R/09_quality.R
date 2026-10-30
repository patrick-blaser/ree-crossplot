# ---------------------------------------------------------------------------
# "Data quality" tab: distribution of the roughness values of the user's data (one curve per
# group) against the reference data
# ---------------------------------------------------------------------------

REF_CLASS_LEVELS <- c("Foraminifera (reference)", "Pore water (reference)", "Seawater, compilations (reference)",
                      "Seawater GEOTRACES ≤ 200 m (reference)", "Seawater GEOTRACES > 200 m (reference)")
REF_CLASS_COLOURS <- setNames(c("#7B2D8E", "#4C9A2A", "#2166AC", "#74A9CF", "#08306B"), REF_CLASS_LEVELS)

#' Class of each reference row for the roughness comparison
ref_class <- function(ref) {
  cl <- ifelse(ref$type == "forams", REF_CLASS_LEVELS[1],
        ifelse(ref$type == "pore water", REF_CLASS_LEVELS[2],
        ifelse(ref$ref == "GEOTRACES2025_shallow", REF_CLASS_LEVELS[4],
        ifelse(ref$ref == "GEOTRACES2025_deep", REF_CLASS_LEVELS[5], REF_CLASS_LEVELS[3]))))
  factor(cl, levels = REF_CLASS_LEVELS)
}

#' Summary table of roughness values per dataset
#' @param ref data.frame(class, v) of reference values (or NULL); @param user data.frame(group, v) (or NULL)
#' @param limit roughness limit used for the share above it
roughness_summary <- function(ref, user, limit = NA) {
  one <- function(name, v, kind) {
    v <- v[is.finite(v)]
    data.frame(Dataset = name, Type = kind, n = length(v),
               Median = if (length(v)) signif(median(v), 3) else NA_real_,
               `90th percentile` = if (length(v)) signif(unname(quantile(v, 0.9)), 3) else NA_real_,
               `Above limit (%)` = if (length(v) && is.finite(limit)) round(100 * mean(v > limit), 1) else NA_real_,
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  out <- list()
  if (!is.null(user) && nrow(user)) for (g in natural_sort(user$group)) out[[length(out) + 1]] <- one(g, user$v[user$group == g], "your data")
  if (!is.null(ref) && nrow(ref)) for (cl in levels(droplevels(ref$class))) out[[length(out) + 1]] <- one(sub(" \\(reference\\)$", "", cl), ref$v[ref$class == cl], "reference")
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

#' Density curves of the roughness
#' @param ref data.frame(class, v); @param user data.frame(group, v)
#' @param opts log (log10 axis), limit (numeric or NA), limit_on (is the filter active?), label (axis title),
#'   group_levels (fixes the colours), legend, base_size, adjust
plot_roughness <- function(ref, user, opts = list()) {
  o <- utils::modifyList(list(log = TRUE, limit = NA_real_, limit_on = FALSE, label = "REE roughness",
                              group_levels = NULL, legend = TRUE, base_size = 14, adjust = 1), opts)
  tr <- if (o$log) log10 else identity
  ok <- function(v) is.finite(v) & (!o$log | v > 0)
  # density curves are computed here (not by geom_density) so that they run out smoothly beyond the data range
  dens_df <- function(d, col) {
    out <- lapply(split(d, d[[col]]), function(g) {
      if (length(unique(g$x)) < 5) return(NULL)               # too few values for a curve (they still show as ticks)
      de <- stats::density(g$x, adjust = o$adjust, cut = 2.5, n = 400)
      data.frame(x = de$x, y = de$y, series = g[[col]][1])
    })
    out <- out[!vapply(out, is.null, logical(1))]
    if (!length(out)) return(NULL)
    r <- do.call(rbind, out); r$series <- factor(r$series, levels = levels(d[[col]])); r
  }
  cols <- character(); lt <- character()
  p <- ggplot2::ggplot()
  if (!is.null(ref) && nrow(ref)) {
    r <- ref[ok(ref$v), , drop = FALSE]; r$x <- tr(r$v); r$series <- droplevels(r$class)
    cols <- c(cols, REF_CLASS_COLOURS[levels(r$series)]); lt <- c(lt, setNames(rep("dashed", nlevels(r$series)), levels(r$series)))
    rd <- dens_df(r, "series")
    if (!is.null(rd)) p <- p + ggplot2::geom_line(data = rd, ggplot2::aes(x = x, y = y, colour = series), linetype = "dashed", linewidth = 0.8)
  }
  if (!is.null(user) && nrow(user)) {
    u <- user[ok(user$v), , drop = FALSE]; u$x <- tr(u$v)
    lv <- if (is.null(o$group_levels)) natural_sort(user$group) else o$group_levels
    u$series <- factor(u$group, levels = lv)
    ucol <- setNames(rep_len(USER_COLOURS, length(lv)), lv)
    present <- levels(droplevels(u$series))
    cols <- c(cols, ucol[present]); lt <- c(lt, setNames(rep("solid", length(present)), present))
    ud <- dens_df(u, "series")
    if (!is.null(ud)) p <- p +
      ggplot2::geom_ribbon(data = ud, ggplot2::aes(x = x, ymin = 0, ymax = y, fill = series), alpha = 0.3, colour = NA, show.legend = FALSE) +
      ggplot2::geom_line(data = ud, ggplot2::aes(x = x, y = y, colour = series), linetype = "solid", linewidth = 1)
    p <- p + ggplot2::geom_rug(data = u, ggplot2::aes(x = x, colour = series), show.legend = TRUE, sides = "b", alpha = 0.6, length = grid::unit(0.03, "npc"))
    p <- p + ggplot2::scale_fill_manual(values = ucol, guide = "none")
  }
  if (is.finite(o$limit) && (!o$log || o$limit > 0)) {
    p <- p + ggplot2::geom_vline(xintercept = tr(o$limit), linetype = "dotted", linewidth = 0.7) +
      ggplot2::annotate("text", x = tr(o$limit), y = Inf, vjust = 1.4, hjust = -0.1, size = 3.6,
                        label = if (o$limit_on) "limit" else "limit (filter off)")
  }
  brk <- c(1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
  p <- p + ggplot2::scale_colour_manual(values = cols, name = NULL, breaks = names(cols))
  if (o$log) {
    p <- p + ggplot2::scale_x_continuous(breaks = log10(brk), labels = format(brk, scientific = FALSE, trim = TRUE, drop0trailing = TRUE),
                                         expand = ggplot2::expansion(mult = 0.03))
  }
  p <- p + ggplot2::labs(x = paste0(o$label, if (o$log) " (log scale)" else ""), y = "Density") +
    ggplot2::theme_classic(base_size = o$base_size) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                   axis.line = ggplot2::element_blank(), legend.position = if (isTRUE(o$legend)) "right" else "none")
  p
}
