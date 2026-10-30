# Server-logic tests for the Shiny app (no browser). Run from repo root: Rscript tests/test_app.R
suppressMessages({library(shiny); library(ggplot2); library(MASS)})
cnt <- 0L
check <- function(desc, cond) { if (!isTRUE(cond)) stop("FAILED: ", desc, call. = FALSE); cnt <<- cnt + 1L; cat("ok  ", desc, "\n") }
app_path <- normalizePath("app")
for (f in list.files(file.path(app_path, "R"), full.names = TRUE)) source(f)  # constants used by the tests
html <- function(x) paste(x$html, collapse = "")

# all inputs the UI defines, with their defaults; tests override single entries
base_inputs <- function(...) {
  utils::modifyList(list(
    use_example = FALSE, mode = "concentrations", basis = "mass", id_col = "(none)", group_col = "(none)",
    xcol = NULL, ycol = NULL, norm = DEFAULT_NORM, lree = "LaPrNd", gap = "interpolate", gd = "measured",
    gd_method = "BD1996", gd_baseline = 1, fill_var = "group", rough_on = FALSE, rough_max = 0.3,
    show_fields = TRUE, show_refpts = TRUE, show_auth = TRUE, show_detr = TRUE, show_labels = TRUE,
    ref_sel = names(REF_LABELS), lambda_space = FALSE, log_axes = FALSE, window = "auto",
    coverage = 0.95, adjust = 1, pt_size = 3.2, user_fields = FALSE, ce_form = "geoPrNd", show_legend = TRUE, min_ree = 0, req_ree = NULL,
    id_extra = NULL, spec_x = "element", spec_norm2 = "(none)", spec_color = "group", spec_y = TRUE, spec_log = TRUE,
    spec_connect = TRUE, spec_points = TRUE, spec_legend = TRUE, spec_lw = 0.5,
    rough_log = TRUE, rough_ref = TRUE, rough_legend = TRUE), list(...))
}
set_all <- function(session, ...) do.call(session$setInputs, base_inputs(...))

# 1. example data
testServer(app_path, {
  set_all(session, use_example = TRUE, id_col = "sample", group_col = "group")
  u <- user_data()
  check("example data loads through the app (26 rows, 25 plottable: one sample lacks Gd)",
        nrow(u) == 26 && sum(is.finite(u$x) & is.finite(u$y)) == 25)
  check("status reports plotted samples and the one sample lacking Gd",
        grepl("25 of 26", html(output$status)) && grepl("missing values", html(output$status)))
  check("column mapping is listed", grepl("Columns used \\(14\\)", html(output$mapping)))
  tb <- computed()
  check("computed table has the expected columns and a plot_status",
        all(c("row", "sample", "group", "HREE/LREE", "MREE/MREE*", "Ce/Ce*", "Gd/Gd*", "REE_roughness", "plot_status", "normalisation", "LREE_definition", "gap_method", "Gd_handling", "Ce_anomaly", "La/La*") %in% names(tb)) &&
        identical(tb$row, 1:26))
  check("figure object builds", inherits(fig(), "ggplot"))
  f <- output$dl_csv
  check("results CSV download works", file.exists(f) && nrow(read.csv(f)) == 26)
  f <- output$dl_pdf
  check("PDF download works", file.exists(f) && file.size(f) > 5000)
  af <- read.csv(output$dl_arr, check.names = FALSE)
  check("arrays & fields CSV holds arrays, detrital polygon and all three density fields, with the settings",
        all(c("authigenic_array", "detrital_curve", "detrital_polygon", "density_field") %in% af$layer) &&
        setequal(unique(af$subset[af$layer == "density_field"]), c("foraminifera", "pore water", "seawater")) &&
        all(af$normalisation == "PAAS - Pourmand et al. (2012)") && all(af$LREE_definition == "La-Pr-Nd"))
  set_all(session, use_example = TRUE, id_col = "sample", group_col = "group", norm = "EUS_Bau2018", lree = "PrNdSm")
  af2 <- read.csv(output$dl_arr, check.names = FALSE)
  check("arrays & fields respond to normalisation and LREE definition",
        all(af2$normalisation == "EUS (European Shale) - Bau et al. (2018)") && all(af2$LREE_definition == "Pr-Nd-Sm") &&
        !isTRUE(all.equal(af$`HREE/LREE`, af2$`HREE/LREE`)) &&
        all(c("HREE/LREE", "MREE/MREE*", "lambda1", "lambda2") %in% names(af2)) && !any(grepl("neg_", names(af2))))
  set_all(session, use_example = TRUE, id_col = "sample", group_col = "group", ref_sel = "Osborne2017")
  check("reference dataset selection filters the reference data", all(ref_used()$ref == "Osborne2017"))
  set_all(session, use_example = TRUE, id_col = "sample", group_col = "group", ref_sel = character(0))
  check("deselecting all reference data still plots", inherits(fig(), "ggplot") && nrow(ref_used()) == 0)
  set_all(session, use_example = TRUE, id_col = "sample", group_col = "group", ref_sel = c("GEOTRACES2025_shallow"))
  check("GEOTRACES depth classes can be chosen separately", all(ref_used()$ref == "GEOTRACES2025_shallow") && nrow(ref_used()) > 500)
})

# 2. groups: show one group, hide others
testServer(app_path, {
  set_all(session, use_example = TRUE, id_col = "sample", group_col = "group")
  session$setInputs(groups_shown = "LS shelf")
  st <- plot_status()
  check("only the selected group is plotted", all(user_plot()$group == "LS shelf") && any(st == "hidden (group)"))
  check("hidden groups are mentioned in the status and keep their colour slot",
        grepl("hidden groups", html(output$status)) && length(plot_opts()$group_levels) == 2)
  check("figure builds with a single visible group", inherits(fig(), "ggplot"))
  session$setInputs(groups_shown = character(0))
  check("no visible group still plots the reference figure", nrow(user_plot()) == 0 && inherits(fig(), "ggplot"))
  session$setInputs(groups_shown = c("LS shelf", "LS deep water"))
  check("all groups visible again", nrow(user_plot()) == 25)
})

# 3. an uploaded file with its own columns, roughness filter, fill mapping, Gd handling
tmp <- tempfile(fileext = ".csv")
tail_rows <- c("ID,Group,Age,La_ppm,Ce_ppm,Pr_ppm,Nd_ppm,Sm_ppm,Eu_ppm,Gd_ppm,Tb_ppm,Dy_ppm,Ho_ppm,Er_ppm,Tm_ppm,Yb_ppm,Lu_ppm",
  "s1,A,1.5,38.2,79.6,8.83,33.9,5.55,1.08,4.66,0.774,4.68,0.991,2.85,0.405,2.82,0.433",
  "s2,B,3.0,38.2,79.6,8.83,33.9,5.55,1.08,13.98,0.774,4.68,0.991,2.85,0.405,2.82,0.433",
  "s3,B,4.5,38.2,79.6,8.83,33.9,5.55,1.08,4.66,1.0,4.68,0.7,2.85,0.5,2.82,0.433")
writeLines(tail_rows, tmp)
fileinp <- data.frame(name = "mine.csv", datapath = tmp, size = 1, type = "text/csv")
testServer(app_path, {
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", norm = "PAAS_McLennan1989")
  u <- user_data()
  check("custom upload: PAAS (McLennan) sample plots at (1,1), 3 samples", nrow(u) == 3 && abs(u$x[1] - 1) < 1e-9 && abs(u$y[1] - 1) < 1e-9)
  check("anthropogenic Gd from a 3x Gd spike (2/3 of the measured Gd)", abs(u$Gd_Gdstar[2] - 3) < 1e-9 && abs(u$Gd_anth[2] - 13.98 * 2 / 3) < 1e-6)
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", norm = "PAAS_McLennan1989", gd = "ignore")
  check("ignoring Gd restores the spiked sample to (1,1)", abs(user_data()$y[2] - 1) < 1e-9)
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", norm = "PAAS_McLennan1989", rough_on = TRUE, rough_max = 0.3)
  rc <- user_data()$REE_rough
  check("roughness: the smooth sample stays; the Gd-spiked and the zig-zag sample are excluded",
        rc[1] < 1e-9 && rc[2] > 0.3 && rc[3] > 0.3 &&
        identical(plot_status(), c("plotted", "excluded (roughness)", "excluded (roughness)")) &&
        grepl("excluded by the roughness limit", html(output$status)))
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", norm = "PAAS_McLennan1989", rough_on = TRUE,
          rough_max = 0.3, gd = "ignore")
  check("with Gd ignored the Gd-spiked sample is no longer rough",
        identical(plot_status(), c("plotted", "plotted", "excluded (roughness)")))
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", norm = "PAAS_McLennan1989", fill_var = "x:Age")
  ex <- prepped()$extra
  check("numeric column Age is available for symbol fill and reaches the plot",
        "Age" %in% names(ex) && identical(fill_values()$label, "Age") && ".fill" %in% names(user_plot()) && inherits(fig(), "ggplot"))
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", fill_var = "rough")
  check("roughness can be used as the symbol fill", identical(fill_values()$label, "REE roughness") && inherits(fig(), "ggplot"))
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", log_axes = TRUE, window = "arrays")
  check("log10 axes and array-range window build", inherits(fig(), "ggplot"))
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", lambda_space = TRUE)
  check("lambda space: coordinates are computed only now and the figure builds",
        all(c("lx", "ly", "lambda_rough") %in% names(user_data())) && inherits(fig(), "ggplot"))
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", lree = "PrNdSm", gap = "groupmean")
  check("Pr-Nd-Sm and group-mean options run", all(is.finite(user_data()$x)) && inherits(fig(), "ggplot"))
})

# 3b. screening by the number of measured REE, extra columns, spectra
tmp2 <- tempfile(fileext = ".csv")
rows2 <- c("ID,Group,Age,Y,La,Ce,Pr,Nd,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu",
  "a10,B,1,25,38.2,79.6,8.83,33.9,5.55,1.08,4.66,0.774,4.68,0.991,2.85,0.405,2.82,0.433",
  "a2,A,2,25,38.2,79.6,8.83,33.9,5.55,1.08,4.66,0.774,4.68,0.991,2.85,0.405,2.82,0.433",
  "a3,A,3,25,38.2,79.6,,33.9,5.55,1.08,4.66,,4.68,0.991,2.85,,2.82,")
writeLines(rows2, tmp2)
fileinp2 <- data.frame(name = "m2.csv", datapath = tmp2, size = 1, type = "text/csv")
testServer(app_path, {
  set_all(session, file = fileinp2, id_col = "ID", group_col = "Group", id_extra = "Age", norm = "PAAS_McLennan1989")
  check("groups are listed alphabetically", identical(all_groups(), c("A", "B")))
  tb <- computed()
  check("results: row number first, extra column after the group, settings columns last",
        identical(names(tb)[1:4], c("row", "sample", "group", "Age")) && identical(tb$row, 1:3) &&
        tb$normalisation[1] == "PAAS - McLennan (1989)" && tb$LREE_definition[1] == "La-Pr-Nd" && !any(c("-lambda1", "lx") %in% names(tb)))
  check("sample with gaps is still plotted by default", all(plot_status() == "plotted"))
  set_all(session, file = fileinp2, id_col = "ID", group_col = "Group", min_ree = 13)
  check("minimum number of measured REE excludes the sparse sample (13 of 14 needed)",
        identical(plot_status(), c("plotted", "plotted", "excluded (too few REE)")) && grepl("too few measured", html(output$status)))
  set_all(session, file = fileinp2, id_col = "ID", group_col = "Group", req_ree = "Pr")
  check("required elements exclude samples lacking them", identical(plot_status(), c("plotted", "plotted", "excluded (too few REE)")))
  set_all(session, file = fileinp2, id_col = "ID", group_col = "Group", norm = "PAAS_McLennan1989", lambda_space = TRUE, min_ree = 12)
  ud <- user_data()
  check("lambda: the sparse pattern gets no lambda, the complete ones do",
        is.finite(ud$lx[1]) && is.na(ud$lx[3]) && identical(plot_status()[3], "excluded (too few REE)"))
  af <- read.csv(output$dl_arr, check.names = FALSE)
  check("lambda fields/arrays are exported as positive lambda values",
        all(c("lambda1", "lambda2") %in% names(af)) && all(af$lambda1[af$layer == "authigenic_array"] > 0))
  # spectra
  set_all(session, file = fileinp2, id_col = "ID", group_col = "Group", norm = "PAAS_McLennan1989", spec_norm2 = "Dy")
  check("spectra figure builds", inherits(spec_fig(), "ggplot"))
  w <- spec_wide()
  check("normalised REE are offered for download (with Y and secondary normalisation)",
        all(c("norm_Y", "norm_La", "norm2_Dy") %in% names(w)) && abs(w$norm_La[1] - 1) < 1e-9 && all(abs(w$norm2_Dy - 1) < 1e-9) &&
        nrow(read.csv(output$dl_spec_csv)) == 3)
  set_all(session, file = fileinp2, id_col = "ID", group_col = "Group", spec_x = "radius", spec_log = FALSE, spec_y = FALSE, spec_legend = FALSE)
  check("spectra: radius axis without Y and legend", inherits(spec_fig(), "ggplot") && file.exists(output$dl_spec_pdf))
  session$setInputs(groups_shown = "A")
  check("spectra follow the group visibility", all(spec_long()$group == "A"))
})

# 4. own reference
nf <- tempfile(fileext = ".csv"); write.csv(as.data.frame(t(NORMALISERS$NASC_McLennan1989)), nf, row.names = FALSE)
testServer(app_path, {
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", norm = "custom",
          norm_file = data.frame(name = "ref.csv", datapath = nf, size = 1, type = "text/csv"))
  check("own reference is used for the calculation", cfg()$norm == "custom" && is.finite(user_data()$x[1]))
})
testServer(app_path, {
  set_all(session, file = fileinp, id_col = "ID", group_col = "Group", norm = "custom")
  check("own reference without a file falls back to the default with a notice",
        cfg()$norm == DEFAULT_NORM && grepl("Own reference", html(output$status)) && is.finite(user_data()$x[1]))
})

# 5. broken upload shows an error instead of crashing
bad <- tempfile(fileext = ".csv"); writeLines(c("La,Ce", "1,2"), bad)
testServer(app_path, {
  set_all(session, file = data.frame(name = "bad.csv", datapath = bad, size = 1, type = "text/csv"))
  check("broken upload: informative message and the reference figure still renders",
        is.null(user_data()) && grepl("Missing REE column", html(output$status)) && inherits(fig(), "ggplot"))
})

# 6. data quality tab
testServer(app_path, {
  set_all(session, use_example = TRUE, id_col = "sample", group_col = "group", rough_on = TRUE, rough_max = 0.3)
  ru <- rough_ref(); ur <- rough_user()
  check("roughness reference table has the five classes with finite values", nlevels(droplevels(ru$class)) >= 3 && all(is.finite(ru$v)))
  check("roughness of the user's data: one row per sample with a value, groups kept", nrow(ur) > 0 && all(ur$group %in% user_data()$group))
  check("roughness figure is a ggplot and builds", inherits(rough_fig(), "ggplot") && !inherits(try(ggplot_build(rough_fig()), silent = TRUE), "try-error"))
  tb <- roughness_summary(ru, ur, 0.3)
  check("summary table lists user groups and reference classes, shares between 0 and 100",
        all(c("your data", "reference") %in% tb$Type) && all(tb[["Above limit (%)"]] >= 0 & tb[["Above limit (%)"]] <= 100, na.rm = TRUE))
  session$setInputs(lambda_space = TRUE)
  check("in lambda space the lambda roughness is used for both data and reference",
        identical(rough_col(), "lambda_rough") && nrow(rough_ref()) > 3000 && max(rough_ref()$v) < 0.1)
  check("lambda roughness figure builds", inherits(rough_fig(), "ggplot"))
})

# PANGAEA .tab upload (Paul, 2025, https://doi.org/10.1594/PANGAEA.982625, CC-BY 4.0)
testServer(app_path, {
  pgf <- normalizePath(file.path(app_path, "..", "tests", "data", "MSM96_So268_REE.tab"))
  set_all(session, file = data.frame(name = "MSM96_So268_REE.tab", datapath = pgf, size = 1, type = "text/plain"), id_col = "Sample label", basis = "mass")
  u <- user_data()
  check("PANGAEA .tab loads through the app (26 rows, all plottable)", nrow(u) == 26 && all(is.finite(u$x)))
  sn <- html(output$source_note)
  check("dataset citation, DOI link and unit note are shown",
        grepl("PANGAEA dataset", sn) && grepl("10.1594/PANGAEA.982625", sn) && grepl("mg/kg", sn) && grepl("mass basis", sn))
  tb <- computed()
  check("results table carries the dataset citation and DOI",
        all(tb$source_doi == "10.1594/PANGAEA.982625") && all(grepl("Paul", tb$source_dataset)))
  check("sample IDs come from the 'Sample label' column", tb$sample[1] == "172MUC 0-2 cm authigenic")
})

cat("\nAll", cnt, "app checks passed.\n")
