# Base-R tests (no testthat needed). Run from the repository root: Rscript tests/run_tests.R
suppressMessages({library(ggplot2); library(MASS)})
owd <- setwd("app"); on.exit(setwd(owd))
for (f in list.files("R", full.names = TRUE)) source(f)

n_pass <- 0L
check <- function(desc, cond) {
  if (!isTRUE(cond)) stop("FAILED: ", desc, call. = FALSE)
  n_pass <<- n_pass + 1L; cat("ok  ", desc, "\n")
}
near <- function(a, b, tol = 1e-9) isTRUE(all(abs(a - b) <= tol))
rcols <- c("HREE_LREE", "MREE_MREEstar", "Ce_Ceanom", "Eu_Euanom")
cfg_fig3 <- default_cfg(norm = "PAAS_McLennan1989")      # the standard used for the published figure

# --- normalisation standards ------------------------------------------------------
for (nm in names(NORMALISERS)) {
  v <- as.data.frame(t(NORMALISERS[[nm]]))
  r <- ree_ratios(v, "mass", default_cfg(norm = nm))
  check(sprintf("%s: the standard itself gives all ratios = 1 (mass)", nm), near(unlist(r[1, rcols]), 1, 1e-12))
}
v <- as.data.frame(t(NORMALISERS$EUS_Bau2018 / ATOMIC_MASS))
check("molar version of a standard gives ratios = 1 (molar basis)",
      near(unlist(ree_ratios(v, "molar", default_cfg(norm = "EUS_Bau2018"))[1, rcols]), 1, 1e-12))
check("molar data normalised as mass would be biased (unit handling matters)",
      abs(ree_ratios(v, "mass", default_cfg(norm = "EUS_Bau2018"))$HREE_LREE - 1) > 0.1)
mass_v <- as.data.frame(t(NORMALISERS$PAAS_Pourmand2012))
check("unit prefix is irrelevant (x1000)", near(ree_ratios(mass_v * 1000, "mass")$HREE_LREE, 1, 1e-12))
check("invalid basis is rejected", tryCatch({ree_ratios(mass_v, "weight"); FALSE}, error = function(e) TRUE))
# mass route and molar route are the same calculation
x <- as.data.frame(t(c(La = 20, Ce = 30, Pr = 4, Nd = 15, Sm = 3, Eu = 0.7, Gd = 3.4, Tb = 0.5, Dy = 3.2, Ho = 0.7, Er = 2, Tm = 0.3, Yb = 1.9, Lu = 0.3)))
xm <- x * 1e-3 * matrix(ATOMIC_MASS[REE], 1)       # the same pattern expressed per mass
check("molar and mass routes agree (data converted to mass first)",
      near(unlist(ree_ratios(x, "molar")[1, rcols]), unlist(ree_ratios(xm, "mass")[1, rcols]), 1e-12))
custom <- default_cfg(norm = "custom", norm_custom = NORMALISERS$NASC_McLennan1989)
check("own reference behaves like the built-in one",
      near(unlist(ree_ratios(x, "mass", custom)[1, rcols]), unlist(ree_ratios(x, "mass", default_cfg(norm = "NASC_McLennan1989"))[1, rcols])))
check("own reference without values is rejected", tryCatch({norm_ppm(default_cfg(norm = "custom")); FALSE}, error = function(e) TRUE))

# --- LREE definition and gap handling -----------------------------------------------
ps <- as.data.frame(t(NORMALISERS$PAAS_Pourmand2012))
check("Pr-Nd-Sm definition: standard gives 1", near(ree_ratios(ps, "mass", default_cfg(lree = "PrNdSm"))$HREE_LREE, 1, 1e-12))
steep <- ps; steep[c("La", "Ce")] <- steep[c("La", "Ce")] * 2          # La enriched: only La-Pr-Nd feels it
check("La anomaly changes HREE/LREE for La-Pr-Nd but not for Pr-Nd-Sm",
      ree_ratios(steep, "mass")$HREE_LREE < 0.9 && near(ree_ratios(steep, "mass", default_cfg(lree = "PrNdSm"))$HREE_LREE, 1, 1e-12))
g <- ps; g[c("Pr", "Tb", "Tm", "Lu")] <- NA
rg <- ree_ratios(g, "mass")
check("missing Pr/Tb/Tm/Lu are estimated and counted", near(rg$HREE_LREE, 1, 1e-12) && rg$n_estimated == 4)
check("no filling when gaps are switched off", is.na(ree_ratios(g, "mass", default_cfg(gap = "none"))$HREE_LREE))
rm_ <- ree_ratios(g, "mass", default_cfg(gap = "groupmean"))
check("group-mean mode uses the measured members of each group (flat PAAS stays at 1)",
      near(rm_$HREE_LREE, 1, 1e-12) && near(rm_$MREE_MREEstar, 1, 1e-12) && rm_$n_estimated == 4)
check("without La, Pr falls back on Nd-Sm so Pr-Nd-Sm still works",
      near(ree_ratios({z <- ps; z[c("La", "Pr")] <- NA; z}, "mass", default_cfg(lree = "PrNdSm"))$HREE_LREE, 1, 1e-12))
check("required columns: Pr-Nd-Sm does not need La",
      length(check_ree_columns(c("Nd", "Sm", "Gd", "Dy", "Yb", "Er"), default_cfg(lree = "PrNdSm"))) == 0 &&
      "La" %in% check_ree_columns(c("Nd", "Sm", "Gd", "Dy", "Yb", "Er"), default_cfg()))
check("required columns: group-mean mode only needs one member per group",
      length(check_ree_columns(c("Nd", "Gd", "Yb"), default_cfg(gap = "groupmean"))) == 0)

# --- gadolinium -------------------------------------------------------------------------
gd3 <- ps; gd3$Gd <- gd3$Gd * 3
r3 <- ree_ratios(gd3, "mass")
check("Gd/Gd* (Bau and Dulski) is 3 for a 3x Gd spike on a flat pattern; Gd_anth = 2/3 of Gd",
      near(r3$Gd_Gdstar, 3, 1e-12) && near(r3$Gd_anth, gd3$Gd * 2 / 3, 1e-9))
ri <- ree_ratios(gd3, "mass", default_cfg(gd = "ignore"))
check("ignoring Gd restores the undisturbed ratios", near(c(ri$MREE_MREEstar, ri$Eu_Euanom), 1, 1e-12) && near(ri$Gd_Gdstar, 3, 1e-12))
check("Gd/Gd* with Sm^2/Nd and Tb-Dy mean also equal 3",
      near(ree_ratios(gd3, "mass", default_cfg(gd_method = "SmNd"))$Gd_Gdstar, 3, 1e-12) &&
      near(ree_ratios(gd3, "mass", default_cfg(gd_method = "TbDy"))$Gd_Gdstar, 3, 1e-12))
check("baseline correction: anomaly equal to the baseline gives no anthropogenic Gd",
      near(ree_ratios(gd3, "mass", default_cfg(gd_baseline = 3))$Gd_anth, 0, 1e-12))

# --- roughness ---------------------------------------------------------------------------
check("smooth pattern (the standard) has zero roughness", near(ree_ratios(ps, "mass")$REE_rough, 0, 1e-12))
rough <- ps; rough[c("Tb", "Ho", "Tm")] <- rough[c("Tb", "Ho", "Tm")] * c(1.3, 0.7, 1.3)
check("a zig-zag heavy tail has clearly positive roughness", ree_roughness(ree_normalise(rough, "mass", NORMALISERS$PAAS_Pourmand2012)) > 0.3)
check("roughness is NA when elements are missing (interpolated gaps are not used)", is.na(ree_ratios(g, "mass")$REE_rough))

# --- lambda (BLambdaR port) --------------------------------------------------------------
ci <- as.data.frame(t(LAMBDA_NORM)); bl <- blambda(ci, "mass")
check("chondrite pattern has lambda1..4 = 0", near(unlist(bl[1, c("lambda1", "lambda2", "lambda3", "lambda4")]), 0, 1e-8))
ref <- load_reference()
fo <- ref[ref$type == "forams" & complete.cases(ref[, REE]), ][c(5, 40, 80), ]
bf <- blambda(fo[, REE], "mass")     # reference values: P. Blaser's Blambda_full() on the same rows
check("lambda1/2/3, tau1 agree with the original BLambdaR function",
      near(bf$lambda1, c(9.689583188, 8.427436775, 12.460652789), 1e-6) &&
      near(bf$lambda2, c(2.80471091, 72.65752292, -14.14013087), 1e-5) &&
      near(bf$lambda3, c(65.01908183, -161.03191216, 48.94159481), 1e-4) &&
      near(bf$tau1, c(-0.0813572240, -0.3452159994, 0.1466490590), 1e-6))
check("Gd/Gd* from the fit and lambda_rough agree with the original",
      near(bf$Gd_MFR, c(0.9983144511, 1.0004495347, 0.9957495710), 1e-6) &&
      near(bf$lambda_rough, c(0.004443096830, 0.001686115143, 0.002441876052), 1e-8))
check("lambda fit tolerates two missing elements",
      is.finite(blambda({z <- fo[1, REE]; z[c("Ho", "Tm")] <- NA; z}, "mass")$lambda1))
check("lambda fit gives NA when fewer than 10 elements remain (underdetermined fits are not reported)",
      is.na(blambda({z <- fo[1, REE]; z[c("Ho", "Tm", "Pr")] <- NA; z}, "mass")$lambda1))
check("Gd anomaly from the polynomial fit is close to 1 for foraminifera", {
  a <- lambda_gd_anomaly(fo[, REE], "mass"); all(abs(a - 1) < 0.15) })
sp <- fo[1, REE]; sp$Gd <- sp$Gd * 3
check("Gd spike is recovered by the polynomial fit (about 3)", abs(lambda_gd_anomaly(sp, "mass") - 3) < 0.5)

# --- arrays vs. the published figure (Fig. 3, PAAS as in the paper) ---------------------
em  <- load_endmembers()
aut <- authigenic_array(em, cfg_fig3)
check("authigenic array starts at the HREE-rich seawater (~6.4, ~0.63)",
      near(c(aut$HREE_LREE[1], aut$MREE_MREEstar[1]), c(6.38, 0.634), c(0.03, 0.01)))
check("authigenic array ends at the nodule (~1.11, ~1.70)",
      near(c(tail(aut$HREE_LREE, 1), tail(aut$MREE_MREEstar, 1)), c(1.114, 1.697), 0.01))
det <- detrital_array(em, cfg_fig3); e <- split(det$endpoints, det$endpoints$end)
check("GAC end-member position", near(c(e$GAC$HREE_LREE, e$GAC$MREE_MREEstar), c(0.513, 1.065), 0.01))
check("CPS end-member position", near(c(e$CPS$HREE_LREE, e$CPS$MREE_MREEstar), c(0.647, 0.895), 0.01))
check("detrital polygon is a closed ring of finite points",
      all(is.finite(det$polygon$HREE_LREE)) && nrow(det$polygon) == 3 * 101)
aut_e <- authigenic_array(em, default_cfg(norm = "EUS_Bau2018", lree = "PrNdSm"))
check("arrays follow the chosen standard and LREE definition", !near(aut_e$HREE_LREE, aut$HREE_LREE, 0.05) && all(is.finite(aut_e$HREE_LREE)))
lam <- read.csv("data/arrays_lambda.csv")
arr <- attach_lambda(list(authigenic = authigenic_array(em), detrital = detrital_array(em)), lam)
check("lambda coordinates are attached to every array point",
      all(is.finite(arr$authigenic$lx)) && all(is.finite(arr$detrital$curves$ly)) && all(is.finite(arr$detrital$polygon$lx)))

# --- reference data ------------------------------------------------------------------------
rr <- reference_ratios(ref, cfg_fig3)
sw <- rr[rr$type == "seawater" & !grepl("GEOTRACES", rr$ref), ]
check("seawater maximum HREE/LREE ~6.4 (as in Fig. 3)", near(max(sw$x, na.rm = TRUE), 6.38, 0.02))
check("ratios recomputed for all reference rows with data", mean(is.finite(rr$x) & is.finite(rr$y)) > 0.99)
fo <- rr[rr$ref == "Osborne2017", ]
check("a foraminiferal point at (2.25, 1.80) as in Fig. 3", any(abs(fo$x - 2.25) < 0.02 & abs(fo$y - 1.80) < 0.02))
gt <- ref[grepl("GEOTRACES", ref$ref), ]
check("GEOTRACES: two depth classes (<= 200 m and > 200 m), all with La and finite coordinates",
      setequal(unique(gt$ref), c("GEOTRACES2025_shallow", "GEOTRACES2025_deep")) && all(is.finite(gt$La)) &&
      all(gt$depth_m[gt$ref == "GEOTRACES2025_shallow"] <= 200) && all(gt$depth_m[gt$ref == "GEOTRACES2025_deep"] > 200) &&
      all(is.finite(gt$lx)) && nrow(gt) > 2000)
check("every reference label has data", all(names(REF_LABELS) %in% ref$ref))
rp <- reference_ratios(ref, default_cfg(lree = "PrNdSm", norm = "EUS_Bau2018"))
fm <- rr$type == "forams"
check("reference ratios respond to the configuration (foraminifera shift by >10% with Pr-Nd-Sm)",
      abs(median(rp$x[fm], na.rm = TRUE) / median(rr$x[fm], na.rm = TRUE) - 1) > 0.1)

# --- density fields -----------------------------------------------------------------------
pip <- function(px, py, poly) { # ray casting
  inside <- FALSE; n <- nrow(poly); j <- n
  for (i in seq_len(n)) {
    if (((poly$y[i] > py) != (poly$y[j] > py)) &&
        (px < (poly$x[j] - poly$x[i]) * (py - poly$y[i]) / (poly$y[j] - poly$y[i]) + poly$x[i])) inside <- !inside
    j <- i
  }
  inside
}
coverage_of <- function(xx, yy, hp) mean(vapply(seq_along(xx), function(i)
  any(vapply(split(hp, hp$piece), function(pc) pip(xx[i], yy[i], pc), TRUE)), TRUE))
set.seed(1); xx <- rnorm(300); yy <- rnorm(300)
for (cv in c(0.5, 0.9)) {
  frac <- coverage_of(xx, yy, hdr_polygons(xx, yy, coverage = cv))
  check(sprintf("HDR polygon at %.0f%% coverage contains ~%.0f%% of the points (got %.0f%%)", 100 * cv, 100 * cv, 100 * frac),
        abs(frac - cv) < 0.04)
}
xb <- rnorm(2500); yb <- rnorm(2500)
frac <- coverage_of(xb[1:600], yb[1:600], hdr_polygons(xb, yb, coverage = 0.9))
check(sprintf("large samples (grid-interpolated density) still honour the coverage (got %.0f%% of a subsample)", 100 * frac),
      abs(frac - 0.9) < 0.05)
check("too few points give no field", is.null(hdr_polygons(1:5, 1:5)))

# --- user input ----------------------------------------------------------------------------
ex <- read_user_file("data/example_user_data.csv")
u <- prepare_user_data(ex, "concentrations", "mass", cfg_fig3, group_col = "group", id_col = "sample")
check("example data: 26 samples, 'basis' column overrides the global setting",
      nrow(u$data) == 26 && all(u$data$basis == "molar"))
check("example data: Labrador Sea deep water falls on the authigenic array side (HREE/LREE 2.5-4.5)",
      all(u$data$x[u$data$group == "LS deep water"] > 2.5 & u$data$x[u$data$group == "LS deep water"] < 4.5))
check("column mapping is reported (14 elements) and extra columns are kept",
      nrow(u$mapping) == 14 && all(c("water_depth_m") %in% names(u$extra)) && nrow(u$extra) == 26)
ul <- prepare_user_data(ex, "concentrations", "mass", cfg_fig3, group_col = "group", lambda = TRUE)
check("lambda option adds lambda coordinates and lambda roughness",
      all(c("lx", "ly", "lambda_rough", "lambda1", "tau1") %in% names(ul$data)) && sum(is.finite(ul$data$lx)) >= 24)

tmp <- tempfile(fileext = ".csv")
writeLines(c("Sample;Site;La (ppb);Ce_ppb;Pr;Nd;Sm;Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu",
             "a;X;38,2;79,6;8,83;33,9;5,55;1,08;4,66;0,774;4,68;0,991;2,85;0,405;2,82;0,433"), tmp)
us <- prepare_user_data(read_user_file(tmp), "concentrations", "mass", cfg_fig3, id_col = "Sample", group_col = "Site")
check("semicolon CSV with decimal commas and 'La (ppb)' headers is read; PAAS plots at (1,1)",
      near(c(us$data$x, us$data$y, us$data$Ce_Ceanom), 1, 1e-9) && us$data$group == "X")
check("mapping shows which header was used for La", us$mapping$column[us$mapping$element == "La"] == "La (ppb)")

bad <- tempfile(fileext = ".csv"); writeLines(c("La,Ce,Nd", "1,2,3"), bad)
err <- tryCatch(prepare_user_data(read_user_file(bad), "concentrations", "mass"), error = function(e) conditionMessage(e))
check("missing REE columns give an informative error", grepl("Missing REE column", err) && grepl("Gd", err))

rt <- tempfile(fileext = ".csv"); writeLines(c("id,HREE/LREE,MREE/MREE*", "s1,1.5,1.2", "s2,3.0,0.9"), rt)
ur <- prepare_user_data(read_user_file(rt), "ratios", xcol = "HREE/LREE", ycol = "MREE/MREE*")
check("pre-computed ratio upload", near(ur$data$x, c(1.5, 3)) && near(ur$data$y, c(1.2, 0.9)))
gr <- guess_ratio_columns(c("id", "HREE/LREE", "MREE/MREE*"))
check("ratio columns are guessed", gr$x == "HREE/LREE" && gr$y == "MREE/MREE*")
nf <- tempfile(fileext = ".csv"); write.csv(as.data.frame(t(NORMALISERS$EUS_Bau2018)), nf, row.names = FALSE)
check("own reference file is read (one row, 14 REE)", near(read_norm_file(nf), NORMALISERS$EUS_Bau2018))
nf2 <- tempfile(fileext = ".csv"); write.csv(rbind(NORMALISERS$EUS_Bau2018, NORMALISERS$EUS_Bau2018), nf2, row.names = FALSE)
check("own reference file with two rows is rejected", tryCatch({read_norm_file(nf2); FALSE}, error = function(e) grepl("exactly one", conditionMessage(e))))

# --- plot builds ------------------------------------------------------------------------------
arrays <- list(authigenic = aut, detrital = det)
tf <- tempfile(fileext = ".pdf")
render <- function(p) { ggsave(tf, p, width = 8, height = 6); file.exists(tf) && file.size(tf) > 5000 }
check("figure with user data renders to PDF", render(plot_martin(rr, arrays, u$data, list(user_fields = TRUE))))
check("figure without user data renders", render(plot_martin(rr, arrays, NULL)))
ud <- u$data; ud$.fill <- ex$water_depth_m
check("continuous fill (depth) renders", render(plot_martin(rr, arrays, ud, list(fill_label = "depth (m)"))))
check("log10 axes render", render(plot_martin(rr, arrays, u$data, list(log = TRUE))))
check("window restricted to the arrays renders and spans the arrays",
      render(plot_martin(rr, arrays, u$data, list(window = "arrays"))) &&
      { l <- plot_limits(rr, arrays, u$data, "ratio", "arrays"); l$xlim[1] < min(aut$HREE_LREE) && l$xlim[2] > max(aut$HREE_LREE) })
arr_l <- attach_lambda(arrays, lam); rl <- ref
check("lambda space renders", render(plot_martin(rl, arr_l, ul$data, list(space = "lambda", user_fields = TRUE))))
check("groups can be fixed to a colour table when some are hidden",
      render(plot_martin(rr, arrays, u$data[u$data$group == "LS shelf", ], list(group_levels = c("LS deep water", "LS shelf")))))
# ---- Y, spectra, natural sort, row ids ----------------------------------------------------------
ex <- read_user_file("data/example_user_data.csv")
check("natural_sort orders numbers numerically and ignores case", identical(natural_sort(c("b", "a10", "a2", "A1")), c("A1", "a2", "a10", "b")))
check("Y column detection", identical(detect_y_column(c("sample", "Yb", "Y_ppm", "La")), "Y_ppm") && is.null(detect_y_column(c("Yb", "y"))))
ex$Y <- 20; names(ex)[names(ex) == "La"] <- "La_ppm"
pp <- prepare_user_data(ex, "concentrations", "mass", default_cfg(), "group", "sample")
check("Y is kept for the spectra and listed among the columns used", all(pp$conc$Y == 20) && "Y" %in% pp$mapping$element && !"Y" %in% names(pp$extra))
check("rows carry the row number of the file", identical(pp$data$row, seq_len(nrow(ex))))
ex2 <- ex; ex2[3, ] <- NA
check("row numbers survive empty rows in the file", identical(prepare_user_data(ex2, "concentrations", "mass", default_cfg())$data$row, setdiff(seq_len(nrow(ex)), 3L)))
check("n_REE_measured counts measured REE", all(pp$data$n_REE_measured %in% 13:14) && is.logical(pp$meas))
nn <- ree_normalise_y(data.frame(Y = 27.31, La = 44.56, Lu = 0.4386, check.names = FALSE)[, c("Y", "La", "Lu")], "mass", default_cfg())
check("Y and La..Lu normalise to 1 for PAAS (Pourmand)", abs(nn[, "Y"] - 1) < 1e-12 && abs(nn[, "La"] - 1) < 1e-12 && is.na(nn[, "Ce"]))
molar <- data.frame(t(NORMALISERS$PAAS_Pourmand2012 * 1e9 / ATOMIC_MASS[REE])); names(molar) <- REE
check("molar data (pmol/kg) normalise to 1 against the molar reference", all(abs(ree_normalise_y(molar, "molar", default_cfg())[, REE] - 1) < 1e-12))
sp <- ree_spectra(pp$conc, "mass", default_cfg(), "Dy")
check("secondary normalisation sets the chosen element to 1", all(abs(sp$n2[, "Dy"] - 1) < 1e-12, na.rm = TRUE))
radii <- load_radii("data/ionic_radii_3plus.csv")
check("Shannon radii table: 16 elements, decreasing from La to Lu", length(radii) == 16 && all(diff(radii[REE]) < 0))
long <- spectra_long(sp$n, data.frame(sid = pp$data$row, row = pp$data$row, sample = pp$data$sample, group = pp$data$group), radii)
check("long spectra table has one row per sample and element", nrow(long) == nrow(pp$conc) * 15)
check("spectra plots build (element axis, radius axis, numeric colour)",
      inherits(plot_spectra(long, list(x = "element")), "ggplot") && inherits(plot_spectra(long, list(x = "radius", log = FALSE)), "ggplot") &&
      inherits(plot_spectra(transform(long, .col = row), list()), "ggplot"))
# window and legend
cfg0 <- default_cfg(); ref0 <- reference_ratios(load_reference("data/reference_data.csv"), cfg0); em0 <- load_endmembers("data/endmembers.csv")
arr0 <- attach_lambda(list(authigenic = authigenic_array(em0, cfg0), detrital = detrital_array(em0, cfg0)), read.csv("data/arrays_lambda.csv"))
us <- data.frame(x = c(rep(2, 300), 50), y = c(rep(1, 300), 40))
l1 <- plot_limits(ref0, arr0, us, "ratio", "auto")
check("a single extreme outlier does not stretch the window (>= 200 points)", l1$xlim[2] < 20 && l1$ylim[2] < 10)
lg <- plot_limits(ref0, arr0, NULL, "ratio", "auto", ext = list(x = c(0.2, 12), y = c(0.4, 2.5)), log = TRUE)
check("log window holds the field outlines", lg$xlim[1] < 0.2 && lg$xlim[2] > 12 && lg$ylim[1] < 0.4)
check("figure without legend carries limits", {g <- plot_martin(ref0, arr0, NULL, list(legend = FALSE, log = TRUE)); is.list(attr(g, "limits"))})

# ---- anomalies and gap filling after Barrat and Bayon (2024) -------------------------------------
lin <- matrix(exp(-0.08 * (REE_Z - 64)), 1, 14, dimnames = list(NULL, REE))   # straight in log-Z (Pm omitted)
nrm <- NORMALISERS$PAAS_Pourmand2012
conc_of <- function(m) as.data.frame(m * rep(nrm, each = nrow(m)))
filled_ok <- function(missing, back) {
  m <- lin; m[, missing] <- NA
  f <- ree_fill(m)$n
  all(abs(f[, back] / lin[, back] - 1) < 1e-12)
}
check("gap filling recovers Pr, Tb, Ho, Tm, Er and Lu on a pattern straight in atomic number (Pm omitted)",
      filled_ok("Pr", "Pr") && filled_ok("Tb", "Tb") && filled_ok("Ho", "Ho") && filled_ok("Tm", "Tm") &&
      filled_ok(c("Tm", "Er"), "Tm") && filled_ok("Lu", "Lu") && filled_ok(c("La", "Pr"), "Pr"))
cx <- lin; cx[, "Ce"] <- cx[, "Ce"] * 0.5
r_geo <- ree_ratios(conc_of(cx), "mass", default_cfg())
r_ar  <- ree_ratios(conc_of(cx), "mass", default_cfg(ce_form = "arithLaPr"))
r_gl  <- ree_ratios(conc_of(cx), "mass", default_cfg(ce_form = "geoLaPr"))
check("Ce/Ce*: both geometric forms return the imposed 0.5 exactly (no La anomaly); the arithmetic form does not",
      near(r_geo$Ce_Ceanom, 0.5, 1e-12) && near(r_gl$Ce_Ceanom, 0.5, 1e-12) && abs(r_ar$Ce_Ceanom - 0.5) > 1e-3 &&
      near(r_ar$Ce_Ceanom, 2 * 0.5 * lin[, "Ce"] / (lin[, "La"] + lin[, "Pr"]), 1e-12))
la <- lin; la[, "La"] <- la[, "La"] * 2
check("with a La anomaly (x2) the La-Pr forms of Ce/Ce* are biased, the Pr-Nd form is not",
      near(ree_ratios(conc_of(la), "mass", default_cfg())$Ce_Ceanom, 1, 1e-12) &&
      abs(ree_ratios(conc_of(la), "mass", default_cfg(ce_form = "geoLaPr"))$Ce_Ceanom - 1) > 0.2)
check("La/La* = La Nd^2 / Pr^3 returns the imposed factor", near(ree_ratios(conc_of(la), "mass", default_cfg())$La_Laanom, 2, 1e-12))
gx <- lin; gx[, "Gd"] <- gx[, "Gd"] * 3
check("Gd/Gd*: geometric methods (Sm-Tb, Nd-Sm, Tb-Dy) return the imposed factor 3 exactly",
      all(vapply(c("BDgeo", "SmNd", "TbDy"), function(m) near(ree_ratios(conc_of(gx), "mass", default_cfg(gd_method = m))$Gd_Gdstar, 3, 1e-12), logical(1))))
check("Gd/Gd* with the arithmetic Bau-Dulski estimate is close to, but not exactly, 3 on the same pattern",
      {v <- ree_ratios(conc_of(gx), "mass", default_cfg(gd_method = "BD1996"))$Gd_Gdstar; abs(v - 3) < 0.02 && abs(v - 3) > 1e-9})

rr <- reference_ratios(load_reference("data/reference_data.csv"), default_cfg())
check("reference data carry REE and lambda roughness", sum(is.finite(rr$REE_rough)) > 3000 && sum(is.finite(rr$lambda_rough)) > 3000)
check("lambda roughness of the reference data matches a fresh fit (first 50 rows)",
      near(blambda(rr[1:50, REE], rr$basis[1:50])$lambda_rough, rr$lambda_rough[1:50], 1e-5))
set.seed(1)
dd <- data.frame(class = factor(rep(REF_CLASS_LEVELS[1:2], each = 40), levels = REF_CLASS_LEVELS), v = rlnorm(80, -1, 0.5))
du <- data.frame(group = rep(c("a", "b", "c"), c(30, 30, 2)), v = rlnorm(62, -0.5, 0.6))
check("plot_roughness builds with log and linear axes, and with tiny groups (rug only)",
      all(vapply(c(TRUE, FALSE), function(lg) inherits(ggplot_build(plot_roughness(dd, du, list(log = lg, limit = 0.3))), "ggplot_built"), logical(1))))
check("roughness_summary share above limit is correct", {
  sm <- roughness_summary(NULL, data.frame(group = "a", v = c(0.1, 0.2, 0.4, 0.5)), 0.3); sm[["Above limit (%)"]] == 50 && sm$n == 4 })

cat("\nAll", n_pass, "checks passed.\n")