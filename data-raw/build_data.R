# Build the files in data/ from the three raw compilation exports.
# Usage: Rscript data-raw/build_data.R <dir with Martin_plot_*_export.csv>
# Pre-computed ratio columns of the exports are deliberately dropped: every ratio is
# recomputed in the app from the raw concentrations (see R/ree_core.R).
args <- commandArgs(trailingOnly = TRUE)
raw <- if (length(args)) args[1] else "data-raw"
for (f in c("01_ree_core.R", "02_lambda.R", "03_arrays.R")) source(file.path("app/R", f))
na <- c("NA", "")
d  <- read.csv(file.path(raw, "Martin_plot_data_export.csv"), na.strings = na, stringsAsFactors = FALSE)
a  <- read.csv(file.path(raw, "Martin_plot_authigenic_array_export.csv"), na.strings = na, stringsAsFactors = FALSE)
dt <- read.csv(file.path(raw, "Martin_plot_detrital_array_export.csv"), na.strings = na, stringsAsFactors = FALSE)

# ---- reference compilation -------------------------------------------------
# Concentration basis per data type: foraminifera are mass-based (ppm-type); seawater and
# pore water are molar (pmol/kg). Verified against the published Fig. 3 (see README).
d$basis <- ifelse(d$type == "forams", "mass", "molar")
d$site  <- ifelse(d$type == "pore water", d$method, NA)
d$method <- ifelse(d$type == "forams", d$method, NA)
d$sample_id <- sprintf("%s_%03d", d$ref, ave(seq_len(nrow(d)), d$ref, FUN = seq_along))
ref <- d[, c("sample_id", "ref", "type", "basis", "lat", "lon", "bottom_depth", "depth",
             "site", "method", "cmtop", "surface", REE)]
names(ref)[names(ref) == "bottom_depth"] <- "bottom_depth_m"
names(ref)[names(ref) == "depth"] <- "depth_m"
names(ref)[names(ref) == "method"] <- "leach_method"
names(ref)[names(ref) == "cmtop"] <- "depth_cm_or_mbsf"

# ---- GEOTRACES IDP2025 seawater (pmol/kg) ---------------------------------------
# Flagged samples are dropped, as are rows that cannot yield both axes under the default settings
# (this removes all rows without La) and rows identical to seawater already in the compilation.
gt <- read.csv(file.path(raw, "GEOTRACES_IDP2025_REE_subset.csv"), na.strings = na, stringsAsFactors = FALSE)
n0 <- nrow(gt)
gt <- gt[gt$flag == "No problem reported (0)" & is.finite(gt$La) & is.finite(gt$depth), ]
r  <- ree_ratios(gt[, REE], "molar", default_cfg(gap = "interpolate"))
gt <- gt[is.finite(r$HREE_LREE) & is.finite(r$MREE_MREEstar), ]
key <- function(x) paste(signif(x$La, 5), signif(x$Nd, 5), signif(x$Yb, 5))
dup <- key(gt) %in% key(ref[ref$type == "seawater", ])
gt <- gt[!dup, ]
gtr <- data.frame(sample_id = sprintf("GEOTRACES_%s_st%s_%d", gt$cruise, gt$station, seq_len(nrow(gt))),
                  ref = ifelse(gt$depth <= 200, "GEOTRACES2025_shallow", "GEOTRACES2025_deep"),
                  type = "seawater", basis = "molar", lat = gt$lat, lon = gt$lon,
                  bottom_depth_m = gt$depth_bottom, depth_m = gt$depth, site = gt$cruise,
                  leach_method = NA, depth_cm_or_mbsf = NA, surface = NA, gt[, REE])
cat(sprintf("GEOTRACES: %d rows read, %d kept (%d without La/flagged/incomplete removed, %d duplicates of the compilation removed)\n",
            n0, nrow(gtr), n0 - nrow(gtr) - sum(dup), sum(dup)))
ref <- rbind(ref, gtr)

# ---- lambda coordinates (-lambda1, -lambda2), BLambdaR default settings; independent of the
# normalisation standard, hence pre-computed here (about a minute) ------------------------------
bl <- blambda(ref[, REE], ref$basis)
ref$lx <- round(-bl$lambda1, 5)
ref$ly <- round(-bl$lambda2, 5)
write.csv(ref, "app/data/reference_data.csv", row.names = FALSE, na = "")

# ---- end-members -----------------------------------------------------------
pick <- function(df, i) as.numeric(df[i, REE])
sw <- d[d$ref == "Hathorne2015", ]
sw <- sw[which.max(ree_ratios(sw[, REE], "molar")$HREE_LREE), ]          # most HREE-rich seawater
nod <- a[nrow(a), ]                                                       # r = 1 row: nodule
gac <- dt[dt$endmember %in% "Greenland exposed Archaean crust", ][1, ]
mor <- dt[dt$endmember %in% "mean_MORB_Gale2013", ][1, ]
cps <- dt[dt$endmember %in% "Canadian Precambrian shield", ]
cps <- cps[is.na(cps$Er), ][1, ]                                          # raw row: Er and Tm missing
em <- rbind(
  data.frame(id = "seawater_HREE", label = "HREE-rich seawater", basis = "molar",
             source = "Hathorne et al. (2015), 47.66S 4.27E, 400 m (pmol/kg)", t(pick(sw, 1))),
  data.frame(id = "nodule_CCFZ", label = "MREE-rich diagenetic Fe-Mn nodule", basis = "mass",
             source = "Bau et al. (2014), Clarion-Clipperton Fracture Zone (ug/g)", t(pick(nod, 1))),
  data.frame(id = "crust_GAC", label = "Greenland Archaean crust (GAC)", basis = "mass",
             source = "Wedepohl et al. (1991)", t(pick(gac, 1))),
  data.frame(id = "MORB_Gale2013", label = "Mean MORB", basis = "mass",
             source = "Gale et al. (2013) (ug/g)", t(pick(mor, 1))),
  data.frame(id = "crust_CPS", label = "Canadian Precambrian shield (CPS)", basis = "mass",
             source = "Shaw et al. (1986)", t(pick(cps, 1)))
)
names(em)[-(1:4)] <- REE
em[REE] <- lapply(em[REE], function(v) signif(v, 7))
write.csv(em, "app/data/endmembers.csv", row.names = FALSE, na = "")
write.csv(transform(arrays_lambda_table(em), lx = round(lx, 5), ly = round(ly, 5)),
          "app/data/arrays_lambda.csv", row.names = FALSE)

# ---- example user data: Labrador Sea seawater (Filippova et al., 2017) -----
fl <- d[d$ref == "Filippova2017", ]
ex <- data.frame(sample = fl$sample_id,
                 group = ifelse(fl$bottom_depth <= 300, "LS shelf", "LS deep water"),
                 water_depth_m = fl$depth, basis = "molar", fl[, REE])
write.csv(ex, "app/data/example_user_data.csv", row.names = FALSE, na = "")
cat("reference rows:", nrow(ref), "| end-members:", nrow(em), "| example rows:", nrow(ex), "\n")
