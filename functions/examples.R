# Examples for REE_functions.R. Run from the repository root: Rscript functions/examples.R
source("functions/REE_functions.R")
library(ggplot2)

dat <- read.csv("app/data/example_user_data.csv")    # seawater, pmol/kg; has a "basis" column

# 1. Ratios, anomalies, anthropogenic Gd and roughness (app defaults: PAAS Pourmand 2012, La-Pr-Nd)
res <- REE_ratios(dat, rough_max = 0.3)
head(res[, c("sample", "HREE_LREE", "MREE_MREEstar", "Ce_Cestar", "Gd_Gdstar", "REE_rough", "screen")])

# the same settings as Fig. 3 of Blaser et al. (2020)
fig3 <- REE_ratios(dat, norm = "PAAS_McLennan1989", ce_form = "arithLaPr")

# Gd contaminated samples: ignore measured Gd, estimate natural Gd geometrically from Sm and Tb
res_gd <- REE_ratios(dat, gd = "ignore", gd_method = "BDgeo")

# 2. Cross-plot
ggplot(res[res$screen == "ok", ], aes(HREE_LREE, MREE_MREEstar, colour = group)) +
  geom_point(size = 2.5) + scale_x_log10() + scale_y_log10() +
  labs(x = "HREE/LREE", y = "MREE/MREE*") + theme_classic()

# 3. Lambda space (experimental): the app plots -lambda1 against -lambda2
lam <- REE_lambda(dat, rough_max = 0.03)
ggplot(lam[lam$screen == "ok", ], aes(-lambda1, -lambda2, colour = group)) +
  geom_point(size = 2.5) + labs(x = expression(-lambda[1]), y = expression(-lambda[2])) + theme_classic()

# 4. Spectra, coloured by water depth, only the screened samples
sp <- REE_spectrum(dat, colour = "water_depth_m", keep = res$screen == "ok")
sp_nd <- REE_spectrum(dat, colour = "group", norm2 = "Nd", x = "radius")
head(sp$normalised)

# 5. Distribution of the REE roughness by group (as in the app's roughness tab)
ggplot(res, aes(REE_rough, colour = group)) + geom_density(na.rm = TRUE) +
  geom_vline(xintercept = 0.3, linetype = 2) + labs(x = "REE roughness") + theme_classic()
