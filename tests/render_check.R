# Render the figure to a PNG for visual checks.
# Usage: Rscript tests/render_check.R out.png [norm] [lree] [space] [log] [example]
#   norm  : id in NORMALISERS (default PAAS_Pourmand2012), lree: LaPrNd|PrNdSm, space: ratio|lambda,
#   log   : TRUE|FALSE, example: TRUE adds the example data
suppressMessages({library(ggplot2); library(MASS)})
a <- commandArgs(TRUE)
out <- normalizePath(a[1], mustWork = FALSE)
get <- function(i, d) if (length(a) >= i && nzchar(a[i])) a[i] else d
owd <- setwd("app"); on.exit(setwd(owd))
for (f in list.files("R", full.names = TRUE)) source(f)
cfg <- default_cfg(norm = get(2, DEFAULT_NORM), lree = get(3, "LaPrNd"))
space <- get(4, "ratio"); lg <- identical(get(5, "FALSE"), "TRUE"); ex <- identical(get(6, "FALSE"), "TRUE")
ref <- reference_ratios(load_reference(), cfg)
em <- load_endmembers()
arr <- attach_lambda(list(authigenic = authigenic_array(em, cfg), detrital = detrital_array(em, cfg)),
                     read.csv("data/arrays_lambda.csv"))
user <- NULL
if (ex) user <- prepare_user_data(read_user_file("data/example_user_data.csv"), "concentrations", "mass", cfg,
                                  group_col = "group", id_col = "sample", lambda = identical(space, "lambda"))$data
p <- plot_martin(ref, arr, user, list(space = space, log = lg))
ggsave(out, p, width = 9, height = 6.6, dpi = 110)
cat("written", out, "\n")
