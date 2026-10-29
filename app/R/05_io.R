# ---------------------------------------------------------------------------
# Reading and checking user uploads
# ---------------------------------------------------------------------------

#' Read a user file (csv / tsv / txt; xlsx if readxl is available)
read_user_file <- function(path, name = basename(path)) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Excel files are not supported in this build - please export your sheet as CSV.")
    }
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }
  first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  first <- sub("^﻿", "", first)
  seps <- c(`,` = lengths(regmatches(first, gregexpr(",", first))),
            `;` = lengths(regmatches(first, gregexpr(";", first))),
            `\t` = lengths(regmatches(first, gregexpr("\t", first))))
  sep <- names(which.max(seps))
  dec <- if (sep == ";") "," else "."
  utils::read.table(path, header = TRUE, sep = sep, dec = dec, quote = "\"", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", "", "n.d.", "nd", "NaN"),
                    fileEncoding = "UTF-8-BOM", strip.white = TRUE)
}

#' Map column names to REE symbols. A column matches when it starts with an element symbol that
#' is not followed by another letter ("La", "La_ppm", "Nd (ppb)", "Yb.nmol/kg", ...).
detect_ree_columns <- function(nms) {
  m <- regmatches(nms, regexpr("^[A-Za-z]{2}(?![A-Za-z])", nms, perl = TRUE))
  sym <- rep(NA_character_, length(nms))
  hit <- grepl("^[A-Za-z]{2}(?![A-Za-z])", nms, perl = TRUE)
  sym[hit] <- m
  idx <- match(tolower(sym), tolower(REE))
  out <- setNames(nms[!is.na(idx)], REE[idx[!is.na(idx)]])
  out[!duplicated(names(out))]
}

#' Name of the yttrium column ("Y", "Y_ppm", "Y (ppb)" ...), or NULL
detect_y_column <- function(nms) {
  hit <- grep("^Y(?![A-Za-z])", nms, perl = TRUE, value = TRUE)
  if (length(hit)) hit[1] else NULL
}

#' Sort labels "naturally" and case-insensitively: sample2 < sample10, a < B
natural_sort <- function(x) {
  x <- unique(as.character(x))
  key <- vapply(x, function(s) {
    parts <- regmatches(s, gregexpr("[0-9]+|[^0-9]+", s))[[1]]
    paste(ifelse(grepl("^[0-9]+$", parts), formatC(suppressWarnings(as.numeric(parts)), width = 15, flag = "0", format = "f", digits = 0),
                 tolower(parts)), collapse = "")
  }, character(1))
  x[order(key, x)]
}

#' Guess the two ratio columns of a table of pre-computed ratios
guess_ratio_columns <- function(nms) {
  x <- grep("HREE|HLREE|H[/_.]?L|slope", nms, ignore.case = TRUE, value = TRUE)
  y <- grep("MREE|bulge", nms, ignore.case = TRUE, value = TRUE)
  list(x = if (length(x)) x[1] else NA_character_, y = if (length(y)) y[1] else NA_character_)
}

to_number <- function(v) {
  if (is.numeric(v)) return(v)
  suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(v)), fixed = TRUE)))
}

#' Read an own normalising reference: a file with exactly one row and the 14 REE (mass basis, ppm)
read_norm_file <- function(path, name = basename(path)) {
  df <- read_user_file(path, name)
  df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]
  if (nrow(df) != 1) stop("The reference file must contain exactly one data row (found ", nrow(df), ").")
  cols <- detect_ree_columns(names(df))
  miss <- setdiff(REE, names(cols))
  if (length(miss)) stop("Own reference: missing column(s) ", paste(miss, collapse = ", "), ".")
  v <- vapply(REE, function(el) to_number(df[[cols[[el]]]]), numeric(1))
  if (any(!is.finite(v) | v <= 0)) stop("Own reference: all values must be positive numbers.")
  v
}

#' Turn an upload into plot-ready data.
#' @param mode "concentrations" or "ratios"
#' @param cfg calculation options, see `default_cfg()`
#' @param lambda also compute the lambda coordinates (-lambda1, -lambda2; slow for large files)
#' @return list(data, messages, mapping, extra, conc, meas): `mapping` lists which file column was used
#'   for which element, `extra` holds all other columns of the file, `conc` the Y + La..Lu table and
#'   `meas` a logical matrix of measured values (all row-aligned with `data`; `data$row` is the row
#'   number in the uploaded file)
prepare_user_data <- function(df, mode = "concentrations", basis = "mass", cfg = default_cfg(),
                              group_col = NULL, id_col = NULL, xcol = NULL, ycol = NULL, lambda = FALSE) {
  msgs <- character()
  keep <- rowSums(!is.na(df)) > 0
  row_id <- which(keep)                       # row number in the file (header not counted)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) stop("The file contains no data rows.")
  sample <- if (!is.null(id_col) && id_col %in% names(df)) as.character(df[[id_col]]) else
    sprintf("row%03d", row_id)
  group <- if (!is.null(group_col) && group_col %in% names(df)) as.character(df[[group_col]]) else
    rep("my data", nrow(df))
  group[is.na(group) | group == ""] <- "(no group)"
  mapping <- data.frame(element = character(), column = character(), stringsAsFactors = FALSE)
  conc <- NULL; meas <- NULL

  if (mode == "ratios") {
    if (is.null(xcol) || is.na(xcol) || !xcol %in% names(df)) stop("Pick the HREE/LREE column.")
    if (is.null(ycol) || is.na(ycol) || !ycol %in% names(df)) stop("Pick the MREE/MREE* column.")
    out <- data.frame(row = row_id, sample = sample, group = group, x = to_number(df[[xcol]]),
                      y = to_number(df[[ycol]]), stringsAsFactors = FALSE)
    mapping <- data.frame(element = c("HREE/LREE", "MREE/MREE*"), column = c(xcol, ycol),
                          stringsAsFactors = FALSE)
    used <- c(xcol, ycol)
    if (isTRUE(lambda)) msgs <- c(msgs, "Lambda space needs REE concentrations; it is not available for pre-computed ratios.")
  } else {
    cols <- detect_ree_columns(names(df))
    miss <- check_ree_columns(names(cols), cfg)
    if (length(miss)) {
      stop("Missing REE column(s): ", paste(miss, collapse = ", "),
           ". Column names must start with the element symbol (e.g. 'La', 'La_ppm').")
    }
    conc <- as.data.frame(matrix(NA_real_, nrow(df), length(REE), dimnames = list(NULL, REE)))
    for (el in names(cols)) conc[[el]] <- to_number(df[[cols[[el]]]])
    ycol_name <- detect_y_column(names(df))
    conc$Y <- if (!is.null(ycol_name)) to_number(df[[ycol_name]]) else NA_real_    # spectra tab only
    cm <- as.matrix(conc[, c("Y", REE)])
    nonpos <- sum(!is.na(cm) & cm <= 0, na.rm = TRUE)
    if (nonpos) msgs <- c(msgs, sprintf("%d zero/negative concentration value(s) were treated as missing.", nonpos))
    meas <- !is.na(cm) & cm > 0; dimnames(meas) <- list(NULL, c("Y", REE))
    b <- if ("basis" %in% tolower(names(df))) df[[names(df)[tolower(names(df)) == "basis"][1]]] else basis
    b <- tolower(trimws(as.character(b)))
    lg <- if (identical(cfg$gd_method, "BLambdaR")) lambda_gd_anomaly(conc, b) else NULL
    r <- ree_ratios(conc, b, cfg, lambda_gd_anom = lg)
    out <- data.frame(row = row_id, sample = sample, group = group, x = r$HREE_LREE, y = r$MREE_MREEstar,
                      Ce_Ceanom = r$Ce_Ceanom, Eu_Euanom = r$Eu_Euanom,
                      Gd_Gdstar = r$Gd_Gdstar, Gd_anth = r$Gd_anth, REE_rough = r$REE_rough,
                      n_estimated = r$n_estimated, n_REE_measured = rowSums(meas[, REE, drop = FALSE]),
                      basis = rep_len(b, nrow(df)), stringsAsFactors = FALSE)
    if (isTRUE(lambda)) {
      bl <- blambda(conc, b, anomalies = c(Gd = identical(cfg$gd, "ignore")))
      out$lx <- -bl$lambda1; out$ly <- -bl$lambda2
      out$lambda_rough <- bl$lambda_rough
      out <- cbind(out, bl[, c("lambda0", "lambda1", "lambda2", "lambda3", "lambda4", "tau1", "tau2", "tau3", "tau4")])
    }
    mapping <- data.frame(element = names(cols), column = unname(cols), stringsAsFactors = FALSE)
    mapping <- mapping[match(REE[REE %in% mapping$element], mapping$element), ]
    used <- unname(cols)
    if (!is.null(ycol_name)) {
      mapping <- rbind(data.frame(element = "Y", column = ycol_name, stringsAsFactors = FALSE), mapping)
      used <- c(used, ycol_name)
    }
    nm <- setdiff(REE, names(cols))
    if (length(nm)) msgs <- c(msgs, sprintf("Not in the file: %s.", paste(nm, collapse = ", ")))
    if (identical(cfg$gap, "interpolate") && any(out$n_estimated > 0, na.rm = TRUE)) {
      msgs <- c(msgs, sprintf("%d sample(s) have estimated Pr/Tb/Tm/Lu (log-linear interpolation of neighbours).",
                              sum(out$n_estimated > 0, na.rm = TRUE)))
    }
    if (identical(cfg$gap, "groupmean") && any(out$n_estimated > 0, na.rm = TRUE)) {
      msgs <- c(msgs, sprintf("%d sample(s) lack some of the elements in the HREE/MREE/LREE groups (group means use the rest).",
                              sum(out$n_estimated > 0, na.rm = TRUE)))
    }
    if (nrow(out) && !any(is.finite(out$REE_rough))) {
      msgs <- c(msgs, "REE roughness could not be calculated (needs Pr, Nd, Sm, Gd, Tb, Dy, Ho, Er, Tm, Yb and Lu measured).")
    }
  }
  extra <- df[, setdiff(names(df), used), drop = FALSE]
  list(data = out, messages = msgs, mapping = mapping, extra = extra, conc = conc, meas = meas)
}
