# MREE/MREE* vs HREE/LREE cross-plot with authigenic and detrital arrays
# Interactive R Shiny app (runs in the browser via shinylive/webR; see README)
library(shiny)
library(ggplot2)
library(MASS)

# R/*.R is sourced automatically by shiny (>= 1.5), in file-name order
ref_all      <- load_reference()
endm         <- load_endmembers()
arrays_lam   <- read.csv("data/arrays_lambda.csv")
ref_sets     <- ref_all[ref_all$ref %in% names(REF_LABELS), ]
example_path <- "data/example_user_data.csv"
radii        <- load_radii()

NORM_CHOICES <- setNames(names(NORMALISER_LABELS), NORMALISER_LABELS)
GD_METHODS <- c("Bau & Dulski (1996): 0.33 Sm + 0.67 Tb" = "BD1996",
                "Geometric Sm–Tb: Sm^(1/3) Tb^(2/3) (Barrat & Bayon 2024)" = "BDgeo",
                "Sm²/Nd (Kulaksız & Bau 2007; Barrat & Bayon 2024)" = "SmNd",
                "Tb²/Dy (Barrat & Bayon 2024)" = "TbDy",
                "Polynomial fit without Gd, Ce, Eu (experimental)" = "BLambdaR")

section <- function(title, ..., open = FALSE) {
  tags$details(class = "sect", open = if (open) NA else NULL, tags$summary(title), div(class = "sect-body", ...))
}

ui <- fluidPage(
  title = "REE cross-plot: MREE/MREE* vs HREE/LREE",
  tags$head(tags$style(HTML("
    body { font-size: 14px; }
    .well { padding: 10px; }
    .small-note { color: #555; font-size: 12px; margin: 2px 0 8px 0; }
    .status-ok { color: #2a6f2a; } .status-warn { color: #9a5b00; } .status-err { color: #a12; font-weight: 600; }
    details.sect { border: 1px solid #d5d9de; border-radius: 5px; margin-bottom: 8px; background: #fff; }
    details.sect > summary { cursor: pointer; font-weight: 600; padding: 7px 10px; background: #f1f4f7;
                             border-radius: 5px; list-style-position: inside; }
    details.sect[open] > summary { border-bottom: 1px solid #d5d9de; border-bottom-left-radius: 0; border-bottom-right-radius: 0; }
    .sect-body { padding: 8px 10px 2px 10px; }
    .sect-body .form-group { margin-bottom: 8px; }
    .sect-body h5 { margin: 10px 0 4px 0; font-weight: 600; font-size: 13px; color: #444; }
    .map-table { font-size: 12px; }
    .dl-row .btn { margin: 0 4px 6px 0; }
    .lam-note { display: none; color: #b3261e; font-size: 12px; margin: -4px 0 8px 0; }
    body.lam .lam-note { display: block; }
    body.lam .lam-off { opacity: 0.45; }
    .spec-opts .form-group { margin-bottom: 6px; }
    .cite-note { background: #f1f4f7; border-left: 4px solid #2166AC; padding: 6px 10px; margin: 0 0 12px 0; font-size: 12.5px; }
  ")),
  # grey out (and annotate) the options that have no effect while lambda space is selected
  tags$script(HTML("$(document).on('shiny:inputchanged', function(e) {
    if (e.name === 'lambda_space') { $('body').toggleClass('lam', !!e.value); } });"))),
  titlePanel("REE cross-plot: MREE/MREE* vs HREE/LREE"),
  div(class = "cite-note",
      HTML("<b>Please cite Blaser et al. (2020)</b> when you use this tool or figures made with it (full reference under \"About &amp; references\"). ",
           "The MREE/MREE* vs HREE/LREE diagram goes back to Martin et al. (2010), Du et al. (2016) added the authigenic and detrital arrays, ",
           "and Blaser et al. (2020) the reference fields of foraminifera, pore water and seawater. Citing the earlier two papers is optional but appreciated.")),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      section("1 · Your data", open = TRUE,
        fileInput("file", NULL, accept = c(".csv", ".tsv", ".txt", ".tab", ".xlsx"),
                  buttonLabel = "Browse…", placeholder = "CSV / TSV / PANGAEA .tab"),
        checkboxInput("use_example", "Use example data (Labrador Sea seawater, Filippova et al. 2017)", FALSE),
        radioButtons("mode", NULL, inline = TRUE,
                     choices = c("REE concentrations" = "concentrations", "Pre-computed ratios" = "ratios")),
        conditionalPanel("input.mode == 'concentrations'",
          radioButtons("basis", "Concentration units are per …", inline = TRUE,
                       choices = c("mass (ppm, ng/g, …)" = "mass", "mole (pmol/kg, …)" = "molar")),
          div(class = "small-note", "A 'basis' column in the file overrides this per row.")),
        conditionalPanel("input.mode == 'ratios'",
          selectInput("xcol", "HREE/LREE column", choices = NULL),
          selectInput("ycol", "MREE/MREE* column", choices = NULL)),
        fluidRow(column(6, selectInput("id_col", "Sample ID column", choices = "(none)")),
                 column(6, selectInput("group_col", "Group column", choices = "(none)"))),
        selectizeInput("id_extra", "Further columns to carry into the results", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "optional")),
        uiOutput("status"),
        uiOutput("source_note"),
        uiOutput("mapping")),
      section("2 · Normalisation & ratios",
        div(class = "lam-off",
          selectInput("norm", "Normalise to", choices = NORM_CHOICES, selected = DEFAULT_NORM),
          conditionalPanel("input.norm == 'custom'",
            fileInput("norm_file", NULL, accept = c(".csv", ".tsv", ".txt"), buttonLabel = "Browse…",
                      placeholder = "one row, La–Lu (and optionally Y), ppm (mass basis)"),
            div(class = "small-note", "Own reference: one data row with the 14 REE (La … Lu, optionally Y) in mass units (ppm).")),
          div(class = "lam-note", "Not used in λ space (λ uses CI chondrite, O'Neill 2016). Still used on the Spectra tab."),
          radioButtons("lree", "LREE in HREE/LREE",
                       choices = c("La – Pr – Nd (Martin et al. 2010)" = "LaPrNd", "Pr – Nd – Sm" = "PrNdSm")),
          div(class = "lam-note", "Not used in λ space."),
          selectInput("gap", "Missing elements (Pr, Tb, Tm, Lu …)",
                      choices = c("Interpolate from neighbours" = "interpolate",
                                  "Mean of measured group members" = "groupmean",
                                  "Not allowed (sample is dropped)" = "none")),
          div(class = "lam-note", "Not used in λ space (only measured elements enter the fit).")),
        selectInput("ce_form", "Ce/Ce* (table only)",
                    choices = c("Geometric from Pr and Nd: Ce Nd / Pr² (Barrat & Bayon 2024)" = "geoPrNd",
                                "Arithmetic La–Pr: 2 Ce / (La + Pr) (Martin et al. 2010)" = "arithLaPr",
                                "Geometric La–Pr: Ce / √(La Pr)" = "geoLaPr"))),
      section("3 · Screening & symbols",
        selectInput("fill_var", "Symbol fill", choices = c("Group (categorical)" = "group")),
        checkboxInput("rough_on", "Exclude samples rougher than …", FALSE),
        sliderInput("rough_max", "REE roughness limit", min = 0, max = 1.5, value = 0.3, step = 0.01),
        div(class = "small-note", "Roughness = visible jaggedness of the REE pattern, a sign of poor data quality (see the Data quality tab). Applies to your data only. Samples without a roughness value (elements missing) are kept."),
        sliderInput("min_ree", "Minimum number of measured REE (of 14)", min = 0, max = 14, value = 0, step = 1),
        selectizeInput("req_ree", "Elements that must be measured", choices = REE, multiple = TRUE,
                       options = list(placeholder = "optional")),
        div(class = "small-note", "Excludes samples with too many gaps. Measured values only (interpolated gaps do not count). In λ space, patterns with fewer than 10 usable elements give no λ at all."),
        uiOutput("groups_ui")),
      section("4 · Plot",
        checkboxInput("lambda_space", "Lambda space (experimental): −λ1 vs −λ2 instead of the ratios", FALSE),
        div(class = "small-note", "λ values after O'Neill (2016), fitted as in Anenburg & Williams (2022). Related web app for exploring and calculating λ values: ",
            tags$a(href = "https://lambdar.rses.anu.edu.au/blambdar/", target = "_blank", rel = "noopener", "BLambdaR"), "."),
        conditionalPanel("!input.lambda_space", checkboxInput("log_axes", "log10 axes (both)", FALSE)),
        checkboxInput("show_legend", "Show legend", TRUE),
        radioButtons("window", "Plot window",
                     choices = c("Default (arrays, fields and the bulk of my data)" = "auto", "Zoom to my data" = "user",
                                 "Only the range of the two arrays" = "arrays")),
        sliderInput("coverage", "Share of data inside field outlines", min = 0.5, max = 0.99, value = 0.95, step = 0.01),
        sliderInput("adjust", "Field smoothness (bandwidth multiplier)", min = 0.5, max = 3, value = 1.0, step = 0.1),
        sliderInput("pt_size", "Point size", min = 1, max = 8, value = 3.2, step = 0.2),
        checkboxInput("user_fields", "Density field for each of my groups (≥ 8 points)", FALSE)),
      section("5 · Reference layers",
        checkboxInput("show_fields", "Density fields (published data)", TRUE),
        checkboxInput("show_refpts", "Individual reference points", TRUE),
        checkboxInput("show_auth", "Authigenic array", TRUE),
        checkboxInput("show_detr", "Detrital array", TRUE),
        checkboxInput("show_labels", "Labels", TRUE),
        checkboxGroupInput("ref_sel", "Reference data included",
                           choices = setNames(names(REF_LABELS), REF_LABELS), selected = names(REF_LABELS))),
      section("6 · Anthropogenic Gd",
        radioButtons("gd", "Gadolinium in the ratios", choices = c("Use measured Gd" = "measured", "Ignore Gd (replace by estimate)" = "ignore")),
        div(class = "lam-off",
          selectInput("gd_method", "Natural Gd (Gd*) from", choices = GD_METHODS),
          numericInput("gd_baseline", "Natural Gd/Gd* baseline (1 = none)", value = 1, min = 0.5, max = 3, step = 0.01),
          div(class = "lam-note", "Not used in λ space (only the choice to ignore Gd, which leaves Gd out of the fit).")),
        div(class = "small-note", "Natural samples typically give 1.1–1.2 with Bau & Dulski. The method and baseline are used for the anthropogenic Gd in the table (Gd/Gd*, Gd_anthropogenic) and, if Gd is ignored, for its replacement."))
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Plot",
          plotOutput("plot", height = "560px", click = "plot_click"),
          div(class = "small-note", "Click a point to see its values."),
          tableOutput("clicked"),
          div(class = "dl-row",
            downloadButton("dl_pdf", "Figure (PDF)"),
            downloadButton("dl_png", "Figure (PNG)"),
            downloadButton("dl_csv", "Results (CSV)"),
            downloadButton("dl_arr", "Arrays & fields (CSV)"))),
        tabPanel("Spectra",
          br(),
          uiOutput("spec_msg"),
          div(class = "spec-opts",
            fluidRow(
              column(4, selectInput("spec_x", "x axis", choices = c("Element" = "element", "Ionic radius (Shannon 1976)" = "radius"))),
              column(4, selectInput("spec_norm2", "Secondary normalisation to …", choices = c("(none)", SPEC_ELEMENTS))),
              column(4, selectInput("spec_color", "Colour by", choices = c("Group (categorical)" = "group")))),
            fluidRow(
              column(8, checkboxInput("spec_y", "Show Y", TRUE), checkboxInput("spec_log", "log10 y axis", TRUE),
                     checkboxInput("spec_connect", "Connect across missing elements", TRUE),
                     checkboxInput("spec_points", "Mark the measured elements", TRUE),
                     checkboxInput("spec_legend", "Show legend", TRUE)),
              column(4, sliderInput("spec_lw", "Line width", min = 0.1, max = 2, value = 0.5, step = 0.1)))),
          plotOutput("spec_plot", height = "480px"),
          div(class = "small-note", "Uses the normalising reference and the unit basis chosen on the left, and the same samples as the Plot tab (hidden groups and screened samples are left out). Measured values only: gaps are not interpolated. Molar data are taken as pmol/kg."),
          div(class = "dl-row",
            downloadButton("dl_spec_pdf", "Spectra (PDF)"),
            downloadButton("dl_spec_csv", "Normalised REE (CSV)"))),
        tabPanel("Data quality",
          br(),
          uiOutput("rough_msg"),
          fluidRow(
            column(4, checkboxInput("rough_log", "Logarithmic axis", TRUE)),
            column(4, checkboxInput("rough_ref", "Show reference data", TRUE)),
            column(4, checkboxInput("rough_legend", "Show legend", TRUE))),
          plotOutput("rough_plot", height = "420px"),
          div(class = "small-note", HTML("Distribution of the roughness of your data (filled, one curve per group; ticks mark the individual samples) compared with the reference data (dashed). ",
              "The dotted line is the limit set under <i>3 · Screening &amp; symbols</i>. Roughness grows with measurement noise, contamination and interferences, but also with falling concentrations: ",
              "seawater is naturally rougher than foraminifera, so compare your data with the reference curve of the same kind of sample. Each curve has unit area; curves need at least 5 different values (smaller groups appear as ticks only).")),
          tableOutput("rough_table"),
          div(class = "dl-row", downloadButton("dl_rough_pdf", "Roughness plot (PDF)"))),
        tabPanel("Computed values", br(), tableOutput("table")),
        tabPanel("About & references", br(), about_ui())
      )
    )
  )
)

server <- function(input, output, session) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # ---- calculation options ---------------------------------------------------------------
  norm_custom <- reactive({
    if (!identical(input$norm, "custom") || is.null(input$norm_file)) return(NULL)
    tryCatch(read_norm_file(input$norm_file$datapath, input$norm_file$name), error = function(e) e)
  })
  norm_msg <- reactive({
    if (!identical(input$norm, "custom")) return(NULL)
    v <- norm_custom()
    if (is.null(v)) "Own reference: no file loaded yet - PAAS (Pourmand) is used meanwhile."
    else if (inherits(v, "error")) paste0(conditionMessage(v), " PAAS (Pourmand) is used meanwhile.") else NULL
  })
  cfg <- reactive({
    v <- norm_custom()
    use_custom <- identical(input$norm, "custom") && is.numeric(v)
    default_cfg(norm = if (identical(input$norm, "custom") && !use_custom) DEFAULT_NORM else (input$norm %||% DEFAULT_NORM),
                norm_custom = if (use_custom) v else NULL,
                lree = input$lree %||% "LaPrNd", gap = input$gap %||% "interpolate",
                gd = input$gd %||% "measured", gd_method = input$gd_method %||% "BD1996",
                gd_baseline = if (is.finite(input$gd_baseline %||% 1)) input$gd_baseline %||% 1 else 1,
                ce_form = input$ce_form %||% "geoPrNd")
  })
  # reference data are always computed with measured Gd and interpolated gaps; only the shale
  # standard and the LREE definition follow the user's choice
  ref_cfg <- reactive({ c <- cfg(); c$gd <- "measured"; c$gd_method <- "BD1996"; c$gap <- "interpolate"; c })
  lambda_on <- reactive(isTRUE(input$lambda_space))
  space <- reactive(if (lambda_on()) "lambda" else "ratio")

  # ---- user data ---------------------------------------------------------------------------
  raw_df <- reactive({
    if (isTRUE(input$use_example)) return(read_user_file(example_path))
    req(input$file)
    read_user_file(input$file$datapath, input$file$name)
  })

  observeEvent(input$file, updateCheckboxInput(session, "use_example", value = FALSE))

  observeEvent(raw_df(), {
    df <- tryCatch(raw_df(), error = function(e) NULL)
    req(df)
    nms <- names(df)
    id_guess <- c("sample", "Sample label", "Sample code/label", "Label", "Sample ID", "Event")
    id_guess <- id_guess[id_guess %in% nms]
    updateSelectInput(session, "id_col", choices = c("(none)", nms),
                      selected = if (length(id_guess)) id_guess[1] else "(none)")
    # unit basis from the column names (e.g. "La [mg/kg]"), unless a 'basis' column decides per row
    if (!"basis" %in% tolower(nms)) {
      ub <- guess_basis_from_units(detect_ree_columns(nms))
      if (ub$basis %in% c("mass", "molar")) updateRadioButtons(session, "basis", selected = ub$basis)
    }
    nonree <- setdiff(nms, c(detect_ree_columns(nms), detect_y_column(nms)))
    updateSelectInput(session, "group_col", choices = c("(none)", nonree),
                      selected = if ("group" %in% nms) "group" else "(none)")
    updateSelectizeInput(session, "id_extra", choices = nonree, selected = character(0))
    g <- guess_ratio_columns(nms)
    updateSelectInput(session, "xcol", choices = nms, selected = g$x)
    updateSelectInput(session, "ycol", choices = nms, selected = g$y)
    enough <- length(check_ree_columns(names(detect_ree_columns(nms)), isolate(cfg()))) == 0
    updateRadioButtons(session, "mode", selected = if (enough || is.na(g$x)) "concentrations" else "ratios")
  })

  prepped <- reactive({
    df <- raw_df()
    none <- function(v) if (is.null(v) || identical(v, "(none)")) NULL else v
    tryCatch(
      prepare_user_data(df, input$mode, input$basis, cfg(), none(input$group_col), none(input$id_col),
                        input$xcol, input$ycol, lambda = lambda_on()),
      error = function(e) list(error = conditionMessage(e)))
  })

  user_all <- reactive({
    p <- tryCatch(prepped(), error = function(e) NULL)
    if (is.null(p) || !is.null(p$error)) return(NULL)
    p$data
  })
  user_data <- user_all   # (kept for the tests)

  all_groups <- reactive({ u <- user_all(); if (is.null(u)) character() else natural_sort(u$group) })

  # group visibility
  output$groups_ui <- renderUI({
    g <- all_groups()
    if (length(g) < 2) return(NULL)
    tagList(
      h5("Groups shown"),
      div(class = "small-note", actionLink("grp_all", "all"), " · ", actionLink("grp_none", "none"),
          " · tick a single box to show one group"),
      checkboxGroupInput("groups_shown", NULL, choices = g, selected = isolate(input$groups_shown) %||% g))
  })
  observeEvent(all_groups(), {
    g <- all_groups()
    if (length(g) >= 2) updateCheckboxGroupInput(session, "groups_shown", choices = g, selected = g)
  })
  observeEvent(input$grp_all, updateCheckboxGroupInput(session, "groups_shown", selected = all_groups()))
  observeEvent(input$grp_none, updateCheckboxGroupInput(session, "groups_shown", selected = character(0)))
  groups_shown <- reactive({
    g <- all_groups()
    if (length(g) < 2) g else (input$groups_shown %||% g)
  })

  # colour/fill choices follow the data and the plotting space
  colour_choices <- reactive({
    ch <- c("Group (categorical)" = "group")
    p <- tryCatch(prepped(), error = function(e) NULL)
    if (is.null(p) || !is.null(p$error)) return(ch)
    rough_ok <- if (lambda_on()) "lambda_rough" %in% names(p$data) && any(is.finite(p$data$lambda_rough)) else
      any(is.finite(p$data$REE_rough))
    if (rough_ok) ch <- c(ch, setNames("rough", if (lambda_on()) "Lambda roughness" else "REE roughness"))
    if ("Gd_Gdstar" %in% names(p$data) && any(is.finite(p$data$Gd_Gdstar))) ch <- c(ch, "Gd/Gd*" = "gdstar")
    ex <- p$extra
    for (nm in names(ex)) {
      v <- to_number(ex[[nm]])
      if (mean(is.finite(v)) >= 0.5 && length(unique(v[is.finite(v)])) > 2) ch <- c(ch, setNames(paste0("x:", nm), nm))
    }
    ch
  })
  observeEvent(colour_choices(), {
    ch <- colour_choices()
    updateSelectInput(session, "fill_var", choices = ch, selected = isolate(input$fill_var))
    updateSelectInput(session, "spec_color", choices = ch, selected = isolate(input$spec_color))
  })
  observeEvent(lambda_on(), {
    if (lambda_on()) {
      updateSliderInput(session, "rough_max", label = "Lambda roughness limit", min = 0, max = 0.1, value = 0.03, step = 0.001)
      updateSliderInput(session, "min_ree", value = 12)    # the 9-parameter fit needs well-measured patterns
      updateRadioButtons(session, "window", selected = "auto")
    } else {
      updateSliderInput(session, "rough_max", label = "REE roughness limit", min = 0, max = 1.5, value = 0.3, step = 0.01)
      updateSliderInput(session, "min_ree", value = 0)
    }
  }, ignoreInit = TRUE)

  rough_col <- reactive(if (lambda_on()) "lambda_rough" else "REE_rough")

  # samples with too few measured REE (NULL if the screening is off or no concentrations were given)
  too_few <- reactive({
    p <- tryCatch(prepped(), error = function(e) NULL)
    if (is.null(p) || !is.null(p$error) || is.null(p$meas)) return(NULL)
    mz <- input$min_ree %||% 0; rq <- input$req_ree
    if (mz <= 0 && length(rq) == 0) return(NULL)
    bad <- rowSums(p$meas[, REE, drop = FALSE]) < mz
    if (length(rq)) bad <- bad | !apply(p$meas[, rq, drop = FALSE], 1, all)
    bad
  })

  # which samples are drawn, and why the others are not
  plot_status <- reactive({
    u <- user_all(); req(u)
    uc <- .space_cols(space())
    st <- rep("plotted", nrow(u))
    st[!is.finite(u[[uc[1]]]) | !is.finite(u[[uc[2]]])] <- "not plottable (missing values)"
    thin <- too_few()
    if (!is.null(thin)) st[thin] <- "excluded (too few REE)"
    rc <- rough_col()
    if (isTRUE(input$rough_on) && rc %in% names(u)) {
      st[st == "plotted" & is.finite(u[[rc]]) & u[[rc]] > (input$rough_max %||% Inf)] <- "excluded (roughness)"
    }
    st[st == "plotted" & !(u$group %in% groups_shown())] <- "hidden (group)"
    st
  })

  colour_values <- function(choice) {
    u <- user_all(); req(u)
    fv <- choice %||% "group"
    if (identical(fv, "rough")) list(v = u[[rough_col()]], label = if (lambda_on()) "lambda roughness" else "REE roughness")
    else if (identical(fv, "gdstar")) list(v = u$Gd_Gdstar, label = "Gd/Gd*")
    else if (startsWith(fv, "x:")) {
      nm <- sub("^x:", "", fv); ex <- tryCatch(prepped()$extra, error = function(e) NULL)
      if (!is.null(ex) && nm %in% names(ex)) list(v = to_number(ex[[nm]]), label = nm) else NULL
    } else NULL
  }
  fill_values <- reactive(colour_values(input$fill_var))

  user_plot <- reactive({
    u <- user_all()
    if (is.null(u)) return(NULL)
    fl <- fill_values()
    if (!is.null(fl)) u$.fill <- fl$v
    u[plot_status() == "plotted", , drop = FALSE]
  })

  output$status <- renderUI({
    if (!isTRUE(input$use_example) && is.null(input$file)) {
      return(tagList(div(class = "small-note", "No data loaded – showing the reference plot only."),
                     lapply(norm_msg(), function(m) div(class = "status-warn", m))))
    }
    p <- tryCatch(prepped(), error = function(e) list(error = conditionMessage(e)))
    if (!is.null(p$error)) return(div(class = "status-err", p$error))
    st <- plot_status()
    n <- length(st)
    lines <- list(div(class = "status-ok", sprintf("%d of %d samples plotted.", sum(st == "plotted"), n)))
    cnt <- function(pat, txt) if (any(st == pat)) div(class = "status-warn", sprintf("%d %s", sum(st == pat), txt))
    lines <- c(lines, list(cnt("not plottable (missing values)", "sample(s) could not be plotted (missing values)."),
                           cnt("excluded (too few REE)", "sample(s) excluded (too few measured REE)."),
                           cnt("excluded (roughness)", "sample(s) excluded by the roughness limit."),
                           cnt("hidden (group)", "sample(s) in hidden groups.")))
    if (sum(st == "plotted") > 0) {
      u <- user_plot(); l <- attr(fig(), "limits"); uc <- .space_cols(space())
      if (!is.null(l)) {
        out <- sum(u[[uc[1]]] < l$xlim[1] | u[[uc[1]]] > l$xlim[2] | u[[uc[2]]] < l$ylim[1] | u[[uc[2]]] > l$ylim[2])
        if (out) lines <- c(lines, list(div(class = "status-warn", sprintf("%d plotted sample(s) lie outside the plot window.", out))))
      }
    }
    tagList(lines, lapply(norm_msg(), function(m) div(class = "status-warn", m)),
            lapply(p$messages, function(m) div(class = "status-warn", m)))
  })

  # provenance and unit messages for the loaded file
  source_info <- reactive({
    df <- tryCatch(raw_df(), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    src <- attr(df, "source")
    nms <- names(df)
    ub <- if (!"basis" %in% tolower(nms)) guess_basis_from_units(detect_ree_columns(nms)) else NULL
    list(src = src, ub = ub)
  })
  output$source_note <- renderUI({
    si <- source_info(); if (is.null(si) || identical(input$mode, "ratios")) return(NULL)
    out <- list()
    if (!is.null(si$src)) {
      doi <- si$src$doi
      out <- c(out, list(div(class = "small-note", tags$b("PANGAEA dataset: "), si$src$citation,
        if (!is.na(doi)) tagList(" (", tags$a(href = paste0("https://doi.org/", doi), target = "_blank", "doi"), ")"),
        if (!is.na(si$src$licence) && grepl("CC-BY|Attribution", si$src$licence))
          " Licence CC-BY: please cite the dataset.")))
    }
    ub <- si$ub
    if (!is.null(ub) && length(ub$units)) {
      if (ub$basis %in% c("mass", "molar")) {
        out <- c(out, list(div(class = "status-ok", sprintf("Units read from the column names: %s (%s basis selected).",
                                                           paste(ub$units, collapse = ", "), ub$basis))))
        if (ub$basis == "molar" && !any(grepl("pmol", ub$units)))
          out <- c(out, list(div(class = "status-warn", "The Spectra tab takes molar data as pmol/kg: its absolute values are scaled accordingly; ratios are unaffected.")))
      } else if (identical(ub$basis, "mixed")) {
        out <- c(out, list(div(class = "status-warn", sprintf("The REE columns have different units (%s): check them, or add a 'basis' column.",
                                                             paste(ub$units, collapse = ", ")))))
      }
    }
    tagList(out)
  })

  output$mapping <- renderUI({
    p <- tryCatch(prepped(), error = function(e) NULL)
    if (is.null(p) || !is.null(p$error) || nrow(p$mapping) == 0) return(NULL)
    m <- p$mapping
    tags$details(tags$summary(sprintf("Columns used (%d)", nrow(m))),
      tags$table(class = "map-table table table-condensed",
        lapply(seq_len(nrow(m)), function(i) tags$tr(tags$td(tags$b(m$element[i])), tags$td("←"), tags$td(m$column[i])))))
  })

  # ---- reference data, arrays ------------------------------------------------------------------
  ref_used <- reactive({
    r <- ref_sets[ref_sets$ref %in% input$ref_sel, ]
    reference_ratios(r, ref_cfg())
  })
  arrays <- reactive({
    a <- list(authigenic = authigenic_array(endm, cfg()), detrital = detrital_array(endm, cfg()))
    attach_lambda(a, arrays_lam)
  })

  plot_opts <- reactive({
    list(space = space(), log = isTRUE(input$log_axes) && !lambda_on(), window = input$window %||% "auto",
         fields = input$show_fields, ref_points = input$show_refpts,
         authigenic = input$show_auth, detrital = input$show_detr,
         labels = input$show_labels, coverage = input$coverage, adjust = input$adjust,
         point_size = input$pt_size, user_fields = input$user_fields,
         group_levels = all_groups(), legend = isTRUE(input$show_legend),
         fill_label = tryCatch(fill_values()$label, error = function(e) NULL))
  })
  fig <- reactive(plot_martin(ref_used(), arrays(), user_plot(), plot_opts()))
  output$plot <- renderPlot(fig(), res = 110)

  output$clicked <- renderTable({
    u <- user_plot(); req(u, input$plot_click, nrow(u) > 0)
    uc <- .space_cols(space())
    u$X <- if (isTRUE(input$log_axes) && !lambda_on()) log10(u[[uc[1]]]) else u[[uc[1]]]
    u$Y <- if (isTRUE(input$log_axes) && !lambda_on()) log10(u[[uc[2]]]) else u[[uc[2]]]
    np <- nearPoints(u, input$plot_click, xvar = "X", yvar = "Y", threshold = 12, maxpoints = 3)
    req(nrow(np) > 0)
    out <- data.frame(sample = np$sample, group = np$group, a = round(np[[uc[1]]], 3), b = round(np[[uc[2]]], 3),
                      check.names = FALSE)
    names(out)[3:4] <- if (lambda_on()) c("-lambda1", "-lambda2") else c("HREE/LREE", "MREE/MREE*")
    out
  })

  settings_cols <- reactive({
    cf <- cfg()
    data.frame(normalisation = norm_label(cf),
               LREE_definition = c(LaPrNd = "La-Pr-Nd", PrNdSm = "Pr-Nd-Sm")[[cf$lree]],
               gap_method = c(interpolate = "interpolate", groupmean = "group mean", none = "none")[[cf$gap]],
               Gd_handling = if (identical(cf$gd, "ignore")) "ignored (replaced by estimate)" else "measured",
               Gd_star_method = unname(names(GD_METHODS)[match(cf$gd_method, GD_METHODS)]),
               Gd_baseline = cf$gd_baseline,
               Ce_anomaly = c(geoPrNd = "Ce Nd / Pr^2", arithLaPr = "2 Ce / (La + Pr)", geoLaPr = "Ce / sqrt(La Pr)")[[cf$ce_form]],
               plot_space = space(), stringsAsFactors = FALSE)
  })

  computed <- reactive({
    p <- prepped(); u <- user_all(); req(u)
    out <- u
    out$plot_status <- plot_status()
    out$lx <- NULL; out$ly <- NULL       # plotted coordinates -lambda1/-lambda2; lambda1/lambda2 are in the table
    num <- vapply(out, is.numeric, logical(1))
    num[intersect(names(num), c("row", "n_estimated", "n_REE_measured"))] <- FALSE
    out[num] <- lapply(out[num], function(v) signif(v, 5))
    nm <- c(x = "HREE/LREE", y = "MREE/MREE*", Ce_Ceanom = "Ce/Ce*", La_Laanom = "La/La*", Eu_Euanom = "Eu/Eu*",
            Gd_Gdstar = "Gd/Gd*", Gd_anth = "Gd_anthropogenic (input units)", REE_rough = "REE_roughness")
    names(out)[match(names(nm), names(out), nomatch = 0)] <- nm[names(nm) %in% names(out)]
    ids <- c("row", "sample", "group")
    ex <- p$extra[, intersect(input$id_extra, names(p$extra)), drop = FALSE]
    ex <- ex[, setdiff(names(ex), c(ids, names(out))), drop = FALSE]
    res <- cbind(out[ids], ex, out[setdiff(names(out), ids)], settings_cols(), row.names = NULL)
    src <- tryCatch(attr(raw_df(), "source"), error = function(e) NULL)
    if (!is.null(src)) { res$source_dataset <- src$citation; res$source_doi <- src$doi }
    res
  })
  output$table <- renderTable(computed(), digits = 3)

  # arrays and fields of the current configuration, in one long table (lambda1/lambda2 are given
  # as lambda values, i.e. without the sign used for the plot axes)
  arrays_fields <- reactive({
    a <- arrays(); cf <- cfg(); sp <- space()
    meta <- data.frame(normalisation = norm_label(cf),
                       LREE_definition = c(LaPrNd = "La-Pr-Nd", PrNdSm = "Pr-Nd-Sm")[[cf$lree]],
                       plot_space = sp, stringsAsFactors = FALSE)
    row <- function(layer, subset, piece, r, hl, mm, l1, l2)
      data.frame(layer = layer, subset = subset, piece = piece, r = r, `HREE/LREE` = hl, `MREE/MREE*` = mm,
                 lambda1 = l1, lambda2 = l2, meta, row.names = NULL, check.names = FALSE)
    au <- a$authigenic; cv <- a$detrital$curves; dp <- a$detrital$polygon
    parts <- list(
      row("authigenic_array", "seawater-nodule mixing", NA, au$r, au$HREE_LREE, au$MREE_MREEstar, -au$lx, -au$ly),
      row("detrital_curve", cv$curve, NA, cv$r, cv$HREE_LREE, cv$MREE_MREEstar, -cv$lx, -cv$ly),
      row("detrital_polygon", "GAC-MORB-CPS", NA, dp$r, dp$HREE_LREE, dp$MREE_MREEstar, -dp$lx, -dp$ly))
    ru <- ref_used(); uc <- .space_cols(sp)
    for (tp in names(FIELD_COLOURS)) {
      s <- ru[ru$type == tp & is.finite(ru[[uc[1]]]) & is.finite(ru[[uc[2]]]), ]
      poly <- hdr_polygons(s[[uc[1]]], s[[uc[2]]], input$coverage, input$adjust)
      if (is.null(poly)) next
      parts[[length(parts) + 1]] <- if (identical(sp, "lambda"))
        row("density_field", FIELD_NAMES[[tp]], poly$piece, NA, NA, NA, -poly$x, -poly$y) else
        row("density_field", FIELD_NAMES[[tp]], poly$piece, NA, poly$x, poly$y, NA, NA)
    }
    out <- do.call(rbind, parts)
    out$field_coverage <- ifelse(out$layer == "density_field", input$coverage, NA)
    out$field_bandwidth_multiplier <- ifelse(out$layer == "density_field", input$adjust, NA)
    out
  })

  # ---- data quality tab --------------------------------------------------------------------------
  rough_label <- reactive(if (lambda_on()) "Lambda roughness" else "REE roughness")
  rough_ref <- reactive({
    if (!isTRUE(input$rough_ref)) return(NULL)
    ru <- ref_used(); rc <- rough_col()
    if (!rc %in% names(ru) || nrow(ru) == 0) return(NULL)
    d <- data.frame(class = ref_class(ru), v = ru[[rc]])
    d[is.finite(d$v), , drop = FALSE]
  })
  rough_user <- reactive({
    u <- user_all(); rc <- rough_col()
    if (is.null(u) || !rc %in% names(u)) return(NULL)
    ok <- u$group %in% groups_shown() & is.finite(u[[rc]])
    if (!any(ok)) return(NULL)
    data.frame(group = u$group[ok], v = u[[rc]][ok], stringsAsFactors = FALSE)
  })
  rough_fig <- reactive({
    plot_roughness(rough_ref(), rough_user(),
                   list(log = isTRUE(input$rough_log), limit = input$rough_max %||% NA_real_, limit_on = isTRUE(input$rough_on),
                        label = rough_label(), group_levels = all_groups(), legend = isTRUE(input$rough_legend)))
  })
  output$rough_plot <- renderPlot({
    validate(need(!is.null(rough_ref()) || !is.null(rough_user()), "Load REE concentrations (or tick the reference data) to see the roughness distribution."))
    suppressWarnings(rough_fig())
  }, res = 110)
  output$rough_table <- renderTable(roughness_summary(rough_ref(), rough_user(), input$rough_max %||% NA_real_), digits = 3, na = "")
  output$rough_msg <- renderUI({
    if (identical(input$mode, "ratios") && (isTRUE(input$use_example) || !is.null(input$file)))
      return(div(class = "status-warn", "Pre-computed ratios carry no roughness: only the reference data are shown."))
    u <- user_all()
    if (is.null(u) || isTRUE(identical(input$mode, "ratios"))) return(NULL)
    rc <- rough_col()
    if (!rc %in% names(u) || !any(is.finite(u[[rc]])))
      return(div(class = "status-warn", "No roughness could be calculated for your data (all of Pr, Nd, Sm, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu must be measured; in λ space at least 10 usable elements)."))
    n0 <- sum(!is.finite(u[[rc]]))
    if (n0) div(class = "status-warn", sprintf("%d sample(s) without roughness value (elements missing) are not shown.", n0))
  })
  output$dl_rough_pdf <- downloadHandler("ree_roughness.pdf", function(file)
    suppressWarnings(ggplot2::ggsave(file, rough_fig(), width = 9, height = 5.5, device = grDevices::pdf)))

  # ---- spectra tab -----------------------------------------------------------------------------
  spec_ref_short <- reactive(sub(" - .*$", "", norm_label(cfg())))
  spec_norm2 <- reactive(if (identical(input$spec_norm2, "(none)") || is.null(input$spec_norm2)) NULL else input$spec_norm2)

  spec_all <- reactive({              # all rows of the upload, normalised
    p <- prepped(); req(is.null(p$error), !is.null(p$conc))
    ree_spectra(p$conc, p$data$basis, cfg(), spec_norm2())
  })
  spec_values <- reactive({ sp <- spec_all(); if (is.null(sp$n2)) sp$n else sp$n2 })
  spec_keep <- reactive(!(plot_status() %in% c("hidden (group)", "excluded (roughness)", "excluded (too few REE)")))

  spec_long <- reactive({
    p <- prepped(); u <- user_all(); m <- spec_values(); keep <- spec_keep()
    cv <- colour_values(input$spec_color)
    meta <- data.frame(sid = u$row, row = u$row, sample = u$sample, group = u$group, stringsAsFactors = FALSE)
    if (!is.null(cv)) meta$.col <- cv$v
    spectra_long(m[keep, , drop = FALSE], meta[keep, , drop = FALSE], radii)
  })
  spec_fig <- reactive({
    cv <- colour_values(input$spec_color)
    ref <- spec_ref_short(); n2 <- spec_norm2()
    ylab <- if (is.null(n2)) sprintf("REE / %s", ref) else sprintf("(REE / %s) / (%s / %s)", ref, n2, ref)
    plot_spectra(spec_long(), list(x = input$spec_x %||% "element", log = isTRUE(input$spec_log), show_y = isTRUE(input$spec_y),
                                   connect = isTRUE(input$spec_connect), points = isTRUE(input$spec_points),
                                   legend = isTRUE(input$spec_legend), line_size = input$spec_lw %||% 0.5, ylab = ylab,
                                   colour_label = if (is.null(cv)) NULL else cv$label, group_levels = all_groups()))
  })
  output$spec_plot <- renderPlot({
    p <- tryCatch(prepped(), error = function(e) NULL)
    validate(need(!is.null(p) && is.null(p$error) && !is.null(p$conc),
                  if (identical(input$mode, "ratios")) "The spectra need REE concentrations, not pre-computed ratios." else
                    "Load REE concentrations to see their patterns."))
    suppressWarnings(spec_fig())
  }, res = 110)
  output$spec_msg <- renderUI({
    p <- tryCatch(prepped(), error = function(e) NULL)
    if (is.null(p) || !is.null(p$error) || is.null(p$conc)) return(NULL)
    msgs <- character()
    if (isTRUE(input$spec_y) && !is.finite(norm_y(cfg())))
      msgs <- c(msgs, "The chosen reference has no Y value: Y is not shown.")
    else if (isTRUE(input$spec_y) && !any(p$meas[, "Y"]))
      msgs <- c(msgs, "No Y column found in the file (a column named Y, Y_ppm, ...).")
    lapply(msgs, function(m) div(class = "status-warn", m))
  })

  spec_wide <- reactive({             # normalised REE of every row of the upload
    p <- prepped(); u <- user_all(); sp <- spec_all(); n2 <- spec_norm2()
    out <- data.frame(row = u$row, sample = u$sample, group = u$group, plot_status = plot_status(), basis = u$basis,
                      normalisation = norm_label(cfg()), secondary_normalisation = if (is.null(n2)) "none" else n2,
                      stringsAsFactors = FALSE)
    cols <- if (any(is.finite(sp$n[, "Y"]))) SPEC_ELEMENTS else REE
    nn <- signif(sp$n[, cols, drop = FALSE], 5); colnames(nn) <- paste0("norm_", cols)
    out <- cbind(out, nn)
    if (!is.null(sp$n2)) { n2m <- signif(sp$n2[, cols, drop = FALSE], 5); colnames(n2m) <- paste0("norm2_", cols); out <- cbind(out, n2m) }
    out
  })

  output$dl_csv <- downloadHandler("ree_results.csv", function(file) utils::write.csv(computed(), file, row.names = FALSE))
  output$dl_arr <- downloadHandler("ree_arrays_and_fields.csv", function(file)
    utils::write.csv(arrays_fields(), file, row.names = FALSE, na = ""))
  output$dl_spec_csv <- downloadHandler("ree_normalised.csv", function(file)
    utils::write.csv(spec_wide(), file, row.names = FALSE, na = ""))
  output$dl_spec_pdf <- downloadHandler("ree_spectra.pdf", function(file)
    suppressWarnings(ggplot2::ggsave(file, spec_fig(), width = 9, height = 6, device = grDevices::pdf)))
  output$dl_pdf <- downloadHandler("ree_crossplot.pdf", function(file) {
    ggplot2::ggsave(file, fig(), width = 9, height = 6.6, device = grDevices::pdf)
  })
  output$dl_png <- downloadHandler("ree_crossplot.png", function(file) {
    ggplot2::ggsave(file, fig(), width = 9, height = 6.6, dpi = 300, device = "png")
  })
}

shinyApp(ui, server)
