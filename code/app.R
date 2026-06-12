# ============================================================================
# auto-indexwc results browser
# ----------------------------------------------------------------------------
# Browses the `output/` tree on the `autogen-results` branch: the coastwide
# index, single diagnostics (mesh / QQ / anisotropy / fixed effects), density
# plots by year, and residual plots by year. One model run at a time.
#
# NOTE: the .rdata / .rds model objects are intentionally NOT on that branch
# (they live in the upstream R-run-models artifacts), so this app shows the
# PNGs + CSV only.
#
# Run:
#   install.packages(c("shiny", "jsonlite"))   # httr only for github + PAT
#   shiny::runApp("app.R")
#
# Data source (set DATA_SOURCE below):
#   "github" - read live from the branch (default; nothing to download).
#              Optionally set GITHUB_PAT to raise the API rate limit.
#   "local"  - read from a checked-out copy; set LOCAL_ROOT to the directory
#              that CONTAINS the `output/` folder. No network used.
# ============================================================================

library(shiny)
library(jsonlite)

# Show real error messages in the UI. Many hosts sanitize Shiny errors to a
# generic "An error has occurred" and don't expose logs, which hides the cause.
options(shiny.sanitize.errors = FALSE)

# ---- Config ----------------------------------------------------------------
DATA_SOURCE <- "github"                       # "github" or "local"
REPO        <- "pfmc-assessments/auto-indexwc"
BRANCH      <- "autogen-results"
LOCAL_ROOT  <- "."                            # used only when DATA_SOURCE == "local"

SECTIONS <- c("index", "diagnostics", "data")

# ---- Path / source helpers -------------------------------------------------
raw_url <- function(path) {
  sprintf("https://raw.githubusercontent.com/%s/%s/%s", REPO, BRANCH, path)
}

img_src <- function(path) {
  if (DATA_SOURCE == "github") raw_url(path) else file.path("localdata", sub("^output/", "", path))
}

read_csv_any <- function(path) {
  src <- if (DATA_SOURCE == "github") raw_url(path) else file.path(LOCAL_ROOT, path)
  tryCatch(utils::read.csv(src, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) NULL)
}

# ---- Discovery -------------------------------------------------------------
list_paths_github <- function() {
  api <- sprintf("https://api.github.com/repos/%s/git/trees/%s?recursive=1", REPO, BRANCH)
  pat <- Sys.getenv("GITHUB_PAT", "")
  if (nzchar(pat) && requireNamespace("httr", quietly = TRUE)) {
    resp   <- httr::GET(api, httr::add_headers(Authorization = paste("token", pat)))
    parsed <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"),
                                 simplifyDataFrame = TRUE)
  } else {
    parsed <- jsonlite::fromJSON(api, simplifyDataFrame = TRUE)
  }
  if (is.null(parsed$tree)) {
    stop("GitHub API returned no tree. ",
         if (!is.null(parsed$message)) parsed$message else "",
         " (Set GITHUB_PAT to avoid rate limits, or use DATA_SOURCE = 'local'.)")
  }
  if (isTRUE(parsed$truncated))
    warning("GitHub tree truncated; some runs may be missing. Use 'local' mode for completeness.")
  parsed$tree$path[parsed$tree$type == "blob"]
}

list_paths_local <- function() {
  root <- file.path(LOCAL_ROOT, "output")
  if (!dir.exists(root)) stop("LOCAL_ROOT/output not found at: ", root)
  file.path("output", list.files(root, recursive = TRUE))
}

build_manifest <- function() {
  paths <- if (DATA_SOURCE == "github") list_paths_github() else list_paths_local()
  paths <- paths[grepl("\\.(png|csv|txt)$", paths, ignore.case = TRUE)]
  
  parts <- strsplit(paths, "/")
  ok <- vapply(parts, function(p) length(p) >= 7 && p[1] == "output" && p[6] %in% SECTIONS,
               logical(1))
  paths <- paths[ok]; parts <- parts[ok]
  if (length(paths) == 0) stop("No browsable files found under output/ on this branch.")
  
  species  <- vapply(parts, `[`, "", 2)
  survey   <- vapply(parts, `[`, "", 3)
  latyears <- vapply(parts, `[`, "", 4)
  family   <- vapply(parts, `[`, "", 5)
  section  <- vapply(parts, `[`, "", 6)
  fname    <- vapply(parts, function(p) p[length(p)], "")
  
  run_id <- paste(species, survey, latyears, family, sep = "|")
  df <- data.frame(run_id, species, survey, latyears, family, section, fname,
                   path = paths, stringsAsFactors = FALSE)
  split(df, df$run_id)
}

# ---- Prettifiers -----------------------------------------------------------
# Species: drop underscores, capitalize first letter only (sentence case).
pretty_species <- function(x) {
  s <- gsub("_", " ", x)
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

pretty_family <- function(x) {
  map <- c(delta_gamma = "Delta-gamma", delta_lognormal = "Delta-lognormal",
           delta_gengamma = "Delta-generalized-gamma", tweedie = "Tweedie")
  ifelse(x %in% names(map), map[x], gsub("_", " ", x))
}

pretty_lat <- function(x) {
  m <- regmatches(x, regexec("lat_([0-9.]+)-([0-9.]+)", x))[[1]]
  if (length(m) == 3) sprintf("%s\u2013%s\u00b0N", m[2], m[3]) else x
}

lat_min <- function(x) {
  m <- regmatches(x, regexec("lat_([0-9.]+)-", x))[[1]]
  if (length(m) == 2) as.numeric(m[2]) else NA_real_
}

# Density/residual plots are paginated 2 years per page (facet_wrap_paginate
# ncol = 2). Map a 1-indexed page to its year range using the run's ACTUAL
# modeled years -- WCGBTS skipped 2020, so don't assume a gapless 2003 start.
page_year_label <- function(years, page, per_page = 2) {
  if (is.null(years) || !length(years)) return(sprintf("Page %02d", page))
  idx <- seq((page - 1) * per_page + 1, page * per_page)
  idx <- idx[idx >= 1 & idx <= length(years)]
  yy <- years[idx]
  if (!length(yy)) return(sprintf("Page %02d", page))
  if (length(yy) == 1) as.character(yy) else paste0(min(yy), "\u2013", max(yy))
}

# Delta families: component 1 = presence-absence, 2 = positive catch.
resid_label <- function(component, family) {
  if (!grepl("delta", family, ignore.case = TRUE)) return("Residuals")
  switch(as.character(component),
         "1" = "Presence-absence",
         "2" = "Positive catch",
         paste0("Component ", component))
}

# ---- Build manifest once at startup ----------------------------------------
# Don't let a failed fetch (e.g. GitHub API rate limit) crash app startup --
# capture the message and show it in the UI instead.
MANIFEST <- tryCatch(build_manifest(),
                     error = function(e) structure(list(), startup_error = conditionMessage(e)))
STARTUP_ERROR <- attr(MANIFEST, "startup_error")
if (DATA_SOURCE == "local") shiny::addResourcePath("localdata", file.path(LOCAL_ROOT, "output"))

species_raw <- if (length(MANIFEST))
  sort(unique(vapply(MANIFEST, function(d) d$species[1], ""))) else character(0)
species_choices <- if (length(species_raw))
  stats::setNames(species_raw, vapply(species_raw, pretty_species, "")) else character(0)

startup_notice <- function() {
  div(style = "color:#b3261e; padding:1rem; border:1px solid #f2b8b5; border-radius:8px;",
      h4("Could not load results"),
      tags$p(STARTUP_ERROR),
      tags$p(tags$b("Likely cause: "), "the GitHub API rate limit (60/hr unauthenticated). ",
             "Set a ", tags$code("GITHUB_PAT"), " environment variable in the deployment ",
             "(a read-only token is enough for a public repo) to raise it to 5,000/hr, ",
             "or switch ", tags$code("DATA_SOURCE"), " to ", tags$code("\"local\""), "."))
}

# ---- UI --------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Inter', system-ui, sans-serif; }
    .run-badge-fail { background:#b3261e; color:#fff; padding:2px 8px; border-radius:10px; font-size:12px; }
    .run-badge-ok   { background:#1b6e3c; color:#fff; padding:2px 8px; border-radius:10px; font-size:12px; }
    .diag-img, .index-img { max-width:100%; border:1px solid #e0e0e0; border-radius:6px; }
    .half-img { max-width:37.5%; border:1px solid #e0e0e0; border-radius:6px; }
    .section-empty { color:#777; font-style:italic; padding:1rem 0; }
    table { font-size: 13px; }
  "))),
  titlePanel("auto-indexwc results browser"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("species", "Species", choices = species_choices),
      uiOutput("lat_picker"),
      uiOutput("family_picker"),
      uiOutput("converge_badge"),
      tags$hr(),
      tags$small(sprintf("Source: %s \u00b7 %s @ %s", DATA_SOURCE, REPO, BRANCH)),
      tags$br(),
      tags$small(sprintf("%d runs across %d species", length(MANIFEST), length(species_raw)))
    ),
    mainPanel(
      width = 9,
      if (!is.null(STARTUP_ERROR)) startup_notice() else tabsetPanel(
        id = "tabs",
        tabPanel("Coastwide Index",
                 br(),
                 uiOutput("index_status"),
                 uiOutput("index_png_cw")),
        tabPanel("Diagnostics",
                 br(),
                 uiOutput("diag_singles")),
        tabPanel("Density",
                 br(),
                 uiOutput("density_picker"),
                 uiOutput("density_image")),
        tabPanel("Residuals",
                 br(),
                 uiOutput("resid_controls"),
                 uiOutput("resid_image")),
        tabPanel("Files",
                 br(),
                 uiOutput("files_list"))
      )
    )
  )
)

# ---- Server ----------------------------------------------------------------
server <- function(input, output, session) {
  
  species_runs <- reactive({
    req(input$species)
    Filter(function(d) d$species[1] == input$species, MANIFEST)
  })
  
  output$lat_picker <- renderUI({
    runs <- species_runs()
    lys  <- unique(vapply(runs, function(d) d$latyears[1], ""))
    lys  <- lys[order(vapply(lys, lat_min, 0))]
    selectInput("latyears", "Latitude range",
                choices = stats::setNames(lys, vapply(lys, pretty_lat, "")))
  })
  
  output$family_picker <- renderUI({
    req(input$latyears)
    runs <- Filter(function(d) d$latyears[1] == input$latyears, species_runs())
    fams <- sort(unique(vapply(runs, function(d) d$family[1], "")))
    selectInput("family", "Model family",
                choices = stats::setNames(fams, vapply(fams, pretty_family, "")))
  })
  
  current <- reactive({
    req(input$species, input$latyears, input$family)
    hit <- Filter(function(d)
      d$species[1] == input$species && d$latyears[1] == input$latyears &&
        d$family[1] == input$family, MANIFEST)
    req(length(hit) >= 1)
    hit[[1]]
  })
  
  sec_files <- function(section) {
    d <- current()
    d[d$section == section, , drop = FALSE]
  }
  has_section <- function(section) nrow(sec_files(section)) > 0
  path_for <- function(section, file) {
    f <- sec_files(section)
    hit <- f$path[f$fname == file]
    if (length(hit)) hit[1] else NA_character_
  }
  
  output$converge_badge <- renderUI({
    req(input$species, input$latyears, input$family)
    if (has_section("index"))
      tags$div(style = "margin-top:8px;", tags$span(class = "run-badge-ok", "Index produced"))
    else
      tags$div(style = "margin-top:8px;",
               tags$span(class = "run-badge-fail", "No index \u2014 likely non-converged"))
  })
  
  # ---- Index ---------------------------------------------------------------
  index_df <- reactive({
    p <- path_for("index", "est_by_area.csv")
    if (is.na(p)) return(NULL)
    df <- read_csv_any(p)
    if (!is.null(df) && "type" %in% names(df)) df <- df[df$type == "index", , drop = FALSE]
    df
  })
  
  output$index_status <- renderUI({
    if (!has_section("index"))
      tags$div(class = "section-empty",
               "This run has no index output \u2014 the model likely failed the sanity check.")
  })
  
  output$index_png_cw <- renderUI({
    p <- path_for("index", "index_coastwide.png")
    if (is.na(p))
      return(tags$div(class = "section-empty",
                      "No index_coastwide.png for this run (model may not have converged)."))
    tags$img(src = img_src(p), class = "half-img")
  })
  
  # ---- Modeled years (for paginated labels) --------------------------------
  model_years <- reactive({
    df <- index_df()
    if (is.null(df) || !"year" %in% names(df)) return(NULL)
    sort(unique(df$year))
  })
  
  # ---- Diagnostics (single plots) ------------------------------------------
  output$diag_singles <- renderUI({
    f <- sec_files("diagnostics")
    defs <- list(c("mesh.png", "Mesh"),
                 c("qq.png", "QQ Plot"),
                 c("anisotropy.png", "Anisotropy"),
                 c("fixed_effects.png", "Fixed effects"))
    present <- Filter(function(d) d[1] %in% f$fname, defs)
    if (!length(present))
      return(tags$div(class = "section-empty", "No single-plot diagnostics for this run."))
    tagList(lapply(present, function(x) {
      p <- path_for("diagnostics", x[1])
      tags$div(style = "margin-bottom:1.2rem;", h4(x[2]),
               tags$img(src = img_src(p), class = "half-img"))
    }))
  })
  
  # ---- Density -------------------------------------------------------------
  density_pages <- reactive({
    f    <- sec_files("diagnostics")
    pngs <- f$fname[grepl("^density_page_\\d+\\.png$", f$fname)]
    if (!length(pngs)) return(NULL)
    pg <- as.integer(sub("^density_page_(\\d+)\\.png$", "\\1", pngs))
    o  <- order(pg)
    years <- model_years()
    data.frame(file = pngs[o],
               label = vapply(pg[o], function(p) page_year_label(years, p), ""),
               stringsAsFactors = FALSE)
  })
  
  output$density_picker <- renderUI({
    d <- density_pages()
    if (is.null(d)) return(tags$div(class = "section-empty", "No density plots for this run."))
    selectInput("density_year", "Years",
                choices = stats::setNames(d$file, d$label), width = "280px")
  })
  
  output$density_image <- renderUI({
    d <- density_pages(); req(!is.null(d))
    fl <- input$density_year; req(fl); req(fl %in% d$file)
    p <- path_for("diagnostics", fl); req(!is.na(p))
    tagList(h4(paste("Density \u2014", d$label[match(fl, d$file)])),
            tags$img(src = img_src(p), class = "diag-img"))
  })
  
  # ---- Residuals -----------------------------------------------------------
  residual_sets <- reactive({
    f    <- sec_files("diagnostics")
    pngs <- f$fname[grepl("^residuals_\\d+_page_\\d+\\.png$", f$fname)]
    if (!length(pngs)) return(NULL)
    comp <- as.integer(sub("^residuals_(\\d+)_page_\\d+\\.png$", "\\1", pngs))
    pg   <- as.integer(sub("^residuals_\\d+_page_(\\d+)\\.png$", "\\1", pngs))
    fam  <- current()$family[1]
    years <- model_years()
    groups <- list()
    for (k in sort(unique(comp))) {
      sel <- comp == k; o <- order(pg[sel])
      groups[[resid_label(k, fam)]] <- data.frame(
        file = pngs[sel][o],
        label = vapply(pg[sel][o], function(p) page_year_label(years, p), ""),
        stringsAsFactors = FALSE)
    }
    groups
  })
  
  output$resid_controls <- renderUI({
    groups <- residual_sets()
    if (is.null(groups)) return(tags$div(class = "section-empty", "No residual plots for this run."))
    tagList(fluidRow(
      column(6, selectInput("resid_set", "Component", choices = names(groups), width = "320px")),
      column(6, uiOutput("resid_year_picker"))
    ))
  })
  
  output$resid_year_picker <- renderUI({
    g <- input$resid_set; req(g)
    groups <- residual_sets(); req(g %in% names(groups))
    grp <- groups[[g]]
    selectInput("resid_year", "Years",
                choices = stats::setNames(grp$file, grp$label), width = "280px")
  })
  
  output$resid_image <- renderUI({
    g <- input$resid_set; req(g)
    groups <- residual_sets(); req(g %in% names(groups))
    fl <- input$resid_year; req(fl); req(fl %in% groups[[g]]$file)
    p <- path_for("diagnostics", fl); req(!is.na(p))
    grp <- groups[[g]]
    tagList(h4(sprintf("Residuals: %s \u2014 %s", g, grp$label[match(fl, grp$file)])),
            tags$img(src = img_src(p), class = "diag-img"))
  })
  
  # ---- Files ---------------------------------------------------------------
  output$files_list <- renderUI({
    d <- current()
    link <- if (DATA_SOURCE == "github") raw_url(d$path) else file.path(LOCAL_ROOT, d$path)
    items <- Map(function(sec, fn, lk)
      tags$li(tags$b(paste0(sec, ": ")), tags$a(href = lk, target = "_blank", fn)),
      d$section, d$fname, link)
    tags$ul(items)
  })
}

shinyApp(ui, server)