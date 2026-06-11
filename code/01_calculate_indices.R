library(nwfscSurvey)
library(indexwc)
library(dplyr)
library(sdmTMB)
library(stringr)
library(tibble)

# Patch save_index_outputs() to fix fs::dir.create() bug in installed package
# (should be fs::dir_create()). Overwrites the package version in this session.
save_index_outputs <- function(
  fit,
  diagnostics,
  indices,
  dir = NULL,
  overwrite = FALSE
) {
  nwfscSurvey::check_dir(dir = dir, verbose = TRUE)
  if (!inherits(fit, "sdmTMB")) {
    cli::cli_abort(c("x" = "{.arg fit} must be an sdmTMB object"))
  }
  if (!is.list(diagnostics)) {
    cli::cli_abort(c("x" = "{.arg diagnostics} must be a list"))
  }
  if (!is.list(indices)) {
    cli::cli_abort(c("x" = "{.arg indices} must be a list"))
  }

  dir_save <- if (is.null(dir)) fit$dir else fs::path(dir, fit$dir)
  fs::dir_create(dir_save, recurse = TRUE)  # fixed: was fs::dir.create()
  if (!file.exists(dir_save)) {
    dir_save <- fs::path(getwd(), fit$dir)
    fs::dir_create(dir_save, recurse = TRUE)
  }
  if (!file.exists(dir_save)) {
    cli::cli_abort("A directory could not be created based upon the dir argument and fit$dir.")
  }
  cli::cli_alert_info("Output will be saved to {dir_save}:")

  dir_data        <- fs::path(dir_save, "data")
  dir_diagnostics <- fs::path(dir_save, "diagnostics")
  dir_index       <- fs::path(dir_save, "index")
  fs::dir_create(dir_data,        recurse = TRUE)
  fs::dir_create(dir_diagnostics, recurse = TRUE)
  fs::dir_create(dir_index,       recurse = TRUE)

  # Check for existing files
  if (!overwrite) {
    existing_files <- c(fs::path(dir_data, "data.rdata"), fs::path(dir_data, "fit.rds"))
    existing <- existing_files[fs::file_exists(existing_files)]
    if (length(existing) > 0) {
      cli::cli_abort(c(
        "x" = "Files already exist in {.path {dir_save}}",
        "i" = "Set {.code overwrite = TRUE} to replace existing files"
      ))
    }
  }

  # Save data and fit
  data_to_save <- fit$data
  save(data_to_save, file = fs::path(dir_data, "data.rdata"))
  saveRDS(fit, file = fs::path(dir_data, "fit.rds"))

  # Save sanity and AIC
  cli::cli_inform(c("*" = "Running sanity check..."))
  utils::write.table(diagnostics$sanity,
    file = fs::path(dir_diagnostics, "sanity_data_frame.csv"),
    append = FALSE, sep = ",", row.names = FALSE)
  write.table(
    rbind(c("AIC", diagnostics$aic), c("NLL", -1 * diagnostics$loglike)),
    file = fs::path(dir_diagnostics, "aic_nll.txt"),
    row.names = FALSE, col.names = FALSE)

  # Save indices
  cli::cli_inform(c("*" = "Saving and plotting indices..."))
  write.csv(indices$indices,
    file = fs::path(dir_index, "est_by_area.csv"), row.names = FALSE)

  if (any(grepl("wide", indices$indices[["area"]], ignore.case = TRUE))) {
    gg_cw <- indexwc::plot_indices(
      data = dplyr::filter(indices$indices, grepl("wide", area, ignore.case = TRUE)),
      save_loc = NULL, file_name = NULL)
    suppressMessages(ggplot2::ggsave(
      filename = fs::path(dir_index, "index_coastwide.png"),
      plot = gg_cw, height = 7, width = 7))
  }
  if (!is.null(indices$plot_indices)) {
    suppressMessages(ggplot2::ggsave(
      filename = fs::path(dir_index, "index_all_areas.png"),
      plot = indices$plot_indices, height = 7, width = 10))
  }

  # Save diagnostic plots and objects
  cli::cli_inform(c("*" = "Plotting diagnostics..."))
  if (!is.null(diagnostics$mesh_plot)) {
    suppressMessages(ggplot2::ggsave(
      filename = fs::path(dir_diagnostics, "mesh.png"),
      plot = diagnostics$mesh_plot, height = 7, width = 7))
  }
  run_diagnostics <- list(
    model = diagnostics$model, formula = diagnostics$formula,
    loglike = diagnostics$loglike, aic = diagnostics$aic,
    effects = diagnostics$effects)
  save(run_diagnostics,
    file = fs::path(dir_diagnostics, "run_diagnostics_and_estimates.rdata"))

  if (!is.null(diagnostics$qq_plot)) {
    suppressMessages(ggplot2::ggsave(
      filename = fs::path(dir_diagnostics, "qq.png"),
      plot = diagnostics$qq_plot, height = 7, width = 7))
  }
  if (!is.null(diagnostics$anisotropy_plot) &&
      inherits(diagnostics$anisotropy_plot, "ggplot")) {
    suppressMessages(ggplot2::ggsave(
      filename = fs::path(dir_diagnostics, "anisotropy.png"),
      plot = diagnostics$anisotropy_plot, height = 7, width = 7))
  }
  if (!is.null(diagnostics$fixed_effects_plot)) {
    suppressMessages(ggplot2::ggsave(
      filename = fs::path(dir_diagnostics, "fixed_effects.png"),
      plot = diagnostics$fixed_effects_plot, height = 7, width = 7))
  }
  if (!is.null(diagnostics$residual_maps_by_year)) {
    for (i in seq_along(diagnostics$residual_maps_by_year)) {
      rp <- diagnostics$residual_maps_by_year[[i]]
      if (!is.null(rp)) {
        n_pages <- ggforce::n_pages(rp)
        for (page in seq_len(n_pages)) {
          suppressMessages(ggplot2::ggsave(
            filename = fs::path(dir_diagnostics,
                                sprintf("residuals_%d_page_%02d.png", i, page)),
            plot = rp + ggforce::facet_wrap_paginate("year", nrow = 1, ncol = 2, page = page),
            height = 5, width = 10))
        }
      }
    }
  }
  if (!is.null(diagnostics$density_plots)) {
    dp <- diagnostics$density_plots[[1]]
    if (!is.null(dp)) {
      n_pages <- ggforce::n_pages(dp)
      for (page in seq_len(n_pages)) {
        suppressMessages(ggplot2::ggsave(
          filename = fs::path(dir_diagnostics,
                              sprintf("density_page_%02d.png", page)),
          plot = dp + ggforce::facet_wrap_paginate("year", nrow = 1, ncol = 2, page = page),
          height = 5, width = 10))
      }
    }
  }

  data_with_residuals <- diagnostics$data_with_residuals
  save(data_with_residuals,
    file = fs::path(dir_diagnostics, "data_with_residuals.rdata"))
  predictions <- diagnostics$predictions
  save(predictions, file = fs::path(dir_diagnostics, "predictions.rdata"))
}

num_batches <- 100
# Read the batch number passed from GitHub Action
args <- commandArgs(trailingOnly = TRUE)
current_batch <- as.numeric(args[1])

# Load configuration
raw_url <- "https://raw.githubusercontent.com/pfmc-assessments/indexwc/main/data/configuration.rda"
temp_file <- tempfile(fileext = ".rda")
download.file(raw_url, temp_file, mode = "wb")
load(temp_file, envir = .GlobalEnv)
config_data <- dplyr::filter(configuration, source == "NWFSC.Combo")
sp_2025 <- c("canary rockfish", "chilipepper", "rougheye rockfish", "sablefish",
             "widow rockfish", "yelloweye rockfish", "yellowtail rockfish")
config_data_2025 <- dplyr::filter(config_data, species %in% sp_2025 & used == TRUE)
config_data_other <- dplyr::filter(config_data, !species %in% sp_2025)
config_data <- rbind(config_data_2025, config_data_other) |>
  dplyr::arrange(tolower(species))
unlink(temp_file)

# Add model index
config_data$index_id <- seq_len(nrow(config_data))
config_data$family_short <- str_replace(config_data$family, "sdmTMB::", "") |>
  str_replace("\\(\\)", "")

# One row per batch. num_batches follows the filtered row count, so
# every row gets processed and batches beyond nrow ignored
num_batches <- nrow(config_data)
cat("Filtered config rows (= number of real batches):", num_batches, "\n")
config_data$batch <- seq_len(nrow(config_data))
config_data <- dplyr::filter(config_data, batch == current_batch)

# Bootstrap biomass-weighted depth (custom, not in indexwc)
bootstrap_year_sample <- function(df, family_name, n_boot = 200) {
  if (family_name == "tweedie") {
    df$biomass <- exp(df$est)
  } else {
    df$biomass <- plogis(df$est1) * exp(df$est2)
  }
  means <- numeric(n_boot)
  for (ii in seq_len(n_boot)) {
    sampled <- df[sample(nrow(df), size = nrow(df), replace = TRUE, prob = df$biomass), ]
    means[ii] <- sum(sampled$depth * sampled$biomass) / sum(sampled$biomass)
  }
  tibble(
    year       = unique(df$year),
    mean_depth = mean(means),
    ci_lower   = quantile(means, 0.025),
    ci_upper   = quantile(means, 0.975)
  )
}

process_species <- function(i) {
  cat("\n=== Processing row", i, ":", config_data$species[i],
      "|", config_data$family[i], "===\n")

  # -- 1. Pull and format data -----------------------------------------------
  my_data <- tryCatch(
    pull_and_format_data(configuration_to_run = config_data[i, ]),
    error = function(e) {
      cat("pull_and_format_data() failed for", config_data$species[i], ":\n",
          conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(my_data)) return(invisible(NULL))

  # -- 2. Fit model ----------------------------------------------------------
  fit <- tryCatch(
    run_sdmtmb(
      data           = my_data$data_filtered[[1]],
      family         = my_data$family,
      formula        = my_data$formula,
      n_knots        = my_data$knots,
      share_range    = my_data$share_range,
      anisotropy     = my_data$anisotropy,
      spatiotemporal = list(my_data$spatiotemporal1, my_data$spatiotemporal2)
    ),
    error = function(e) {
      cat("run_sdmtmb() failed for", config_data$species[i], ":\n",
          conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(fit)) return(invisible(NULL))

  # -- 3. Build prediction grid with effort = 1 (offset = 0) ---------------
  # diagnose() and calc_index_areas() call predict() internally; the model was
  # fit with offset = log(effort), so the prediction grid needs an effort
  # column set to 1 (log(1) = 0) to avoid an offset length mismatch error.
  pred_grid <- tryCatch({
    g <- lookup_grid(
      x             = fit$data[["survey_name"]][1],
      max_latitude  = fit$ranges$latitude_max,
      min_latitude  = fit$ranges$latitude_min,
      max_longitude = fit$ranges$longitude_max,
      min_longitude = fit$ranges$longitude_min,
      max_depth     = abs(fit$ranges$depth_max),
      years         = sort(unique(fit$data$year)),
      data          = california_current_grid
    )
    g$effort <- 1
    g
  }, error = function(e) {
    cat("lookup_grid() failed for", config_data$species[i], ":\n",
        conditionMessage(e), "\n")
    NULL
  })
  if (is.null(pred_grid)) return(invisible(NULL))

  # -- 4. Diagnostics --------------------------------------------------------
  diag <- tryCatch(
    diagnose(fit = fit, prediction_grid = pred_grid),
    error = function(e) {
      cat("diagnose() failed for", config_data$species[i], ":\n",
          conditionMessage(e), "\n")
      return(NULL)
    }
  )

  # Stop here if hessian not ok -- manually save sanity CSV, skip index
  # diag$sanity is a data frame with columns: names, logical, text
  # (returned by sanity_data() which wraps sdmTMB::sanity())
  hessian_ok <- !is.null(diag) &&
    isTRUE(diag$sanity$logical[diag$sanity$names == "hessian_ok"])
  if (!hessian_ok) {
    cat("Hessian not ok or diagnostics failed for", config_data$species[i],
        "-- skipping index\n")
    if (!is.null(diag) && !is.null(fit$dir)) {
      tryCatch({
        diag_dir <- file.path("output", fit$dir, "diagnostics")
        dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
        write.csv(diag$sanity,
                  file      = file.path(diag_dir, "sanity_data_frame.csv"),
                  row.names = FALSE)
        write.table(
          rbind(c("AIC", diag$aic), c("NLL", -1 * diag$loglike)),
          file      = file.path(diag_dir, "aic_nll.txt"),
          row.names = FALSE, col.names = FALSE
        )
      }, error = function(e) cat("Manual diagnostics save failed:",
                                 conditionMessage(e), "\n"))
    }
    return(invisible(NULL))
  }

  # -- 5. Calculate indices --------------------------------------------------
  index <- tryCatch(
    calc_index_areas(
      data            = fit$data,
      fit             = fit,
      prediction_grid = pred_grid,
      boundaries      = c("Coastwide")#c("Coastwide", "CA", "OR", "WA")
    ),
    error = function(e) {
      cat("calc_index_areas() failed for", config_data$species[i], ":\n",
          conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(index)) return(invisible(NULL))

  # -- 6. Biomass-weighted depth (custom) ------------------------------------
  # Uses the full-grid prediction object stored in the Coastwide results
  pred_data <- tryCatch({
    pred_df    <- index$results[["Coastwide"]]$pred$data
    mean_depth <- pred_df |>
      dplyr::group_by(year) |>
      dplyr::group_split() |>
      lapply(bootstrap_year_sample, family_name = config_data$family_short[i]) |>
      dplyr::bind_rows()
    mean_depth
  }, error = function(e) {
    cat("Biomass-weighted depth failed for", config_data$species[i], ":\n",
        conditionMessage(e), "\n")
    NULL
  })

  # -- 7. Save all outputs ---------------------------------------------------
  tryCatch(
    save_index_outputs(
      fit         = fit,
      diagnostics = diag,
      indices     = index,
      dir         = "output"
    ),
    error = function(e) cat("save_index_outputs() failed:", conditionMessage(e), "\n")
  )

  # Save biomass-weighted depth alongside the index files
  if (!is.null(pred_data)) {
    out_dir <- fs::path("output", fit$dir, "index")
    fs::dir_create(out_dir, recurse = TRUE)
    write.csv(
      pred_data,
      file      = fs::path(out_dir, "biomass_weighted_depth.csv"),
      row.names = FALSE
    )
  }
}

for (spp in seq_len(nrow(config_data))) {
  process_species(spp)
}
