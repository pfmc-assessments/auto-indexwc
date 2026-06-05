library(nwfscSurvey)
library(indexwc)
library(dplyr)
library(sdmTMB)
library(stringr)
library(tibble)

num_batches <- 34
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

# Add model index and clean family names
config_data$index_id <- seq_len(nrow(config_data))
config_data$family <- str_replace(config_data$family, "sdmTMB::", "")
config_data$family <- str_replace(config_data$family, "\\(\\)", "")

# Assign batch numbers in a round-robin fashion and filter to focal batch
config_data$batch <- rep(1:num_batches, length.out = nrow(config_data))
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

  # -- 3. Diagnostics --------------------------------------------------------
  diag <- tryCatch(
    diagnose(fit = fit),
    error = function(e) {
      cat("diagnose() failed for", config_data$species[i], ":\n",
          conditionMessage(e), "\n")
      return(NULL)
    }
  )

  # Stop here if hessian not ok -- save diagnostics if available, skip index
  if (is.null(diag) || !isTRUE(diag$sanity$hessian_ok)) {
    cat("Hessian not ok or diagnostics failed for", config_data$species[i],
        "-- skipping index\n")
    if (!is.null(diag)) {
      tryCatch(
        save_index_outputs(
          fit         = fit,
          diagnostics = diag,
          indices     = list(indices = data.frame(), plot_indices = NULL),
          dir         = "output"
        ),
        error = function(e) cat("save_index_outputs() failed (no index):",
                                conditionMessage(e), "\n")
      )
    }
    return(invisible(NULL))
  }

  # -- 4. Calculate indices --------------------------------------------------
  index <- tryCatch(
    calc_index_areas(
      data       = fit$data,
      fit        = fit,
      boundaries = c("Coastwide", "CA", "OR", "WA")
    ),
    error = function(e) {
      cat("calc_index_areas() failed for", config_data$species[i], ":\n",
          conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(index)) return(invisible(NULL))

  # -- 5. Biomass-weighted depth (custom) ------------------------------------
  # Uses the full-grid prediction object stored in the Coastwide results
  pred_data <- tryCatch({
    pred_df    <- index$results[["Coastwide"]]$pred$data
    mean_depth <- pred_df |>
      dplyr::group_by(year) |>
      dplyr::group_split() |>
      lapply(bootstrap_year_sample, family_name = config_data$family[i]) |>
      dplyr::bind_rows()
    mean_depth
  }, error = function(e) {
    cat("Biomass-weighted depth failed for", config_data$species[i], ":\n",
        conditionMessage(e), "\n")
    NULL
  })

  # -- 6. Save all outputs ---------------------------------------------------
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
