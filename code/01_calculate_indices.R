library(nwfscSurvey)
library(indexwc)
#library(lubridate)
library(dplyr)
library(sdmTMB)
library(stringr)
library(tibble)
# library(future)
# library(future.apply)

num_batches <- 34
# Read the batch number passed from GitHub Action
args <- commandArgs(trailingOnly = TRUE)
current_batch <- as.numeric(args[1]) # This gets the batch number

# This replaces the csv code
raw_url <- "https://raw.githubusercontent.com/pfmc-assessments/indexwc/main/data/configuration.rda"
temp_file <- tempfile(fileext = ".rda")
download.file(raw_url, temp_file, mode = "wb")
load(temp_file, envir = .GlobalEnv)
config_data <- dplyr::filter(configuration, source == "NWFSC.Combo")
sp_2025 <- c("canary rockfish", "chilipepper", "rougheye rockfish", "sablefish", "widow rockfish", "yelloweye rockfish", "yellowtail rockfish")
config_data_2025 <- dplyr::filter(config_data, species %in% sp_2025 & used == TRUE)
config_data_other <- dplyr::filter(config_data, !species %in% sp_2025)
config_data <- rbind(config_data_2025, config_data_other) |>
  dplyr::arrange(tolower(species))
unlink(temp_file)

# add model index
config_data$index_id <- seq_len(nrow(config_data))

# replace longnames in family
config_data$family <- str_replace(config_data$family, "sdmTMB::", "")
config_data$family <- str_replace(config_data$family, "\\(\\)", "")

# Use nwfscSurvey to get the data
haul <- nwfscSurvey::pull_haul(survey = "NWFSC.Combo")
dat <- nwfscSurvey::pull_catch(survey = "NWFSC.Combo",
                               common_name = config_data$species)
names(dat) <- tolower(names(dat))
dat <- dplyr::left_join(dat, haul[,c("trawl_id","area_swept_ha_der")])
# convert date string to doy
#dat$yday <- lubridate::yday(as.Date(dat$date))
# filter out a few bad locations
dat <- dplyr::filter(dat, !is.na(longitude_dd),
                !is.na(latitude_dd))

# Set default CRS for estimation / prediction
crs_out <- 32610

# add X, Y
dat <- sdmTMB::add_utm_columns(dat,
                               ll_names = c("longitude_dd","latitude_dd"),
                               utm_crs = crs_out)
# temp fix for rougheye - blackspotted
dat$common_name[which(dat$common_name == "rougheye and blackspotted rockfish")] <- "rougheye rockfish"

# Assign batch numbers in a round-robin fashion
config_data$batch <- rep(1:num_batches, length.out = nrow(config_data))
# Filter out only focal batch
config_data <- dplyr::filter(config_data, batch == current_batch)
# Plan for parallelization (adjust number of workers as needed)
#plan(multisession, workers = parallel::detectCores() - 1)

process_species <- function(i) {
  sub <- dplyr::filter(dat, common_name == config_data$species[i])
  #sub <- dplyr::mutate(sub, zday = (yday - mean(sub$yday)) / sd(sub$yday))
  sub$pass_scaled <- sub$pass - mean(range(sub$pass)) # -0.5, 0.5
  # make sure depth is negative, like config file
  sub$depth_m <- -sub$depth_m
  # apply the year, latitude, and depth filters if used
  sub <- dplyr::filter(sub,
                       latitude_dd >= config_data$min_latitude[i],
                       latitude_dd < config_data$max_latitude[i],
                       year >= config_data$min_year[i],
                       year <= config_data$max_year[i],
                       depth_m <= config_data$min_depth[i],
                       depth_m > config_data$max_depth[i]) |>
    dplyr::rename(catch_weight = total_catch_wt_kg)

  # make a mesh based on settings in config
  mesh <- sdmTMB::make_mesh(sub, xy_cols = c("X","Y"),
                            n_knots = config_data$knots[i])
  sub$fyear <- as.factor(sub$year) # year as factor
  sub$catch_weight = sub$catch_weight * 0.001 # convert to mt, matching indexwc & SS
  sub$area_km2 <- sub$area_swept_ha_der * 0.01 # convert to km2
  
  # include depth_scaled for shortspine and chilipepper
  sub$neg_depth <- -sub$depth_m
  mean_neg_depth <- mean(sub$neg_depth)
  sd_neg_depth <- sd(sub$neg_depth)
  sub$depth_scaled <- - (sub$neg_depth - mean_neg_depth) / sd_neg_depth
  sub$depth_scaled_squared <- sub$depth_scaled^2

  # include split_mendocino for yellowtail
  sub$split_mendocino <- ifelse(sub$latitude_dd > 40.1666667, "N", "S")
  
  # this is to help with printing, if done below
  st <- if(config_data$family[i] == "tweedie") {
    config_data$spatiotemporal1[i]
  } else {
    list(config_data$spatiotemporal1[i],
         config_data$spatiotemporal2[i])
  }

  # fit the model using arguments in configuration file
  # initialize to NULL and wrap in try() to avoid
  # 'system is computationally singular' error
  
  # This is from the indexwc vignette, https://github.com/pfmc-assessments/indexwc/blob/main/vignettes/a3_multiple_area_indices.Rmd
  # fit simple linear model to construct the coefficient mapping
  lm <- lm(
    formula = as.formula(config_data$formula[i]),
    data = sub
  )
  coef_names <- names(coef(lm))
  pres_not_identifiable <- names(which(is.na(coef(lm))))
  lm_pos <- lm(
    formula = as.formula(config_data$formula[i]),
    data = dplyr::filter(sub, catch_weight > 0)
  )
  pos_not_identifiable <- names(which(is.na(coef(lm_pos))))
  
  # create mapping for covariates not identifiable
  map_pres <- coef_names
  map_pres[coef_names %in% pres_not_identifiable] <- NA
  map_pres <- factor(map_pres)
  map_pos <- coef_names
  map_pos[coef_names %in% pos_not_identifiable] <- NA
  map_pos <- factor(map_pos)
  
  # initial values for these
  start_pres <- rep(0, length(coef_names))
  start_pres[coef_names %in% pres_not_identifiable] <- -20
  start_pos <- rep(0, length(coef_names))
  start_pos[coef_names %in% pos_not_identifiable] <- -20
  
  # Only modify control list if there are parameters that aren't identifiable
  if(length(pres_not_identifiable) + length(pos_not_identifiable) == 0) {
    sdmTMBcontrol <- sdmTMB::sdmTMBcontrol(newton_loops = 3)
  } else {
    sdmTMBcontrol <- sdmTMB::sdmTMBcontrol(
      map = list(b_j = map_pres,
                 b_j2 = map_pos),
      start = list(b_j = start_pres, b_j2 = start_pos),
      newton_loops = 2)
  }
  
  fit <- NULL
  fit <- suppressWarnings(try(sdmTMB(formula = as.formula(config_data$formula[i]),
                time = "year",
                offset = log(sub$area_km2),
                mesh = mesh,
                data = sub,
                spatial="on",
                spatiotemporal=st,
                anisotropy = config_data$anisotropy[i],
                family = get(config_data$family[i])(),
                share_range = config_data$share_range[i],
                control = sdmTMBcontrol)
                ,
             silent = TRUE))

  # create output directory if it doesn't exist
  if (!dir.exists("diagnostics")) {
    dir.create("diagnostics")
  }
  # Check write access
  file.access("diagnostics", mode = 2)
  if(inherits(fit, "try-error")) {
    san <- list(all_ok = FALSE)
    cat("Fitting failed for ", config_data$species[i], "\n")
  } else {
    san <- sanity(fit, silent = TRUE)
  }
  write.csv(san, file=paste0("diagnostics/sanity_",
                             stringr::str_replace_all(tolower(sub$common_name[1]),
                                                      "[^a-z0-9]+", "_"),
                             ".csv"), row.names=FALSE)

  if(class(fit) == "sdmTMB" & san$hessian_ok == TRUE) {
      # make predictions
      wcgbts_grid <- indexwc::california_current_grid

      # first filter the grid like with the data
      wcgbts_grid <- dplyr::filter(wcgbts_grid,
                                   latitude >= config_data$min_latitude[i],
                                   latitude < config_data$max_latitude[i],
                                   depth <= config_data$min_depth[i],
                                   depth > config_data$max_depth[i],
                                   area_km2_WCGBTS > 0)
      
      wcgbts_grid$neg_depth <- - wcgbts_grid$depth
      wcgbts_grid$depth_scaled <- - (wcgbts_grid$neg_depth - mean_neg_depth) / sd_neg_depth
      wcgbts_grid$depth_scaled_squared <- wcgbts_grid$depth_scaled^2
      #print(nrow(wcgbts_grid))
      # Add calendar date -- predicting to jul 1
      #wcgbts_grid$zday <- (182 - mean(sub$yday)) / sd(sub$yday)
      wcgbts_grid$pass_scaled <- 0
      # add X-Y
      wcgbts_grid <- sdmTMB::add_utm_columns(wcgbts_grid,
                                             ll_names = c("longitude","latitude"),
                                             utm_crs = crs_out)

      # replicate grid
      wcgbts_grid <- replicate_df(wcgbts_grid, time_name = "year",
                                  time_values = unique(sub$year))
      wcgbts_grid$fyear <- as.factor(wcgbts_grid$year)

      # Make coastwide index
      cat("Making predictions for ", config_data$species[i], "\n")
      pred_all <- predict(fit, wcgbts_grid, return_tmb_object = TRUE)
      cat("Generating coastwide index for ", config_data$species[i], "\n")
      index_all <- get_index(pred_all,
                             area = wcgbts_grid$area_km2_WCGBTS,
                             bias_correct = TRUE)
      index_all$index <- "Coastwide"

      # bootstrapping function for biomass wrighted depth
      bootstrap_year_sample <- function(df, n_boot = 200) {
        if(config_data$family[i] == "tweedie") {
          df$biomass <- exp(df$est)
        } else {
          df$biomass <- plogis(df$est1) * exp(df$est2)
        }
        means <- 0
        for(ii in 1:n_boot) {
          sampled <- df[sample(1:nrow(df), size = nrow(df), replace = TRUE, prob = df$biomass), ]
          means[ii] <- sum(sampled$depth * sampled$biomass) / sum(sampled$biomass)
        }

        tibble(
          year = unique(df$year),
          mean_depth = mean(means),
          ci_lower = quantile(means, 0.025),
          ci_upper = quantile(means, 0.975)
        )
      }
      # calculate biomass weighted depth
      cat("Bootstrapping depth for ", config_data$species[i], "\n")
      mean_depth <- pred_all$data |>
        group_by(year) |>
        group_split() |>
        lapply(bootstrap_year_sample) |>
        bind_rows()

      process_region <- function(region_code, region_name) {
        sub_grid <- dplyr::filter(as.data.frame(wcgbts_grid), split_state == region_code)

        if (nrow(sub_grid) > 0) {
          pred <- predict(fit, sub_grid, return_tmb_object = TRUE)
          index <- get_index(pred, area = sub_grid$area_km2_WCGBTS, bias_correct = TRUE)
          index$index <- region_name
        } else {
          index <- data.frame(
            year = sort(unique(wcgbts_grid$year)),
            est = NA,
            lwr = NA,
            upr = NA,
            log_est = NA,
            se = NA,
            se_natural = NA,
            type = "index",
            index = region_name
          )
        }

        return(index)
      }

      cat("Generating California index for ", config_data$species[i], "\n")
      index_CA <- process_region(region_code = "C",
                                 region_name = "California")
      cat("Generating Oregon index for ", config_data$species[i], "\n")
      index_OR <- process_region(region_code = "O",
                                 region_name = "Oregon")
      cat("Generating Washington index for ", config_data$species[i], "\n")
      index_WA <- process_region(region_code = "W",
                                 region_name = "Washington")

      indices <- rbind(index_all, index_CA, index_OR, index_WA)
      indices$index_id <- config_data$index[i]
      indices$common_name <- sub$common_name[1]
      indices$family <- config_data$family[i]
      indices$formula <- config_data$formula[i]
      indices$min_depth <- config_data$min_depth[i]
      indices$max_depth <- config_data$max_depth[i]
      indices$min_latitude <- config_data$min_latitude[i]
      indices$max_latitude <- config_data$max_latitude[i]
      indices$anisotropy <- config_data$anisotropy[i]
      indices$knots <- config_data$knots[i]
      indices$spatiotemporal1 <- config_data$spatiotemporal1[i]
      indices$spatiotemporal2 <- config_data$spatiotemporal2[i]
      indices$share_range <- config_data$share_range[i]

      # append date as attribute
      attr(indices, "date") <- Sys.Date()

      # create output directory if it doesn't exist
      if (!dir.exists("output")) {
        dir.create("output")
      }
      # Check write access
      file.access("output", mode = 2)

      write.csv(indices,
              paste0("output/",
                     stringr::str_replace_all(tolower(sub$common_name[1]),
                        "[^a-z0-9]+", "_"),
                        "_index.csv"), row.names=FALSE)
      write.csv(mean_depth,
              paste0("output/",
                    "biomass_weighted_depth_",
                    stringr::str_replace_all(tolower(sub$common_name[1]),
                        "[^a-z0-9]+", "_"),
                        ".csv"), row.names=FALSE)
  }
}

# Apply process_species in parallel
#future_lapply(1:nrow(config_data), process_species, future.seed = TRUE)
for(spp in 1:nrow(config_data)) {
  process_species(spp)
}
