library(surveyjoin)
library(lubridate)
library(dplyr)
library(sdmTMB)
library(stringr)
library(dbplyr)
# library(future)
# library(future.apply)

num_batches <- 8
# Read the batch number passed from GitHub Action
args <- commandArgs(trailingOnly = TRUE)
current_batch <- as.numeric(args[1]) # This gets the batch number

url <- "https://raw.githubusercontent.com/pfmc-assessments/indexwc/main/data-raw/configuration.csv"
config_data <- read.csv(url, stringsAsFactors = FALSE)
config_data <- dplyr::filter(config_data,
                             source == "NWFSC.Combo")
# add model index
config_data$index_id <- seq_len(nrow(config_data))

# switch signs on the depth -- negative but pos in data
config_data$max_depth <- -config_data$max_depth
config_data$min_depth <- -config_data$min_depth
# drop out pass_scaled and put in yday instead
config_data$formula <- str_replace(config_data$formula,
                                   "pass_scaled",
                                   "zday + I(zday^2)")
# replace longnames in family
config_data$family <- str_replace(config_data$family, "sdmTMB::", "")
config_data$family <- str_replace(config_data$family, "\\(\\)", "")

# handful of rockfishes aren't in surveyjoin:
#[1] "aurora rockfish"       "blackgill rockfish"    "chilipepper"
#[4] "greenspotted rockfish" "longspine thornyhead"  "rougheye rockfish"
#[7] "stripetail rockfish"   "yelloweye rockfish"    "yellowtail rockfish"

# Use the surveyjoin data
surveyjoin::cache_data()
surveyjoin::load_sql_data()
dat <- surveyjoin::get_data()

# cut down data for only species in the config file
dat <- dplyr::filter(dat,
                     common_name %in% tolower(config_data$species))

# for illustrative purposes, focus initially on WCGBTS
dat <- dplyr::filter(dat, survey_name == "NWFSC.Combo")

# convert date string to doy
dat$yday <- lubridate::yday(lubridate::ymd(dat$date))

dat <- dplyr::rename(dat, latitude_dd = lat_start,
                     longitude_dd = lon_start) |>
  dplyr::filter(!is.na(longitude_dd),
                !is.na(latitude_dd))

# filter fields for smaller file size
dat <- dplyr::select(dat,
                     #event_id,
                     common_name,
                     year,
                     yday,
                     depth_m,
                     effort,
                     catch_weight,
                     scientific_name,
                     longitude_dd,
                     latitude_dd
)

crs_out <- 32610

# Load data
#dat <- readRDS("data/wcgbts.rds")
# add X, Y
dat <- sdmTMB::add_utm_columns(dat,
                               ll_names = c("longitude_dd","latitude_dd"),
                               utm_crs = crs_out)

# Filter config data down based on available species
spp <- unique(dat$common_name)
spp <- spp[which(spp%in%config_data$species)]
config_data <- dplyr::filter(config_data, tolower(species) %in% spp)

# Assign batch numbers in a round-robin fashion
config_data$batch <- rep(1:num_batches, length.out = nrow(config_data))
# Filter out only focal batch
config_data <- dplyr::filter(config_data, batch == current_batch)

# Plan for parallelization (adjust number of workers as needed)
#plan(multisession, workers = parallel::detectCores() - 1)

process_species <- function(i) {
  sub <- dplyr::filter(dat, common_name == unique(tolower(config_data$species[i])))
  sub <- dplyr::mutate(sub, zday = (yday - mean(sub$yday)) / sd(sub$yday))
  # apply the year, latitude, and depth filters if used
  sub <- dplyr::filter(sub,
                       latitude_dd >= config_data$min_latitude[i],
                       latitude_dd < config_data$max_latitude[i],
                       year >= config_data$min_year[i],
                       year <= config_data$max_year[i],
                       depth_m >= config_data$min_depth[i],
                       depth_m <= config_data$max_depth[i])

  # make a mesh based on settings in config
  mesh <- sdmTMB::make_mesh(sub, xy_cols = c("X","Y"),
                            n_knots = config_data$knots[i])
  sub$fyear <- as.factor(sub$year) # year as factor

  # fit the model using arguments in configuration file
  # initialize to NULL and wrap in try() to avoid
  # 'system is computationally singular' errors
  fit <- NULL
  fit <- try(sdmTMB(formula = as.formula(config_data$formula[i]),
                time = "year",
                offset = log(sub$effort),
                mesh = mesh,
                data = sub,
                spatial="on",
                spatiotemporal=list(config_data$spatiotemporal1[i],
                                    config_data$spatiotemporal2[i]),
                anisotropy = config_data$anisotropy[i],
                family = get(config_data$family[i])(),
                share_range = config_data$share_range[i]), silent = TRUE)

  if(class(fit) == "sdmTMB") {
      # create output directory if it doesn't exist
      if (!dir.exists("diagnostics")) {
        dir.create("diagnostics")
      }
      # Check write access
      file.access("diagnostics", mode = 2)

      san <- sanity(fit, silent=TRUE)
      saveRDS(san, file=paste0("diagnostics/sanity_",
                               config_data$index[i], ".rds"))

      # make predictions
      wcgbts_grid <- surveyjoin::nwfsc_grid
      # first filter the grid like with the data
      wcgbts_grid <- dplyr::filter(wcgbts_grid,
                                   lat >= config_data$min_latitude[i],
                                   lat < config_data$max_latitude[i],
                                   depth_m >= config_data$min_depth[i],
                                   depth_m < config_data$max_depth[i],
                                   area > 0)
      # Add calendar date -- predicting to jul 1
      wcgbts_grid$zday <- (182 - mean(sub$yday)) / sd(sub$yday)
      # add X-Y
      wcgbts_grid <- sdmTMB::add_utm_columns(wcgbts_grid,
                                             ll_names = c("lon","lat"),
                                             utm_crs = crs_out)

      # replicate grid
      wcgbts_grid <- replicate_df(wcgbts_grid, time_name = "year",
                                  time_values = unique(sub$year))
      wcgbts_grid$fyear <- as.factor(wcgbts_grid$year)

      # Make coastwide index
      pred_all <- predict(fit, wcgbts_grid, return_tmb_object = TRUE)
      index_all <- get_index(pred_all,
                             area = wcgbts_grid$area,
                             bias_correct = TRUE)
      index_all$index <- "Coastwide"

      process_region <- function(region_code, region_name) {
        sub_grid <- dplyr::filter(as.data.frame(wcgbts_grid), split_state == region_code)

        if (nrow(sub_grid) > 0) {
          pred <- predict(fit, sub_grid, return_tmb_object = TRUE)
          index <- get_index(pred, area = sub_grid$area, bias_correct = TRUE)
          index$index <- region_name
        } else {
          index <- data.frame(
            year = sort(unique(wcgbts_grid$year)),
            est = NA,
            lwr = NA,
            upr = NA,
            log_est = NA,
            se = NA,
            type = "index",
            index = region_name
          )
        }

        return(index)
      }

      index_CA <- process_region(region_code = "C",
                                 region_name = "California")
      index_OR <- process_region(region_code = "O",
                                 region_name = "Oregon")
      index_WA <- process_region(region_code = "W",
                                 region_name = "Washington")

      indices <- rbind(index_all, index_CA, index_OR, index_WA)
      indices$index_id <- config_data$index[i]
      indices$common_name <- sub$common_name[1]

      # append date as attribute
      attr(indices, "date") <- Sys.Date()

      # create output directory if it doesn't exist
      if (!dir.exists("output")) {
        dir.create("output")
      }
      # Check write access
      file.access("output", mode = 2)

      saveRDS(indices,
              paste0("output/",
                     sub$common_name[1],"_",
                     config_data$index_id[i],".rds"))
  }
}

# Apply process_species in parallel
#future_lapply(1:nrow(config_data), process_species, future.seed = TRUE)
for(i in 1:nrow(config_data)) {
   process_species(i)
}
