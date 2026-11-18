# ======================================================================================================
#  That's Odd!  -  outlier detection for Hollandse Luchten PM2.5 sensors
#  Hollandse Luchten Hackathon 2026 - Challenge 1 "Spot the outliers"
#  Team: Piet Hoogkamer, Alexander (Sasha) Kuznetsov, Shirin Janghorban
# ======================================================================================================
#
#  WHAT THIS SCRIPT DOES
#  It follows the team's workflow (document "Workflow", 15 Nov) step by step. The step numbers in this file
#  are the step numbers in that document:
#
#    Step 1  Start: load the data and clean it (one clock, missing locations, impossible values)
#    Step 2  Select an area (a cluster from Sasha's file, or k-medoids clusters on location)
#    Step 3  Select a larger time window (2-3 months)
#    Step 4  Find weird SENSORS in that larger window
#              a. too much missing data (> 25 %)
#              b. PCA: sensors far from the others of the same sensor type
#              c. other rules: frozen, does not follow the others, reads structurally too high/low
#              d. mark them on a map
#              e. a table with the reason for every marked sensor
#    Step 5  Select an individual sensor (CFG$selected_sensor; all sensors are analysed anyway)
#    Step 6  Select a smaller time window of 2 weeks, sliding over the larger window
#    Step 7  Find weird TIME VALUES inside each 2-week window
#              a. STL of the sensor's own series -> remainder
#              b. PCA with nearby sensors -> residual
#              c. difference with the median of nearby sensors -> residual, plus the rules
#                 (frozen, out of range, impossible, drift)
#              Each residual is judged with LOCAL limits: the moving median of the residual
#              +/- 4 local robust standard deviations, in a moving window of 13 hours.
#              The three detectors vote. A value is only called an error when the comparison with
#              nearby sensors supports it ("unusual is not the same as wrong").
#              d. mark suspect values in time-series plots, with a symbol per cause
#              e. a table with the reason for every marked value
#    Step 8  Return: change the settings below (area, period, sensor) and run again
#    Extra   Validation: plant known faults in healthy sensors and count how many we find
#
#  HOW TO RUN
#    1. Put the data files in a folder "data" next to this script (see INPUT FILES below)
#    2. In R / RStudio:   source("thats_odd.R");  main()
#    3. Results appear in the folder "output"
#
#  PACKAGES
#    data.table and cluster (always), MASS (comes with R), and arrow OR nanoparquet to read .parquet:
#    install.packages(c("data.table", "cluster", "arrow"))
#    Note: arrow / nanoparquet are not attached with library(); they are only called when a .parquet
#    file is read. So in the RStudio "Packages" tab their box stays unticked - that is normal.
#    The versions used are written to output/session_info.txt on every run.
#
#  INPUT FILES (folder "data"; .parquet, .csv or .csv.gz all work)
#    data_sensors              kit_id, timestamp_from, pm25_kal, pm25_raw, rh, temp, latitude, longitude
#    metadata_sensors          kit_id, pm25_sensor (sensor type)
#    reference_measurements    official stations (fallback when a sensor has too few neighbours)
#    reference_station_metadata.csv
#    knmi_weather_data         humidity (U), temperature (T), wind (DD, FH)
#    Sasha's cluster file      optional, see CFG$cluster_file
# ======================================================================================================

suppressPackageStartupMessages({ library(data.table); library(cluster) })

# ======================================================================================================
#  SETTINGS - every number used in the analysis is here, with the reason for it.
#  The document "Method, validation and examples" explains each one in more detail.
# ======================================================================================================
CFG <- list(
  data_dir       = "data",
  out_dir        = "output",

  # --- Step 2: area -------------------------------------------------------------------------------
  cluster_file   = NA,             # path to Sasha's cluster file. NA = make clusters with k-medoids.
  n_clusters     = 8,              # only used without Sasha's file (8 = what Sasha used)
  areas          = c("Haarlem", "Zaandam"),     # which clusters to analyse. NULL = all clusters.
  train_area     = "Haarlem",      # validation: reported separately for a train area ...
  test_area      = "Zaandam",      # ... and a test area

  # --- Step 3: larger time window -------------------------------------------------------------------
  period_start   = "2025-10-01",   # UTC. The sensor archive ends on 2025-12-31 23:00.
  period_end     = "2026-01-01",   # exclusive. 3 months (the workflow asks for 2-3 months).

  # --- Step 4: weird sensors ------------------------------------------------------------------------
  min_coverage   = 0.75,           # 4a. more than 25 % missing = set apart. ONE rule for all plots.
  frozen_share   = 0.10,           # 4c. frozen in >= 10 % of the hours (about 9 days in 3 months)
  min_corr       = 0.60,           # 4c. correlation with the median of the other sensors below 0.6
  ratio_limit    = 1.75,           # 4c. median ratio to the others >= 1.75x or <= 1/1.75x
  pca_min_group  = 5,              # 4b. PCA is done per sensor type; a type needs >= 5 sensors in the area
  pca_factor     = 3,              # 4b. weird if its distance to the median of the OTHER sensors is
                                   #     more than 3x the typical (median) distance in its group

  # --- Step 5: one sensor to look at in detail ----------------------------------------------------
  selected_sensor = NA,            # e.g. "HLL_hl_device_249". NA = the sensor with the most flagged hours
  selected_window = NA,            # window number for the detail plots. NA = the window with the most flags

  # --- Step 6: smaller, sliding window -------------------------------------------------------------
  window_days    = 14,             # 2 weeks
  step_days      = 7,              # slide one week at a time (set 1 to slide day by day)

  # --- Step 7: weird time values --------------------------------------------------------------------
  local_half     = 6,              # local limits: window of 6 hours before + the hour + 6 hours after = 13 hours
  limit_sd       = 4,              # limit in local robust standard deviations (3 is the textbook value;
  compare_limits = c(3, 4),        #   the validation runs 3 and 4 side by side)
  min_abs        = 5,              # a flagged hour must also differ >= 5 ug/m3 from what was expected
  use_log        = TRUE,           # step 7 residuals on the log scale (log(1 + x)) or on the original scale
  compare_log    = c(TRUE, FALSE), # the validation runs both
  sd_floor_log   = 0.05,           # smallest local robust sd on the log scale (avoids dividing by ~0) ...
  sd_floor_ug    = 0.5,            # ... and on the original scale, in ug/m3
  stl_max_gap    = 0.25,           # STL needs a full series: skip STL if > 25 % of the window is missing
  nb_radius_km   = 1,              # 7b. nearby sensors within 1 km ...
  nb_radius_max  = 3,              # ... widened step by step up to 3 km if fewer than nb_min are found
  nb_min         = 3,              # at least 3 neighbours, otherwise the cluster, then the official station
  pca_components = 1,              # 7b. PCs removed in the local PCA (PC1 = what all sensors do together)
  stuck_h        = 6,              # frozen: same raw value >= 6 hours in a row
  range_min      = 0,              # out of range: negative ...
  range_max      = 200,            # ... or above 200 ug/m3 (QA list)
  implaus        = 100,            # impossible: > 100 ug/m3 AND > 3x expected
  implaus_ratio  = 3,
  drift_h        = 48,             # drift: rolling 48-hour median of the ratio to the neighbours ...
  drift_factor   = 1.75,           # ... >= 1.75x (or <= 1/1.75x) ...
  drift_min_abs  = 4,              # ... and >= 4 ug/m3
  event_factor   = 2,              # area-wide event: the area median is >= 2x its window median ...
  event_min      = 20,             # ... and >= 20 ug/m3
  humid_rh       = 90,             # moisture: sensor RH >= 90 % or KNMI RH >= 95 % while reading too high.
  humid_knmi     = 95,             #   The sensor's own RH is only used when its temperature is right (see temp_diff)
  temp_diff      = 10,             # sensor temperature > 10 C off KNMI = its RH is wrong too ...
  temp_note_share = 0.20,          # ... and if that happens in >= 20 % of hours: mounting note on the report card
  big_abs        = 20,             # "hard evidence" for red: >= 20 ug/m3 AND >= 2.5x off
  big_ratio      = 2.5,
  yellow_from    = 30,
  red_from       = 70,
  seed           = 42
)

# Towns used to give clusters a readable name (the town where most of a cluster's sensors are).
TOWNS <- data.table(
  name = c("Haarlem","Amsterdam","Zaandam","IJmuiden","Beverwijk","Hoofddorp","Hilversum","Alkmaar",
           "Purmerend","Hoorn","Heemskerk","Velsen","Laren","Huizen","Weesp","Amstelveen","Blaricum",
           "Bussum","Naarden","Eemnes","Wijk aan Zee","Velserbroek","Heemstede","Bloemendaal"),
  lat  = c(52.381,52.372,52.438,52.458,52.487,52.303,52.223,52.632,52.505,52.643,52.511,52.462,52.257,
           52.300,52.307,52.303,52.287,52.276,52.295,52.256,52.494,52.428,52.355,52.403),
  lon  = c(4.636,4.894,4.826,4.605,4.657,4.690,5.176,4.750,4.959,5.060,4.671,4.640,5.226,5.240,5.042,
           4.859,5.243,5.162,5.161,5.257,4.602,4.646,4.623,4.610))

# ======================================================================================================
#  Small helpers
# ======================================================================================================
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")

short_id <- function(sensor_id) sub("HLL_hl_device_", "#", sensor_id)

read_table <- function(name, columns = NULL) {
  # Reads data/<name>.parquet, .csv.gz or .csv - whichever exists.
  base_name <- file.path(CFG$data_dir, name)
  for (extension in c(".parquet", ".csv.gz", ".csv")) {
    file_name <- paste0(base_name, extension)
    if (!file.exists(file_name)) next
    if (extension == ".parquet") {
      if (requireNamespace("arrow", quietly = TRUE)) return(as.data.table(arrow::read_parquet(file_name, col_select = columns)))
      if (requireNamespace("nanoparquet", quietly = TRUE)) {
        table <- as.data.table(nanoparquet::read_parquet(file_name))
        return(if (is.null(columns)) table else table[, ..columns])
      }
      stop("Install 'arrow' or 'nanoparquet' to read ", file_name)
    }
    return(fread(file_name, select = columns))
  }
  stop("File not found: ", base_name, " (.parquet / .csv.gz / .csv)")
}

as_utc <- function(time_values) {
  # The sensor and KNMI timestamps are UTC without a time zone label.
  if (inherits(time_values, "POSIXct")) return(as.POSIXct(format(time_values, "%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = "UTC"))
  as.POSIXct(substr(as.character(time_values), 1, 19), tz = "UTC",
             tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"))
}

km_between <- function(lat1, lon1, lat2, lon2) {
  # Distance in km between two points on the earth (haversine formula).
  degrees_to_radians <- pi / 180
  earth_radius_km <- 6371
  haversine <- sin((lat2 - lat1) * degrees_to_radians / 2)^2 +
    cos(lat1 * degrees_to_radians) * cos(lat2 * degrees_to_radians) * sin((lon2 - lon1) * degrees_to_radians / 2)^2
  central_angle <- 2 * asin(sqrt(pmin(1, haversine)))       # in radians
  earth_radius_km * central_angle
}

to_matrix <- function(long_table, value_column, hours, sensor_ids) {
  # Turns a long table into a matrix: one row per hour, one column per sensor. Missing stays NA, never 0.
  result <- matrix(NA_real_, length(hours), length(sensor_ids), dimnames = list(NULL, sensor_ids))
  row_nr <- match(as.numeric(long_table$t), as.numeric(hours))
  col_nr <- match(long_table$kit_id, sensor_ids)
  ok <- !is.na(row_nr) & !is.na(col_nr)
  result[cbind(row_nr[ok], col_nr[ok])] <- long_table[[value_column]][ok]
  result
}

run_length_same <- function(values) {
  # For every hour: how many hours in a row the value has been exactly the same (0 where missing).
  runs <- rle(values)
  run_length <- rep(runs$lengths, runs$lengths)
  run_length[is.na(values)] <- 0L
  run_length
}

row_median <- function(value_matrix) {
  apply(value_matrix, 1, function(row) if (all(is.na(row))) NA_real_ else median(row, na.rm = TRUE))
}

rolling <- function(values, half_width, summary_function) {
  # Centred moving window of 2*half_width+1 values; NA when fewer than half the window has data.
  n <- length(values); offsets <- -half_width:half_width
  shifted <- sapply(offsets, function(offset) {
    index <- seq_len(n) + offset; out <- rep(NA_real_, n); inside <- index >= 1 & index <= n
    out[inside] <- values[index[inside]]; out
  })
  apply(shifted, 1, function(window) {
    window <- window[!is.na(window)]
    if (length(window) < max(3, length(offsets) %/% 2)) NA_real_ else summary_function(window)
  })
}

local_limits <- function(residual, sd_floor) {
  # Local limits for a residual series, hour by hour:
  #   centre = median of the 13 residuals around the hour (6 before, the hour, 6 after)
  #   spread = 1.4826 x MAD of those 13 residuals (equals the standard deviation for normal data,
  #            but is not pulled up by the outliers themselves)
  #   z      = (residual - centre) / spread      -> |z| >= CFG$limit_sd is unusual
  # Because centre and spread move with the window, the limits are different for every hour:
  # wide when the data are restless, narrow when they are calm.
  centre <- rolling(residual, CFG$local_half, median)
  spread <- pmax(1.4826 * rolling(residual, CFG$local_half, function(window) median(abs(window - median(window)))), sd_floor)
  list(z = (residual - centre) / spread, centre = centre, spread = spread)
}

to_scale <- function(values, use_log) {
  # The scale on which residuals are computed: log(1 + x) or the original ug/m3.
  # log(1 + x) turns "2x too high" into the same step at 5 and at 50 ug/m3; the +1 handles zero.
  if (use_log) log1p(pmax(values, 0)) else values
}

from_scale <- function(values, use_log) if (use_log) expm1(values) else values

# ======================================================================================================
#  STEP 1 - Start: load and clean
# ======================================================================================================
load_data <- function() {
  period_start <- as.POSIXct(CFG$period_start, tz = "UTC")
  period_end   <- as.POSIXct(CFG$period_end, tz = "UTC")
  hours        <- seq(period_start, period_end - 3600, by = "hour")

  say("Step 1 - loading sensor data")
  sensor_data <- read_table("data_sensors", c("kit_id","timestamp_from","pm25_raw","pm25_kal","rh","temp","latitude","longitude"))
  sensor_data[, t := as_utc(timestamp_from)][, timestamp_from := NULL]
  sensor_data <- sensor_data[t >= period_start & t < period_end]
  sensor_data[abs(latitude) < 1 | abs(longitude) < 1, c("latitude","longitude") := NA]   # (0,0) = no GPS fix
  sensor_data[rh <= 0 | rh > 100, rh := NA]                                              # impossible humidity
  n_duplicates <- nrow(sensor_data) - uniqueN(sensor_data, by = c("kit_id","t"))
  sensor_data <- unique(sensor_data, by = c("kit_id","t"))
  say("Step 1 - duplicate sensor-hours removed:", n_duplicates)

  sensors <- sensor_data[!is.na(latitude), .(lat = median(latitude), lon = median(longitude)), by = kit_id]   # one spot per sensor
  sensor_data <- merge(sensor_data, sensors, by = "kit_id")
  moved <- sensor_data[!is.na(latitude), .(moved_share = mean(km_between(latitude, longitude, lat, lon) > 0.2)), by = kit_id]
  sensors <- merge(sensors, moved, by = "kit_id", all.x = TRUE)
  sensor_types <- tryCatch(read_table("metadata_sensors", c("kit_id","pm25_sensor")), error = function(e) NULL)
  if (!is.null(sensor_types)) sensors <- merge(sensors, sensor_types, by = "kit_id", all.x = TRUE) else sensors[, pm25_sensor := "unknown"]
  sensors[is.na(pm25_sensor) | pm25_sensor == "", pm25_sensor := "unknown"]

  sensor_ids  <- sensors$kit_id
  pm25        <- to_matrix(sensor_data, "pm25_kal", hours, sensor_ids)
  pm25_raw    <- to_matrix(sensor_data, "pm25_raw", hours, sensor_ids)
  rh_sensor   <- to_matrix(sensor_data, "rh", hours, sensor_ids)
  temp_sensor <- to_matrix(sensor_data, "temp", hours, sensor_ids)

  say("Step 1 - official stations (always UTC+1, so minus 1 hour)")
  reference_data <- read_table("reference_measurements", c("meetlocatie_id","component","begindatumtijd","waarde"))
  reference_data <- reference_data[component == "PM2.5"]
  reference_data[, t := as.POSIXct(substr(begindatumtijd, 1, 19), tz = "UTC", format = "%Y-%m-%dT%H:%M:%S") - 3600]
  reference_data <- reference_data[t >= period_start & t < period_end]
  reference_stations <- fread(file.path(CFG$data_dir, "reference_station_metadata.csv"))[meetlocatie_id %in% unique(reference_data$meetlocatie_id)]
  reference <- to_matrix(reference_data[, .(kit_id = meetlocatie_id, t, waarde)], "waarde", hours, reference_stations$meetlocatie_id)

  say("Step 1 - KNMI weather (time = end of the hour, so minus 1 hour)")
  knmi_data <- read_table("knmi_weather_data", c("time","wsi","lat","lon","T","U","DD","FH"))
  knmi_data[, t := as_utc(time) - 3600]
  knmi_data <- knmi_data[t >= period_start & t < period_end]
  knmi_stations <- knmi_data[, .(lat = median(lat), lon = median(lon), n = sum(!is.na(U))), by = wsi][n > 0]
  knmi <- lapply(c(U = "U", T = "T", DD = "DD", FH = "FH"),
                 function(variable) to_matrix(knmi_data[, .(kit_id = wsi, t, value = get(variable))], "value", hours, knmi_stations$wsi))
  list(hours = hours, sensors = sensors, pm25 = pm25, pm25_raw = pm25_raw, rh_sensor = rh_sensor, temp_sensor = temp_sensor,
       reference = reference, reference_stations = reference_stations, knmi = knmi, knmi_stations = knmi_stations)
}

# ======================================================================================================
#  STEP 2 - Select an area (cluster)
# ======================================================================================================
nearest_town <- function(lat, lon) {
  TOWNS$name[sapply(seq_along(lat), function(i) which.min(km_between(lat[i], lon[i], TOWNS$lat, TOWNS$lon)))]
}

read_cluster_file <- function(path) {
  # Sasha's file: a sensor id column plus one 0/1 column per cluster. All zeros = the last cluster.
  # A single column with the cluster name or number also works.
  cluster_table <- fread(path)
  id_column <- intersect(c("kit_id","kit","sensor","sensor_id","id"), names(cluster_table))[1]
  if (is.na(id_column)) id_column <- names(cluster_table)[1]
  setnames(cluster_table, id_column, "kit_id")
  cluster_columns <- setdiff(names(cluster_table), "kit_id")
  if (length(cluster_columns) == 1) return(cluster_table[, .(kit_id, cluster_id = as.character(get(cluster_columns)))])
  one_hot <- as.matrix(cluster_table[, ..cluster_columns])
  label <- ifelse(rowSums(one_hot, na.rm = TRUE) == 0, "last cluster (all zeros)",
                  cluster_columns[max.col(one_hot, ties.method = "first")])
  cluster_table[, .(kit_id, cluster_id = label)]
}

make_areas <- function(dat) {
  sensors <- copy(dat$sensors)
  if (!is.na(CFG$cluster_file)) {
    say("Step 2 - clusters from Sasha's file:", CFG$cluster_file)
    sensors <- merge(sensors, read_cluster_file(CFG$cluster_file), by = "kit_id")
  } else {
    say("Step 2 - no cluster file: k-medoids (PAM) on location, k =", CFG$n_clusters)
    set.seed(CFG$seed)
    xy_km <- cbind(sensors$lat * 111, sensors$lon * 111 * cos(52.4 * pi / 180))     # degrees -> km
    sensors[, cluster_id := as.character(pam(xy_km, k = CFG$n_clusters)$clustering)]
  }
  # Readable name: the town where MOST of the cluster's sensors are (not the town nearest its centre).
  sensors[, town := nearest_town(lat, lon)]
  cluster_names <- sensors[, .N, by = .(cluster_id, town)][order(cluster_id, -N)][, .SD[1], by = cluster_id]
  cluster_names[, area := make.unique(town, sep = "_")]
  sensors <- merge(sensors, cluster_names[, .(cluster_id, area)], by = "cluster_id")
  dat$sensors <- sensors
  say("Step 2 - areas:", paste(sensors[, .N, by = area][order(-N), paste0(area, " (", N, ")")], collapse = ", "))
  dat
}

# ======================================================================================================
#  STEP 4 - Find weird sensors in the larger window
# ======================================================================================================
weird_sensors <- function(dat, area_name) {
  sensor_ids <- dat$sensors[area == area_name, kit_id]
  pm25 <- dat$pm25[, sensor_ids, drop = FALSE]; pm25_raw <- dat$pm25_raw[, sensor_ids, drop = FALSE]
  sensor_table <- data.table(kit_id = sensor_ids, area = area_name,
                             sensor_type = dat$sensors[match(sensor_ids, kit_id), pm25_sensor],
                             coverage = colMeans(!is.na(pm25)),
                             frozen_share = apply(pm25_raw, 2, function(values) mean(run_length_same(values) >= CFG$stuck_h)))
  # 4a. too much missing data - reported separately: these sensors are not odd, they just have no data
  sensor_table[, group := fifelse(coverage < CFG$min_coverage, "too little data", "")]
  sensor_table[group == "" & frozen_share >= CFG$frozen_share, group := "weird"]
  sensor_table[, reason := fifelse(coverage < CFG$min_coverage, sprintf("only %.0f%% of hours have data (limit %.0f%%)", 100 * coverage, 100 * CFG$min_coverage),
                           fifelse(frozen_share >= CFG$frozen_share, sprintf("frozen in %.0f%% of hours", 100 * frozen_share), ""))]
  sensor_table[, `:=`(corr = NA_real_, slope = NA_real_, ratio = NA_real_, pc1 = NA_real_, pc2 = NA_real_,
                      pca_dist = NA_real_, pca_typical = NA_real_)]
  candidates <- sensor_table[group == "", kit_id]
  if (length(candidates) >= 4) {
    log_pm25 <- log1p(pmax(pm25[, candidates, drop = FALSE], 0))
    # 4c. compare each sensor with the median of the OTHER sensors: correlation, slope and ratio
    for (sensor in candidates) {
      others <- row_median(log_pm25[, setdiff(candidates, sensor), drop = FALSE])
      usable <- !is.na(log_pm25[, sensor]) & !is.na(others) & pm25[, sensor] < CFG$range_max
      sensor_table[kit_id == sensor, `:=`(corr = cor(log_pm25[usable, sensor], others[usable]),
                                          slope = coef(lm(log_pm25[usable, sensor] ~ others[usable]))[2],
                                          ratio = exp(median(log_pm25[usable, sensor] - others[usable])))]
    }
    # 4b. PCA, done separately for each sensor type (the two types form two groups).
    #     Every sensor is one point, every hour one variable (log scale, gaps filled with the group median).
    #     Distance of a sensor = Euclidean distance in the PC1/PC2 plane to the median point of the
    #     OTHER sensors of its group (the sensor itself is left out of that median).
    #     Weird if that distance is more than CFG$pca_factor x the typical (median) distance in the group.
    explained <- list()
    for (type in unique(sensor_table[kit_id %in% candidates, sensor_type])) {
      group_ids <- sensor_table[kit_id %in% candidates & sensor_type == type, kit_id]
      if (length(group_ids) < CFG$pca_min_group) next
      group_log <- log_pm25[, group_ids, drop = FALSE]
      group_median <- row_median(group_log)
      for (col in seq_len(ncol(group_log))) { missing <- is.na(group_log[, col]); group_log[missing, col] <- group_median[missing] }
      group_log[is.na(group_log)] <- median(group_log, na.rm = TRUE)
      pca <- prcomp(t(group_log), center = TRUE, scale. = FALSE)
      scores <- pca$x[, 1:2, drop = FALSE]
      distance <- sapply(seq_along(group_ids), function(i) {
        others_centre <- apply(scores[-i, , drop = FALSE], 2, median)
        sqrt(sum((scores[i, ] - others_centre)^2))
      })
      sensor_table[match(group_ids, kit_id), `:=`(pc1 = scores[, 1], pc2 = scores[, 2], pca_dist = distance,
                                                   pca_typical = median(distance))]
      explained[[type]] <- round(100 * pca$sdev[1:2]^2 / sum(pca$sdev^2), 1)
    }
    attr(sensor_table, "pca_var") <- explained
    sensor_table[kit_id %in% candidates & corr < CFG$min_corr,
                 `:=`(group = "weird", reason = sprintf("does not follow the others (correlation %.2f < %.2f)", corr, CFG$min_corr))]
    sensor_table[kit_id %in% candidates & group == "" & (ratio >= CFG$ratio_limit | ratio <= 1 / CFG$ratio_limit),
                 `:=`(group = "weird", reason = sprintf("reads structurally %s (%.2fx the others)", ifelse(ratio > 1, "too high", "too low"), ratio))]
    sensor_table[kit_id %in% candidates & group == "" & !is.na(pca_dist) & pca_dist > CFG$pca_factor * pca_typical,
                 `:=`(group = "weird", reason = sprintf("far from the other %s sensors in the PCA plot (distance %.1f, typical %.1f, limit %.0fx)",
                                                        sensor_type, pca_dist, pca_typical, CFG$pca_factor))]
  }
  sensor_table[group == "", `:=`(group = "kept", reason = "ok")]
  sensor_table[]
}

plot_pca <- function(sensor_table, area_name, file_name) {
  # One panel per sensor type. The cross marks the median point of the group; triangles are weird sensors.
  explained <- attr(sensor_table, "pca_var")
  types <- names(explained)
  if (!length(types)) return(invisible())
  png(file_name, 900 * length(types), 760, type = "cairo"); par(mfrow = c(1, length(types)))
  for (type in types) {
    group <- sensor_table[sensor_type == type & !is.na(pc1)]
    plot(group$pc1, group$pc2, pch = ifelse(group$group == "weird", 17, 19), cex = 1.6,
         col = ifelse(group$group == "weird", "#d64541", "#1f5fbf"),
         xlab = sprintf("PC1 (%s%% of variance)", explained[[type]][1]), ylab = sprintf("PC2 (%s%% of variance)", explained[[type]][2]),
         main = sprintf("Step 4b - %s, %s (%d sensors)", area_name, type, nrow(group)))
    points(median(group$pc1), median(group$pc2), pch = 3, cex = 2.5, lwd = 2)
    text(group$pc1, group$pc2, short_id(group$kit_id), pos = 3, cex = 0.75)
    legend("bottomleft", c("kept", "weird", "median of the group"), col = c("#1f5fbf", "#d64541", "black"), pch = c(19, 17, 3), bty = "n")
  }
  dev.off()
}

plot_cross_tables <- function(dat, sensor_table, area_name, file_name) {
  # Correlation and slope of every sensor against every other sensor (same sensors as the PCA).
  sensor_ids <- sensor_table[group != "too little data", kit_id]
  if (length(sensor_ids) < 3) return(invisible())
  pm25 <- dat$pm25[, sensor_ids, drop = FALSE]; n <- length(sensor_ids)
  correlation <- matrix(NA_real_, n, n, dimnames = list(sensor_ids, sensor_ids)); slope <- correlation
  for (row_sensor in seq_len(n)) for (col_sensor in seq_len(n)) {
    usable <- !is.na(pm25[, row_sensor]) & !is.na(pm25[, col_sensor]) &
              pm25[, row_sensor] < CFG$range_max & pm25[, col_sensor] < CFG$range_max
    if (sum(usable) > 48 && sd(pm25[usable, row_sensor]) > 0 && sd(pm25[usable, col_sensor]) > 0) {   # a frozen sensor has no correlation
      correlation[row_sensor, col_sensor] <- cor(pm25[usable, row_sensor], pm25[usable, col_sensor])
      slope[row_sensor, col_sensor] <- coef(lm(pm25[usable, row_sensor] ~ pm25[usable, col_sensor]))[2]
    }
  }
  labels <- short_id(sensor_ids)
  fwrite(data.table(sensor = labels, round(correlation, 2)), sub(".png$", "_correlation.csv", file_name))
  fwrite(data.table(sensor = labels, round(slope, 2)), sub(".png$", "_slope.csv", file_name))
  png(file_name, 1500, 760, type = "cairo"); par(mfrow = c(1, 2), mar = c(6, 6, 3, 1))
  image(1:n, 1:n, t(correlation)[, n:1], col = colorRampPalette(c("#d64541","#f3d36b","#2e9e5b"))(20), zlim = c(0.4, 1),
        axes = FALSE, xlab = "", ylab = "", main = paste(area_name, "- correlation (red = below 0.8)"))
  axis(1, 1:n, labels, las = 2, cex.axis = 0.7); axis(2, 1:n, rev(labels), las = 2, cex.axis = 0.7)
  slope_clipped <- slope; slope_clipped[] <- pmin(2, pmax(0, slope))
  image(1:n, 1:n, t(slope_clipped)[, n:1], col = colorRampPalette(c("#d64541","#2e9e5b","#d64541"))(21), zlim = c(0, 2),
        axes = FALSE, xlab = "", ylab = "", main = paste(area_name, "- slope of row sensor on column sensor (green = 1)"))
  axis(1, 1:n, labels, las = 2, cex.axis = 0.7); axis(2, 1:n, rev(labels), las = 2, cex.axis = 0.7)
  dev.off()
}

plot_scatter_vs_others <- function(dat, sensor_table, area_name, out_dir) {
  # Scatter plot: a weird sensor against the average of the kept sensors.
  # The line is fitted ROBUSTLY (MASS::rlm) on normal values only (no frozen hours, nothing above 200),
  # so a few extreme values outside the picture cannot pull the line away from the cloud.
  kept <- sensor_table[group == "kept", kit_id]
  for (sensor in sensor_table[group == "weird", kit_id]) {
    sensor_values <- dat$pm25[, sensor]
    others_mean <- rowMeans(dat$pm25[, kept, drop = FALSE], na.rm = TRUE)
    frozen <- run_length_same(dat$pm25_raw[, sensor]) >= CFG$stuck_h
    both_present <- !is.na(sensor_values) & !is.na(others_mean)
    used <- both_present & !frozen & sensor_values < CFG$range_max & others_mean < CFG$range_max
    if (sum(used) < 48) next
    fit <- tryCatch(MASS::rlm(sensor_values[used] ~ others_mean[used], maxit = 50),
                    error = function(e) lm(sensor_values[used] ~ others_mean[used]))
    axis_limits <- c(0, quantile(c(sensor_values[used], others_mean[used]), 0.995))
    png(file.path(out_dir, paste0("step4_scatter_", sub("HLL_hl_device_", "", sensor), ".png")), 760, 760, type = "cairo")
    plot(others_mean[used], sensor_values[used], pch = 16, col = rgb(0.12, 0.37, 0.75, 0.25), xlim = axis_limits, ylim = axis_limits,
         xlab = "Average of the other (kept) sensors in the area (ug/m3)", ylab = paste("Sensor", short_id(sensor), "(ug/m3)"),
         main = paste0(short_id(sensor), ": ", sensor_table[kit_id == sensor, reason]))
    abline(0, 1, lty = 2, col = "grey40"); abline(coef(fit), col = "#d64541", lwd = 2)
    legend("topleft", bty = "n", lty = c(2, 1), col = c("grey40", "#d64541"),
           c("perfect agreement (slope 1)", sprintf("robust fit: slope %.2f, r = %.2f", coef(fit)[2], cor(sensor_values[used], others_mean[used]))))
    mtext(sprintf("Hours used in the fit: %d.   Hours left out of the fit: %d (frozen, or above 200 ug/m3).",
                  sum(used), sum(both_present) - sum(used)), side = 1, line = 4, cex = 0.8)
    dev.off()
  }
}

plot_map <- function(dat, all_sensor_tables, file_name) {
  # 4d. map: kept = green, weird = red triangle, too little data = grey
  map_points <- merge(all_sensor_tables[, .(kit_id, group)], dat$sensors[, .(kit_id, lat, lon)], by = "kit_id")
  group_colours <- c(kept = "#2e9e5b", weird = "#d64541", `too little data` = "#9aa3b2")
  png(file_name, 1000, 900, type = "cairo")
  plot(map_points$lon, map_points$lat, asp = 1 / cos(52.4 * pi / 180), pch = ifelse(map_points$group == "weird", 17, 19),
       col = group_colours[map_points$group], cex = 1.4, xlab = "longitude", ylab = "latitude", main = "Step 4d - weird sensors on the map")
  towns_shown <- TOWNS[lat %between% range(map_points$lat) & lon %between% range(map_points$lon)]
  text(towns_shown$lon, towns_shown$lat, towns_shown$name, col = "grey40", cex = 0.85)
  legend("topleft", names(group_colours), col = group_colours, pch = c(19, 17, 19), bty = "n")
  dev.off()
}

# ======================================================================================================
#  STEP 7 - detectors for weird time values (used inside every 2-week window)
# ======================================================================================================
neighbours_of <- function(dat, sensor, pool) {
  # 7b. Nearby sensors: within 1 km, widened in 1 km steps up to 3 km until at least 3 are found.
  this_sensor <- dat$sensors[kit_id == sensor]; candidates <- dat$sensors[kit_id %in% setdiff(pool, sensor)]
  distance_km <- km_between(this_sensor$lat, this_sensor$lon, candidates$lat, candidates$lon)
  for (radius in seq(CFG$nb_radius_km, CFG$nb_radius_max, by = 1)) {
    neighbours <- candidates$kit_id[distance_km <= radius]
    if (length(neighbours) >= CFG$nb_min) return(list(ids = neighbours, radius = radius))
  }
  list(ids = setdiff(pool, sensor), radius = NA_real_)   # too few nearby: use the whole area
}

expected_value <- function(dat, sensor, neighbours, window_rows) {
  # Median of the nearby sensors at the same hour; fallback: nearest official station within 25 km.
  neighbour_values <- dat$pm25[window_rows, neighbours, drop = FALSE]
  n_present <- rowSums(!is.na(neighbour_values))
  expected <- row_median(neighbour_values); expected[n_present < CFG$nb_min] <- NA
  source_ <- ifelse(is.na(expected), "official station", "nearby sensors")
  this_sensor <- dat$sensors[kit_id == sensor]
  distance_km <- km_between(this_sensor$lat, this_sensor$lon, dat$reference_stations$breedtegraad, dat$reference_stations$lengtegraad)
  for (station in dat$reference_stations$meetlocatie_id[order(distance_km)][sort(distance_km) <= 25][1:3]) {
    if (is.na(station)) next
    fill <- is.na(expected) & !is.na(dat$reference[window_rows, station])
    expected[fill] <- dat$reference[window_rows[fill], station]
  }
  source_[is.na(expected)] <- "none"
  list(value = expected, source = source_, n_present = n_present)
}

stl_remainder <- function(pm25, use_log) {
  # 7a. STL of the sensor's own series: trend + daily pattern (24 h) + remainder.
  # Only the hours of this 2-week window are used. Gaps are filled by linear interpolation
  # (STL needs a complete series); the filled hours are never judged themselves.
  series <- to_scale(pm25, use_log)
  if (mean(is.na(series)) > CFG$stl_max_gap || sum(!is.na(series)) < 72) return(list(remainder = rep(NA_real_, length(series))))
  series_filled <- approx(seq_along(series), series, xout = seq_along(series), rule = 2)$y
  fit <- stl(ts(series_filled, frequency = 24), s.window = "periodic", robust = TRUE)
  parts <- fit$time.series
  remainder <- as.numeric(parts[, "remainder"]); remainder[is.na(series)] <- NA
  list(remainder = remainder,
       fitted    = from_scale(series_filled - as.numeric(parts[, "remainder"]), use_log),   # trend + daily pattern, in ug/m3
       trend     = as.numeric(parts[, "trend"]), seasonal = as.numeric(parts[, "seasonal"]), observed = series)
}

local_pca_residual <- function(dat, sensor, neighbours, window_rows, use_log) {
  # 7b. PCA with the nearby sensors: remove what ALL of them do together (PC1); what is left for
  # the sensor is its own deviation. Each sensor is centred by its median and scaled by its MAD.
  group_ids <- c(sensor, neighbours)
  values <- to_scale(dat$pm25[window_rows, group_ids, drop = FALSE], use_log)
  centre <- apply(values, 2, median, na.rm = TRUE); scale_ <- pmax(apply(values, 2, mad, na.rm = TRUE), 1e-3)
  standardised <- sweep(sweep(values, 2, centre), 2, scale_, "/")
  available <- !is.na(standardised); standardised[!available] <- 0
  complete_rows <- rowMeans(available) >= 0.8
  if (sum(complete_rows) < 48 || length(group_ids) < 3) return(list(residual = rep(NA_real_, length(window_rows)), explained = NA))
  pca <- prcomp(standardised[complete_rows, , drop = FALSE], center = FALSE)
  loadings <- pca$rotation[, seq_len(CFG$pca_components), drop = FALSE]
  residual <- (standardised - standardised %*% loadings %*% t(loadings))[, 1]
  residual[!available[, 1]] <- NA
  list(residual = residual, explained = round(100 * sum(pca$sdev[seq_len(CFG$pca_components)]^2) / sum(pca$sdev^2), 1))
}

analyse_window <- function(dat, area_name, kept, window_rows, limit = CFG$limit_sd, use_log = CFG$use_log,
                           pm25 = dat$pm25, pm25_raw = dat$pm25_raw) {
  # Runs every detector for every kept sensor of the area inside ONE 2-week window.
  # Everything below only sees the hours of this window (window_rows).
  dat_used <- dat; dat_used$pm25 <- pm25; dat_used$pm25_raw <- pm25_raw
  sd_floor <- if (use_log) CFG$sd_floor_log else CFG$sd_floor_ug
  area_median <- row_median(pm25[window_rows, kept, drop = FALSE])
  event <- !is.na(area_median) & area_median >= CFG$event_min &
           area_median >= CFG$event_factor * median(area_median, na.rm = TRUE)
  results <- vector("list", length(kept))
  for (i in seq_along(kept)) {
    sensor <- kept[i]
    neighbours <- neighbours_of(dat_used, sensor, kept)
    expected <- expected_value(dat_used, sensor, neighbours$ids, window_rows)
    value <- pm25[window_rows, sensor]; raw_value <- pm25_raw[window_rows, sensor]
    # detector 1 (7a): STL remainder of the sensor's own series, judged with local limits
    stl_part <- stl_remainder(value, use_log); lim_stl <- local_limits(stl_part$remainder, sd_floor)
    stl_fitted <- if (is.null(stl_part$fitted)) rep(NA_real_, length(value)) else stl_part$fitted
    # detector 2 (7c): difference with the median of the nearby sensors, judged with local limits
    resid_nb <- to_scale(value, use_log) - to_scale(expected$value, use_log); lim_nb <- local_limits(resid_nb, sd_floor)
    # detector 3 (7b): local PCA residual (in standardised units), judged with local limits
    local_pca <- local_pca_residual(dat_used, sensor, neighbours$ids, window_rows, use_log)
    lim_pca <- local_limits(local_pca$residual, CFG$sd_floor_log)
    abs_diff <- value - expected$value
    # the ratio to the neighbours (log scale) is used by the drift rule, the 'hard evidence' rule and to
    # decide between spike and drop - the same in both settings of use_log
    log_ratio <- log1p(pmax(value, 0)) - log1p(pmax(expected$value, 0))
    big_enough <- !is.na(abs_diff) & abs(abs_diff) >= CFG$min_abs
    f_stl <- !is.na(lim_stl$z) & abs(lim_stl$z) >= limit & !is.na(stl_fitted) & abs(value - stl_fitted) >= CFG$min_abs
    f_nb  <- !is.na(lim_nb$z)  & abs(lim_nb$z)  >= limit & big_enough
    f_pca <- !is.na(lim_pca$z) & abs(lim_pca$z) >= limit & big_enough
    # rules (7c)
    frozen     <- run_length_same(raw_value) >= CFG$stuck_h & !is.na(value)
    out_rng    <- !is.na(value) & (value < CFG$range_min | value > CFG$range_max)
    impossible <- !is.na(value) & value > CFG$implaus & (is.na(expected$value) | value > CFG$implaus_ratio * expected$value)
    # drift: the sensor slowly moves away from its neighbours - the 48-hour median of the ratio is far from 1
    drift_ratio_log <- rolling(log_ratio, CFG$drift_h %/% 2, median)
    drift_abs       <- rolling(abs(abs_diff), CFG$drift_h %/% 2, median)
    drift <- !is.na(drift_ratio_log) & abs(drift_ratio_log) >= log(CFG$drift_factor) & !is.na(drift_abs) & drift_abs >= CFG$drift_min_abs
    results[[i]] <- data.table(t = dat$hours[window_rows], kit_id = sensor, area = area_name, pm25 = value,
                               expected = expected$value, exp_src = expected$source, n_neigh = expected$n_present,
                               nb_radius_km = neighbours$radius, stl_fitted,
                               resid_stl = stl_part$remainder, resid_nb, resid_pca = local_pca$residual,
                               z_stl = lim_stl$z, z_nb = lim_nb$z, z_pca = lim_pca$z,
                               centre_stl = lim_stl$centre, spread_stl = lim_stl$spread,
                               centre_nb = lim_nb$centre, spread_nb = lim_nb$spread,
                               centre_pca = lim_pca$centre, spread_pca = lim_pca$spread,
                               f_stl, f_nb, f_pca, frozen, out_rng, impossible, drift, drift_factor = exp(drift_ratio_log),
                               r_nb = log_ratio, absd = abs_diff, event, local_pca_var = local_pca$explained,
                               rh = dat$rh_sensor[window_rows, sensor], temp = dat$temp_sensor[window_rows, sensor])
  }
  rbindlist(results)
}

classify <- function(hours_table) {
  # 7e. Ensemble: how many of the three detectors agree, plus the hard rules -> cause, score, confidence.
  h <- hours_table
  h[, votes := f_stl + f_nb + f_pca]
  h[, compared := f_nb | f_pca]           # supported by the comparison with nearby sensors
  # Unusual is not wrong: during an area-wide event only frozen / out-of-range still count as errors.
  h[event == TRUE, `:=`(f_stl = FALSE, f_nb = FALSE, f_pca = FALSE, impossible = FALSE, drift = FALSE, compared = FALSE)]
  h[event == TRUE, votes := 0L]
  # Moisture. RH is the amount of water relative to the maximum for the temperature, so when the sensor's
  # temperature is more than 10 C off KNMI, its RH is wrong too and only the KNMI humidity is used.
  h[, temp_wrong := !is.na(temp) & !is.na(knmi_T) & abs(temp - knmi_T) > CFG$temp_diff]
  h[, humid := (!temp_wrong & !is.na(rh) & rh >= CFG$humid_rh) | (!is.na(knmi_U) & knmi_U >= CFG$humid_knmi)]
  h[, big := !is.na(absd) & abs(absd) >= CFG$big_abs & !is.na(r_nb) & abs(r_nb) >= log(CFG$big_ratio)]
  h[, fac := fifelse(is.na(drift_factor), 1, pmax(drift_factor, 1 / drift_factor))]
  h[, zmax := pmax(fifelse(f_nb, abs(z_nb), 0), fifelse(f_pca, abs(z_pca), 0), na.rm = TRUE)]
  h[, score := 0]
  h[compared == TRUE, score := pmin(85, 30 + 10 * (zmax - CFG$limit_sd) + 10 * (votes - 1))]
  h[compared == FALSE & f_stl == TRUE, score := 20]               # only unusual for itself: information, not an error
  h[drift == TRUE, score := pmax(score, fifelse(fac >= CFG$big_ratio, 60, 45))]
  h[frozen | out_rng | impossible, score := 95]
  h[, hard := frozen | out_rng | impossible | big | (drift & fac >= CFG$big_ratio)]
  h[hard == FALSE, score := pmin(score, CFG$red_from - 5)]        # no red without hard evidence
  h[is.na(pm25), score := NA]
  h[, cause := fifelse(is.na(pm25), "missing",
               fifelse(frozen, "frozen",
               fifelse(out_rng, "out of range",
               fifelse(impossible, "impossible value",
               fifelse((compared | drift) & humid & !is.na(r_nb) & r_nb > 0, "moisture",
               fifelse(drift, "drift",
               fifelse(compared & !is.na(r_nb) & r_nb < 0, "sudden drop",
               fifelse(compared, "spike",
               fifelse(f_stl, "unusual for this sensor only",
               fifelse(event, "area-wide event", "ok"))))))))))]
  h[, status := fifelse(is.na(pm25), "grey", fifelse(score >= CFG$red_from, "red", fifelse(score >= CFG$yellow_from, "yellow",
                fifelse(cause == "area-wide event", "blue", "green"))))]
  # confidence: 3 detectors agree = high, 2 = medium, 1 = low; hard rules always high
  h[, confidence := fifelse(frozen | out_rng | impossible, "High",
                    fifelse(status %in% c("green","blue","grey"), NA_character_,
                    fifelse(votes >= 3, "High", fifelse(votes == 2, "Medium", "Low"))))]
  h[]
}

# ======================================================================================================
#  STEP 6 - slide the 2-week window over the larger window and combine the results
# ======================================================================================================
make_windows <- function(hours) {
  n_hours <- length(hours); window_hours <- CFG$window_days * 24; step_hours <- CFG$step_days * 24
  starts <- seq(1, max(1, n_hours - window_hours + 1), by = step_hours)
  if (tail(starts, 1) + window_hours - 1 < n_hours) starts <- c(starts, n_hours - window_hours + 1)   # last window ends at the end of the period
  lapply(seq_along(starts), function(i) list(id = i, rows = starts[i]:(starts[i] + window_hours - 1)))
}

run_area <- function(dat, area_name, kept, limit = CFG$limit_sd, use_log = CFG$use_log,
                     pm25 = dat$pm25, pm25_raw = dat$pm25_raw, quiet = FALSE) {
  windows <- make_windows(dat$hours)
  per_window <- vector("list", length(windows))
  for (window in windows) {
    window_table <- analyse_window(dat, area_name, kept, window$rows, limit, use_log, pm25, pm25_raw)
    window_table[, `:=`(window = window$id, pos = seq_len(.N) - 1L), by = kit_id]
    # how central an hour is in this window (edges of a window are judged less well)
    window_table[, centrality := pmin(pos, CFG$window_days * 24 - 1 - pos)]
    per_window[[window$id]] <- window_table
  }
  all_windows <- rbindlist(per_window)
  # weather: nearest KNMI station (second nearest when the nearest has no value)
  weather <- rbindlist(lapply(kept, function(sensor) {
    this_sensor <- dat$sensors[kit_id == sensor]
    distance_km <- km_between(this_sensor$lat, this_sensor$lon, dat$knmi_stations$lat, dat$knmi_stations$lon)
    nearest_two <- dat$knmi_stations$wsi[order(distance_km)][1:2]
    pick <- function(knmi_matrix) { v <- knmi_matrix[, nearest_two[1]]; v[is.na(v)] <- knmi_matrix[is.na(v), nearest_two[2]]; v }
    data.table(t = dat$hours, kit_id = sensor, knmi_U = pick(dat$knmi$U), knmi_T = pick(dat$knmi$T),
               knmi_FH = pick(dat$knmi$FH), knmi_DD = pick(dat$knmi$DD))
  }))
  all_windows <- merge(all_windows, weather, by = c("t","kit_id"))
  all_windows <- classify(all_windows)
  # stability: in how many of the windows that contain an hour was it flagged?
  all_windows[, flagged := status %in% c("yellow","red")]
  stability <- all_windows[, .(windows = .N, windows_flagged = sum(flagged)), by = .(t, kit_id)]
  # final verdict per hour = the window in which the hour lies most centrally (least edge effect)
  final <- all_windows[order(-centrality)][, .SD[1], by = .(t, kit_id)]
  final <- merge(final, stability, by = c("t","kit_id"))
  if (!quiet) say(sprintf("Step 6/7 - %s: %d windows of %d days, sliding %d days", area_name, length(windows), CFG$window_days, CFG$step_days))
  list(final = final[order(kit_id, t)], all_windows = all_windows)
}

# ======================================================================================================
#  Plots and tables for step 7
# ======================================================================================================
CAUSE_SYMBOL <- c(`frozen` = 15, `out of range` = 4, `impossible value` = 8, `moisture` = 1, `drift` = 13,
                  `sudden drop` = 25, `spike` = 17)
STATUS_COLOUR <- c(red = "#d64541", yellow = "#e3a21a")
EVENT_FILL <- rgb(0.23, 0.48, 0.84, 0.15)

shade_events <- function(times, event) {
  # light blue background for hours with an area-wide event (real pollution, not a sensor error)
  hit <- which(event %in% TRUE)
  if (length(hit)) rect(times[hit], par("usr")[3], times[hit] + 3600, par("usr")[4], col = EVENT_FILL, border = NA)
}

plot_overview <- function(final, sensor, file_name) {
  # 7d. The whole larger window for one sensor: measured, expected, and a symbol per flagged hour.
  d <- final[kit_id == sensor][order(t)]
  if (!nrow(d)) return(invisible())
  png(file_name, 1500, 600, type = "cairo"); par(mar = c(3, 4.5, 3, 1))
  y_max <- min(150, max(c(d$pm25, d$expected), na.rm = TRUE) * 1.1)
  plot(d$t, d$pm25, type = "n", ylim = c(0, y_max), xlab = "", ylab = "PM2.5 (ug/m3)",
       main = sprintf("Step 7d - sensor %s (%s), %s to %s", short_id(sensor), d$area[1], format(min(d$t), "%d %b %Y"), format(max(d$t), "%d %b %Y")))
  shade_events(d$t, d$event)
  lines(d$t, d$expected, col = "grey20", lwd = 1.2)
  lines(d$t, d$pm25, col = "#1f5fbf")
  marked <- d[status %in% names(STATUS_COLOUR) & cause %in% names(CAUSE_SYMBOL)]
  if (nrow(marked)) points(marked$t, pmin(marked$pm25, y_max), pch = CAUSE_SYMBOL[marked$cause], cex = 1.4, lwd = 2,
                           col = STATUS_COLOUR[marked$status], bg = STATUS_COLOUR[marked$status])
  causes_shown <- intersect(names(CAUSE_SYMBOL), unique(marked$cause))
  legend("topleft", bty = "n", cex = 0.85, ncol = 3,
         legend = c("sensor", "expected (median of nearby sensors)", "blue band = area-wide event (real, not an error)",
                    "red = certainly wrong", "orange = doubtful", causes_shown),
         col = c("#1f5fbf", "grey20", rgb(0.23, 0.48, 0.84, 0.45), STATUS_COLOUR, rep("grey30", length(causes_shown))),
         lty = c(1, 1, NA, NA, NA, rep(NA, length(causes_shown))),
         pch = c(NA, NA, 15, 15, 15, CAUSE_SYMBOL[causes_shown]),
         pt.cex = c(1, 1, 2, 1.5, 1.5, rep(1.2, length(causes_shown))))
  dev.off()
}

plot_residuals <- function(window_table, sensor, file_name) {
  # 7d. One 2-week window, one panel per detector: the residual (grey line), the local limits
  # (red dashed: moving median +/- 4 local robust sd over 13 hours) and the hours flagged by that detector.
  d <- window_table[kit_id == sensor][order(t)]
  if (!nrow(d)) return(invisible())
  panels <- list(
    list("resid_stl", "centre_stl", "spread_stl", "f_stl", "7a - STL remainder of the sensor's own series"),
    list("resid_pca", "centre_pca", "spread_pca", "f_pca", "7b - local PCA residual (what is left after removing what all nearby sensors do together)"),
    list("resid_nb",  "centre_nb",  "spread_nb",  "f_nb",  "7c - difference with the median of the nearby sensors"))
  png(file_name, 1500, 1200, type = "cairo"); par(mfrow = c(3, 1), mar = c(3, 4.5, 3, 1))
  for (panel in panels) {
    residual <- d[[panel[[1]]]]; upper <- d[[panel[[2]]]] + CFG$limit_sd * d[[panel[[3]]]]; lower <- d[[panel[[2]]]] - CFG$limit_sd * d[[panel[[3]]]]
    y_range <- range(c(residual, upper, lower), na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(y_range))) { plot.new(); title(paste(panel[[5]], "- not available in this window")); next }
    plot(d$t, residual, type = "n", ylim = y_range, xlab = "", ylab = "residual",
         main = sprintf("%s\nred dashed = local limits: moving median +/- %s local robust sd (13-hour window)", panel[[5]], CFG$limit_sd))
    shade_events(d$t, d$event)
    lines(d$t, residual, col = "grey35")
    lines(d$t, upper, col = "#d64541", lty = 2); lines(d$t, lower, col = "#d64541", lty = 2)
    hit <- d[[panel[[4]]]] %in% TRUE & d$status %in% names(STATUS_COLOUR)
    if (any(hit)) points(d$t[hit], residual[hit], pch = 19, cex = 1.3, col = STATUS_COLOUR[d$status[hit]])
  }
  dev.off()
}

plot_stl_decomposition <- function(dat, sensor, window_rows, use_log, file_name) {
  # Step 5/7a. The STL model of one sensor in one 2-week window: data = trend + daily pattern + remainder.
  stl_part <- stl_remainder(dat$pm25[window_rows, sensor], use_log)
  if (is.null(stl_part$trend)) return(invisible())
  times <- dat$hours[window_rows]; scale_label <- if (use_log) "log(1 + PM2.5)" else "PM2.5 (ug/m3)"
  png(file_name, 1400, 1100, type = "cairo"); par(mfrow = c(4, 1), mar = c(2.5, 4.5, 2.5, 1))
  plot(times, stl_part$observed, type = "l", col = "#1f5fbf", xlab = "", ylab = scale_label,
       main = sprintf("STL of sensor %s, %s to %s - measured", short_id(sensor), format(min(times), "%d %b"), format(max(times), "%d %b %Y")))
  plot(times, stl_part$trend, type = "l", xlab = "", ylab = "trend", main = "trend (slow changes over days)")
  plot(times, stl_part$seasonal, type = "l", xlab = "", ylab = "daily pattern", main = "daily pattern (repeats every 24 hours)")
  plot(times, stl_part$remainder, type = "h", col = "grey35", xlab = "", ylab = "remainder", main = "remainder = measured - trend - daily pattern")
  abline(h = 0, col = "grey60")
  dev.off()
}

plot_log_comparison <- function(dat, sensor, window_rows, file_name) {
  # The effect of the log transformation. The same sensor and window, STL + local limits,
  # left on the log scale, right on the original scale.
  times <- dat$hours[window_rows]; value <- dat$pm25[window_rows, sensor]
  png(file_name, 1600, 900, type = "cairo"); par(mfcol = c(2, 2), mar = c(3, 4.5, 3.5, 1))
  for (use_log in c(TRUE, FALSE)) {
    stl_part <- stl_remainder(value, use_log)
    if (is.null(stl_part$fitted)) { plot.new(); plot.new(); next }
    limits <- local_limits(stl_part$remainder, if (use_log) CFG$sd_floor_log else CFG$sd_floor_ug)
    flagged <- !is.na(limits$z) & abs(limits$z) >= CFG$limit_sd & abs(value - stl_part$fitted) >= CFG$min_abs
    scale_name <- if (use_log) "LOG scale: log(1 + PM2.5)" else "ORIGINAL scale: PM2.5 in ug/m3"
    plot(times, value, type = "l", col = "#1f5fbf", xlab = "", ylab = "PM2.5 (ug/m3)",
         main = sprintf("%s - measured (blue) and STL model (black)", scale_name))
    lines(times, stl_part$fitted, col = "grey15")
    points(times[flagged], value[flagged], pch = 19, col = "#d64541", cex = 1.3)
    upper <- limits$centre + CFG$limit_sd * limits$spread; lower <- limits$centre - CFG$limit_sd * limits$spread
    plot(times, stl_part$remainder, type = "l", col = "grey35", ylim = range(c(stl_part$remainder, upper, lower), na.rm = TRUE),
         xlab = "", ylab = "STL remainder",
         main = sprintf("STL remainder with local limits (+/- %s local robust sd) - %d hours flagged", CFG$limit_sd, sum(flagged)))
    lines(times, upper, col = "#d64541", lty = 2); lines(times, lower, col = "#d64541", lty = 2)
    points(times[flagged], stl_part$remainder[flagged], pch = 19, col = "#d64541", cex = 1.3)
  }
  dev.off()
}

episodes <- function(final) {
  marked <- final[status %in% c("yellow","red","blue")][order(kit_id, t)]
  if (!nrow(marked)) return(data.table())
  marked[, episode := cumsum(c(TRUE, diff(as.numeric(t)) > 7200 | kit_id[-1] != kit_id[-.N] | cause[-1] != cause[-.N]))]
  marked[, .(kit_id = kit_id[1], area = area[1], cause = cause[1], start = min(t), end = max(t) + 3600, hours = .N,
             score = max(score, na.rm = TRUE), status = status[which.max(score)], confidence = confidence[which.max(score)],
             detectors_agreeing = max(votes), peak = max(pm25, na.rm = TRUE), expected = median(expected, na.rm = TRUE),
             flagged_in_windows = paste0(max(windows_flagged), "/", max(windows))), by = episode][, episode := NULL][]
}

# ======================================================================================================
#  EXTRA - validation: planted faults in healthy sensors
# ======================================================================================================
plant_faults <- function(dat, kept, seed = CFG$seed) {
  # Known faults in at most 1 in 4 sensors (so the neighbours stay honest), at random times:
  # 4 spikes (+25..+60), 2 impossible values (300..900), 12 frozen hours, 4 days of x2.2 drift.
  set.seed(seed)
  pm25 <- dat$pm25; pm25_raw <- dat$pm25_raw
  victims <- sample(kept, max(1, min(8, floor(length(kept) / 4))))
  truth <- list()
  for (sensor in victims) {
    candidate_hours <- which(!is.na(pm25[, sensor]))
    candidate_hours <- candidate_hours[candidate_hours > 72 & candidate_hours < nrow(pm25) - 120]
    label <- rep("", nrow(pm25))
    for (hour in sample(candidate_hours, 4)) {
      added <- runif(1, 25, 60)
      pm25[hour, sensor] <- pm25[hour, sensor] + added; pm25_raw[hour, sensor] <- pm25_raw[hour, sensor] + added; label[hour] <- "spike"
    }
    for (hour in sample(candidate_hours, 2)) if (label[hour] == "") {
      value <- runif(1, 300, 900)
      pm25[hour, sensor] <- value; pm25_raw[hour, sensor] <- value; label[hour] <- "impossible"
    }
    frozen_start <- sample(candidate_hours, 1); frozen_hours <- frozen_start:(frozen_start + 11); frozen_hours <- frozen_hours[label[frozen_hours] == ""]
    pm25[frozen_hours, sensor] <- pm25[frozen_start, sensor]; pm25_raw[frozen_hours, sensor] <- pm25_raw[frozen_start, sensor]; label[frozen_hours] <- "frozen"
    drift_start <- sample(candidate_hours, 1); drift_hours <- drift_start:(drift_start + 95); drift_hours <- drift_hours[label[drift_hours] == ""]
    pm25[drift_hours, sensor] <- pm25[drift_hours, sensor] * 2.2; pm25_raw[drift_hours, sensor] <- pm25_raw[drift_hours, sensor] * 2.2; label[drift_hours] <- "drift"
    truth[[sensor]] <- data.table(t = dat$hours, kit_id = sensor, planted = label)
  }
  list(pm25 = pm25, pm25_raw = pm25_raw, truth = rbindlist(truth))
}

validate <- function(dat, area_name, kept, role) {
  # Runs the whole step 6/7 on data with planted faults, for every combination of
  # limit (3 or 4 local sd) and scale (log or original), and counts what was found.
  planted <- plant_faults(dat, kept)
  rows <- list()
  for (use_log in CFG$compare_log) for (limit in CFG$compare_limits) {
    final <- run_area(dat, area_name, kept, limit, use_log, planted$pm25, planted$pm25_raw, quiet = TRUE)$final
    check <- merge(final, planted$truth, by = c("t","kit_id"))[order(kit_id, t)]
    # hours within 24 hours of a planted fault are not counted as false alarms
    check[, near := frollapply(as.integer(planted != ""), 49, max, align = "center", fill = 0), by = kit_id]
    rates <- function(flag) {
      check[, flag_ := flag]
      found <- check[planted != "" & !is.na(pm25), .(rate = mean(flag_, na.rm = TRUE)), by = planted]
      found <- setNames(found$rate, found$planted)
      c(found[c("spike","impossible","frozen","drift")],
        false_alarm = check[planted == "" & near == 0 & !is.na(pm25), mean(flag_, na.rm = TRUE)])
    }
    add_row <- function(method, flag, lim = limit, scale = if (use_log) "log" else "original") {
      rows[[length(rows) + 1]] <<- data.table(area = area_name, role = role, scale = scale, limit = lim, method = method,
                                              t(round(100 * rates(flag), 1)))
    }
    add_row("STL own series (7a)", check$f_stl)
    add_row("vs nearby sensors (7c)", check$f_nb)
    add_row("local PCA (7b)", check$f_pca)
    add_row("rules only", check$frozen | check$out_rng | check$impossible | check$drift)
    add_row("FINAL (ensemble + rules)", check$status %in% c("yellow","red"))
    rows[[length(rows) + 1]] <- data.table(area = area_name, role = role, scale = if (use_log) "log" else "original", limit = limit,
                                           method = "FINAL - false RED",
                                           false_alarm = round(100 * check[planted == "" & near == 0 & !is.na(pm25), mean(status == "red")], 2))
    if (use_log == CFG$compare_log[1] && limit == CFG$compare_limits[1]) {
      check[, naive3 := { m <- mean(pm25, na.rm = TRUE); s <- sd(pm25, na.rm = TRUE); abs(pm25 - m) > 3 * s }, by = kit_id]
      add_row("naive mean +/- 3 sd (whole period)", check$naive3, lim = NA, scale = "original")
    }
  }
  rbindlist(rows, fill = TRUE)
}

# ======================================================================================================
#  MAIN
# ======================================================================================================
main <- function() {
  dir.create(CFG$out_dir, showWarnings = FALSE)
  dat <- make_areas(load_data())
  areas <- if (is.null(CFG$areas)) unique(dat$sensors$area) else intersect(CFG$areas, unique(dat$sensors$area))
  if (!length(areas)) stop("None of CFG$areas found. Available: ", paste(unique(dat$sensors$area), collapse = ", "))
  sensor_tables <- list(); final_tables <- list(); window_tables <- list(); validation_tables <- list()
  for (area_name in areas) {
    # ---- step 4: weird sensors in the larger window
    sensor_table <- weird_sensors(dat, area_name); sensor_tables[[area_name]] <- sensor_table
    say(sprintf("Step 4 - %s: %d sensors | kept %d | weird %d | too little data %d", area_name, nrow(sensor_table),
                sum(sensor_table$group == "kept"), sum(sensor_table$group == "weird"), sum(sensor_table$group == "too little data")))
    plot_pca(sensor_table, area_name, file.path(CFG$out_dir, paste0("step4b_pca_", area_name, ".png")))
    plot_cross_tables(dat, sensor_table, area_name, file.path(CFG$out_dir, paste0("step4c_cross_", area_name, ".png")))
    plot_scatter_vs_others(dat, sensor_table, area_name, CFG$out_dir)
    kept <- sensor_table[group == "kept", kit_id]
    if (length(kept) < CFG$nb_min + 1) { say("   too few kept sensors for step 7"); next }
    # ---- steps 6 and 7: sliding 2-week windows, weird time values
    area_result <- run_area(dat, area_name, kept)
    final_tables[[area_name]] <- area_result$final
    window_tables[[area_name]] <- area_result$all_windows
    fwrite(area_result$all_windows[, .(window, t, kit_id, status, cause, score)], file.path(CFG$out_dir, paste0("step6_all_windows_", area_name, ".csv.gz")))
    # ---- extra: validation with planted faults
    role <- if (area_name == CFG$train_area) "train" else if (area_name == CFG$test_area) "test" else "other"
    if (role != "other") validation_tables[[area_name]] <- validate(dat, area_name, kept, role)
  }
  all_sensors <- rbindlist(sensor_tables, fill = TRUE); final <- rbindlist(final_tables, fill = TRUE)
  plot_map(dat, all_sensors, file.path(CFG$out_dir, "step4d_map.png"))
  fwrite(all_sensors, file.path(CFG$out_dir, "step4_weird_sensors.csv"))
  fwrite(final[, .(t, kit_id, area, pm25, expected, exp_src, nb_radius_km, z_stl = round(z_stl, 2), z_nb = round(z_nb, 2),
                   z_pca = round(z_pca, 2), f_stl, f_nb, f_pca, votes, cause, score, confidence, status, windows, windows_flagged,
                   rh, temp, knmi_U, knmi_T)],
         file.path(CFG$out_dir, "step7_sensor_hours.csv.gz"))
  episode_table <- episodes(final); fwrite(episode_table, file.path(CFG$out_dir, "step7_episodes.csv"))
  validation <- rbindlist(validation_tables, fill = TRUE); fwrite(validation, file.path(CFG$out_dir, "validation.csv"))

  # ---- step 7d: overview plot per area for the kept sensor with the most flagged hours
  for (area_name in names(final_tables)) {
    sensor <- final_tables[[area_name]][, .(n = sum(status %in% c("yellow","red"))), by = kit_id][order(-n)]$kit_id[1]
    plot_overview(final_tables[[area_name]], sensor, file.path(CFG$out_dir, paste0("step7_overview_", area_name, ".png")))
  }

  # ---- step 5: one selected sensor, looked at in detail in one 2-week window
  selected <- CFG$selected_sensor
  if (is.na(selected)) selected <- final[, .(n = sum(status %in% c("yellow","red"))), by = kit_id][order(-n)]$kit_id[1]
  if (selected %in% final$kit_id) {
    area_name <- final[kit_id == selected, area[1]]
    windows <- make_windows(dat$hours)
    window_id <- CFG$selected_window
    if (is.na(window_id)) {
      per_window <- window_tables[[area_name]][kit_id == selected, .(n = sum(status %in% c("yellow","red"))), by = window]
      window_id <- per_window[order(-n)]$window[1]
    }
    rows <- windows[[window_id]]$rows
    tag <- paste0(sub("HLL_hl_device_", "", selected), "_window", window_id)
    say(sprintf("Step 5 - selected sensor %s, window %d (%s to %s)", short_id(selected), window_id,
                format(dat$hours[min(rows)], "%d %b"), format(dat$hours[max(rows)], "%d %b %Y")))
    plot_overview(final, selected, file.path(CFG$out_dir, paste0("step5_selected_", sub("HLL_hl_device_", "", selected), "_overview.png")))
    plot_residuals(window_tables[[area_name]][window == window_id], selected, file.path(CFG$out_dir, paste0("step7_residuals_", tag, ".png")))
    plot_stl_decomposition(dat, selected, rows, CFG$use_log, file.path(CFG$out_dir, paste0("step7a_stl_decomposition_", tag, ".png")))
    plot_log_comparison(dat, selected, rows, file.path(CFG$out_dir, paste0("step7a_log_vs_original_", tag, ".png")))
  }

  # ---- report card per sensor
  card <- final[, .(hours = .N, green = round(mean(status %in% c("green","blue")), 3), yellow = round(mean(status == "yellow"), 3),
                    red = round(mean(status == "red"), 3), unstable_hours = sum(windows_flagged > 0 & windows_flagged < windows),
                    temp_wrong_share = round(mean(temp_wrong, na.rm = TRUE), 3),
                    main_issue = { issues <- cause[!cause %in% c("ok","missing","area-wide event","unusual for this sensor only")]
                                   if (length(issues)) names(sort(-table(issues)))[1] else "none" }),
                by = .(kit_id, area)]
  card[, note := fifelse(!is.na(temp_wrong_share) & temp_wrong_share >= CFG$temp_note_share,
                         "temperature often > 10 C off KNMI: check mounting (sun / warm wall); its own RH is not used", "")]
  card <- merge(all_sensors[, .(kit_id, area, sensor_type, coverage = round(coverage, 3), group, reason)], card, by = c("kit_id","area"), all.x = TRUE)
  fwrite(card, file.path(CFG$out_dir, "report_card.csv"))
  say("RESULTS"); print(all_sensors[, .N, by = .(area, group)]); print(round(prop.table(table(final$status)), 4)); print(validation)
  writeLines(capture.output(sessionInfo()), file.path(CFG$out_dir, "session_info.txt"))   # R and package versions used
  say("done ->", CFG$out_dir)
  invisible(list(dat = dat, sensors = all_sensors, final = final, episodes = episode_table, validation = validation))
}

if (sys.nframe() == 0) main()
