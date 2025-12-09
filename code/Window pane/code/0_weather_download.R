#####################################
####### Download weather data ####### 
#####################################

####### Author: Vinicius Garnica
####### Date: Apr, 2025

### Load packages ------------------------------------------------------------------------------------ 
pacman::p_load(data.table,
               furrr,       # Consider removing if not used
               openmeteo,
               tidyverse)

# Clear the workspace
rm(list = ls())

### Load data ------------------------------------------------------------------------------------ 
data = read.csv("data/locations.csv")

### Modeling the anthesis day as anchor point ------------------------------------------------------------------------------
data$start_date = as.Date(paste0(data$year - 1, "-10-01"))
data$end_date = as.Date(paste0(data$year, "-06-01"))

### Download function ----------------------------------------------------------------------------
download_weather_data = function(lat, lon, start_date, end_date, retries = 2) {
  for (attempt in 1:retries) {
    result = tryCatch(
      {
        weather_history(
          c(lat, lon),
          start = start_date,
          end = end_date,
          hourly = list("temperature_2m", "relative_humidity_2m", "precipitation", "et0_fao_evapotranspiration"),
          timezone = "auto"
        )
      },
      error = function(e) {
        warning(paste("Attempt", attempt, "failed for location:", lat, lon, "-", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(result)) return(result)
    Sys.sleep(2)  
  }
  return(NULL)
}

### Download weather data for each location -----------------------------------------------------
weather_list = list()
failed_sites = list()

for (i in seq_len(nrow(data))) {
  site = as.character(data$site[i])
  
  lat = data$lat[i]
  lon = data$lon[i]
  start_date = data$start_date[i]
  end_date = data$end_date[i]
  
  # Download weather data
  weather = download_weather_data(lat, lon, start_date, end_date, retries = 3)
  
  if (!is.null(weather)) {
    weather$site = site
    weather_list[[site]] = weather  # Save to cache
  } else {
    failed_sites[[site]] = site
  }
}

### Combine all weather data into one dataframe -------------------------------------------------
weather_hour = bind_rows(weather_list)
failed = bind_rows(failed_sites)
length(unique(weather_hour$site))


### Save the weather data -----------------------------------------------------------------------------------
save(weather_hour, file = "data/weather_hour.Rdata")
