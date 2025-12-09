
###########################################
########## Processing weather data ######## 
###########################################

####### Author: Vinicius Garnica
####### Date: Apr, 2025

### Load packages -----------------------------------------------------------------------------------------------
pacman::p_load(lubridate,
               ggplot2,
               suncalc,
               data.table,
               stringr,
               dplyr)

### Load and wrangle data sets -----------------------------------------------------------------------------------
rm(list = ls())
load("data/weather_hour.Rdata")

weather_hour = weather_hour %>% 
  rename('TIME' = 'datetime',
         'TMP'='hourly_temperature_2m',
         'HMD' ='hourly_relative_humidity_2m',
         'RNF' = 'hourly_precipitation',
         'ETO' = 'hourly_et0_fao_evapotranspiration',
         'SITE' = 'site') %>%
  mutate(DATE = as.Date(TIME),
         YEAR = year(DATE),
         DOY =yday(DATE)) %>%
  filter(month(DATE) %in% c(1,2,3,4,5,6))


### Functions -----------------------------------------------------------------------------------

create_site_list = function(year, sites) {
  start_date = as.Date(paste0(year, "-01-01"))
  end_date = as.Date(paste0(year, "-06-01"))
  
  site_list = expand.grid(date = seq.Date(start_date, end_date, by = 1), 
                           sites = sites) %>%
    mutate(month = month(date),
           doy = yday(date)) %>%
    filter(month %in% c(1, 2, 3, 4, 5, 6))
  
  return(site_list)
}


### Function to add dawn periods
add_dawn = function(df) {
  index_dusk = which(df$PERIOD2 == "dawn")
  for (index in index_dusk) {
    start_row = max(1, index - 4)
    end_row = min(nrow(df), index + 3)
    df$PERIOD2[start_row:end_row] = "dawn"
  }
  return(df)
}

### Function to add dusk periods
add_dusk = function(df) {
  index_dusk = which(df$PERIOD2 == "dusk")
  for (index in index_dusk) {
    start_row = max(1, index - 3)
    end_row = min(nrow(df), index +4)
    df$PERIOD2[start_row:end_row] = "dusk"
  }
  return(df)
}
### Creating intra-day periods -----------------------------------------------------------------------------------

### To get site-specific sunset and sunrise hours, we will be using the suncalc R package and site coordinates.

site_info <- data.frame(
  site = unique(weather_hour$SITE),
  year = paste0("20",str_sub(unique(weather_hour$SITE), -2, -1))
)

site_lists <- site_info %>%
  group_by(year) %>%
  summarize(sites = list(site)) %>%
  rowwise() %>%
  mutate(site_list = list(create_site_list(year, sites))) %>%
  select(year, site_list)

sites = do.call(rbind,site_lists$site_list)

dt=read.csv("data/locations.csv") %>% select(-region)

sun_dt = left_join(sites,dt,by=c('sites' = 'site'))

### Obtain sunrise and sunset times for each site and date ----------------------------------------------------------
###This will be merged into the dataset to calculate nighttime and daytime hours.
sun = getSunlightTimes(data =sun_dt, keep = c("sunrise", "sunset"), tz = "EST")

### Add DOY column to the sun data
sun$doy = yday(sun$date)

### Join sun data with dt_list and remove unnecessary columns
sun = left_join(sun, sun_dt, by = c("date", "lat", "lon", "doy")) %>% select(-date, -lat, -lon)

### Annotate daytime and nighttime hours at each location
weather_hour= left_join(weather_hour,sun,by= c("DOY"="doy","SITE"="sites")) %>%
  na.omit() %>%
  group_by(YEAR,SITE) %>%
  mutate(PERIOD1=ifelse(TIME>=sunrise & TIME<=sunset ,"daytime","nighttime"), # period 1 will be using to get daytime and nighttime hours at each location 
         PERIOD3=rleid(PERIOD1)) %>% # we want to keep nighttime hours consecutive when grouping by. If group by DATE, we end up with nonconsecutive nighttime hours, instead of the same day.
  select(-sunrise,-sunset,-month) %>%
  group_by(YEAR,SITE) %>%
  mutate(PERIOD2 = case_when(PERIOD1=="nighttime" & lag(PERIOD1)=="daytime"~"dusk",
                             PERIOD1=="daytime" & lag(PERIOD1)=="nighttime"~"dawn")) %>%
  ungroup()  %>%
  select(-year,)


### Group by the "SITE" column and apply the function within each group
weather_hour = weather_hour %>%
  group_by(SITE, YEAR) %>%
  group_modify(~ add_dusk(.x)) %>%
  group_modify(~ add_dawn(.x))

View(weather_hour) # view results

### Check for missing data -----------------------------------------------------------------------------------  
anyNA(weather_hour$RNF)
anyNA(weather_hour$TMP)
anyNA(weather_hour$HMD)
anyNA(weather_hour$ETO)

table(weather_hour$SITE)


### Save -----------------------------------------------------------------------------------
save(weather_hour, file = "code/Window pane/data/weather_data.RData")

