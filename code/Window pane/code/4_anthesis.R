#############################################################################
######## Modeling wheat anthesis to set anchoring point of window pane ######
#############################################################################

####### Authors: Vinicius Garnica
####### Date: Dec, 2024

### Load packages -----------------------------------------------------------------------------------------------
pacman::p_load(data.table,
               furrr,
               meteor,
               tidyverse)

# /share/snb2023/vcastel
### Load data sets
rm(list = ls())

### Load data ------------------------------------------------------------------------------------ 
load("data/weather_hour.Rdata") 
load("code/Window pane/data/loadings.RData") 
load("code/Window pane/data/fixed_window.RData") 

### Formulas ------------------------------------------------------------------------------

veff = function(Temp) {
  result = numeric(length(Temp))
  
  condition1 = Temp < -4 | Temp > 17
  condition2 = Temp >= -4 & Temp < 3
  condition3 = Temp >= 3 & Temp < 10
  condition4 = Temp >= 10 & Temp <= 17
  
  result[condition1] = 0
  result[condition2] = (Temp[condition2] - (-4)) / (3 - (-4))
  result[condition3] = 1
  result[condition4] = (17 - Temp[condition4]) / (17 - 10)
  
  return(result)
}

fv = function(vdd, Vbase, Vsat) {
  result = numeric(length(vdd))
  
  condition1 = vdd < Vbase
  condition2 = vdd >= Vbase & vdd <= Vsat
  condition3 = vdd > Vsat
  
  result[condition1] = 0
  result[condition2] = (vdd[condition2] - Vbase) / (Vsat - Vbase)
  result[condition3] = 1
  
  return(result)
}

tt = function(Temp, Tbase, sigma) {
  result = numeric(length(Temp))
  
  condition1 = Temp <= Tbase | Temp > 37
  condition2 = Temp > Tbase & Temp <= 26
  condition3 = Temp > 26 & Temp <= 37
  
  result[condition1] = 0
  result[condition2] = 26 * (exp(-((Temp[condition2] - 26) / (2 * sigma))^2))
  result[condition3] = 26 * (1 - ((Temp[condition3] - 26) / (37 - 26))^2)
  
  return(result)
}


fp = function(ph,pbase,popt){
  result = numeric(length(ph))
  
  condition1 = ph < pbase
  condition2 = ph >= pbase & ph <= popt
  condition3 = ph > popt
  
  result[condition1] = 0
  result[condition2] = (ph[condition2] - pbase) / (popt - pbase)
  result[condition3] = 1
  
  return(result)
}

ts = function(Temp, Tbase) {
  return(sin((pi/2) * (Temp - Tbase) / (26 - Tbase)))
}


### Data wrangling ------------------------------------------------------------------------------
outcome_list = lapply(c("sev"), function(element) {
  df = as.data.frame(loadings[[element]])
  df = tibble::rownames_to_column(df, "SITE")
  df$VAR = element
  return(df)
})

print(outcome_list)

fa_loadings=do.call(rbind, outcome_list) %>% 
  mutate_at(vars(SITE,VAR),factor)

### Modeling the anthesis day as anchor point ------------------------------------------------------------------------------

dt = read.csv("data/locations.csv")

daily_weather = weather_hour %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(site,date) %>%
  reframe(T_air=mean(hourly_temperature_2m)) %>%
  drop_na() %>%
  left_join(.,dt %>% select(site,lat,region))


### Observed anthesis date ------------------------------------------------------------------------------

anthesis_dates = data.frame(
  site = c("RW22", "CL22", "MR22", "KS22", "KS23", "KS24", "PY22", "PY23", "SB22", "SB23", "SB24", "LB23", "UN23", "OX23", "AL24", "BE24", "RO24", "LB24"),
  observed_DATE = as.Date(c("2022-04-13", "2022-04-24", "2022-04-22", "2022-04-17", "2023-04-10", "2024-04-13", "2022-04-20", "2023-04-22", "2022-04-23", "2023-04-14", "2024-04-15", "2023-04-18", "2023-04-20", "2023-04-22", "2024-04-15", "2024-04-13", "2024-04-27", "2024-04-18"))
)

### Predicted anthesis date ------------------------------------------------------------------------------

anthesis = daily_weather %>%
  mutate(
    photo = photoperiod(date,lat),
  filtered = case_when(
    region == "Piedmont" & !(month(date) %in% c(9)) & (!(month(date)==10 & day(date) < 10)) ~ TRUE,
    region == "Southeastern Plains" & !(month(date) %in% c(9)) & (!(month(date)==10 & day(date) < 20)) ~ TRUE,
    region == "Middle Atlantic Coastal Plain" & !(month(date) %in% c(9)) & (!(month(date)==10 & day(date) < 30)) ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  filter(filtered == TRUE) %>%
  mutate(
    tt= tt(T_air,1.5,7.6),
    ts  = ts(T_air,1.5),
    veff = veff(T_air),
    fpa= fp(photo, 5,20)
  ) %>%
  group_by(site) %>%
  mutate(
    vdd= cumsum(veff),
    fv = fv(vdd,30,80),
    PVT = (tt*fv*fpa*ts),
    aPVT = 148+cumsum(PVT)) %>%
  filter(aPVT >= 500) %>%
  slice(1) %>%
  select(site, date) %>% 
  left_join(anthesis_dates, by = "site") %>%
  rename(predicted_DATE = date) %>% 
  mutate(predicted_DOY = yday(predicted_DATE),
         observed_DOY = yday(observed_DATE),
         diff=predicted_DOY-observed_DOY) %>%
  rename(SITE = site) %>%
  ungroup();anthesis



fa_loadings = left_join(fa_loadings,anthesis) 


write.csv(anthesis,"data/anthesis.csv")



### Merging and calcuating LAG ------------------------------------------------------------------------------

### Day of anthesis as anchor point
window_anthesis = 
  lapply(lapply(fixed_window, left_join, fa_loadings %>% select(SITE,fa1,fa2,fa3,VAR,predicted_DATE)), function(x) {
    x %>% 
      mutate(LAG= as.numeric(predicted_DATE-DATE)) %>% 
      select(-predicted_DATE) %>%
      na.omit()
  })


### Combine data sets with the same window length -------------------------------------------------------------------------------------
### This is the same as combining data sets with different intra-day periods.

pattern = c("SITE","DATE","DOY","VAR","fa1","fa2","fa3","LAG")

merge_datasets = function(data,lags,pattern) {
  result = data[[lags[1]]]
  for (i in lags[-1]) {
    result = left_join(result, data[[i]], by = pattern)
  }
  return(result)
}

window_sizes = list(c(1, 5, 9, 13, 17), c(2, 6, 10, 14, 18), c(3, 7, 11, 15, 19),c(4, 8, 12, 16, 20))

rolling_windows =lapply(window_sizes, function(window_size) {merge_datasets(window_anthesis, window_size, pattern)})

### Save --------------------------------------------------------------------------------
save(rolling_windows,file = "code/Window pane/data/rolling_windows.RData")
