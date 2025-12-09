
############################################
####### Weather variable engineering ####### 
############################################

####### Authors: Vinicius Garnica
####### Date: Apr, 2025

### Load packages -----------------------------------------------------------------------------------------------
pacman::p_load(data.table,
               lubridate,
               tidyverse)

### Load data sets
rm(list = ls())
load("code/Window pane/data/weather_data.RData")

### Transform
weather_hour$YEAR = as.factor(weather_hour$YEAR)
weather_hour$SITE = as.factor(weather_hour$SITE)

### Glimpse
glimpse(weather_hour)

### Indicator variables -------------------------------------------------------------------------------------------

# Robust peak detection algorithm (using z-scores) https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data
ThresholdingAlgo = function(y,lag,threshold,influence) {
  signals = rep(0,length(y))
  filteredY = y[0:lag]
  avgFilter = NULL
  stdFilter = NULL
  avgFilter[lag] = mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] = sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] = 1;
      } else {
        signals[i] = -1;
      }
      filteredY[i] = influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] = 0
      filteredY[i] = y[i]
    }
    avgFilter[i] = mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] = sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

### Feature engineering of hourly data -------------------------------------------------------------------------------------------


data_hour = 
  weather_hour %>%
  group_by(SITE) %>%
  mutate(
    # Relative humidity
    RH.L40 = ifelse(HMD <= 40, 1, 0), # HMD lower than 40% (1,0)
    RH.G90 = ifelse(HMD >= 90, 1, 0), # HMD greater than 90% (1,0)

    # Temperature 
    T.7T10 = ifelse(TMP >= 7 & TMP <= 10, 1, 0),
    T.10T13 = ifelse(TMP >= 10 & TMP <= 13, 1, 0), 
    T.13T16 = ifelse(TMP >= 13 & TMP <= 16, 1, 0), 
    T.16T19 = ifelse(TMP >= 16 & TMP <= 19, 1, 0), 
    T.19T22 = ifelse(TMP >= 19 & TMP <= 22, 1, 0), 
    T.22T25 = ifelse(TMP >= 22 & TMP <= 25, 1, 0),
    T.25T28 = ifelse(TMP >= 25 & TMP <= 28, 1, 0),
    T.G28 = ifelse(TMP >= 28, 1, 0),
    
    # Rainfall
    R.EPISODE = ifelse(RNF > 0, 1, 0), 
    
    # ETo
    ETo.G0.4 = ifelse(ETO > 0.4, 1, 0), 

    # Combination of variables
    TRH.16T19nRH.L40 = ifelse(TMP >= 16 & TMP <= 19 & HMD <= 35, 1, 0), 
    TRH.19T22nRH.L40 = ifelse(TMP >= 19 & TMP <= 22 & HMD <= 35, 1, 0), 
    TRH.22T25nRH.L40 = ifelse(TMP >= 22 & TMP <= 25 & HMD <= 35, 1, 0), 
    TRH.25T28nRH.L40 = ifelse(TMP >= 25 & TMP <= 28 & HMD <= 35, 1, 0),
    
    TRH.7T10nRH.G85 = ifelse(TMP >= 7 & TMP <= 10 & HMD >= 85, 1, 0),      
    TRH.10T13nRH.G85 = ifelse(TMP >= 10 & TMP <= 13 & HMD >= 85, 1, 0),   
    TRH.13T16nRH.G85 = ifelse(TMP >= 13 & TMP <= 16 & HMD >= 85, 1, 0), 
    TRH.16T19nRH.G85 = ifelse(TMP >= 16 & TMP <= 19 & HMD >= 85, 1, 0), 
    TRH.19T22nRH.G85 = ifelse(TMP >= 19 & TMP <= 22 & HMD >= 85, 1, 0), 
    
    TR.3T7nR.G.0.2 = ifelse(TMP >= 3 & TMP <= 7 & RNF > 0.2, 1, 0), 
    TR.7T10nR.G.0.2 = ifelse(TMP >= 7 & TMP <= 10 & RNF > 0.2, 1, 0), 
    TR.10T13nR.G.0.2 = ifelse(TMP >= 10 & TMP <= 13 & RNF > 0.2, 1, 0), 
    TR.13T16nR.G.0.2 = ifelse(TMP >= 13 & TMP <= 16 & RNF >0.2, 1, 0), 
    TR.16T19nR.G.0.2 = ifelse(TMP >= 16 & TMP <= 19 & RNF >0.2, 1, 0),
    TR.19T22nR.G.0.2 = ifelse(TMP >= 19 & TMP <= 22 & RNF >0.2, 1, 0),

    RH8.peak4 = ThresholdingAlgo(HMD,8,4,1)$signals,
    T8.peak4 = ThresholdingAlgo(TMP,8,4,1)$signals

  ) %>%
  left_join(., weather_hour %>% select(SITE, DATE, YEAR, TIME, PERIOD1, PERIOD2, PERIOD3, DOY))


### Supporting functions for variables -------------------------------------------------------------------------------------------

### Define a function to calculate the maximum consecutive hours above a certain level
MAX_rl = function(var, level) {
  x = ifelse(var >= level, 1, 0)
  z = rle(x) 
  ifelse(length(z$lengths) == 1,
         ifelse(z$values > 0, z$lengths, 0),
         max((z$lengths[z$values > 0])))
}



### Define a function to count events where a variable is above or equal to a level for n hours
COUNT_rl_above = function(var, level, n) {
  x = ifelse(var >= level, 1, 0)
  z = rle(x)
  sum(z$values[z$lengths >= n])
}

### Define a function to count events where a variable is below or equal to a level for n hours
COUNT_rl_below = function(var, level, n) {
  x = ifelse(var <= level, 1, 0)
  z = rle(x)
  sum(z$values[z$lengths >= n])
}



### Creating a function to summarize variables in daytime and nighttime -------------------------------------------------------------------------------------
### Daytime function. Here we we don't need to use PERIOD from weather_hour because both 24h, daytime, dusk and dusk are subsets of the same calendar day.
### However, nighttime period collect data from two calendar days. Think about it.

FUNCTION_day = function(dataset, period = "24h") {
  dd = case_when(period == "24h" ~ "24h",
                 period == "daytime" ~ "daytime",
                 period == "dawn" ~ "dawn",
                 period == "dusk" ~ "dusk")
  
  dataset %>% group_by(SITE, DATE) %>%
    summarise(
      nhours._   = n(),
      RH.L40._   = sum(RH.L40),
      RH.G90._   = sum(RH.G90),
      RH.90.rl.count4._ = COUNT_rl_above(var = HMD, level = 90, n=4),
      RH.90.rl.count6._ = COUNT_rl_above(var = HMD, level = 90, n=6),
      RH.30.rl.count4._ = COUNT_rl_below(var = HMD, level = 30, n=4),
      RH.30.rl.count6._ = COUNT_rl_below(var = HMD, level = 30, n=6),
      
      T.7T10._     = sum(T.7T10),      
      T.10T13._  = sum(T.10T13),
      T.13T16._  = sum(T.13T16),  
      T.16T19._  = sum(T.16T19),
      T.19T22._  = sum(T.19T22),
      T.22T25._  = sum(T.22T25),
      T.25T28._  = sum(T.25T28),
      T.G28._  = sum(T.G28),
      
      TRH.16T19nRH.L40._  = sum(TRH.16T19nRH.L40),
      TRH.19T22nRH.L40._  = sum(TRH.19T22nRH.L40),
      TRH.22T25nRH.L40._  = sum(TRH.22T25nRH.L40),
      TRH.25T28nRH.L40._ = sum(TRH.25T28nRH.L40),
      
      TRH.7T10nRH.G85._  = sum(TRH.7T10nRH.G85),
      TRH.10T13nRH.G85._  = sum(TRH.10T13nRH.G85),
      TRH.13T16nRH.G85._ = sum(TRH.13T16nRH.G85),
      TRH.16T19nRH.G85._ = sum(TRH.16T19nRH.G85),
      TRH.19T22nRH.G85._ = sum(TRH.19T22nRH.G85),
      
      TR.3T7nR.G.0.2._  = sum(TR.3T7nR.G.0.2),
      TR.7T10nR.G.0.2._  = sum(TR.7T10nR.G.0.2),
      TR.10T13nR.G.0.2._  = sum(TR.10T13nR.G.0.2),
      TR.13T16nR.G.0.2._  = sum(TR.13T16nR.G.0.2), 
      TR.16T19nR.G.0.2._  = sum(TR.16T19nR.G.0.2),
      TR.19T22nR.G.0.2._  = sum(TR.19T22nR.G.0.2),
      
      R.S._      = sum(RNF),
      R.AH._     = sum(R.EPISODE),
      R.0.2.rl.count6._ = COUNT_rl_above(var = RNF, level = 0.2, n=6),
      R.0.4.rl.count4._ = COUNT_rl_above(var = RNF, level = 0.4, n=4),
      R.0.6.rl.count2._ = COUNT_rl_above(var = RNF, level = 0.6, n=2),
      
      ETo.G0.4._ = sum(ETo.G0.4),
      
      RH8.peak4._ = sum(rle(RH8.peak4)$values!=0),
      T8.peak4._ = sum(rle(T8.peak4)$values!=0)
    ) %>%
    setNames(gsub("_", dd, names(.))) # Change variable names to reflect the time of interest
}


### Define a function to summarize and create new variables for nighttime 

FUNCTION_nighttime = function(dataset, period = "nighttime"){  # The "period" argument is use to change the variables names, so we can combine then we can combine 24h and overnighttime data with no conflict   
  dd = case_when(period == "nighttime"~ "nighttime")
  
  dataset %>% group_by(SITE, PERIOD3) %>% # group by includes period3 instead of day because we want to keep nighttime hours consecutive. If group by DATE, we end up with nonconsecutive nighttime hours.
      summarise(
        nhours._   = n(),
        RH.L40._   = sum(RH.L40),
        RH.G90._   = sum(RH.G90),
        RH.90.rl.count4._ = COUNT_rl_above(var = HMD, level = 90, n=4),
        RH.90.rl.count6._ = COUNT_rl_above(var = HMD, level = 90, n=6),
        RH.30.rl.count4._ = COUNT_rl_below(var = HMD, level = 30, n=4),
        RH.30.rl.count6._ = COUNT_rl_below(var = HMD, level = 30, n=6),
        
        T.7T10._     = sum(T.7T10),      
        T.10T13._  = sum(T.10T13),
        T.13T16._  = sum(T.13T16),  
        T.16T19._  = sum(T.16T19),
        T.19T22._  = sum(T.19T22),
        T.22T25._  = sum(T.22T25),
        T.25T28._  = sum(T.25T28),
        T.G28._  = sum(T.G28),
        
        TRH.16T19nRH.L40._  = sum(TRH.16T19nRH.L40),
        TRH.19T22nRH.L40._  = sum(TRH.19T22nRH.L40),
        TRH.22T25nRH.L40._  = sum(TRH.22T25nRH.L40),
        TRH.25T28nRH.L40._ = sum(TRH.25T28nRH.L40),
        
        TRH.7T10nRH.G85._  = sum(TRH.7T10nRH.G85),
        TRH.10T13nRH.G85._  = sum(TRH.10T13nRH.G85),
        TRH.13T16nRH.G85._ = sum(TRH.13T16nRH.G85),
        TRH.16T19nRH.G85._ = sum(TRH.16T19nRH.G85),
        TRH.19T22nRH.G85._ = sum(TRH.19T22nRH.G85),
        
        TR.3T7nR.G.0.2._  = sum(TR.3T7nR.G.0.2),
        TR.7T10nR.G.0.2._  = sum(TR.7T10nR.G.0.2),
        TR.10T13nR.G.0.2._  = sum(TR.10T13nR.G.0.2),
        TR.13T16nR.G.0.2._  = sum(TR.13T16nR.G.0.2), 
        TR.16T19nR.G.0.2._  = sum(TR.16T19nR.G.0.2),
        TR.19T22nR.G.0.2._  = sum(TR.19T22nR.G.0.2),
        
        R.S._      = sum(RNF),
        R.AH._     = sum(R.EPISODE),
        R.0.2.rl.count6._ = COUNT_rl_above(var = RNF, level = 0.2, n=6),
        R.0.4.rl.count4._ = COUNT_rl_above(var = RNF, level = 0.4, n=4),
        R.0.6.rl.count2._ = COUNT_rl_above(var = RNF, level = 0.6, n=2),
        
        ETo.G0.4._ = sum(ETo.G0.4),
        
        RH8.peak4._ = sum(rle(RH8.peak4)$values!=0),
        T8.peak4._ = sum(rle(T8.peak4)$values!=0)
      )  %>%
    setNames(gsub("_", dd, names(.))) # Change variable names to reflect the time of interest
}



### Applying formula to different intra-day periods -------------------------------------------------------------------------------------------

data_24h  = FUNCTION_day(data_hour, period = "24h") %>% mutate(DOY =yday(DATE)) %>% relocate(DOY,.after=DATE)# applying function above to entire day or 24h
data_day  = data_hour %>% filter(PERIOD1=="daytime") %>% FUNCTION_day(., period = "daytime") %>% mutate(DOY =yday(DATE)) %>% relocate(DOY,.after=DATE)# filtering and applying functions above to hours with sunlight
data_dawn  = data_hour %>% filter(PERIOD2=="dawn") %>% FUNCTION_day(., period = "dawn")%>% mutate(DOY =yday(DATE)) %>% relocate(DOY,.after=DATE) # filtering and applying functions around dawn (3 before and 3 after)
data_dusk  = data_hour %>% filter(PERIOD2=="dusk") %>% FUNCTION_day(., period = "dusk")%>% mutate(DOY =yday(DATE)) %>% relocate(DOY,.after=DATE) # filtering and applying functions around dusk (3 before and 3 after)
data_nighttime = data_hour %>% filter(PERIOD1=="nighttime") %>% FUNCTION_nighttime(., period = "nighttime") # filtering and applying functions above to hours without sunlight

### Merging date to the nighttime set. Period 3 will drop

data_nighttime= left_join(data_nighttime %>% filter(!nhours.nighttime<10),unique(data_hour %>% 
                                 filter(PERIOD1=="nighttime") %>% 
                                 select(PERIOD3,DATE,SITE)) %>%
               distinct(PERIOD3,.keep_all = TRUE), by=c("SITE","PERIOD3")) %>% # period consists of consecutive nighttime hours
  relocate(DATE,.before = PERIOD3) %>%
  select(-PERIOD3) %>%  # relocating merged column and removing period3
  distinct(SITE,DATE,.keep_all = TRUE) %>% 
  mutate(DOY =yday(DATE)) %>% 
  relocate(DOY,.after=DATE)


### Save -------------------------------------------------------------------------------------------
### create a list of data frames
engineered_variables = list("d24h"= data_24h,
                         "daytime"= data_day,
                         "nighttime" = data_nighttime,
                         "dawn" = data_dawn,
                         "dusk" = data_dusk) 

table(data_dusk$SITE)
table(data_24h$SITE)
table(data_dawn$SITE)
table(data_nighttime$SITE)

save(engineered_variables, file = "code/Window pane/data/engineered_variables.RData")


