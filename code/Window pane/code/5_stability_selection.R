
### Load packages -----------------------------------------------------------------------------------------------
library(data.table) # using library() instead of packman::p_loa because of HPC
library(furrr)
library(tidyverse)
library(future)
library(sharp)

### Load data sets
rm(list = ls())
load("code/Window pane/data/rolling_windows.RData")

### Subsetting locations that are going to be used in the CV0 scheme (prediction in unknown environments)
dt= lapply(rolling_windows, function(x) {x %>% filter(!SITE %in% c("OX23","LB24","SB24","CL22"))})


### Set seed
set.seed(123)

# Set up parallel backend using furrr
plan(multisession)


### Variable selection ------------------------------------------------------------------------------------------- 


selection_lasso = function(data,response){
  y_reg = data %>% select({{response}}) 
  x_reg = data %>% select(-c(SITE, DATE, DOY, matches("fa")))
  
  stab_reg = VariableSelection(
    xdata = x_reg,
    ydata = y_reg,
    K = 1000,
    seed = 123,
    resampling = "bootstrap",
    verbose = FALSE)
  
  if(any(is.na(Stable(stab_reg)))){return(NA)}
  else {
    covars= names(Stable(stab_reg))[Stable(stab_reg) != 0]
    dataset= bind_cols(y_reg,x_reg[covars])
    return(dataset)
  }
}

### Stability selection -------------------------------------------------------------------------------------------

compute_stable = function(data) {
  data %>%
    na.omit() %>%
    nest(data = -c(VAR, LAG)) %>%
    filter(map(data, nrow) >= 12) %>%   # at least 12 environments for correlation 
    mutate(
      fa1 = map(data, ~ selection_lasso(.x, fa1)),
      fa2 = map(data, ~ selection_lasso(.x, fa2)),
      fa3 = map(data, ~ selection_lasso(.x, fa3))
    ) %>%
    select(-data)
}

stable_7 = compute_stable(dt[[1]])
stable_14 = compute_stable(dt[[2]])
stable_21 = compute_stable(dt[[3]])
stable_28 = compute_stable(dt[[4]])


plan(sequential) 

### Save -------------------------------------------------------------------------------------------
stable_results = list(stable_7, stable_14, stable_21, stable_28)
lags = c(7, 14, 21, 28)

for (i in seq_along(stable_results)) {
  result = stable_results[[i]]
  lag = lags[i]
  save(result, file = paste0("code/Window pane/data/stable/stable_", lag, ".RData"))
}

