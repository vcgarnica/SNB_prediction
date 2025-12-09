################################
######## Window pane  ##########
################################

####### Authors: Vinicius Garnica
####### Date: Apr, 2025

### Load packages -----------------------------------------------------------------------------------------------
pacman::p_load(data.table,
               tidyverse,
               slider,
               doParallel,
               purrr)

### Load data sets
rm(list = ls())
load("code/Window pane/data/engineered_variables.RData")

# Register the parallel backend -----------------------------------------------------------------------------------------------
cl = makeCluster(detectCores()-1)
registerDoParallel(cl)


### Rolling average and sum functions -----------------------------------------------------------------------------------------------
# Function to compute rolling sum for specified columns and window size
slide_sum = function(data, wind){
  data %>%
    mutate(across(-matches("SITE|DATE|DOY|nhours"),
                  ~ slide_dbl(.x, ~sum(.x), 
                              .before = wind - 1, 
                              .complete = TRUE),  
                  .names = "{.col}.sum_")) %>%
    select(DATE, DOY, matches("sum_")) %>%
    rename_with(~ paste0(., wind, sep = ""), -c("DATE","DOY"))
}

### Preparing data sets--------------------------------------------------------------------------------

# Function to compute rolling mean and sum for different window sizes
wind_function = function(DF, wind) {
  res = DF %>%
    ungroup() %>%
    group_by(SITE) %>%
    nest() %>%
    summarise(BLOCK = map(data, ~slide_sum(.x, wind))) %>%
    unnest(data = ., cols = BLOCK)
  return(res)
}


### Length of windows from 7 to 28. Users should define this windows based on expert knowledge about the system under investigation.
wind_size = seq(7, 28, 7)

### Loops for each variable and window size --------------------------------------------------------------------------------
output = vector('list')


for (i in names(engineered_variables)) {
  output[[i]] = map(wind_size, ~ wind_function(engineered_variables[[i]], .x))
}

### Stop the parallel backend
stopCluster(cl)

### Since the output is a list, we clean and omit NAs. Why are there NAs?
### Note that the first 4 values resulting from the computation of the rolling average with a 5 day window will be NAs. That makes sense.
### Similarly, the first 9 values of the rolling average (or sum) for a 10 day window will be NAs. 

remove_nas = function(df) {
  df[complete.cases(df), ]
}

fixed_window = map(output , ~ map(.x, remove_nas)) %>% 
  purrr::flatten() 

### Quality check --------------------------------------------------------------------------------

#fixed_window[[20]] %>% # change data sets by changing the number from 1 to 30. There should be one only observation per day
#  group_by(DATE, SITE) %>% 
#  summarise(n=n()) %>%
#  ggplot(aes(x = day(DATE), y = month(DATE), fill = n)) +
#  geom_tile() +
#  labs(x = NULL) +
#  facet_wrap(~SITE) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

### Save --------------------------------------------------------------------------------
save(fixed_window, file = "code/Window pane/data/fixed_window.RData")
