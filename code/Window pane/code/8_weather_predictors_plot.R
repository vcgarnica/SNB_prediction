######################################
######### Association plots ########## 
######################################

####### Author: Vinicius Garnica
####### Date: Sep, 2024


### Load packages ------------------------------------------------------------------------------
pacman::p_load(purrr,
               tidyverse,
               gridExtra)


theme_set(theme_bw())
### Load data sets ------------------------------------------------------------------------------ 
rm(list = ls())
load("data/pre_anthesis_variables.RData")
load("data/SNB.RData")

### Left join function ------------------------------------------------------------------------------ 

SNB_data = SNB


# Function to create plots ------------------------------------------------------------------------------ 
plot_predictor = function(predictor) {
  a=SNB_data %>% group_by(site) %>% summarise(sev=mean(sev)) %>%
    left_join(.,pre_anthesis_variables , by = c("site" = "SITE")) %>% 
    filter(!site %in% c("OX23","LB24","SB24","CL22")) %>% select(-site) 
  
  ggplot(a ,aes_string(x = predictor, y = "sev")) +
    geom_point() +
    geom_smooth(method="lm", se = T) +
    labs(title = predictor,
         x = "Cumulative hours or events over optimal period",
         y = "Average SNB severity (%)")+
    theme(plot.title = element_text(face = "italic",size = 11))
}

# Function to save plots in batches
save_plots_in_batches = function(predictors, batch_size = 20) {
  num_batches = ceiling(length(predictors) / batch_size)
  for (i in 1:num_batches) {
    start_index = (i - 1) * batch_size + 1
    end_index = min(i * batch_size, length(predictors))
    predictors_batch = predictors[start_index:end_index]
    
    plots = map(predictors_batch, plot_predictor)
    
    filename = paste0("code/Window pane/figures/weather_predictors_batch_", i, ".png")
    g = arrangeGrob(grobs = plots, ncol = 4)
    ggsave(filename, g, width = 16, height = 20)
  }
}

### Filter numeric predictors ------------------------------------------------------------------------------ 
numeric_predictors = pre_anthesis_variables %>%
  select(where(is.numeric)) %>%
  names()


# Save the plots ------------------------------------------------------------------------------ 
save_plots_in_batches(numeric_predictors)
