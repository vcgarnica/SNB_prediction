#############################
######### Heatmaps ########## 
#############################

####### Author: Vinicius Garnica
####### Date: Sep, 2024

### Load packages ------------------------------------------------------------------------------
pacman::p_load(data.table,
               tidyverse,
               ggthemes,
               ggpmisc,
               cowplot,
               officer,
               flextable,
               gridExtra,
               patchwork)


theme_set(theme_bw())



### Load data sets ------------------------------------------------------------------------------ 
rm(list = ls())
load("code/Window pane/results/result_7.RData")
load("code/Window pane/results/result_14.RData")
load("code/Window pane/results/result_21.RData")
load("code/Window pane/results/result_28.RData")

    
### Filtering data sets by loading factor ------------------------------------------------------------------------------ 

data.fa1 = rbind(result_7,result_14,result_21,result_28) %>% ungroup() %>%
  filter(var1=="fa1" & VAR=="sev")     

data.fa2 = rbind(result_7,result_14,result_21,result_28) %>%ungroup() %>%
  filter(var1=="fa2" & VAR=="sev" )     

data.fa3 = rbind(result_7,result_14,result_21,result_28) %>%ungroup() %>%
  filter(var1=="fa3" & VAR=="sev")       

### Functions ------------------------------------------------------------------------------

process_data = function(data) {
  data = data %>%
    mutate(
      significant = case_when(
        mean < 0 & upper_95 < 0 ~ TRUE,
        mean > 0 & lower_95 > 0 ~ TRUE,
        mean >= 0 & lower_95 <= 0 ~ FALSE,
        mean <= 0 & upper_95 >= 0 ~ FALSE
      ),
      correlation = case_when(
        mean <= -0.6 & significant == TRUE ~ "Strongly negative",
        mean > -0.6 & mean <= -0.2  ~ "Moderately negative",
        mean > 0.2 & mean < 0.6  ~ "Moderately positive",
        mean >= 0.6 & significant == TRUE ~ "Strongly positive",
        TRUE ~ "Non-significant or weak"
      ),
      Support = ifelse(
        correlation != "Non-significant or weak",
        "Drop",
        "Stay"
      )
    )
  
  # Number of distinct variables in var2 (assuming var2 is a column in the data)
  n_distinct_var2 = n_distinct(na.omit(data$var2))
  
  # Max and min correlation values
  correlation_summary = data %>%
    ungroup() %>%
    summarise(
      max_mean = max(mean, na.rm = TRUE),
      min_mean = min(mean, na.rm = TRUE)
    )
  
  # Reorder the correlation levels
  data$correlation = factor(
    data$correlation,
    levels = rev(c("Strongly negative", "Moderately negative", "Non-significant or weak", "Moderately positive", "Strongly positive"))
  )
  
  # Returning the modified data and the additional information as a list
  list(
    modified_data = data,
    n_distinct_var2 = n_distinct_var2,
    correlation_summary = correlation_summary
  )
}

### Filtering all weather elements that display consecutive relationship with outcomes 

filter_data = function(data,var1) {
  
  results = data %>%
    complete(var2, LAG = seq(-40, 100), fill = list(var1 = {{var1}}, significant = FALSE)) %>% # fill NA for LAG observations that are non-significant
    arrange(var2,LAG) %>%
    group_by(var2) %>%
    mutate(consecutive = rleid(Support)) %>%
    group_by(var1, var2,consecutive) %>%
    arrange(var2,LAG) %>%
    mutate(c_consecutive = n()) %>%
    filter(c_consecutive>=7 & correlation %in% c("Moderately negative","Strongly negative", "Strongly positive","Moderately positive")) %>%
    ungroup() %>%
    select(-consecutive) %>%
    distinct(var1, var2)
  
  ### Combine
  filtered_data = left_join(results,data)
  
  # Number of distinct variables in var2 (assuming var2 is a column in the data)
  n_distinct_var2 = n_distinct(na.omit(filtered_data$var2))
  
  # Max and min correlation values
  correlation_summary = filtered_data %>%
    ungroup() %>%
    summarise(
      max_mean = max(mean, na.rm = TRUE),
      min_mean = min(mean, na.rm = TRUE)
    )
  
  # Returning the modified data and the additional information as a list
  list(
    filtered_data = filtered_data,
    n_distinct_var2 = n_distinct_var2,
    correlation_summary = correlation_summary
  )
}

### Heatmap function

heatmap_plot = function(data, var1, palette) {
  
  # Returning the time series boundaries
  lag_max = max(data$LAG)
  lag_min = min(data$LAG)
  
  # Data wrangling
  a = data %>%
    filter(var1=={{var1}}) %>%
    ungroup() %>%
    complete(var2, LAG = seq(lag_min-4,lag_max+5, by = 1), var1, fill = list(correlation = "Non-significant or weak")) %>%  # filling the dataset that was removed during stability selection with non-significant values
    group_by(var2) %>%
    mutate(consecutive_count = max(1, cumsum(correlation %in% c("Strongly positive", "Moderately positive")) + 1),
           max_consecutive = max(consecutive_count)*LAG) %>%
    ungroup() %>%
    mutate(
      var1 = case_when(var1=="fa1"~as.character("hat(λ)[1]"),
                       var1=="fa2"~as.character("hat(λ)[2]"),
                       var1=="fa3"~as.character("hat(λ)[3]"))) 
  
  p = ggplot(a, aes(x = LAG, y=reorder(var2,max_consecutive), fill = mean)) +
    geom_tile(color = "white", linewidth = 0.2) +
    labs(y = NULL,x = "Days before anthesis",fill="Mean correlation")+
    coord_equal(xlim = c(lag_max, lag_min)) +
    facet_grid(col = vars(var1), labeller = label_parsed) + 
    scale_fill_viridis_c(option = "plasma",direction = -1,na.value = palette)+
    scale_x_reverse(breaks = seq(100, -40, by = -10)) +
    theme(legend.position = "bottom", 
          axis.text.y = element_text(size=8),
          strip.text = element_text(hjust = 0),
          axis.title.x = element_text(size=10,margin = margin(t = 10)),
          strip.text.x = element_text(size=16),
          plot.title = element_text(size=19),
          strip.background = element_blank())
  
  return(p)
}

fill_gaps <- function(column) {
  for (i in 2:(length(column) - 1)) {
    if (is.na(column[i]) && !is.na(column[i - 1]) && !is.na(column[i + 1])) {
      column[i] <- 0.5
    } else if (is.na(column[i]) && is.na(column[i + 1]) && !is.na(column[i + 2]) && !is.na(column[i - 1])) {
      column[i] <- 0.5
      column[i + 1] <- 0.5
    } else if (is.na(column[i]) && is.na(column[i + 1]) && is.na(column[i + 2]) && !is.na(column[i + 3]) && !is.na(column[i - 1])) {
      column[i] <- 0.5
      column[i + 1] <- 0.5
      column[i + 2] <- 0.5
    }
  }
  return(column)
}

### Calculate area function

calculate_area = function(subset_data, sev_data, fa) {
  var_sym = sym(unique(subset_data$var2))
  lag_input = subset_data$LAG
  
  column_input = sev_data %>%
    select(SITE, DATE, LAG, !!var_sym)
  
  row_input = column_input %>%
    filter(LAG %in% lag_input)
  
  row_input %>%
    group_by(SITE) %>%
    summarise(area = sum(!!var_sym)) %>%
    rename(!!paste0(fa, max(subset_data$LAG), '_', min(subset_data$LAG), '.', var_sym) := area)
}


### Summarize function

summ_function = function(data){
  res = data %>%
    dplyr::summarise(across(where(is.numeric), list(
      min = ~round(min(., na.rm = TRUE),1),
      q1 = ~round(quantile(., 0.25, na.rm = TRUE),1),
      mean = ~round(mean(., na.rm = TRUE),1),
      q3 = ~round(quantile(., 0.75, na.rm = TRUE),1),
      max = ~round(max(., na.rm = TRUE),1)
    ))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "variable_statistic",
      values_to = "value"
    ) %>%
    mutate(variable_statistic = sub("_(?!.*_)", "_DELIMITER_", variable_statistic, perl = TRUE)) %>%
    separate(variable_statistic, into = c("variable", "statistic"), sep = "_DELIMITER_") %>%
    pivot_wider(names_from = statistic, values_from = value)
  return(res)
}


### Summarize function

export_word_table = function(data, file_path) {
  
  doc = read_docx()
  
  doc = doc %>%
    body_add_flextable(value = flextable(data))
  
  print(doc, target = file_path)
}



### Descriptive analysis  ------------------------------------------------------------------------------

process_data(rbind(result_7,result_14,result_21,result_28))


#### FA1  ---------------------------------------------------------------------------------------

# Data wrangling
fa1 = process_data(data.fa1)
fa1.filtered = filter_data(fa1$modified_data,"fa1")

### Summaries
fa1.filtered$correlation_summary

### Plots
heatmap_plot(fa1.filtered$filtered_data, "fa1","#21918c") 
ggsave("code/Window pane/figures/figure2.tiff", width =20,height =9,units = "in",limitsize = FALSE)

#### FA2  ---------------------------------------------------------------------------------------

### Data wrangling
fa2 = process_data(data.fa2)
fa2.filtered = filter_data(fa2$modified_data,"fa2")

### Summaries
fa2.filtered$correlation_summary

### Plots
heatmap_plot(fa2.filtered$filtered_data, "fa2","#bc3754") 
ggsave("code/Window pane/figures/figure3.tiff", width =20,height =9,units = "in",limitsize = FALSE)

#### FA3  ---------------------------------------------------------------------------------------

### Data wrangling
fa3 = process_data(data.fa3)
fa3.filtered = filter_data(fa3$modified_data,"fa3")

### Summaries
fa3.filtered$correlation_summary

### Plots
heatmap_plot(fa3.filtered$filtered_data, "fa3","deeppink3") 
ggsave("code/Window pane/figures/figure4.tiff", width =20,height =9,units = "in",limitsize = FALSE)


######################################
######### Creating predictors ######## 
######################################


### Load window pane data ---------------------------------------------------------------------------------------

load('code/Window pane/data/rolling_windows.RData')

sev_data=rolling_windows %>%
  reduce(function(x, y) full_join(x, y, by = c('SITE','DATE','DOY', 'LAG',"fa1","fa2","fa3")))  


site_df = data.frame(SITE = factor(unique(sev_data$SITE)))

### FA1 second-level feature engineering ---------------------------------------------------------------------------------------

# Step 1: Filter and process consistent variables
consistent_fa1_vars=
data.fa1 %>%
  complete(var2, LAG = seq(-40, 100), fill = list(var1 = "fa1", significant = FALSE)) %>% # fill NA for LAG observations that are non-significant
  mutate(
    correlation = case_when(
      mean <= -0.2 ~ "Negative",
      mean >= 0.2 ~ "Positive"                    
    )
  ) %>%
  arrange(var2,LAG) %>%
  filter(LAG>=0) %>%
  mutate(Support = ifelse(correlation %in% c("Positive","Negative"),"Stay","Drop"),
         consecutive = rleid(Support)) %>%
  filter(Support!="Drop") %>%
  group_by(var1, var2,consecutive) %>%
  mutate(c_consecutive = n()) %>%
  filter(c_consecutive>=7) %>%
  arrange(var2,LAG) %>%
  ungroup() %>%
  distinct(var2, LAG)

# Step 2: Split data based on LAG differences
split_data_fa1 =consistent_fa1_vars %>%
  group_by(var2) %>%
  group_split() %>%
  map(function(df) {
    if (any(abs(diff(df$LAG)) > 4)) {
      split_points = which(abs(diff(df$LAG)) >4)
      split1 = df[1:split_points[1], ]
      split2 = df[(split_points[1] + 1):nrow(df), ]
      list(split1, split2)
    } else {
      list(df)
    }
  }) %>%
  flatten() 


# Step 3: Calculate area and Reduce dataset
fa1_df = reduce(map(split_data_fa1, ~ calculate_area(.x, sev_data,"fa1.")), left_join, by = "SITE") %>%
  left_join(site_df, by = "SITE")  %>%  dplyr::select(where(~ all(!is.na(.)))); fa1_df


fa1_library = split_data_fa1
save(fa1_library, file = "data/fa1_library.RData")

### FA2 second-level feature engineering ---------------------------------------------------------------------------------------

# Step 1: Filter and process consistent variables
consistent_fa2_vars=
  data.fa2 %>%
  complete(var2, LAG = seq(-40, 100), fill = list(var1 = "fa2", significant = FALSE)) %>% # fill NA for LAG observations that are non-significant
  mutate(
    correlation = case_when(
        mean <= -0.2 ~ "Negative",
        mean >= 0.2 ~ "Positive"                    
      )
    ) %>%
      arrange(var2,LAG) %>%
      filter(LAG>=0) %>%
      mutate(Support = ifelse(correlation %in% c("Positive","Negative"),"Stay","Drop"),
             consecutive = rleid(Support)) %>%
      filter(Support!="Drop") %>%
      group_by(var1, var2,consecutive) %>%
      mutate(c_consecutive = n()) %>%
      filter(c_consecutive>=7) %>%
  arrange(var2,LAG) %>%
  ungroup() %>%
  distinct(var2, LAG)




# Step 2: Split data based on LAG differences
split_data_fa2 =consistent_fa2_vars %>%
  group_by(var2) %>%
  group_split() %>%
  map(function(df) {
    if (any(abs(diff(df$LAG)) > 4)) {
      split_points = which(abs(diff(df$LAG)) >4)
      split1 = df[1:split_points[1], ]
      split2 = df[(split_points[1] + 1):nrow(df), ]
      list(split1, split2)
    } else {
      list(df)
    }
  }) %>%
  flatten() 

fa2_library = split_data_fa2
save(fa2_library, file = "data/fa2_library.RData")

# Step 3: Calculate area and Reduce dataset
fa2_df = reduce(map(split_data_fa2, ~ calculate_area(.x, sev_data,"fa2.")), left_join, by = "SITE") %>%
  left_join(site_df, by = "SITE") %>%  dplyr::select(where(~ all(!is.na(.))));fa2_df



### FA3 second-level feature engineering ---------------------------------------------------------------------------------------


# Step 1: Filter and process consistent variables
consistent_fa3_vars=
  data.fa3 %>%
  complete(var2, LAG = seq(-40, 100), fill = list(var1 = "fa3", significant = FALSE)) %>% # fill NA for LAG observations that are non-significant
  mutate(
    correlation = case_when(
      mean <= -0.2 ~ "Negative",
      mean >= 0.2 ~ "Positive"                    
    )
  ) %>%
  arrange(var2,LAG) %>%
  filter(LAG>=0) %>%
  mutate(Support = ifelse(correlation %in% c("Positive","Negative"),"Stay","Drop"),
         consecutive = rleid(Support)) %>%
  filter(Support!="Drop") %>%
  group_by(var1, var2,consecutive) %>%
  mutate(c_consecutive = n()) %>%
  filter(c_consecutive>=7) %>%
  arrange(var2,LAG) %>%
  ungroup() %>%
  distinct(var2, LAG)



# Step 2: Split data based on LAG differences
split_data_fa3 =consistent_fa3_vars %>%
  group_by(var2) %>%
  group_split() %>%
  map(function(df) {
    if (any(abs(diff(df$LAG)) > 4)) {
      split_points = which(abs(diff(df$LAG)) >4)
      split1 = df[1:split_points[1], ]
      split2 = df[(split_points[1] + 1):nrow(df), ]
      list(split1, split2)
    } else {
      list(df)
    }
  }) %>%
  flatten() 

fa3_library = split_data_fa3
save(fa3_library, file = "data/fa3_library.RData")

# Step 3: Calculate area and Reduce dataset
fa3_df = reduce(map(split_data_fa3, ~ calculate_area(.x, sev_data,"fa3.")), left_join, by = "SITE") %>%
  left_join(site_df, by = "SITE") %>%  dplyr::select(where(~ all(!is.na(.))));fa3_df


## Data wrangling -----------------------------------------------------------------------

# List of data frames

# Reduce dataset
pre_anthesis_variables = reduce(list(fa1_df,
                                     fa2_df,
                                     fa3_df
                                     ), left_join, by = "SITE");pre_anthesis_variables 

names(pre_anthesis_variables)

save(pre_anthesis_variables, file = "data/pre_anthesis_variables.RData")


## Export -----------------------------------------------------------------------

fa1_result = split_data_fa1 %>%
  map(~ .x %>%
        group_by(var2) %>%
        summarise(max_LAG = max(LAG), 
                  min_LAG = min(LAG),
                  Duration =max_LAG-min_LAG +1)) %>%
  do.call(rbind,.) %>%
  mutate(second_level = paste0("fa1.", max_LAG, "_", min_LAG, ".", var2)) %>%
  left_join(summ_function(fa1_df), by = c("second_level" = "variable")) %>%
#  select(-second_level) %>%
  arrange(desc(max_LAG)) %>% 
  na.omit()


export_word_table(fa1_result, "code/Window pane/results/fa1_variables.docx")


fa2_result = split_data_fa2 %>%
  map(~ .x %>%
        group_by(var2) %>%
        summarise(max_LAG = max(LAG), 
                  min_LAG = min(LAG),
                  Duration =max_LAG-min_LAG +1)) %>%
  do.call(rbind,.) %>%
  mutate(second_level = paste0("fa2.", max_LAG, "_", min_LAG, ".", var2)) %>%
  left_join(summ_function(fa2_df), by = c("second_level" = "variable")) %>%
#  select(-second_level) %>%
  arrange(desc(max_LAG))%>% 
  na.omit()


export_word_table(fa2_result, "code/Window pane/results/fa2_variables.docx")

fa3_result = split_data_fa3 %>%
  map(~ .x %>%
        group_by(var2) %>%
        summarise(max_LAG = max(LAG), 
                  min_LAG = min(LAG),
                  Duration =max_LAG-min_LAG +1)) %>%
  do.call(rbind,.) %>%
  mutate(second_level = paste0("fa3.", max_LAG, "_", min_LAG, ".", var2)) %>%
  left_join(summ_function(fa3_df), by = c("second_level" = "variable")) %>%
#  select(-second_level) %>%
  arrange(desc(max_LAG))%>% 
  na.omit()


export_word_table(fa3_result, "code/Window pane/results/fa3_variables.docx")
