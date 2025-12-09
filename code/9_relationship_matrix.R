#####################################################
####### Matrices of environmental similarity ########
#####################################################

####### Author: Vinicius Garnica
####### Date: Dec, 2024

### Load packages ------------------------------------------------------------------------------------ 
pacman::p_load(data.table,
               Matrix,
               matrixNormal,
               tidyverse)

### Set the path for workflow
rm(list = ls())

### Functions ------------------------------------------------------------------------------------ 

filter_correlated_vars = function(df, threshold = 0.5) {
  cor_matrix = cor(df %>%select_if(is.numeric), use = "pairwise.complete.obs")
  
  correlated_pairs = which(abs(cor_matrix) > threshold & upper.tri(cor_matrix, diag = FALSE), arr.ind = TRUE)
  exclude_cols = unique(c(correlated_pairs[, 1], correlated_pairs[, 2]))

  df_filtered = df %>%
    dplyr::select(-one_of(names(df)[exclude_cols]))
  
  return(df_filtered)
}

heatmap_plot = function(correlation_matrix) {
  heatmap_df = as.data.frame(correlation_matrix)
  heatmap_df$Row = rownames(heatmap_df)
  heatmap_df_long = reshape2::melt(heatmap_df, id.vars = "Row", variable.name = "Column", value.name = "Value")
  
  eclid=hclust(dist(1-correlation_matrix),method = "complete")
  ordered_labels=eclid$labels[eclid$order]
  
  heatmap_df_long=heatmap_df_long %>%
    mutate(Row = factor(Row, levels = ordered_labels),
           Column = factor(Column, levels = ordered_labels))
  
  ggplot(heatmap_df_long, aes(x = Row, y = Column, fill = Value, label = round(Value, 2)))+
    geom_tile(data = subset(heatmap_df_long, Row == Column), fill = "#FF0000") +
    viridis::scale_fill_viridis(option = "A", limits = c(0, 1),direction = -1) +
    geom_tile(color = "white") +
    geom_text(color = "white", size = 2.5,fontface = "bold") +
    scale_y_discrete(limits = ordered_labels)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 11,angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11,angle = 0, hjust = 1),
          legend.position = "bottom",
          strip.background = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size=15),
          strip.text = element_text(size = 15, hjust = 0)) +
    labs(fill = "Correlation")
}


### Load data ------------------------------------------------------------------------------------ 

### All weather variables
load('data/pre_anthesis_variables.RData')

# We are going to build relatedness matrices using a subset of pre-anthesis predictors
# We first subset a few pre-anthesis predictors that we believe have a contribution to SBN dynamics. These predictors are going to be used directly in the model

selected_pre_anthesis_predictors =  c("fa1.77_71.RH.G90.dusk.sum_14",
                                      "fa1.66_58.R.0.6.rl.count2.daytime.sum_28",
                                      "fa1.61_55.T.16T19.daytime.sum_28",  
                                      "fa1.47_38.TRH.7T10nRH.G85.dawn.sum_14", 
                                      "fa1.20_6.TR.19T22nR.G.0.2.daytime.sum_28", 
                                      "fa1.14_7.TR.13T16nR.G.0.2.dawn.sum_21") 

pre_anthesis_variables %>% 
  filter(!SITE %in% c("OX23","LB24","SB24","CL22")) %>%
  dplyr::select(selected_pre_anthesis_predictors)

#The rest of FA1 meteorological predictors will be used for relatedness matrix
weather_df = pre_anthesis_variables %>% 
  filter(!SITE %in% c("OX23","LB24","SB24","CL22")) %>%
  dplyr::select(-any_of(selected_pre_anthesis_predictors)) %>%
  dplyr::select(starts_with(c("fa3","fa2","fa1")), SITE)   

### Apply a filter to remove highly correlated variables and too many zeros ------------------------------------------
zero_proportions = weather_df %>%
  summarise(across(where(is.numeric), ~ mean(. == 0,na.rm=T)))

matrix_data = weather_df %>%
  dplyr::select(SITE,intersect(names(zero_proportions), 
                              names(zero_proportions)[unlist(zero_proportions) <= 0.5])) %>%
  filter_correlated_vars(threshold = 0.72) 


# save predictors list
matrix_predictors = names(matrix_data);matrix_predictors
write.csv(matrix_predictors,"data/matrix_predictors.csv")

# Centers the data, and scales by standard deviation
matrix_data_scaled=scale(matrix_data[,-1], center=TRUE, scale=TRUE)


### Matrix of Euclidean distance ------------------------------------------------------------------------------------ 
# Following Thomson et al. (2018) 
### Scales so that the values are between 0 and 1
K_matrix = as.matrix(dist(matrix_data_scaled, method = "euclidean", diag = TRUE, upper = TRUE))
euclimat = 1- K_matrix/max(K_matrix, na.rm=TRUE)
rownames(euclimat) = colnames(euclimat) = unique(weather_df$SITE)
class(euclimat)
is.symmetric.matrix(euclimat)
is.positive.definite(euclimat)

heatmap_plot(euclimat)


### Save
save(euclimat, file = "data/euclimat.RData")


### Variance-covariance matrix ------------------------------------------------------------------------------------ 
# Similar to Jarquín et al. (2014) 
#varmat = as.matrix(var(t(matrix_data_scaled), na.rm=TRUE))
#rownames(varmat) = colnames(varmat) = unique(weather_df$SITE)
#heatmap(varmat)
#is.symmetric.matrix(varmat)
#is.positive.definite(varmat)

# Make sure matrix positive definite
#varmat = nearPD(varmat)$mat

### Save
#save(varmat, file = "data/varmat.RData")

matrix_predictors

### Reference 
#Jarquín, Diego, José Crossa, Xavier Lacaze, Philippe Du Cheyron, Joëlle Daucourt, Josiane Lorgeou, François Piraux, et al. 2014. 
# A Reaction Norm Model for Genomic Selection Using High-Dimensional Genomic and Environmental Data.” Theoretical and Applied Genetics 127 (3): 595–607. doi:10.1007/s00122-013-2243-1.

#Thomson, Caroline Elizabeth, Isabel Sophie Winney, Oceane Salles, and Benoit Pujol. 2018. 
# A Guide to Using a Multiple-Matrix Animal Model to Disentangle Genetic and Nongenetic Causes of Phenotypic Variance. bioRxiv, May, 318451. doi:10.1101/318451.