#######################
####### Priors ########
#######################

####### Author: Vinicius Garnica
####### Date: July, 2024

### Load packages ------------------------------------------------------------------------------------ 
pacman::p_load(metafor,
               data.table,
               tidyverse)

### Set the path for workflow
rm(list = ls())


### Wheat Residue Effect at 20% soil coverage ----------------------------------------------------------
# L. K. Mehra, C. Cowger, R. Weisz, and P. S. Ojiambo
# Quantifying the Effects of Wheat Residue on Severity of Stagonospora nodorum Blotch and Yield in Winter Wheat. 
# Phytopathology, 2015 105:11, 1417-1426
# Data from tables 5 and 4

trial = c(1:8)
untreated=c(29.6,15,23.4,25.4,19.2,19.8,8,19.4)
residue20=c(39,28,26,29.4,24.2,29.2,24.2,24.2)
se= c(1.61,1.38,0.45,1.31,0.90,0.96,1.37,0.69)
n= rep(5,8)
data=tibble(trial,untreated,residue20,se,n)
data$sd = se*sqrt(n)

dat = escalc(measure="MD",
             m1i = residue20,
             m2i = untreated,
             sd1i = sd,
             sd2i = sd,
             n1i = n,
             n2i = n,
             data=data)

res1 = rma(yi, vi, data=dat)

### Size effect and standard deviation
c(res1$b,1.6787*sqrt(8))


# Kaur 2023

### Wheat Residue Effect at 20% soil coverage ----------------------------------------------------------
# Navjot Kaur, Hillary L. Mehl, David Langston, and David C. Haak
# Evaluation of Stagonospora Nodorum Blotch Severity and Parastagonospora nodorum Population Structure and Genetic Diversity Across Multiple Locations and Wheat Varieties in Virginia
# Phytopathology 2024 114:1, 258-268
# Data from table 1


bullet=c(20,5.3,4.3,1,6.3,1.8,30,8.7,7,4.5,37.5,13.3,36.7,3.5)
sem_bullet=c(7.6,1.6,2.9,0.8,0.8,0.7,2.9,1.3,1.5,1,0,1.7,4.4,0.8)
hiliard = c(3.5,1.2,1.3,0.5,6.2,0.4,20,7.3,8.5,1.3,5.3,8.7,25,0.5)
sem_hiliard=c(1.5,0.7,0.9,0,0.9,0.2,5.8,1.3,1.5,0.8,2.2,1.3,2.9,0)
MAS_61 = c(19.3,1.2,13.3,0.2,3.7,0.5,19.2,6,5.5,1.2,18.3,13.3,30,0.5)
sem_MAS_61 =c(7.2,0.7,8.5,0.2,2,0,3.6,0,0,0.7,7.1,1.7,7.6,0)
shirley=c(2.7,0.9,0.5,0.4,2,0.5,13.3,7.3,5.5,1.2,5,8.7,19.3,1.2)
sem_shirley = c(1.3,0.8,0,0.2,0.8,0,1.7,1.3,0,0.7,0,1.3,1.7,0.7)
USG3316 = c(9.5,1.2,5.8,1.3,1.3,1.2,25,7.3,8.6,3,29.2,10,40,0.5)
sem_USG3316 = c(1,0.7,1.7,0.8,0.8,0.7,0,1.3,3.2,1.4,4.4,0,7.6,0)
n= rep(3,14)

data=tibble(bullet,sem_bullet,hiliard,sem_hiliard,MAS_61,sem_MAS_61,shirley, sem_shirley,USG3316,sem_USG3316,n)


calculate_combined_stats = function(row) {

  values = as.numeric(row[c("bullet", "hiliard", "MAS_61", "shirley")])
  sems = as.numeric(row[c("sem_bullet", "sem_hiliard", "sem_MAS_61", "sem_shirley")])
  sample_sizes = rep(3, length(values))
  
  combined_mean = mean(values)
  
  variances = (sems^2) * sample_sizes
  within_group_var = sum((sample_sizes - 1) * variances)
  between_group_var = sum(sample_sizes * (values - combined_mean)^2)
  combined_variance = (within_group_var + between_group_var) / (sum(sample_sizes) - 1)
  
  combined_sem = sqrt(combined_variance / sum(sample_sizes))
  
  combined_sd = combined_sem * sqrt(sum(sample_sizes))
  
  return(c(combined_mean, combined_sem, combined_sd))
}


results_MR = t(apply(data, 1, calculate_combined_stats))
colnames(results_MR) = c("MR_mean", "MR_sem", "MR_sd")
data = bind_cols(data, as_tibble(results_MR))
data$MS_mean = data$USG3316
data$MS_sem = data$sem_USG3316
data$MS_sd = data$MS_sem*sqrt(3)

mean(data$MR_mean)
mean(data$MS_mean)

mean(data$MR_sd)
mean(data$MS_sd)
