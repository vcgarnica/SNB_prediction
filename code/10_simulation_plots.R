##########################
######### Plots ########## 
##########################

####### Author: Vinicius Garnica
####### Date: Mar, 2025

### Load packages ------------------------------------------------------------------------------
pacman::p_load(data.table,
               tidyverse,
               ggthemes,
               cowplot,
               gridExtra,
               patchwork)


theme_set(theme_minimal())
rm(list = ls())

### Load data sets ------------------------------------------------------------------------------
load("simulation/sev_c.RData")
load("simulation/sev_u.RData")


### Data wrangling ------------------------------------------------------------------------------
### Check frequency of inverse signs
sev_u %>%
  mutate(sign_diff = ifelse(sign(mean) == sign(true), "Equal", "Different")) %>%
  group_by(variable,m) %>%
  summarize(equal_sign_freq = sum(sign_diff == "Equal") / n(),
            diff_sign_freq = sum(sign_diff == "Different") / n())

sev_c %>%
  mutate(sign_diff = ifelse(sign(mean) == sign(true), "Equal", "Different")) %>%
  group_by(variable,m) %>%
  summarize(equal_sign_freq = sum(sign_diff == "Equal") / n(),
            diff_sign_freq = sum(sign_diff == "Different") / n())

unique(sev_c$variable)


### Line plots  ------------------------------------------------------------------------------

data = rbind(sev_c, sev_u) %>%
  mutate(
    cov = case_when(
      cov == "uncorrelated" ~ as.character("I[j]^'*'"),
      cov == "correlated" ~ as.character("K[j[obs]]^'*'")
    ),
    variable = case_when(
      variable == "b_Intercept" ~ as.character("β[Intercept]^'*'"),
      variable == "b_residueY" ~ as.character("β[RESIDUE[Y]]^'*'"),
      variable == "b_resistanceMS" ~ as.character("β[RESISTANCE[MS]]^'*'"),
      variable == "b_resistanceS" ~ as.character("β[RESISTANCE[S]]^'*'"),
      variable == "b_onsetPRE" ~ as.character("β[ONSET[PRE]]^'*'"),
      variable == "b_w1" ~ as.character("κ[fa1.77_71.RH.G90.dusk.sum_14]^'*'"),
      variable == "b_w2" ~ as.character("κ[fa1.66_58.R.0.6.rl.count2.daytime.sum_28]^'*'"),
      variable == "b_w3" ~ as.character("κ[fa1.61_55.T.16T19.daytime.sum_28]^'*'"),
      variable == "b_w4" ~ as.character("κ[fa1.47_38.TRH.7T10nRH.G85.dawn.sum_14]^'*'"),
      variable == "b_w5" ~ as.character("κ[fa1.20_6.TR.19T22nR.G.0.2.daytime.sum_28]^'*'"),
      variable == "b_w6" ~ as.character("κ[fa1.14_7.TR.13T16nR.G.0.2.dawn.sum_21]^'*'"),
      variable == "phi" ~ as.character("ϕ^'*'"),
      variable == "sd_site__Intercept" ~ as.character("σ[α]^'*'")
    ),
    k= case_when(
      k == 0.3 ~ as.character("σ[α] == 0.3"),
      k == 0.9 ~ as.character("σ[α] == 0.9")
    )
  ) 

data %>% 
  filter(str_detect(variable, "β") | str_detect(variable, "σ" ) | str_detect(variable,"ϕ")) %>%
  ggplot(aes(x = m, y = mean, ymin = q5, ymax = q95, color = factor(cov))) +
  geom_line(linewidth =1,alpha=0.7) +
  geom_point(position = position_dodge(width = 0),size =2,alpha=0.5) +
  geom_hline(aes(yintercept = true), linetype = "dashed", color = "black",alpha=0.7) +
  geom_errorbar(position = position_dodge(width = 0), width = 0.5,alpha=0.5) +
  facet_grid(variable~k, labeller = label_parsed,scales = "free") +
  scale_color_viridis_d(begin = 0.7, end = 0.1, labels = c(expression(italic(I[j]^"*")), expression(italic(K[j[obs]]^"*"))))+
    labs(x = expression("Number of grouping levels (" * italic(j) * ")"), 
         y = "Parameter estimate",
         color = "Relatedness matrix")+
    scale_x_continuous(breaks = c(10, 25, 50, 100, 150,250)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size=6),
          strip.text.y.right =  element_text(size = 5.5,face = "italic")) 
  

ggsave("results/figures/figS4.tiff",width =5,height = 9,units = "in",limitsize = FALSE)


data %>%
  filter(!(str_detect(variable, "β") | str_detect(variable, "σ" ) | str_detect(variable,"ϕ"))) %>%
  ggplot(aes(x = m, y = mean, ymin = q5, ymax = q95, color = factor(cov))) +
  geom_line(linewidth =1,alpha=0.7) +
  geom_point(position = position_dodge(width = 0),size =2,alpha=0.5) +
  geom_hline(aes(yintercept = true), linetype = "dashed", color = "black",alpha=0.7) +
  geom_errorbar(position = position_dodge(width = 0), width = 0.5,alpha=0.5) +
  facet_grid(variable~k, labeller = label_parsed,scales = "free") +
  scale_color_viridis_d(begin = 0.7, end = 0.1, labels = c(expression(italic(I[j]^"*")), expression(italic(K[j[obs]]^"*"))))+
  labs(x = expression("Number of grouping levels (" * italic(j) * ")"), 
       y = "Parameter estimate",
       color = "Relatedness matrix")+
  scale_x_continuous(breaks = c(10, 25, 50, 100, 150,250)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=6),
        strip.text.y.right =  element_text(size = 5.5,face = "italic")) 


ggsave("results/figures/figS5.tiff",width =5,height = 9,units = "in",limitsize = FALSE)


data %>%
  filter(cov == "I[j]^'*'") %>%
  mutate(sign_diff = ifelse(sign(mean) == sign(true), "Equal", "Different")) %>%
  group_by(variable) %>%
  summarize(equal_sign_freq = sum(sign_diff == "Equal") / n(),
            diff_sign_freq = sum(sign_diff == "Different") / n()) %>%
  filter(diff_sign_freq>0)

data %>%
  filter(cov == "K[j[obs]]^'*'") %>%
  mutate(sign_diff = ifelse(sign(mean) == sign(true), "Equal", "Different")) %>%
  group_by(variable) %>%
  summarize(equal_sign_freq = sum(sign_diff == "Equal") / n(),
            diff_sign_freq = sum(sign_diff == "Different") / n()) %>%
  filter(diff_sign_freq>0)


data %>%
  filter(cov == "I[j]^'*'") %>%
  mutate(sign_diff = ifelse(sign(mean) == sign(true), "Equal", "Different")) %>%
  group_by(variable,m) %>%
  summarize(equal_sign_freq = sum(sign_diff == "Equal") / n(),
            diff_sign_freq = sum(sign_diff == "Different") / n()) %>%
  filter(diff_sign_freq>0) %>%
  arrange(m)

data %>%
  filter(cov == "K[j[obs]]^'*'") %>%
  mutate(sign_diff = ifelse(sign(mean) == sign(true), "Equal", "Different")) %>%
  group_by(variable,m) %>%
  summarize(equal_sign_freq = sum(sign_diff == "Equal") / n(),
            diff_sign_freq = sum(sign_diff == "Different") / n()) %>%
  filter(diff_sign_freq>0) %>%
  arrange(m)

