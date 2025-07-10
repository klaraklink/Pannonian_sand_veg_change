library(janitor) # v. 2.2.0
library(readxl) # v. 1.4.3
library(weimea) # v. 0.1.4
library(tidyverse) # v. 2.0.0
# dplyr v. 1.1.4
# forcats v. 1.0.0
# readr v. 2.1.5
# stringr v. 1.5.0
# tibble v. 3.2.1
# tidyr v. 1.3.1

# functions ---------------------------------------------------------------

combine.cover <- function(x){
  #x = spe_cov_mod$cover_perc
  while (length(x)>1){
    x[2] <- x[1]+(100-x[1])*x[2]/100
    x <- x[-1]
  }
  return(x)
}

# load data ---------------------------------------------------------------

head <- read_csv('data/head_resurvey_paired.csv') |> 
  left_join(read_csv('data/classification.csv')) |> 
  group_by(rs_plot) |> 
  fill(cluster, .direction = 'updown') |> 
  ungroup() |> 
  mutate(time = fct_relevel(time, 'old', 'new'), 
         cluster = fct_relevel(cluster, 'Corynephorion', 'Armerion', 'Festucion valesiacae')) |> 
  filter(!is.na(rs_plot)) |> 
  arrange(releve_nr)

annuals <- read_csv('data/annuals.csv') |> mutate(valid_name = str_trim(valid_name))
vernal_geo <- tibble(valid_name = c('Gagea pratensis agg.', 'Muscari neglectum', 
                                    'Ornithogalum umbellatum agg.', 'Poa bulbosa', 
                                    'Neotinea ustulata'))

spe <- read_csv('data/spe_resurvey_paired.csv') |> 
  semi_join(head) |> 
  arrange(releve_nr) |> 
  filter(is.na(exclude)) |> 
  summarise(cover_perc = combine.cover(to_br_bl_old), .by = c('releve_nr', 'valid_name')) |>  # merge species covers across layers
  anti_join(annuals) |> 
  anti_join(vernal_geo)

eiv <- read_xlsx('data/Ellenberg_disturbance.xlsx') |> 
  clean_names() |> 
  mutate(across(light:sd_soil_disturbance, ~as.numeric(.x)), 
         across(c(taxon_name, species_level_name), ~str_replace(.x, 'aggr.', 'agg.'))) |> 
  left_join(read_csv('data/taxonomy_eu_convert.csv')) |> 
  mutate(valid_name = coalesce(valid_name, taxon_name)) |> 
  select(valid_name, light:nutrients, disturbance_severity:soil_disturbance)

glimpse(eiv)
glimpse(spe)
glimpse(head)

spe_eiv <- spe |> 
  left_join(eiv)

spe_eiv_cor <- spe_eiv |> 
  distinct(valid_name, .keep_all = T) 

# correlations of eivs
cor(spe_eiv_cor |> select(light:soil_disturbance, -ends_with('herblayer')), use = 'pairwise.complete.obs') 
png('plots/correlation_eivs.png', width = 1500, height = 1100, res = 200)
corrplot::corrplot(cor(spe_eiv_cor |> select(light:soil_disturbance,  -ends_with('herblayer')), use = 'pairwise.complete.obs'), 
         method = 'number', type = 'upper', tl.col = 'black', tl.srt = 45, diag = F)
dev.off()


# calculate community unweighted means --------------------------

head_cwm <- head |> 
  select(releve_nr, cluster, year, time, rs_site, rs_plot) |> 
  left_join(spe_eiv |> 
              summarise(across(light:soil_disturbance, ~mean(.x, na.rm = T)), 
                        .by = 'releve_nr')) |> 
  pivot_longer(cols = light:soil_disturbance, names_to = 'variable') |> 
  mutate(variable = str_remove(variable, '_cum') |> 
           fct_relevel('light', 'temperature', 'moisture', 'reaction', 'nutrients', 
                       'disturbance_severity', 'disturbance_severity_herblayer', 
                       'disturbance_frequency', 'disturbance_frequency_herblayer', 
                       'mowing_frequency', 'grazing_pressure', 'soil_disturbance')) |> 
  filter(str_detect(variable, 'herblayer', negate = T))


glimpse(head_cwm)

head_cwm |> 
  ggplot(aes(time, value, fill = time))+
  geom_boxplot()+
  facet_wrap(~variable, scales = 'free_y', nrow = 2, 
             labeller = labeller(variable = c('light' = '(a) Light', 
                                              'temperature' = '(b) Temperature', 
                                              'moisture' = '(c) Moisture', 
                                              'reaction' = '(d) Reaction', 
                                              'nutrients' = '(e) Nutrients',
                                              'disturbance_severity' = '(f) Disturbance severity', 
                                              'disturbance_frequency' = '(g) Disturbance frequency' , 
                                              'mowing_frequency' = '(h) Mowing frequency', 
                                              'grazing_pressure' = '(i) Grazing pressure', 
                                              'soil_disturbance' = '(j) Soil disturbance')))+
  theme_bw()+
  theme(axis.title = element_blank(), legend.position = 'none', 
        strip.background = element_blank())

ggsave('plots/boxplots_eivs_unweighted.png', width = 10, height = 6)

head_cwm |> 
  ggplot(aes(cluster, value, fill = time))+
  geom_boxplot()+
  facet_wrap(~variable, scales = 'free_y')+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.title = element_blank(), legend.position = 'bottom')

ggsave('plots/boxplots_eivs_unweighted_clusters.png', width = 12, height = 8)

# testing -----------------------------------------------------------------

# species data to wide format for weimea, pa data for unweighted mean
spe_pa <- spe |> 
  select(releve_nr, valid_name, cover_perc) |> 
  mutate(pres_abs = 1) |> 
  select(-cover_perc) |> 
  arrange(valid_name) |> 
  pivot_wider(names_from = valid_name, values_from = pres_abs, values_fill = 0) |> 
  arrange(releve_nr) |> 
  as.data.frame() |> 
  column_to_rownames('releve_nr')

eiv2 <- spe_eiv |> 
  distinct(valid_name, .keep_all = T) |> 
  select(-releve_nr, -cover_perc) |> 
  arrange(valid_name) |> 
  as.data.frame() |> 
  column_to_rownames('valid_name')

# Calculate EIV using cwm function from weimea package
eiv_cwm <- cwm(com = spe_pa, traits = list(light = eiv2$light, 
                                           temperature = eiv2$temperature, 
                                           moisture = eiv2$moisture, 
                                           reaction = eiv2$reaction, 
                                           nutrients = eiv2$nutrients, 
                                           disturbance_severity = eiv2$disturbance_severity, 
                                           disturbance_frequency = eiv2$disturbance_frequency, 
                                           mowing_frequency = eiv2$mowing_frequency, 
                                           grazing_pressure = eiv2$grazing_pressure, 
                                           soil_disturbance = eiv2$soil_disturbance))


# Standard one-way ANOVA, parametric test, F-value = 16.24
summary(aov(eiv_cwm$soil_disturbance ~ head$time)) 

# Repeated measure ANOVA (between pairs of plots sampled in different time), parametric test (we care about "Error: Within" part)
summary(aov(eiv_cwm$light ~ head$time + Error(head$rs_plot))) 

# set number of permutations
no_perm <- 999

# 1) obtain observed F-value
ANOVA_obs <- summary(aov(eiv_cwm$temperature ~ head$time + Error(head$rs_plot)))
F_obs <- ANOVA_obs$`Error: Within`[[1]]$`F value`[1]

# 2) row-based permutation test (permute time)
time_perm <- replicate(no_perm, expr = sample(head$time))
F_row <- sapply(as.data.frame(time_perm), FUN = function (time_rand)
  summary(aov(eiv_cwm$temperature ~ time_rand + Error (head$rs_plot)))$`Error: Within`[[1]]$`F value`[1]
)

P_row <- sum(c(F_row, F_obs) >= F_obs)/(no_perm+1)

# 3) col-based test (permute eiv_cwm before calculating CWM)
rand_cwm <- randomize(eiv_cwm, permutations = no_perm, parallel = 8)
F_col <- sapply(rand_cwm, FUN = function(cwm_rand) 
  summary(aov(as.matrix(cwm_rand) ~ head$time + Error(head$rs_plot)))$`Error: Within`[[1]]$`F value`[1]
)

P_col <- sum(c(F_col, F_obs) >= F_obs)/(no_perm+1)

# 4) combine P_row and P_col into P_max
P_max <- max(P_row, P_col)
P_max

### whole dataset
# light ns aov
# temperature ns aov, Pmax 0.089 no annuals decreasing +
# moisture P_row 0.03, P_col 0.102, 0.05 no annuals increasing +
# reaction ns aov
# nutrients P_row 0.021, P_col 0.081 increasing +, 0.03 no annuals
# disturbance_severity P_row 0.032, P_col 0.12
# disturbance_frequency ns aov
# mowing_frequency ns aov
# grazing_pressure P_max 0.049, decreasing *, 0.013 no annuals
# soil disturbance P_row 0.038, P_col 0.146

# for vegetation types separately -------------------------------------------

head2 <- head |> 
  filter(cluster == 'Festucion valesiacae') # change accordign to vegetation type

# species data to wide format for weimea, pa data for unweighted mean
spe_pa <- spe |> 
  semi_join(head2) |> 
  select(releve_nr, valid_name, cover_perc) |> 
  mutate(pres_abs = 1) |> 
  select(-cover_perc) |> 
  arrange(valid_name) |> 
  pivot_wider(names_from = valid_name, values_from = pres_abs, values_fill = 0) |> 
  arrange(releve_nr) |> 
  as.data.frame() |> 
  column_to_rownames('releve_nr')

eiv2 <- spe_eiv |> 
  semi_join(head2) |> 
  distinct(valid_name, .keep_all = T) |> 
  select(-releve_nr, -cover_perc) |> 
  arrange(valid_name) |> 
  as.data.frame() |> 
  column_to_rownames('valid_name')

# Calculate EIV using cwm function from weimea package
eiv_cwm <- cwm(com = spe_pa, traits = list(light = eiv2$light, 
                                           temperature = eiv2$temperature, 
                                           moisture = eiv2$moisture, 
                                           reaction = eiv2$reaction, 
                                           nutrients = eiv2$nutrients, 
                                           disturbance_severity = eiv2$disturbance_severity, 
                                           disturbance_frequency = eiv2$disturbance_frequency, 
                                           mowing_frequency = eiv2$mowing_frequency, 
                                           grazing_pressure = eiv2$grazing_pressure, 
                                           soil_disturbance = eiv2$soil_disturbance))

# Standard one-way ANOVA, parametric test, F-value = 16.24
summary(aov(eiv_cwm$soil_disturbance ~ head2$time)) 

# Repeated measure ANOVA (between pairs of plots sampled in different time), parametric test (we care about "Error: Within" part)
summary(aov(eiv_cwm$soil_disturbance ~ head2$time + Error(head2$rs_plot))) 

# set number of permutations
no_perm <- 999

# 1) obtain observed F-value
ANOVA_obs <- summary(aov(eiv_cwm$soil_disturbance ~ head2$time + Error(head2$rs_plot)))
F_obs <- ANOVA_obs$`Error: Within`[[1]]$`F value`[1]

# 2) row-based permutation test (permute time)
time_perm <- replicate(no_perm, expr = sample(head2$time))
F_row <- sapply(as.data.frame(time_perm), FUN = function (time_rand)
  summary(aov(eiv_cwm$soil_disturbance ~ time_rand + Error (head2$rs_plot)))$`Error: Within`[[1]]$`F value`[1]
)

P_row <- sum(c(F_row, F_obs) >= F_obs)/(no_perm+1)

# 3) col-based test (permute eiv_cwm before calculating CWM)
rand_cwm <- randomize(eiv_cwm, permutations = no_perm, parallel = 8)
F_col <- sapply(rand_cwm, FUN = function(cwm_rand) 
  summary(aov(as.matrix(cwm_rand) ~ head2$time + Error(head2$rs_plot)))$`Error: Within`[[1]]$`F value`[1]
)

P_col <- sum(c(F_col, F_obs) >= F_obs)/(no_perm+1)

# 4) combine P_row and P_col into P_max
P_max <- max(P_row, P_col)
P_max

### Corynephorion
# light ns aov
# temperature P_max 0.019 decreasing *, 0.024 no annuals
# moisture ns aov
# reaction ns aov
# nutrients P_row 0.033, P_col 0.14
# disturbance_severity aov
# disturbance_frequency ns aov
# mowing_frequency ns aov
# grazing_pressure ns aov, 0.088 no annuals, decreasing +
# soil disturbance ns aov

### Armerion
# light ns aov
# temperature ns aov
# moisture ns aov, 0.04 no annuals, increasing *
# reaction ns aov
# nutrients ns_aov, 0.038 no annuals, increasing *
# disturbance_severity ns_aov, 0.046 no annuals, increasing *
# disturbance_frequency ns aov
# mowing_frequency ns aov
# grazing_pressure ns_aov, 0.026 no annuals, decreasing *
# soil disturbance ns aov

### Festucion valesiacae
# no annuals - no significant results, only increasing soil disturbance + 0.071




