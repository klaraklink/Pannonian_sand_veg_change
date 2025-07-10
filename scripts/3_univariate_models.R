library(lmerTest) # v. 3.1-3
library(broom) # v. 1.0.5
library(geepack) # v. 1.3.12
library(broom.mixed) # v. 0.2.9.6
library(patchwork) # v. 1.3.0
library(readxl) # v. 1.4.3
library(ggeffects) # v. 1.7.2
library(tidyverse) # v. 2.0.0
# dplyr v. 1.1.4
# forcats v. 1.0.0
# readr v. 2.1.5
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
  filter(!is.na(rs_plot)) |> 
  left_join(read_csv('data/classification.csv')) |> 
  group_by(rs_plot) |> 
  fill(cluster, .direction = 'updown') |> 
  ungroup() |> 
  mutate(time = fct_relevel(time, 'old', 'new'), 
         cluster = fct_relevel(cluster, 'Corynephorion', 'Armerion', 'Festucion valesiacae')) |> 
  arrange(releve_nr)

annuals <- read_csv('data/annuals.csv') 
vernal_geo <- tibble(valid_name = c('Gagea pratensis agg.', 'Muscari neglectum', 
                'Ornithogalum umbellatum agg.', 'Poa bulbosa', 
                'Neotinea ustulata'))

spe <- read_csv('data/spe_resurvey_paired.csv') |> 
  semi_join(head) |> 
  filter(is.na(exclude)) |> 
  summarise(to_br_bl_old = combine.cover(to_br_bl_old), .by = c('releve_nr', 'valid_name')) |> 
  arrange(releve_nr) |> 
  anti_join(annuals) |> 
  anti_join(vernal_geo)

glimpse(spe)
glimpse(head)

diag <- read_csv('data/diag_czskat.csv') 

aliens <- read_xlsx('data/species_sands_czskat-2025-01-16-IA.xlsx') 

redlist <- read_csv('data/red_list_czskat.csv') |> 
  mutate(RegionFloraVeg = case_when(RegionFloraVeg == 'Czechia' ~ 'CZ', 
                                    RegionFloraVeg == 'Slovakia' ~ 'SK',
                                    RegionFloraVeg == 'Austria' ~ 'AT'), 
         Redlisted = ifelse(Status %in% c('CR', 'EN', 'VU', 'G', 'NE', 'EX', 'CR(PE)'), 1, 0)) |> 
  rename('RL_status' = 'Status')

woody <- read_csv('data/woody.csv') |> 
  mutate(woody = if_else(woody == T, 1, 0, 0))

# calculate species richness, nr of specialists etc. ----------------------------------------------
spe_trait <- spe |> 
  left_join(diag, by = join_by('valid_name' == 'Taxon')) |> 
  left_join(head |> select(releve_nr, country)) |> 
  left_join(aliens) |> 
  left_join(redlist, by = join_by('valid_name' == 'TaxonKaplan', 'country' == 'RegionFloraVeg')) |> 
  left_join(woody)

head2 <- head |> 
  left_join(spe |>
              group_by(releve_nr) |> 
              count(name = 'richness')) |> 
  left_join(spe_trait |>
              filter(habitat == 'sand') |> 
              distinct(valid_name, releve_nr) |> 
              group_by(releve_nr) |> 
              count(name = 'nr_sand_spec')) |> 
  left_join(spe_trait |>
              filter(habitat == 'dry grassland') |> 
              distinct(valid_name, releve_nr) |> 
              group_by(releve_nr) |> 
              count(name = 'nr_grass_spec')) |> 
  left_join(spe_trait |>
              filter(habitat == 'ruderal') |> 
              distinct(valid_name, releve_nr) |> 
              group_by(releve_nr) |> 
              count(name = 'nr_ruderal_spec')) |> 
  left_join(spe_trait |>
              filter(status != 'native') |> 
              distinct(valid_name, releve_nr) |> 
              group_by(releve_nr) |> 
              count(name = 'nr_alien')) |> 
  left_join(spe_trait |>
              filter(status == 'neo') |> 
              distinct(valid_name, releve_nr) |> 
              group_by(releve_nr) |> 
              count(name = 'nr_neo')) |> 
  left_join(spe_trait |>
              filter(status == 'arch') |> 
              distinct(valid_name, releve_nr) |> 
              group_by(releve_nr) |> 
              count(name = 'nr_arch')) |> 
  left_join(spe_trait |>
              filter(Redlisted == 1) |> 
              distinct(valid_name, releve_nr) |> 
              group_by(releve_nr) |> 
              count(name = 'nr_redlist')) |> 
  left_join(spe_trait |>
              filter(woody == 1) |> 
              distinct(valid_name, releve_nr) |> 
              group_by(releve_nr) |> 
              count(name = 'nr_woody')) |> 
  left_join(spe_trait |>
              filter(woody == 1) |> 
              distinct(valid_name, releve_nr, to_br_bl_old) |> 
              group_by(releve_nr) |> 
              summarize(cov_woody = sum(to_br_bl_old))) |> 
  mutate(prop_sand_spec = nr_sand_spec / richness, 
         prop_grass_spec = nr_grass_spec / richness, 
         prop_ruderal_spec = nr_ruderal_spec / richness,
         prop_alien = nr_alien / richness, 
         prop_neo = nr_neo / richness,
         prop_redlist = nr_redlist / richness, 
         prop_woody = nr_woody / richness, 
         across(c(starts_with('prop'), cov_woody), ~ifelse(is.na(.x), 0, .x))) |> 
  arrange(rs_plot, year) 


# build SAR model
sar <- lm(log(richness) ~ log(surf_area), data = head2)
summary(sar) 
coef(sar)

head3 <- head2 |> 
  group_by(rs_plot) |> 
  mutate(rich_sar = exp(log(richness) + (log(25) - log(surf_area))*coef(sar)[2]), 
         year_diff = year - lag(year), 
         year_diff = ifelse(is.na(year_diff), lead(year_diff), year_diff),
         rich_diff = rich_sar - lag(rich_sar),
         nr_sand_diff = nr_sand_spec - lag(nr_sand_spec), 
         prop_sand_diff = prop_sand_spec - lag(prop_sand_spec), 
         nr_alien_diff = nr_alien - lag(nr_alien), 
         prop_alien_diff = prop_alien - lag(prop_alien),
         nr_grass_diff = nr_grass_spec - lag(nr_grass_spec),
         prop_grass_diff = prop_grass_spec - lag(prop_grass_spec),
         nr_ruderal_diff = nr_ruderal_spec - lag(nr_ruderal_spec),
         prop_ruderal_diff = prop_ruderal_spec - lag(prop_ruderal_spec),
         nr_redlist_diff = nr_redlist - lag(nr_redlist),
         prop_redlist_diff = prop_redlist - lag(prop_redlist),
         prop_woody_diff = prop_woody - lag(prop_woody),  
         rs_plot = as.factor(rs_plot)) |> 
  ungroup()

# species richness --------------------------------------------------------

# boxplots for comparison new vs old by country
a <- head3 |> 
  ggplot(aes(cluster, rich_sar, fill = time))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.title = element_blank())+
  labs(y = bquote('Predicted number of species per 25'~m^2))

# model species richness
m1 <- lmer(rich_sar ~ time*cluster + (1|rs_plot), data = head3)

summary(m1)
anova(m1)
tidy(m1, conf.int = T) -> mm
plot(m1)

# significant change in time, increasing number of species over time 
# in Arm grasslands significantly more species per plot than in Cory
# in Fes.val, species richness is significantly more increasing over time than in Cory and Arm


# sand specialists -------------------------------------------------------------

# boxplot sand specialists

b <- head3 |> 
  ggplot(aes(cluster, prop_sand_spec, fill = time))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.title = element_blank())+
  scale_y_continuous(limits = c(0, 1))+
  labs(y = 'Proportion of sand specialists')


# model proportion of sand specialists
m2 <- geeglm(cbind(nr_sand_spec, richness-nr_sand_spec) ~ time * cluster, id = rs_plot, 
             family = 'binomial', corstr = 'exchangeable', data = head3)

summary(m2)
anova(m2)
tidy(m2, conf.int = T) |> mutate(estimate = exp(estimate)/(1+exp(estimate)))
plot(m2)

# significant decline over time
# in Cory significantly more sand specialists than in Arm and Fesval and significantly more decreasing

# dry grassland specialists -------------------------------------------------------------

# boxplot dry grassland specialists

c <- head3 |> 
  ggplot(aes(cluster, prop_grass_spec, fill = time))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.title = element_blank())+
  scale_y_continuous(limits = c(0, 1))+
  labs(y = 'Proportion of dry grassland specialists')


# model proportion of dry grassland specialists
m3 <- geeglm(cbind(nr_grass_spec, richness-nr_grass_spec) ~ time + cluster, id = rs_plot, 
             family = 'binomial', corstr = 'exchangeable', data = head3)

summary(m3)
anova(m3)
tidy(m3, conf.int = T) |> mutate(estimate = exp(estimate)/(1+exp(estimate)))
plot(m3)

# significant decrease over time
# significantly more dry grassland specialists in Arm and Fes.val than in Cory

# ruderal species -------------------------------------------------------------

# boxplot ruderal species

d <- head3 |> 
  ggplot(aes(cluster, prop_ruderal_spec, fill = time))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.title = element_blank())+
  scale_y_continuous(limits = c(0, 1))+
  labs(y = 'Proportion of ruderal species')

# model proportion of ruderal species
m4 <- geeglm(cbind(nr_ruderal_spec, richness-nr_ruderal_spec) ~ time + cluster, id = rs_plot, 
             family = 'binomial', corstr = 'exchangeable', data = head3)

summary(m4)
anova(m4)
tidy(m4, conf.int = T) |> mutate(estimate = exp(estimate)/(1+exp(estimate)))
plot(m4)

# significant increase of ruderal species over time
# significantly higher proportion of ruderal species in Cory than in Arm and Fes.val

# Red-List species -----------------------------------------------------------

# boxplot red list
e <- head3 |> 
  ggplot(aes(cluster, prop_redlist, fill = time))+
  geom_boxplot()+
  theme_bw()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.title = element_blank())+
  scale_y_continuous(limits = c(0, 1))+
  labs(y = 'Proportion of threatened species')


# proportion of redlist species
m5 <- geeglm(cbind(nr_redlist, richness - nr_redlist) ~ time + cluster, id = rs_plot, 
             family = 'binomial', corstr = 'exchangeable', data = head3)
summary(m5)
anova(m5)
tidy(m5, conf.int = T) 
plot(m5)

predict_response(m6,  terms = c('cluster', 'time')) |> plot()


# alien species -----------------------------------------------------------

# boxplot alien
f <- head3 |> 
  ggplot(aes(cluster, prop_alien, fill = time))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.title = element_blank())+
  scale_y_continuous(limits = c(0, 1))+
  labs(y = 'Proportion of alien species')

# proportion of alien species
m6 <- geeglm(cbind(nr_alien, richness - nr_alien) ~ time + cluster, id = rs_plot, 
             family = 'binomial', corstr = 'exchangeable', data = head3)
summary(m6)
anova(m6)
tidy(m6, conf.int = T) 
plot(m6)

# significantly higher proportion of alien species in Cory than in Arm and Fes.val

# woody species -----------------------------------------------------------

# boxplot woody
g <- head3 |> 
  ggplot(aes(cluster, cov_woody, fill = time))+
  geom_boxplot()+
  theme_bw()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.title = element_blank())+
  #scale_y_continuous(limits = c(0, 1))+
  labs(y = 'Proportion of woody species')


# proportion of woody species
m7 <- geeglm(cbind(nr_woody, richness - nr_woody) ~ time + cluster, id = rs_plot, 
             family = 'binomial', corstr = 'exchangeable', data = head3)
summary(m7)
anova(m7)
tidy(m7, conf.int = T) 
plot(m7)


# combine all plots and save
a + b + c + d + e + f +
  plot_layout(nrow = 3, axes = 'collect_x', guides = 'collect') &
  theme(legend.position = 'bottom', legend.title = element_blank(), 
        axis.title.x = element_blank()) & 
  plot_annotation(tag_levels = 'A')

ggsave('plots/boxplots_clusters.png', width = 8, height = 10)


#  trend magnitude vs interval length -----------------------------------------
# makes sense onlu for Corynephorion

### richness
head5 <- head3 |> 
  filter(!is.na(rich_diff) & cluster == 'Corynephorion')

head5$rich_diff |> hist()


# difference in species richness vs time between surveys
m8 <- lm(rich_diff ~ year_diff, data = head5)

summary(m8)
anova(m8)
tidy(m8, conf.int = T)

# not significant

### specialists
head5$prop_sand_diff |> hist()

m9 <- lm(prop_sand_diff ~ year_diff, data = head5)

summary(m9)
anova(m9)
tidy(m9, conf.int = T)

### dry grassland specialists
m10 <- lm(prop_grass_diff ~ year_diff, data = head5)

summary(m10)
anova(m10)
tidy(m10, conf.int = T)

### ruderal species
m11 <- lm(prop_ruderal_diff ~ year_diff, data = head5)

summary(m11)
anova(m11)
tidy(m11, conf.int = T)

### red list species
m12 <- lm(prop_redlist_diff ~ year_diff, data = head5)

summary(m12)
anova(m12)
tidy(m12, conf.int = T)

### alien species
head5$prop_alien_diff |> hist()

m13 <- lm(prop_alien_diff ~ year_diff, data = head5)

summary(m13)
anova(m13)
tidy(m13, conf.int = T)

### woody species
m14 <- lm(prop_woody_diff ~ year_diff, data = head5)

summary(m14)
anova(m14)
tidy(m14, conf.int = T)


# check which alien species are present in old and new plot observ --------

spe_alien <- spe_trait |> 
  filter(status != 'native') |> 
    left_join(head3, by = 'releve_nr') |> 
    distinct(releve_nr, valid_name, .keep_all = T) |> 
  group_by(valid_name, time, status) |>  
  count() |> 
  ungroup()

spe_alien_unique <- spe_alien |> 
  mutate(n_times = n_distinct(time), .by = 'valid_name') |> 
  filter(n_times == 1)

spe_alien_both <- spe_alien |> 
  mutate(n_times = n_distinct(time), .by = 'valid_name') |> 
  filter(n_times == 2)

# in how many plots do the alien species occur
spe_trait |> 
  filter(status != 'native') |> 
  left_join(head3, by = 'releve_nr') |> 
  distinct(releve_nr, valid_name, .keep_all = T) |> 
  summarize(n_plots = n_distinct(releve_nr), n = n(), .by  = c('time', 'status'))

spe_trait |> 
  filter(status != 'native') |> 
  left_join(head3, by = 'releve_nr') |> 
  distinct(releve_nr, valid_name, .keep_all = T) |> 
  summarize(n_alien = n_distinct(valid_name), n = n(), .by  = c('time', 'releve_nr')) |> 
  arrange(desc(n))
