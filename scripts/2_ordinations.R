library(vegan) # v. 2.6-8
library(patchwork) # v. 1.3.0
library(ggrepel) # v. 0.9.6
library(broom) # v. 1.0.5
library(tidyverse) # v. 2.0.0
# dplyr v. 1.1.4
# forcats v. 1.0.0
# purrr v. 1.0.2
# readr v. 2.1.5
# stringr v. 1.5.0
# tibble v. 3.2.1
# tidyr v. 1.3.1

# functions ---------------------------------------------------------------

# merge same species, random overlaps
combine.cover <- function(x){
  #x = spe_cov_mod$cover_perc
  while (length(x)>1){
    x[2] <- x[1]+(100-x[1])*x[2]/100
    x <- x[-1]
  }
  return(x)
}

# data to wide format with square-root transformation of species covers
to_wide <- function(data_long){
  data_long |> 
    select(releve_nr, valid_name, cover_perc = to_br_bl_old) |> 
    pivot_wider(names_from = valid_name, values_from = cover_perc, values_fill = 0) |> 
    as.data.frame() |> 
    column_to_rownames('releve_nr') |> 
    sqrt()
}

# load data ---------------------------------------------------------------

head <- read_csv('data/head_resurvey_paired.csv') |> 
  left_join(read_csv('data/classification.csv')) |> 
  group_by(rs_plot) |> 
  fill(cluster, .direction = 'updown') |> 
  mutate(year0 = year - lag(year), 
         year0 = ifelse(is.na(year0), 0, year0)) |> 
  ungroup() |> 
  mutate(time = fct_relevel(time, 'old', 'new')) |> 
  filter(!is.na(rs_plot)) |> 
  arrange(releve_nr)

annuals <- read_csv('data/annuals.csv') |> mutate(valid_name = str_trim(valid_name))

vernal_geo <- tibble(valid_name = c('Gagea pratensis agg.', 'Muscari neglectum', 
                                    'Ornithogalum umbellatum agg.', 'Poa bulbosa', 
                                    'Neotinea ustulata'))

spe <- read_csv('data/spe_resurvey_paired.csv') |> 
  semi_join(head) |> 
  mutate(n = n(), .by = 'valid_name') |> 
  filter(is.na(exclude) & n > 1) |> 
  summarise(to_br_bl_old = combine.cover(to_br_bl_old), .by = c('releve_nr', 'valid_name')) |> 
  arrange(releve_nr) |> 
  anti_join(annuals) |> 
  anti_join(vernal_geo)
  
glimpse(spe)
glimpse(head)

# data preparation
spe_wide <- spe |> 
  to_wide()

# Corynephorion
head_cory <- head |> 
  filter(cluster == 'Corynephorion')

spe_wide_cory <- spe |> 
  semi_join(head_cory) |> 
  to_wide()

# Armerion
head_arm <- head |> 
  filter(cluster == 'Armerion')

spe_wide_arm <- spe |> 
  semi_join(head_arm) |> 
  to_wide()

# Festucion valesiacae
head_fesval <- head |> 
  filter(cluster == 'Festucion valesiacae')

spe_wide_fesval <- spe |> 
  semi_join(head_fesval) |> 
  to_wide()


# ordinations --------------------------------------------------------------

# unconstrained ordination
pcoa <- capscale(spe_wide ~ 1, distance = "bray", sqrt.dist = T)
screeplot(pcoa, bstick = T)
labs <- paste0('PCo', 1:9, ' (', round(eigenvals(pcoa)/sum(eigenvals(pcoa)) * 100, 2), '%)')

head_ord <- bind_cols(head, as_tibble(scores(pcoa, choices = 1:3)$sites)) 

# ordination plot 
pcoa_si <- head_ord %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(fill = cluster, shape = time),
             size = 3) +
  scale_shape_manual(values = 21:22,
                     labels = c("old", "new")) +
  geom_path(aes(MDS1, MDS2, group = rs_plot), arrow = arrow(length = unit(0.3, "cm"), ends = 'first')) +
  theme_bw() +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0, 0),
        legend.justification = c(0, 0), 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        plot.title = element_text(size = 11))+
  guides(fill = guide_legend(override.aes = list(pch = 21)))+
  labs(x = labs[1], y = labs[2])

# ordination plot 
pcoa_si_time <- head_ord %>%
  arrange(rs_plot, year) |> 
  mutate(year_diff = year - lag(year), .by = (rs_plot)) |> 
  fill(year_diff, .direction = 'updown') |>
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(fill = year_diff, shape = time),
             size = 3) +
  scale_shape_manual(values = 21:22,
                     labels = c("old", "new")) +
  geom_path(aes(MDS1, MDS2, group = rs_plot), arrow = arrow(length = unit(0.3, "cm"), ends = 'last')) +
  theme_bw() +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0, 0),
        legend.justification = c(0, 0), 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        plot.title = element_text(size = 11))+
  guides(fill = guide_legend(override.aes = list(pch = 21)))+
  labs(x = labs[1], y = labs[2])

spe_fit_pcoa <- envfit(pcoa, spe_wide, display = 'si', choices = 1:2, permutations = 0) 

spe_pcoa <- as_tibble(scores(pcoa, display = 'species', scaling='si', correlation=T), 
                      rownames = 'species_name') %>% 
  mutate(r = spe_fit_pcoa$vectors$r, # goodness of fit
         short_name = str_c(str_split_i(species_name, '\\s', 1) %>% str_sub(., 1, 3), # make short species names Genus species -> Gen.spe
                            str_split_i(species_name, '\\s', 2) %>% str_sub(., 1, 3), sep = '.'))

# ordination plot with species
pcoa_spe <- spe_pcoa %>% 
  slice_max(r, n = 25) %>% # select 40 best-fitting species
  ggplot() + 
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), alpha = 0.6, 
               arrow = arrow(length = unit(0.2, 'cm')), show.legend = F)+ # arrows for species
  geom_text(aes(MDS1 + 0.02 * sign(MDS1), MDS2 + 0.02 * sign(MDS2), label = short_name), # label arrows
            alpha = 0.6, cex = 3.5)+
  theme_bw()+
  labs(x = labs[1], y = labs[2])

# bind together ordination plots for sites and species
pcoa_si + pcoa_spe + plot_annotation(title = 'PCoA, All vegetation types (n = 86)', tag_levels = 'A')
ggsave('plots/pcoa_cluster.png', width = 12, height = 6)

# clusters separately ----------------------------------------------------

# Corynephorion
pcoa_cory <- capscale(spe_wide_cory ~ 1, distance = "bray", sqrt.dist = T)
screeplot(pcoa_cory, bstick = T)
labs <- paste0('PCo', 1:9, ' (', round(eigenvals(pcoa_cory)/sum(eigenvals(pcoa_cory)) * 100, 2), '%)')

head_ord_cory <- bind_cols(head_cory, as_tibble(scores(pcoa_cory)$sites)) 

# ordination plot 
pcoa_si_cory <- head_ord_cory %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(fill = time, shape = time),
             size = 3) +
  # geom_text(aes(color = time, label = releve_nr),
  #            size = 3) +
  scale_shape_manual(values = 21:22,
                     labels = c("old", "new")) +
  geom_path(aes(MDS1, MDS2, group = rs_plot), arrow = arrow(length = unit(0.3, "cm"), ends = 'first'))+
  theme_bw() +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0, 0),
        legend.justification = c(0, 0), 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        plot.title = element_text(size = 11))+
  labs(x = labs[1], y = labs[2])


spe_fit_pcoa_cory <- envfit(pcoa_cory, spe_wide_cory, display = 'si', choices = 1:2, permutations = 0) 

spe_pcoa_cory <- as_tibble(scores(pcoa_cory, display = 'species', scaling='si', correlation=T), 
                      rownames = 'species_name') %>% 
  mutate(r = spe_fit_pcoa_cory$vectors$r, # goodness of fit
         short_name = str_c(str_split_i(species_name, '\\s', 1) %>% str_sub(., 1, 3), # make short species names Genus species -> Gen.spe
                            str_split_i(species_name, '\\s', 2) %>% str_sub(., 1, 3), sep = '.'))

# ordination plot with species
pcoa_spe_cory <- spe_pcoa_cory %>% 
  slice_max(r, n = 25) %>% # select 40 best-fitting species
  ggplot() + 
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), alpha = 0.6, 
               arrow = arrow(length = unit(0.2, 'cm')), show.legend = F)+ # arrows for species
  geom_text(aes(MDS1+0.01*sign(MDS1), MDS2+0.01*sign(MDS2), label = short_name), # label arrows
            alpha = 0.6, cex = 3.5)+
  theme_bw()+
  labs(x = labs[1], y = labs[2])

# bind together ordination plots for sites and species
pcoa_si_cory + pcoa_spe_cory + plot_annotation(title = 'PCoA, Corynephorion (n = 29)', tag_levels = 'A') &
  scale_y_reverse() & scale_x_reverse()
ggsave('plots/pcoa_cory.png', width = 12, height = 6)


# Armerion
pcoa_arm <- capscale(spe_wide_arm ~ 1, distance = "bray", sqrt.dist = T)
screeplot(pcoa_arm, bstick = T)
labs <- paste0('PCo', 1:9, ' (', round(eigenvals(pcoa_arm)/sum(eigenvals(pcoa_arm)) * 100, 2), '%)')

head_ord_arm <- bind_cols(head_arm, as_tibble(scores(pcoa_arm)$sites)) 

# ordination plot 
pcoa_si_arm <- head_ord_arm %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(fill = time, shape = time),
             size = 3) +
  scale_shape_manual(values = 21:22,
                     labels = c("old", "new")) +
  geom_path(aes(MDS1, MDS2, group = rs_plot), arrow = arrow(length = unit(0.3, "cm"), ends = 'first')) +
  theme_bw() +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0, 0),
        legend.justification = c(0, 0), 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        plot.title = element_text(size = 11))+
  labs(x = labs[1], y = labs[2])


spe_fit_pcoa_arm <- envfit(pcoa_arm, spe_wide_arm, display = 'si', choices = 1:2, permutations = 0) 

spe_pcoa_arm <- as_tibble(scores(pcoa_arm, display = 'species', scaling='si', correlation=T), 
                         rownames = 'species_name') %>% 
  mutate(r = spe_fit_pcoa_arm$vectors$r, # goodness of fit
         short_name = str_c(str_split_i(species_name, '\\s', 1) %>% str_sub(., 1, 3), # make short species names Genus species -> Gen.spe
                            str_split_i(species_name, '\\s', 2) %>% str_sub(., 1, 3), sep = '.'))

# ordination plot with species
pcoa_spe_arm <- spe_pcoa_arm %>% 
  slice_max(r, n = 25) %>% # select 40 best-fitting species
  ggplot() + 
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), alpha = 0.6, 
               arrow = arrow(length = unit(0.2, 'cm')), show.legend = F)+ # arrows for species
  geom_text(aes(MDS1+0.02*sign(MDS1), MDS2+0.02*sign(MDS2), label = short_name), # label arrows
            alpha = 0.6, cex = 3.5)+
  theme_bw()+
  labs(x = labs[1], y = labs[2])

# bind together ordination plots for sites and species
pcoa_si_arm + pcoa_spe_arm + plot_annotation(title = 'PCoA, Armerion (n = 40)',  tag_levels = 'A') #&
  #scale_y_reverse() 
ggsave('plots/pcoa_arm.png', width = 12, height = 6)

# Festucion valesiacae
pcoa_fesval <- capscale(spe_wide_fesval ~ 1, distance = "bray", sqrt.dist = T)
screeplot(pcoa_fesval, bstick = T)
labs <- paste0('PCo', 1:9, ' (', round(eigenvals(pcoa_fesval)/sum(eigenvals(pcoa_fesval)) * 100, 2), '%)')

head_ord_fesval <- bind_cols(head_fesval, as_tibble(scores(pcoa_fesval)$sites)) 

# ordination plot 
pcoa_si_fesval <- head_ord_fesval %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(fill = time, shape = time),
             size = 3) +
  scale_shape_manual(values = 21:22,
                     labels = c("old", "new")) +
  geom_path(aes(MDS1, MDS2, group = rs_plot), arrow = arrow(length = unit(0.3, "cm"), ends = 'first')) +
  theme_bw() +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0, 0),
        legend.justification = c(0, 0), 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        plot.title = element_text(size = 11))+
  labs(x = labs[1], y = labs[2])


spe_fit_pcoa_fesval <- envfit(pcoa_fesval, spe_wide_fesval, display = 'si', choices = 1:2, permutations = 0) 

spe_pcoa_fesval <- as_tibble(scores(pcoa_fesval, display = 'species', scaling='si', correlation=T), 
                         rownames = 'species_name') %>% 
  mutate(r = spe_fit_pcoa_fesval$vectors$r, # goodness of fit
         short_name = str_c(str_split_i(species_name, '\\s', 1) %>% str_sub(., 1, 3), # make short species names Genus species -> Gen.spe
                            str_split_i(species_name, '\\s', 2) %>% str_sub(., 1, 3), sep = '.'))

# ordination plot with species
pcoa_spe_fesval <- spe_pcoa_fesval %>% 
  slice_max(r, n = 25) %>% # select 40 best-fitting species
  ggplot() + 
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), alpha = 0.6, 
               arrow = arrow(length = unit(0.2, 'cm')), show.legend = F)+ # arrows for species
  geom_text(aes(MDS1+0.008*sign(MDS1), MDS2+0.008*sign(MDS2), label = short_name), # label arrows
            alpha = 0.6, cex = 3.5)+
  theme_bw()+
  labs(x = labs[1], y = labs[2])

# bind together ordination plots for sites and species
pcoa_si_fesval + pcoa_spe_fesval + plot_annotation(title = 'PCoA, Festucion valesiacae (n = 17)',  tag_levels = 'A')
ggsave('plots/pcoa_fesval.png', width = 12, height = 6)

# constrained ordination --------------------------------------------------

# effect of time
dbrda <- capscale(spe_wide ~ head$year0 + Condition(as.factor(head$rs_plot)), 
                 distance = "bray", sqrt.dist = T)

eig <- paste0(names(eigenvals(dbrda)), ' (', round(eigenvals(dbrda)/sum(eigenvals(dbrda))*100, 1), '%)')

# significance test
anova(dbrda, permutations = how(blocks = as.factor(head$rs_plot), nperm = 999))
# significant change in time

# significantly increasing and decreasing species
dbrda_spe <- envfit(dbrda, spe_wide, choices = 1, 
                   permutations = how(blocks = as.factor(head$rs_plot), nperm = 999), display = "lc")

# filter species significantly correlated with the constrained axis
sig_spe <- as_tibble(dbrda_spe$vectors$arrows, rownames = "SpecName_Layer") %>% 
  bind_cols(p = dbrda_spe$vectors$pvals, r = dbrda_spe$vectors$r)  %>%  
  filter(p < 0.05) %>% 
  arrange(CAP1, p) |> 
  write_csv('results/increasing_decreasing_species_all.csv')

# effect of time, clusters separately
# Festucion valesiacae
dbrda_fesval <- capscale(spe_wide_fesval ~ head_fesval$year0 + Condition(as.factor(head_fesval$rs_plot)), 
                  distance = "bray", sqrt.dist = T)

eig_fesval <- paste0(names(eigenvals(dbrda_fesval)), ' (', round(eigenvals(dbrda_fesval)/sum(eigenvals(dbrda_fesval))*100, 1), '%)')

# significance test
anova(dbrda_fesval, permutations = how(blocks = as.factor(head_fesval$rs_plot), nperm = 999))
# significant change in time

# significantly increasing and decreasing species
dbrda_spe_fesval <- envfit(dbrda_fesval, spe_wide_fesval, choices = 1, 
                    permutations = how(blocks = as.factor(head_fesval$rs_plot), nperm = 999), display = "lc")

# filter species significantly correlated with the constrained axis
sig_spe_fesval <- as_tibble(dbrda_spe_fesval$vectors$arrows, rownames = "SpecName_Layer") %>% 
  bind_cols(p = dbrda_spe_fesval$vectors$pvals, r = dbrda_spe_fesval$vectors$r)  %>%  
  filter(p < 0.05) %>% 
  arrange(CAP1, p)|> 
  write_csv('results/increasing_decreasing_species_fesval.csv')


# Armerion
dbrda_arm <- capscale(spe_wide_arm ~ head_arm$year0 + Condition(as.factor(head_arm$rs_plot)), 
                     distance = "bray", sqrt.dist = T)

eig_arm <- paste0(names(eigenvals(dbrda_arm)), ' (', round(eigenvals(dbrda_arm)/sum(eigenvals(dbrda_arm))*100, 1), '%)')

# significance test
anova(dbrda_arm, permutations = how(blocks = as.factor(head_arm$rs_plot), nperm = 999))
# significant change in time

# significantly increasing and decreasing species
dbrda_spe_arm <- envfit(dbrda_arm, spe_wide_arm, choices = 1, 
                       permutations = how(blocks = as.factor(head_arm$rs_plot), nperm = 999), display = "lc")

# filter species significantly correlated with the constrained axis
sig_spe_arm <- as_tibble(dbrda_spe_arm$vectors$arrows, rownames = "SpecName_Layer") %>% 
  bind_cols(p = dbrda_spe_arm$vectors$pvals, r = dbrda_spe_arm$vectors$r)  %>%  
  filter(p < 0.05) %>% 
  arrange(CAP1, p) |> 
  write_csv('results/increasing_decreasing_species_arm.csv')

# Corynephorion
dbrda_cory <- capscale(spe_wide_cory ~ head_cory$year0 + Condition(as.factor(head_cory$rs_plot)), 
                     distance = "bray", sqrt.dist = T)

eig_cory <- paste0(names(eigenvals(dbrda_cory)), ' (', round(eigenvals(dbrda_cory)/sum(eigenvals(dbrda_cory))*100, 1), '%)')

# significance test
anova(dbrda_cory, permutations = how(blocks = as.factor(head_cory$rs_plot), nperm = 999))
# significant change in time

# significantly increasing and decreasing species
dbrda_spe_cory <- envfit(dbrda_cory, spe_wide_cory, choices = 1, 
                       permutations = how(blocks = as.factor(head_cory$rs_plot), nperm = 999), display = "lc")

# filter species significantly correlated with the constrained axis
sig_spe_cory <- as_tibble(dbrda_spe_cory$vectors$arrows, rownames = "SpecName_Layer") %>% 
  bind_cols(p = dbrda_spe_cory$vectors$pvals, r = dbrda_spe_cory$vectors$r)  %>%  
  filter(p < 0.05) %>% 
  arrange(CAP1, p) |> 
  write_csv('results/increasing_decreasing_species_cory.csv')

dbrda_all + dbrda_cory + dbrda_arm + dbrda_fesval
ggsave('plots/dbrda_time_clusters.png', width = 12, height = 12) 


# shortcuts of species names
spe_pcoa |> 
  slice_max(r, n = 25) |> 
  bind_rows(spe_pcoa_cory |> 
    slice_max(r, n = 25)) |> 
  bind_rows(spe_pcoa_arm |> 
    slice_max(r, n = 25)) |> 
  bind_rows(spe_pcoa_fesval |> 
    slice_max(r, n = 25)) |> 
  arrange(short_name) |> 
  distinct(short_name, species_name) |> 
  write_csv('results/species_shortcuts.csv')


# table with frequencies of individual species ----------------------------

spe_freq <- spe |> 
  left_join(head) |> 
  group_by(valid_name, time) |> 
  summarise(n = n(), .groups = 'drop') |> 
  arrange(valid_name, desc(n)) |> 
  pivot_wider(names_from = time, values_from = n, values_fill = 0) |>
  #mutate(bin_test = map2(new, old, ~binom.test(.x, .x + .y, p = 0.5, alternative = 'two.sided')))
  write_csv('results/spe_freq.csv')

spe_freq_cl <- spe |> 
  left_join(head) |> 
  group_by(valid_name, time, cluster) |> 
  summarise(n = n(), .groups = 'drop') |> 
  arrange(valid_name, desc(n)) |> 
  pivot_wider(names_from = time, values_from = n, values_fill = 0) |>
  #mutate(bin_test = map2(new, old, ~binom.test(.x, .x + .y, p = 0.5, alternative = 'two.sided')))
  write_csv('results/spe_freq_cluster.csv')


# dissimilarity between the old and new plots vs interval length ----------

dissim <- head |> 
  select(rs_plot, year, time) |> 
  distinct() |> 
  pivot_wider(names_from = time, values_from = year) |> 
  mutate(time_diff = new - old) |> 
  select(-old, -new) 

# Join metadata to species data
data_combined <- head %>%
  select(releve_nr, rs_plot) |> 
  bind_cols(spe_wide)

# Nest by pair_id to create sub-tibbles of 2 plots each
results <- data_combined %>%
  group_by(rs_plot) %>%
  nest() %>%
  mutate(
    bray_curtis = map_dbl(data, ~ {
      sp_data <- select(.x, -rs_plot, -releve_nr)
      dist(vegdist(sp_data, method = "bray")) # only 1 distance in a pair
    }),
    plot1 = map_chr(data, ~ .x$releve_nr[1]),
    plot2 = map_chr(data, ~ .x$releve_nr[2])
  ) %>%
  select(rs_plot, plot1, plot2, bray_curtis) |> 
  left_join(head |> 
              select(rs_plot, year, time) |> 
              distinct() |> 
              pivot_wider(names_from = time, values_from = year) |> 
              mutate(time_diff = new - old) |> 
              select(-old, -new)) |> 
  left_join(head |> distinct(rs_plot, cluster))

results |> 
  ggplot(aes(time_diff, bray_curtis, fill = cluster)) +
  geom_point(pch = 21) +
  geom_smooth(method = 'lm') +
  labs(x = 'Time difference (years)', y = 'Bray-Curtis dissimilarity') +
  theme_bw() +
  theme(legend.position = 'none')

results2 <- results |> filter(cluster == 'Corynephorion')

m1 <- lm(bray_curtis ~ time_diff, data = results2)

summary(m1)
anova(m1)
tidy(m2, conf.int = T)


