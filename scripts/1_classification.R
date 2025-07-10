library(isopam) # v. 2.0
library(vegan) # v. 2.6-8
library(indicspecies) # v. 1.7.15
library(factoextra) # v. 1.0.7
library(tidyverse) # v. 2.0.0
# dplyr v. 1.1.4
# readr v. 2.1.5
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
    select(releve_nr, valid_name, layer, cover_perc = to_br_bl_old) |> 
    summarise(cover_perc = combine.cover(cover_perc), .by = c('releve_nr', 'valid_name')) |> 
    pivot_wider(names_from = valid_name, values_from = cover_perc, values_fill = 0) |> 
    as.data.frame() |> 
    column_to_rownames('releve_nr') |> 
    mutate(across(everything(), ~as.numeric(.x))) 
}

# load data ---------------------------------------------------------------

# only old plots for classification
head <- read_csv('data/head_resurvey_paired.csv') |> 
  filter(!is.na(rs_plot) & time == 'old') |> 
  arrange(releve_nr)

spe <- read_csv('data/spe_resurvey_paired.csv') |> 
  semi_join(head) |> 
  group_by(valid_name) |> 
  mutate(n = n()) |> 
  filter(is.na(exclude) & n > 1) |> 
  ungroup() |> 
  arrange(releve_nr)

glimpse(head)

# data preparation
spe_wide <- spe |> 
  to_wide() |> 
  mutate(across(everything(), ~as.numeric(.x)))

# classification ------------------------------------------------
spe_sqrt <- spe_wide |> 
  sqrt() # square-root transformation

# Isopam analysis with max 5 groups per division and 3 hierarchical levels
iso.1 <- isopam(spe_sqrt, c.max = 5, l.max = 3, sieve = F) 

# Silhouette plot for the resulting classification
silhouette(iso.1$flat$lev.1, vegdist(spe_sqrt, method = 'bray')) |> 
  fviz_silhouette() 

aa <- isotab(iso.1) # diagnostic species summary
aa # display of the summary

#Isopam classification based on diagnostic species and their weighing
iso.2 <- isopam(spe_sqrt, c.max = 5, l.max = 3, sieve = T) 

# Silhouette plot for the resulting classification
silhouette(iso.2$flat$lev.1, vegdist(spe_sqrt, method = 'bray')) |> 
  fviz_silhouette()

bb <- isotab(iso.2, level = 1)# diagnostic species summary
bb  # display of the summary

# PCoA to display the clustering of both Isopam analyses with different setting of the sieve parameter in ordination space
pcoa.1 <- capscale(spe_sqrt~1, dist = "bray", sqrt.dist = F)

# extract scores for sites + add cluster identity
spe_pcoa <- scores(pcoa.1, choices = 1:2, display = "sites") |> 
  as_tibble(rownames = 'releve_nr') |> 
  mutate(cluster = iso.1$flat$lev.1, releve_nr = as.numeric(releve_nr)) |> 
  left_join(head)
         #, iso2 = iso.2$flat$level.1) |> 
 # pivot_longer(iso1:iso2, values_to = 'cluster', names_to = 'clustering') # to long format to plot both plots at once

# select points to draw convex hulls around
spe_hull <- spe_pcoa |> 
  group_by(cluster) |> 
  slice(chull(MDS1, MDS2)) # select only points at the border of the convex hull

# ggplot
spe_pcoa |> 
  ggplot(aes(MDS1, MDS2, color = as.factor(cluster), fill = as.factor(cluster)))+
  geom_text(aes(label = releve_nr), pch = 21)+
  geom_polygon(data = spe_hull, alpha = 0.2)+
  theme_bw()+
  theme(legend.position = 'bottom')

# save final classification -----------------------------------------------

iso.1$flat$level.1 |> 
  as_tibble(rownames = 'releve_nr') |> 
  rename('cluster' ='value') |> 
  mutate(cluster = case_when(cluster == 1 ~ 'Armerion', 
                             cluster == 2 ~ 'Corynephorion', 
                             cluster == 3 ~ 'Festucion valesiacae')) |> 
  write_csv('data/classification.csv')

