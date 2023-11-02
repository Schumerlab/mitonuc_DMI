library(tidyverse)
dists <- read.table("input_files/xbir_31M_dist.txt", sep = " ", col.names = c("id", "AA1", "AA2", "dist")) %>%
 mutate(replicate = sapply(id, FUN = function(x) str_split(x, pattern = fixed("."))[[1]][2]),
        template = sapply(id, FUN = function(x) str_split(str_split(x, pattern = fixed("."))[[1]][1], pattern = "-")[[1]][2]))
dists$replicate[is.na(dists$replicate)] <- "RaptorX"
dists$template[is.na(dists$template)] <- "RaptorX"

out <- group_by(dists, AA1, AA2, template) %>%
  summarize(mean_dist = mean(dist),
            sd_dist = sd(dist),
            se_dist = sd(dist)/sqrt(n()))

### New rounds of distance calculations

Y79dists <- read.table("input_files/xbir_79Y_dist.txt", sep = " ", col.names = c("id", "AA1", "AA2", "dist")) %>%
  mutate(replicate = sapply(id, FUN = function(x) str_split(x, pattern = fixed("."))[[1]][2]),
         template = sapply(id, FUN = function(x) str_split(str_split(x, pattern = fixed("."))[[1]][1], pattern = "-")[[1]][2]))
Y79dists$replicate[is.na(Y79dists$replicate)] <- "RaptorX"
Y79dists$template[is.na(Y79dists$template)] <- "RaptorX"

y79out <- group_by(Y79dists, AA1, AA2, template) %>%
  summarize(mean_dist = mean(dist),
            sd_dist = sd(dist),
            se_dist = sd(dist)/sqrt(n()))

### New rounds of distance calculations

L140dists <- read.table("input_files/xbir_140L_dist.txt", sep = " ", col.names = c("id", "AA1", "AA2", "dist")) %>%
  mutate(replicate = sapply(id, FUN = function(x) str_split(x, pattern = fixed("."))[[1]][2]),
         template = sapply(id, FUN = function(x) str_split(str_split(x, pattern = fixed("."))[[1]][1], pattern = "-")[[1]][2]))
L140dists$replicate[is.na(L140dists$replicate)] <- "RaptorX"
L140dists$template[is.na(L140dists$template)] <- "RaptorX"

L140out <- group_by(L140dists, AA1, AA2, template) %>%
  summarize(mean_dist = mean(dist),
            sd_dist = sd(dist),
            se_dist = sd(dist)/sqrt(n()))

bind_rows(out, y79out, L140out) %>% View()
