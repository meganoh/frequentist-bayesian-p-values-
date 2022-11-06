library(brms)
library(rstanarm)
library(bayestestR)
library(BayesFactor)
library(cmdstanr)
library(emmeans)
library(parallel)
library(bridgesampling)
library(tidyverse)
library(extraDistr)
#read in results 
group_2_n_100 <- bind_rows(results_2groups, .id = "id")

alpha.p <- 0.05 #set threshold for significance 
iter = 1000

group_2_n_20 <- group_2_n_20 %>% 
  mutate(significance = if_else(pval<alpha.p, "sig", "nonsig")) #p<.05 --> rej null TRUE, do not rej FALSE
group_2_n_20_typeIerr <- group_2_n_20 %>% 
  group_by(test) %>% 
  count(significance)
group_2_n_20_typeIerr <- group_2_n_20_typeIerr %>% pivot_wider(names_from = "significance", values_from = "n") %>% 
  mutate(rate_typeIerr = sig/iter, n = 20)

group_2_n_30 <- group_2_n_30 %>% 
  mutate(significance = if_else(pval<alpha.p, "sig", "nonsig")) #p<.05 --> rej null TRUE, do not rej FALSE
group_2_n_30_typeIerr <- group_2_n_30 %>% 
  group_by(test) %>% 
  count(significance)
group_2_n_30_typeIerr <- group_2_n_30_typeIerr %>% pivot_wider(names_from = "significance", values_from = "n") %>% 
  mutate(rate_typeIerr = sig/iter, n = 30)

group_2_n_50 <- group_2_n_50 %>% 
  mutate(significance = if_else(pval<alpha.p, "sig", "nonsig")) #p<.05 --> rej null TRUE, do not rej FALSE
group_2_n_50_typeIerr <- group_2_n_50 %>% 
  group_by(test) %>% 
  count(significance)
group_2_n_50_typeIerr <- group_2_n_50_typeIerr %>% pivot_wider(names_from = "significance", values_from = "n") %>% 
  mutate(rate_typeIerr = sig/iter, n = 50)

group_2_n_100 <- group_2_n_100 %>% 
  mutate(significance = if_else(pval<alpha.p, "sig", "nonsig")) #p<.05 --> rej null TRUE, do not rej FALSE
group_2_n_100_typeIerr <- group_2_n_100 %>% 
  group_by(test) %>% 
  count(significance)
group_2_n_100_typeIerr <- group_2_n_100_typeIerr %>% pivot_wider(names_from = "significance", values_from = "n") %>% 
  mutate(rate_typeIerr = sig/iter, n = 100)

typeIerr_2groups <- bind_rows(group_2_n_20_typeIerr, group_2_n_30_typeIerr, group_2_n_50_typeIerr, group_2_n_100_typeIerr, .id = "id")

save(typeIerr_2groups, file = "typeIerr_2groups.rda")


#3groups
group_3_n_100 <- bind_rows(results_3groups, .id = "id")

alpha.p <- 0.05 #set threshold for significance 
iter = 1000

group_3_n_20 <- group_3_n_20 %>% 
  mutate(significance = if_else(pval<alpha.p, "sig", "nonsig")) #p<.05 --> rej null TRUE, do not rej FALSE
group_3_n_20_typeIerr <- group_3_n_20 %>% 
  group_by(test) %>% 
  count(significance)
group_3_n_20_typeIerr <- group_3_n_20_typeIerr %>% pivot_wider(names_from = "significance", values_from = "n") %>% 
  mutate(rate_typeIerr = sig/iter, n = 20)

group_3_n_30 <- group_3_n_30 %>% 
  mutate(significance = if_else(pval<alpha.p, "sig", "nonsig")) #p<.05 --> rej null TRUE, do not rej FALSE
group_3_n_30_typeIerr <- group_3_n_30 %>% 
  group_by(test) %>% 
  count(significance)
group_3_n_30_typeIerr <- group_3_n_30_typeIerr %>% pivot_wider(names_from = "significance", values_from = "n") %>% 
  mutate(rate_typeIerr = sig/iter, n = 30)

group_3_n_50 <- group_3_n_50 %>% 
  mutate(significance = if_else(pval<alpha.p, "sig", "nonsig")) #p<.05 --> rej null TRUE, do not rej FALSE
group_3_n_50_typeIerr <- group_3_n_50 %>% 
  group_by(test) %>% 
  count(significance)
group_3_n_50_typeIerr <- group_3_n_50_typeIerr %>% pivot_wider(names_from = "significance", values_from = "n") %>% 
  mutate(rate_typeIerr = sig/iter, n = 50)

group_3_n_100 <- group_3_n_100 %>% 
  mutate(significance = if_else(pval<alpha.p, "sig", "nonsig")) #p<.05 --> rej null TRUE, do not rej FALSE
group_3_n_100_typeIerr <- group_3_n_100 %>% 
  group_by(test) %>% 
  count(significance)
group_3_n_100_typeIerr <- group_3_n_100_typeIerr %>% pivot_wider(names_from = "significance", values_from = "n") %>% 
  mutate(rate_typeIerr = sig/iter, n = 100)

typeIerr_3groups <- bind_rows(group_3_n_20_typeIerr, group_3_n_30_typeIerr, group_3_n_50_typeIerr, group_3_n_100_typeIerr, .id = "id")

save(typeIerr_3groups, file = "typeIerr_3groups.rda")
