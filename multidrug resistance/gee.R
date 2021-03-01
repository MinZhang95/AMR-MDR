library(dplyr)
library(ggplot2)
library(gee)

gee_obs <- obs %>% 
  mutate(id = 1:nrow(obs)) %>% 
  select(id, contains(".concl")) %>% 
  tidyr::gather(., key = "drug", value = "concl", var1.concl:var2.concl) %>% 
  arrange(id)

mod <- gee(concl ~ drug, family = binomial(link = "logit"), 
           data = gee_obs, id = id, corstr = "unstructured")
summary(mod)

cor(obs$var1.concl, obs$var2.concl, method = "pearson")
cor(obs$var1.concl, obs$var2.concl, method = "spearman")
cohen.kappa(obs %>% select(contains(".concl")))


0.0974 / (sqrt(0.1145*(1-0.1145))*sqrt(0.113*(1-0.113)))
