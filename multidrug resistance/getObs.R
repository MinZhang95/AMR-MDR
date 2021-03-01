library(dplyr)
library(stringr)
# Intermediate would be considered as susceptible !!!
# ----- human -----
population = "human"
genus = "Salmonella"
# humanDat <- read.csv("/work/STAT/minz/mdr/data/humanDatWithMonth.csv")

serotype = "Newport"
antibiotics = c("AMC", "CEP")

obs <- humanDat %>% 
  filter(Genus == genus, Serotype == serotype) %>% 
  dplyr::select(contains(antibiotics), -contains("Pred")) %>% 
  rename_at(vars(starts_with(antibiotics[1])), funs(str_replace(., antibiotics[1], "var1"))) %>% 
  rename_at(vars(starts_with(antibiotics[2])), funs(str_replace(., antibiotics[2], "var2"))) %>% 
  tidyr::drop_na() %>% 
  mutate(var1.y = round(log2(var1), 0), var2.y = round(log2(var2), 0)) %>% 
  mutate(var1.concl = ifelse(var1.Concl == "R", 1, 0), 
         var2.concl = ifelse(var2.Concl == "R", 1, 0)) %>% 
  mutate(var1.ct = ifelse(var1.Sign == "<=", -1, ifelse(var1.Sign == "=", 0, 1)), 
         var2.ct = ifelse(var2.Sign == "<=", -1, ifelse(var2.Sign == "=", 0, 1))) %>% 
  mutate(var1.lb = ifelse(var1.ct == -1, -999, ifelse(var1.ct == 1, var1.y, var1.y-1)), 
         var1.ub = ifelse(var1.ct == 1, 999, var1.y), 
         var2.lb = ifelse(var2.ct == -1, -999, ifelse(var2.ct == 1, var2.y, var2.y-1)), 
         var2.ub = ifelse(var2.ct == 1, 999, var2.y)) %>% 
  dplyr::select(var1.ct, var1.y, var1.concl, var2.ct, var2.y, var2.concl, 
                var1.lb, var1.ub, var2.lb, var2.ub)

serotype <- ifelse(serotype == "I 4,[5],12:i:-", "I4512i", serotype)
write.csv(obs, paste0("data/candidate eg/", population, "_", serotype, "_", antibiotics[1], "_", antibiotics[2], ".csv"), row.names = F)

# ----- animal and retail -----
population = "animal"
genus = "Salmonella"
serotype = "Montevideo"
antibiotics = c("CEP", "CHL")
(logcutoffs = antibiotics %>% 
  lapply(., function(x) read.csv("abb_full_cutoff.csv") 
         %>% filter(abb == x) %>% pull(resistant) %>% log2()) %>% unlist())

obs <- get(paste0(population, "Dat")) %>%
  filter(Genus == genus, Serotype == serotype) %>%
  dplyr::select(contains(antibiotics)) %>%
  rename_at(vars(starts_with(antibiotics[1])), funs(str_replace(., antibiotics[1], "var1"))) %>% 
  rename_at(vars(starts_with(antibiotics[2])), funs(str_replace(., antibiotics[2], "var2"))) %>% 
  tidyr::drop_na() %>% 
  mutate(var1.y = round(log2(var1), 0), var2.y = round(log2(var2), 0)) %>% 
  mutate(var1.concl = ifelse(var1.y >= logcutoffs[1], 1, 0), 
         var2.concl = ifelse(var2.y >= logcutoffs[2], 1, 0)) %>% 
  mutate(var1.concl = ifelse(var1.Sign == "<=", 0, var1.concl), 
         var2.concl = ifelse(var2.Sign == "<=", 0, var2.concl)) %>% 
  mutate(var1.concl = ifelse(var1.Sign == ">",  1, var1.concl), 
         var2.concl = ifelse(var2.Sign == ">",  1, var2.concl)) %>% 
  mutate(var1.ct = ifelse(var1.Sign == "<=", -1, ifelse(var1.Sign == "=", 0, 1)), 
         var2.ct = ifelse(var2.Sign == "<=", -1, ifelse(var2.Sign == "=", 0, 1))) %>% 
  mutate(var1.lb = ifelse(var1.ct == -1, -999, ifelse(var1.ct == 1, var1.y, var1.y-1)), 
         var1.ub = ifelse(var1.ct == 1, 999, var1.y), 
         var2.lb = ifelse(var2.ct == -1, -999, ifelse(var2.ct == 1, var2.y, var2.y-1)), 
         var2.ub = ifelse(var2.ct == 1, 999, var2.y)) %>% 
  dplyr::select(var1.ct, var1.y, var1.concl, var2.ct, var2.y, var2.concl, 
                var1.lb, var1.ub, var2.lb, var2.ub)

# serotype <- ifelse(serotype == "I 4,[5],12:i:-", "I4512i", serotype)
write.csv(obs, paste0("data/candidate eg/temp/", population, "_", serotype, "_", antibiotics[1], "_", antibiotics[2], ".csv"), row.names = F)
