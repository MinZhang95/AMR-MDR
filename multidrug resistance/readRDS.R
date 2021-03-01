library(dplyr)
library(rstan)
obs <- read.csv("./data/human_typhi_KAN_SMX.csv")
(p1 <- mean(obs$var1.concl))
(p2 <- mean(obs$var2.concl))
(x <- nrow(obs[obs$var1.concl==1 & obs$var2.concl==1,]) / nrow(obs))
obs %>% count(var1.concl, var2.concl)

res <- readRDS("~/Downloads/human_typhi_KAN_SMX_mod3.rds")
(resTab <- res %>% summary() %>% .[["summary"]] %>% round(., 4))
traceplot(res, inc_warmup=T, c("p1", "p2"))
write.csv(resTab, "~/Downloads/resTab.csv")

cor(obs$var1.y, obs$var2.y)
cor(obs$var1.y, obs$var2.y, method = "spearman")

p1mean = resTab[1, 1]; p2mean = resTab[2, 1]
(deltaLower = -min(p1mean*p2mean, (1-p1mean)*(1-p2mean)))
(deltaUpper = min(p1mean,p2mean)-p1mean*p2mean)








obs <- read.csv("./data/human_I4512i_STR_TET.csv")
(p1 <- mean(obs$var1.concl))
(p2 <- mean(obs$var2.concl))
(x <- nrow(obs[obs$var1.concl==1 & obs$var2.concl==1,]) / nrow(obs))
obs %>% count(var1.concl, var2.concl)

res <- readRDS("~/Downloads/human_I4512i_STR_TET.rds")
res %>% summary() %>% .[["summary"]] %>% round(., 4)
traceplot(res, inc_warmup=T, c("R[1,2]", "delta"))


res <- readRDS("~/Downloads/sim33.rds")
traceplot(res, c("cor_bin"))
