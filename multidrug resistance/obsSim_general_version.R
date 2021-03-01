library(MASS)
library(dplyr)
library(tidyr)

pair <- "human_I4512i_AMP_NAL"
# (cutoff1 <- read.csv("abb_full_cutoff.csv") %>% filter(abb == "AMP") %>% pull(resistant) %>% log2())
# (cutoff2 <- read.csv("abb_full_cutoff.csv") %>% filter(abb == "NAL") %>% pull(resistant) %>% log2())
cutoff1 <- 5; cutoff2 <- 5

dat <- read.csv(paste0("/work/STAT/minz/mdr/data/", pair, ".csv"))
res <- read.csv(paste0("/work/STAT/minz/mdr/application/res/", pair, ".csv"))

n <- dat %>% nrow()
var1.y.floor <- dat %>% pull(var1.y) %>% min()
var1.y.ceiling <- dat %>% pull(var1.y) %>% max()
var2.y.floor <- dat %>% pull(var2.y) %>% min()
var2.y.ceiling <- dat %>% pull(var2.y) %>% max()

Beta1 <- res %>% filter(grepl("Beta1", X)) %>% pull(mean)
Beta2 <- res %>% filter(grepl("Beta2", X)) %>% pull(mean)
Sigma <- res %>% filter(grepl("Sigma", X)) %>% pull(mean)
rho <- res %>% filter(grepl("R", X)) %>% pull(mean) %>% .[2]
CovMat <- matrix(c(Sigma[1]^2, rho * Sigma[1] * Sigma[2], 
                   rho * Sigma[1] * Sigma[2], Sigma[2]^2), nrow = 2)
p <- res %>% filter(grepl("p", X)) %>% pull(mean)
p1 <- p[1]; p2 <- p[2]
delta <- res %>% filter(grepl("delta", X)) %>% pull(mean)
p_nn <- (1-p1) * (1-p2) + delta
p_nr <- p2 * (1-p1) - delta
p_rn <- p1 * (1-p2) - delta
p_rr <- p1 * p2 + delta


#----- Generate bivariate points -----
if (!dir.exists(paste0("/work/STAT/minz/mdr/simulation/obsSim/", pair, "/"))) {
  dir.create(file.path(paste0("/work/STAT/minz/mdr/simulation/obsSim/", pair, "/")))
}
set.seed(920924)
for (i in 1:100) {
  n_cell <- rmultinom(n, size = 1, prob = c(p_nn, p_nr, p_rn, p_rr)) %>% apply(., 1, sum)
  n_nn <- n_cell[1]
  n_nr <- n_cell[2]
  n_rn <- n_cell[3]
  n_rr <- n_cell[4]
  NA_mat <- matrix(c(NA, NA), ncol = 2)
  obsSim <- rbind(
    ifelse(n_nn == 0, list(NA_mat), list(mvrnorm(n_nn, mu = c(Beta1[1], Beta2[1]), Sigma = CovMat))) %>% do.call(rbind, .), 
    ifelse(n_nr == 0, list(NA_mat), list(mvrnorm(n_nr, mu = c(Beta1[1], Beta2[2]), Sigma = CovMat))) %>% do.call(rbind, .), 
    ifelse(n_rn == 0, list(NA_mat), list(mvrnorm(n_rn, mu = c(Beta1[2], Beta2[1]), Sigma = CovMat))) %>% do.call(rbind, .), 
    ifelse(n_rr == 0, list(NA_mat), list(mvrnorm(n_rr, mu = c(Beta1[2], Beta2[2]), Sigma = CovMat))) %>% do.call(rbind, .)
  ) %>% as.data.frame() %>% drop_na()
  colnames(obsSim) <- c("var1.lt", "var2.lt")
  obsSim <- obsSim %>% 
    mutate(
      var1.ct = ifelse(var1.lt <= var1.y.floor, -1, ifelse(var1.lt > var1.y.ceiling, 1, 0)),
      var2.ct = ifelse(var2.lt <= var2.y.floor, -1, ifelse(var2.lt > var2.y.ceiling, 1, 0)), 
      var1.y = ifelse(var1.lt <= var1.y.floor, var1.y.floor, ifelse(var1.lt > var1.y.ceiling, var1.y.ceiling, ceiling(var1.lt))), 
      var2.y = ifelse(var2.lt <= var2.y.floor, var2.y.floor, ifelse(var2.lt > var2.y.ceiling, var2.y.ceiling, ceiling(var2.lt)))
    ) %>% 
    mutate(
      var1.concl = ifelse(var1.y >= cutoff1, 1, 0),
      var2.concl = ifelse(var2.y >= cutoff2, 1, 0)
    ) %>% 
    mutate(
      var1.lb = ifelse(var1.ct == -1, -999, ifelse(var1.ct == 1, var1.y, var1.y-1)), 
      var1.ub = ifelse(var1.ct == 1, 999, var1.y), 
      var2.lb = ifelse(var2.ct == -1, -999, ifelse(var2.ct == 1, var2.y, var2.y-1)), 
      var2.ub = ifelse(var2.ct == 1, 999, var2.y)
    )
  write.csv(obsSim, paste0("/work/STAT/minz/mdr/simulation/obsSim/", pair, "/sim", i, ".csv"), row.names = F)
}

