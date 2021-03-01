# conduct simulation using human data Sal.Heidelberg tested by AMC and CEP

#----- Set known parameters -----
n <- 500

Beta1 <- c(-0.5273, 3.8115)
Beta2 <- c(0.5793, 4.1723)
centerDf <- data.frame(
  x = c(Beta1[1], Beta1[2], Beta1[1], Beta1[2]), 
  y = c(Beta2[1], Beta2[1], Beta2[2], Beta2[2])
)
Sigma <- c(0.5012, 0.6694)
rho <- 0.6306
CovMat <- matrix(c(Sigma[1]^2, rho * Sigma[1] * Sigma[2], 
                   rho * Sigma[1] * Sigma[2], Sigma[2]^2), nrow = 2)

p1 <- 0.1145
p2 <- 0.113
# delta <- 0.0974
delta <- 0.05040436
# delta <- 0
p_nn <- (1-p1) * (1-p2) + delta
p_nr <- p2 * (1-p1) - delta
p_rn <- p1 * (1-p2) - delta
p_rr <- p1 * p2 + delta

#----- Generate bivariate points -----
library(MASS)
library(dplyr)
library(tidyr)
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
      var1.ct = ifelse(var1.lt <= -1, -1, ifelse(var1.lt > 5, 1, 0)),
      var2.ct = ifelse(var2.lt <= -1, -1, ifelse(var2.lt > 5, 1, 0)), 
      var1.y = ifelse(var1.lt <= -1, -1, ifelse(var1.lt > 5, 5, ceiling(var1.lt))), 
      var2.y = ifelse(var2.lt <= -1, -1, ifelse(var2.lt > 5, 5, ceiling(var2.lt)))
    ) %>% 
    mutate(
      var1.concl = ifelse(var1.y >= 5, 1, 0),
      var2.concl = ifelse(var2.y >= 5, 1, 0)
    ) %>% 
    mutate(
      var1.lb = ifelse(var1.ct == -1, -999, ifelse(var1.ct == 1, var1.y, var1.y-1)), 
      var1.ub = ifelse(var1.ct == 1, 999, var1.y), 
      var2.lb = ifelse(var2.ct == -1, -999, ifelse(var2.ct == 1, var2.y, var2.y-1)), 
      var2.ub = ifelse(var2.ct == 1, 999, var2.y)
    )
  write.csv(obsSim, paste0("/work/STAT/minz/mdr/simulation/obsSim/obsSim_delta_", delta, "/sim", i, ".csv"), row.names = F)
}

#------ Plot simulated obs -----
library(gridExtra)
library(ggplot2)
p1 <- obsSim %>% ggplot(., aes(x = var1.lt, y = var2.lt)) + 
  geom_point() + 
  geom_point(data = centerDf, aes(x = x, y = y), color = "red", cex = 3, shape = 17) +
  geom_hline(yintercept = Beta2, linetype = "dashed", color = "red", alpha = 0.3) +
  geom_vline(xintercept = Beta1, linetype = "dashed", color = "red", alpha = 0.3) +
  annotate("text", x = Beta1[1], y = -2, label = "paste(beta[1], \" [1]\")", parse = T, col = "red", size = 5) +
  annotate("text", x = Beta1[2], y = -2, label = "paste(beta[1], \" [2]\")", parse = T, col = "red", size = 5) +
  annotate("text", y = Beta2[1], x = -2, label = "paste(beta[2], \" [1]\")", parse = T, col = "red", size = 5) +
  annotate("text", y = Beta2[2], x = -2, label = "paste(beta[2], \" [2]\")", parse = T, col = "red", size = 5) +
  geom_errorbar(data = centerDf[1, ] + 1.5, aes(x = x, y = y, ymin = -1.5, ymax = 2.7), width = 0.1, col = "blue") +
  geom_errorbar(data = centerDf[1, ] + 2.3, aes(x = x, y = y, xmin = -1.8, xmax = 0.8), width = 0.1, col = "blue") +
  annotate("text", x = Beta1[1], y = 3.1, label = "sigma[1]", parse = T, col = "blue", size = 5) +
  annotate("text", y = Beta2[1], x = 1.3, label = "sigma[2]", parse = T, col = "blue", size = 5) +
  scale_x_continuous(name="Drug 1 latent MIC", breaks = -2:6, limits = c(-2, 6)) +
  scale_y_continuous(name="Drug 2 latent MIC", breaks = -2:6, limits = c(-2, 6)) +
  theme_bw() +
  theme(text = element_text(size = 15)) 

p2 <- obsSim %>% ggplot(., aes(x = var1.y, y = var2.y)) +
  geom_jitter(width = 0.2, height = 0.2) +
  scale_x_continuous(name="Drug 1 censored MIC", breaks = -2:6, limits = c(-2, 6)) +
  scale_y_continuous(name="Drug 1 censored MIC", breaks = -2:6, limits = c(-2, 6)) +
  theme_bw() +
  theme(text = element_text(size = 15)) 

grid.arrange(p1, p2, ncol = 2)

