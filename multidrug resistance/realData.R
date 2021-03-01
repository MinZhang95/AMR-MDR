library(dplyr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

datls <- list(N = nrow(obs),
              lower1 = obs$var1.lb, upper1 = obs$var1.ub, 
              lower2 = obs$var2.lb, upper2 = obs$var2.ub)
## init
p1 <- mean(obs$var1.concl)
p2 <- mean(obs$var2.concl)
delta <- nrow(obs[obs$var1.concl==1 & obs$var2.concl==1,]) / nrow(obs) - p1*p2
Beta1 <- c(mean(obs$var1.y[obs$var1.concl==0]), 
           mean(obs$var1.y[obs$var1.concl==1]))
Beta2 <- c(mean(obs$var2.y[obs$var2.concl==0]), 
           mean(obs$var2.y[obs$var2.concl==1]))
Sigma <- c(ifelse(sd(obs$var1.y)==0, 0.00001, sd(obs$var1.y)), 
           ifelse(sd(obs$var2.y)==0, 0.00001, sd(obs$var2.y)))
rho <- cor(obs$var1.y, obs$var2.y, method = "spearman")
rho <- ifelse(is.na(rho), 0.01, rho)
R <- matrix(c(1, rho, rho, 1), nrow = 2)
initf <- function() {
  list(
    p1 = p1, p2 = p2, delta = delta, 
    Beta1 = Beta1, Beta2 = Beta2, 
    Sigma = Sigma, R = R
  )
}
## run
ptm <- proc.time()
res <- stan(
  file = "./mod3.stan", 
  data = datls, 
  init = initf, 
  chains = 3, 
  iter = 10000, 
  pars = c("p1", "p2", "delta", 
           "Beta1", "Beta2", 
           "Sigma", "R", 
           "p_nn", "p_rn", "p_nr", "p_rr")
)
proc.time() - ptm

print(res)
