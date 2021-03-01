#----- Set known parameters -----
library(dplyr)
pair <- "human_Heidelberg_AMP_CEP"
antibiotics <- c(strsplit(pair, "_")[[1]][3], strsplit(pair, "_")[[1]][4])
cutoffs <- antibiotics %>% 
  lapply(., function(x) read.csv("/work/STAT/minz/mdr/abb_full_cutoff.csv") 
         %>% filter(abb == x) %>% pull(resistant) %>% log2()) %>% unlist()
n <- read.csv(paste0("/work/STAT/minz/mdr/data/", pair, ".csv")) %>% nrow()

res <- read.csv(paste0("/work/STAT/minz/mdr/application/res/", pair, ".csv"))

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
cor_bin <- res %>% filter(grepl("cor_bin", X)) %>% pull(mean)

parGiven <- as.data.frame(
  matrix(c(p1, p2, delta, 
           Beta1[1], Beta1[2], Beta2[1], Beta2[2], Sigma[1], Sigma[2], 
           1, rho, rho, 1, 
           p_nn, p_rn, p_nr, p_rr, 
           cor_bin), 
         ncol = 1)
)

#----- Get metrics -----
library("abind")
res.ls <- list()
id <- 1
for (i in 1:100) {
  fname <- paste0("/work/STAT/minz/mdr/simulation/resSim/", pair, "/sim", i, ".csv")
  if (file.exists(fname)) {
    res.df <- read.csv(fname)
    if (res.df[res.df$X == "cor_bin", "Rhat"] < 1.03) {
      res.ls[[id]] <- res.df
      id <- id + 1
    }
  }
}

# res.matrix <- abind(res.ls, along=3)
# apply(res.matrix, c(1,2), mean)
# apply(res.matrix, c(1,2), sd)
res.ls.mean <- lapply(res.ls, function(x) as.data.frame(x[1:18, "mean"]))
res.ls.mean.stack <- abind(res.ls.mean, along = 3)
mean <- apply(res.ls.mean.stack, c(1,2), mean)
sd <- apply(res.ls.mean.stack, c(1,2), sd)

res.ls.mean.diff <- lapply(res.ls, function(x) abs(as.data.frame(x[1:18, "mean"]) - parGiven))
res.ls.mean.diff.stack <- abind(res.ls.mean.diff, along = 3)
mae <- apply(res.ls.mean.diff.stack, c(1,2), mean) # mean abs error

res.ls.mean.diff2 <- lapply(res.ls, function(x) (as.data.frame(x[1:18, "mean"]) - parGiven)^2)
res.ls.mean.diff2.stack <- abind(res.ls.mean.diff2, along = 3)
rmse <- sqrt(apply(res.ls.mean.diff2.stack, c(1,2), mean)) # root of mean squared error

summ_all <- data.frame(
  mean = mean, sd = sd, mae = mae, rmse = rmse
)

row.names(summ_all) <- res.ls[[1]][1:18, "X"]
colnames(summ_all) <- c("mean", "sd", "mae", "rmse")

# correlation in binary conclusion

phi.df <- data.frame(
  cor_bin = double(),
  phiBayes = double(), 
  phiGee = double(), 
  rhoNaive = double()
)

id <- 1
for (i in 1:100) {
  resname <- paste0("/work/STAT/minz/mdr/simulation/resSim/", pair, "/sim", i, ".csv")
  obsname <- paste0("/work/STAT/minz/mdr/simulation/obsSim/", pair, "/sim", i, ".csv")
  if (file.exists(resname)) {
    res.df <- read.csv(resname)
    obs.df <- read.csv(obsname)
    
    if (res.df[res.df$X == "cor_bin", "Rhat"] < 1.03) {
      obs.p1 <- mean(obs.df$var1.concl)
      obs.p2 <- mean(obs.df$var2.concl)
      obs.p1p2 <- nrow(obs.df[obs.df$var1.concl==1 & obs.df$var2.concl==1,]) / nrow(obs.df)
      phiGee <- (obs.p1p2 - obs.p1 * obs.p2) / (sqrt(obs.p1 * (1 - obs.p1)) * sqrt(obs.p2 * (1 - obs.p2)))
      
      phiBayes <- res.df[res.df$X == "cor_bin", "mean"]
      
      a <- obs.df %>% filter(var1.concl == 0, var2.concl == 0) %>% pull(var1.y)
      b <- obs.df %>% filter(var1.concl == 0, var2.concl == 0) %>% pull(var2.y)
      rhoNaive <- cor(a, b)
      
      phi.df[id, ] <- c(cor_bin, phiBayes, phiGee, rhoNaive)
      id <- id + 1
    }
  }
}

summ_phi <- data.frame(
  phiBayes = double(), 
  phiGee = double(), 
  rhoNaive = double()
)

summ_phi[1, ] <- phi.df %>% apply(., 2, mean) %>% .[c("phiBayes", "phiGee", "rhoNaive")]
summ_phi[2, ] <- phi.df %>% apply(., 2, sd) %>% .[c("phiBayes", "phiGee", "rhoNaive")]

phi.df <- phi.df %>% 
  mutate(phiBayes_bias = phiBayes - cor_bin, 
         phiGee_bias = phiGee - cor_bin, 
         rhoNaive_bias = rhoNaive - rho)

summ_phi[3, ] <- c(mean(abs(phi.df$phiBayes_bias)), mean(abs(phi.df$phiGee_bias)), mean(abs(phi.df$rhoNaive_bias)))
summ_phi[4, ] <- c(sqrt(mean((phi.df$phiBayes_bias)^2)), sqrt(mean((phi.df$phiGee_bias)^2)), sqrt(mean((phi.df$rhoNaive_bias)^2)))

row.names(summ_phi) <- c("mean", "sd", "mae", "rmse")

summ_combined <- rbind(summ_all[-18, ], t(summ_phi)) %>% 
  mutate(truth = c(parGiven$V1, cor_bin, rho), .before = mean) 
rownames(summ_combined) <- c(rownames(summ_all)[-18], rownames(t(summ_phi)))
summ_combined
