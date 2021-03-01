library(reshape2)
library(gridExtra)
library(ggplot2)
library(stringr)
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat) 
}

pop <- "human"

popSeroRange <- get(paste0(pop, "Dat")) %>% 
  count(Genus, Serotype) %>% arrange(-n) %>% 
  filter(Genus=="Salmonella") %>% as.data.frame() %>% 
  filter(if (pop=="retail") n>500 else n>1000)

for (sero in 1:nrow(popSeroRange)) {
  
  popSub <- get(paste0(pop, "Dat")) %>% filter(Genus == "Salmonella", Serotype == popSeroRange$Serotype[sero])
  
  cormat_all <- popSub %>% select(get(paste0(pop, "Anti"))) %>% 
    select_if(~sum(!is.na(.)) > 1) %>% 
    as.matrix() %>% log2() %>% round(., 0) %>% 
    Hmisc::rcorr(type = "spearman")
  
  upper.cor <- get_upper_tri(cormat_all$r)
  
  # heatPlot <- ggplot(data = melt(upper.cor, na.rm = T), aes(Var2, Var1, fill = value)) +
  #   geom_tile(color = "white") + 
  #   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
  #                        midpoint = 0, limit = c(-1, 1), space = "Lab", 
  #                        name = "Spearman\nCorrelation") + 
  #   theme_minimal() + 
  #   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 270, vjust=0.25), 
  #         axis.title.y = element_blank()) + 
  #   scale_y_discrete(position = "right") +
  #   coord_fixed()
  # # heatPlot
  
  n <- melt(get_upper_tri(cormat_all$n)) %>% rename(n = value) %>% 
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
  p <- melt(get_upper_tri(cormat_all$P)) %>% rename(p = value) %>% 
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
  
  cor.rank <- melt(upper.cor) %>% filter(Var1 != Var2, !is.na(value)) %>% arrange(-value) %>% 
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>% 
    left_join(., n, c("Var1" = "Var1", "Var2" = "Var2")) %>% 
    left_join(., p, c("Var1" = "Var1", "Var2" = "Var2")) %>% 
    mutate(p = round(p, 3)) %>% 
    left_join(., abb2Full, c("Var1" = "abb")) %>% 
    left_join(., abb2Full, c("Var2" = "abb")) %>% 
    rename(Name1 = full.x, Name2 = full.y) %>% 
    filter(n>=100) %>% 
    # to make sure two drugs in two classes
    left_join(., read.csv("abb_full_cutoff.csv") %>% select(abb, class), c("Var1" = "abb")) %>% 
    left_join(., read.csv("abb_full_cutoff.csv") %>% select(abb, class), c("Var2" = "abb")) %>%
    filter(class.x != class.y)
  
  if (nrow(cor.rank) == 0) {next}
  
  scatterPlot <- list()
  id <- 1
  for (i in 1:nrow(cor.rank)) {
    pair <- cor.rank[i, ]
    
    # to make sure the data set not oversized
    if (pair$n > 3000) {next}
    
    antibiotics <- c(pair$Var1, pair$Var2)
    cutoffs <- antibiotics %>% 
      lapply(., function(x) read.csv("abb_full_cutoff.csv") 
             %>% filter(abb == x) %>% pull(resistant) %>% log2()) %>% unlist()
    
    scatterDat <- popSub %>% 
      select(Year, contains(antibiotics, ignore.case = F), -contains("Pred"), -contains("Sign")) %>% 
      rename_all(funs(str_replace(., antibiotics[1], "Var1"))) %>% 
      rename_all(funs(str_replace(., antibiotics[2], "Var2"))) %>% 
      # {if (pop != "human") select(., Year, Var1, Var2)} %>% 
      mutate(Var1 = round(log2(Var1), 0), Var2 = round(log2(Var2), 0)) 
    
    if (pop != "human") {
      scatterDat$Var1.Concl = ifelse(scatterDat$Var1 >= cutoffs[1], 1, 0)
      scatterDat$Var2.Concl = ifelse(scatterDat$Var2 >= cutoffs[2], 1, 0)
    } else {
      scatterDat$Var1.Concl = ifelse(scatterDat$Var1.Concl == "R", 1, 0)
      scatterDat$Var2.Concl = ifelse(scatterDat$Var2.Concl == "R", 1, 0)
    }
    
    scatterDat <- scatterDat %>% 
      tidyr::drop_na() %>% 
      mutate(Sum.Concl = Var1.Concl + Var2.Concl) %>% 
      mutate(MDR = ifelse(Sum.Concl == 2, "MDR", ifelse(Sum.Concl == 1, "SDR", "NR"))) 
    
    phi.test <- cor.test(scatterDat$Var1.Concl, scatterDat$Var2.Concl)
    scatterDat_nn <- scatterDat %>% filter(MDR == "NR")
    if (nrow(scatterDat_nn) > 1) {
      rho.test <- cor.test(scatterDat_nn$Var1, scatterDat_nn$Var2, method = "spearman", exact = F)
    } else {
      next
    }
    
    if (is.na(phi.test$estimate)) {next}
    
    scatterPlot[[id]] <- scatterDat %>% 
      ggplot(., aes(x=Var1, y=Var2, color = MDR)) +
      geom_jitter(width = 0.15, height = 0.15) +
      scale_color_manual(breaks = c("MDR", "SDR", "NR"),
                         values=c("red", "orange", "black")) + 
      labs(title = paste0(popSeroRange$Serotype[sero], " in ", pop, " ", 
                          min(scatterDat$Year), "~", max(scatterDat$Year), ": n=", 
                          pair$n, "\nrho=", round(rho.test$estimate, 2), " (p=", round(rho.test$p.value, 2), 
                          "), phi=", round(phi.test$estimate, 2), " (p=", round(phi.test$p.value, 2), ")"),
           x = paste0(pair$Name1, " (", pair$Var1, ")"), 
           y = paste0(pair$Name2, " (", pair$Var2, ")")) +
      theme(legend.position = "none", 
            text = element_text(size = 17))
    id <- id + 1
  }  
  
  seroFileName <- ifelse(popSeroRange$Serotype[sero]=="I 4,[5],12:i:-", "I4512i", popSeroRange$Serotype[sero])
  pdf(paste0("./scatterPlots/", pop, "_", seroFileName, ".pdf"), onefile = TRUE)
  for (id in 1:length(scatterPlot)) {
    print(scatterPlot[[id]] )
  }
  dev.off()
}
