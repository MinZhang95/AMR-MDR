my_packages <- c("readxl", "openxlsx", "dplyr", "ggplot2", "stringr", "Hmisc", "gridExtra", "openxlsx")
lapply(my_packages, library, character.only = TRUE)

#----- Read and clean human CDC data -----
humanDat <- read.csv("./data/humanDatWithMonth.csv", stringsAsFactors = FALSE)
humanDat <- humanDat %>% 
  rename(Year = Data.Year, Month = Month.of.collection) %>% 
  rename_at(vars(ends_with(".Equiv")), funs(str_replace(., ".Equiv", ".Sign"))) %>% 
  rename_at(vars(ends_with(".Rslt")), funs(str_replace(., ".Rslt", ""))) %>% 
  arrange(Year, Month)
## To match the drug abbr in animal and retail datasets
humanDat <- humanDat %>% 
  rename_at(vars(starts_with("AUG")), funs(str_replace(., "AUG", "AMC"))) %>% 
  rename_at(vars(starts_with("AZM")), funs(str_replace(., "AZM", "AZI")))

#----- Read and clean food animal USDA data -----
animalDat1 <- read.csv("../corAmr/data/NARMS-Bacterial-Isolate-Level-Data-Food-Producing-Animals-USDA-ARS-1997-2005.csv", stringsAsFactors = F)
animalDat2 <- read.csv("../corAmr/data/NARMS-Bacterial-Isolate-Level-Data-Food-Producing-Animals-USDA-ARS-2006-2013.csv", stringsAsFactors = F)
animalDat <- rbind(animalDat1, animalDat2)
names(animalDat)[1:17] <- stringr::str_to_title(names(animalDat)[1:17])
animalDat <- animalDat %>% 
  rename_at(vars(ends_with(" Sign")), funs(str_replace(., " Sign", ".Sign"))) %>% 
  rename_at(vars(starts_with("Str")), funs(str_replace(., "Str", "STR"))) %>% 
  mutate(Genus = ifelse(Genus=="S", "Salmonella", Genus)) %>% 
  arrange(Year, Month)

#----- Read and clean retail FDA data -----
retailDat <- read.csv("../corAmr/data/NARMS-Bacterial-Isolate-Level-Data-Retail-Meats.csv", stringsAsFactors = F)
retailDat <- retailDat %>% 
  rename(Genus = GENUS, Species = SPECIES, Serotype = SEROTYPE) %>% 
  rename_at(vars(ends_with(" Sign")), funs(str_replace(., " Sign", ".Sign"))) %>% 
  rename_at(vars(starts_with("Str")), funs(str_replace(., "Str", "STR"))) %>% 
  mutate(Genus = ifelse(Genus=="S", "Salmonella", Genus)) %>% 
  arrange(Year, Month)

#----- 
humanAnti <- humanDat %>% dplyr::select(ends_with("Sign")) %>% colnames() %>% gsub(".Sign", "", .)
retailAnti <- retailDat %>% dplyr::select(ends_with("Sign")) %>% colnames() %>% gsub(".Sign", "", .)
animalAnti <- animalDat %>% dplyr::select(ends_with("Sign")) %>% colnames() %>% gsub(".Sign", "", .)

#----- 
abb2Full <- readxl::read_xlsx("../corAmr/data/cutoff3.xlsx") %>%
  select(`CDC NARMS Ab`, `CDC NARMS description`) %>% 
  rename(abb = `CDC NARMS Ab`, full = `CDC NARMS description`) %>% 
  filter(!is.na(abb)) %>% as.data.frame()

abbAll <- read.csv("abb_full_cutoff.csv") %>% 
  filter(class != "") %>%
  pull(abb) %>% as.character()
varsIn2Classes <- expand.grid(var1 = abbAll, var2 = abbAll) %>% 
  filter(var1 != var2) %>% 
  left_join(., read.csv("abb_full_cutoff.csv") %>% select(abb, class), c("var1" = "abb")) %>% 
  left_join(., read.csv("abb_full_cutoff.csv") %>% select(abb, class), c("var2" = "abb")) %>%
  filter(class.x != class.y)
