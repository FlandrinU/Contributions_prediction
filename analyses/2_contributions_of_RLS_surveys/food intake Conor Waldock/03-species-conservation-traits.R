### script to get species addition information from fishbase including
### IUCN status, vulnerability, resiliance, toxicity, etc

library(tidyverse)
library(plyr)
library(rfishbase)
library(rredlist)

# read in master list of species
nuts <- read.csv('processed-data/nutrient-content/Spp_NutrientPred_REEF_FUTURESJune2021.csv')

sp_name <- nuts %>% dplyr::select(valid_name_FishBase)

list_fields()$table %>% unique
list_fields() %>% filter(table == 'refrens') %>% View
list_fields() %>% filter(table == 'species') %>% View
list_fields() %>% filter(columns == 'Conservation') # doesn't return anything

# search fields
all_fields <- list_fields()
all_columns <- list_fields() %>% .$columns
all_columns[grepl('cig', all_columns, ignore.case = T)]
all_fields[grepl('cig', all_fields$columns, ignore.case = T),]

table(references(species(sp_name$valid_name_FishBase, 'SpecCode')[[1]], 'Ciguatera'))

# fishbase data----

# fish vulnerability 
vulnerability <- species(sp_name$valid_name_FishBase, fields = c('Species','Vulnerability')) %>% plyr::rename(., c('Species' = 'valid_name_FishBase'))
              
# fish resilience
resilience <- stocks(sp_name$valid_name_FishBase,  fields = c('Species', 'Resilience')) %>% plyr::rename(., c('Species' = 'valid_name_FishBase'))

resilience_filtered <- resilience %>% unique() %>% na.omit() %>% group_by(valid_name_FishBase) %>% do(data.frame(., count = length(.$valid_name_FishBase))) %>% filter(count != 2) %>% dplyr::select(-count)

resilience_filtered <- left_join(sp_name, resilience_filtered)

table(resilience_filtered$Resilience)

# model fish common size
all_lengths    <- species(sp_name$valid_name_FishBase, fields = c('Species','Length', 'CommonLength', 'FamCode', 'GenCode')) %>% plyr::rename(., c('Species' = 'valid_name_FishBase'))
lapply(all_lengths, function(x) table(is.na(x)))

library(lme4)
library(lmtest)
library(DHARMa)
plot(log(all_lengths$Length)~ log(all_lengths$CommonLength))
# plot(lm(log(all_lengths$Length)~ log(all_lengths$CommonLength)))
model1 <- lmer(log(CommonLength) ~ log(Length) + (log(Length)|FamCode/GenCode), data = all_lengths)
model2 <- lmer(log(CommonLength) ~ log(Length) + (1|FamCode/GenCode), data = all_lengths)
model3 <- lmer(log(CommonLength) ~ log(Length) + (1|FamCode), data = all_lengths)
model4 <- lm(log(Length) ~ log(CommonLength), data = all_lengths)
lrtest(model1, model2) # model including genus level coding is better

# check residuals
model1_check <- simulateResiduals(model1)
plot(model1_check) # not perfect but adequate

# make predictions for new common lengths 
all_lengths$CommonLength_model <- all_lengths$CommonLength
all_lengths$CommonLength_model[is.na(all_lengths$CommonLength)] <- exp(predict(model1, newdata = all_lengths[is.na(all_lengths$CommonLength),], allow.new.levels=T))
plot(log(all_lengths$CommonLength_model)[order(is.na(all_lengths$CommonLength))]~ log(all_lengths$Length)[order(is.na(all_lengths$CommonLength))], 
     col = ifelse(is.na(all_lengths$CommonLength)[order(is.na(all_lengths$CommonLength))], 'red', 'blue'))


# check data properties before joining
nrow(resilience_filtered); nrow(vulnerability); nrow(all_lengths); nrow(sp_name)

# fish toxicity? 
stocks(sp_name$valid_name_FishBase) %>% names()


# iucn data ----

# IUCN categories
categories <- c("DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX", "LRlc", "LRnt", "LRcd")
IUCN_species <- lapply(categories, function(x) rl_sp_category(x, key = 'a81affb206f29aec7b78341c5f1f9011090be036629c463bcd4bf525245a6ac4'))

# loop through and get the species of interest from the table and assign column of category
IUCN_categories <- lapply(IUCN_species, function(x){
  y <- data.frame(x$result[which(x$result$scientific_name %in% sp_name$valid_name_FishBase), ])
  y$IUCN_category <- rep(x$category, nrow(y))
  return(y)
    })

# check heads
lapply(IUCN_categories, head)
lapply(IUCN_categories, nrow)

# bind outputs together
IUCN_all <- bind_rows(IUCN_categories[which(lapply(IUCN_categories, nrow)!=0)])

# clean categories
unique(IUCN_all$IUCN_category)
IUCN_all$IUCN_category[which(IUCN_all$IUCN_category == 'LR/lc')] <- 'LC'
IUCN_all$IUCN_category[which(IUCN_all$IUCN_category == 'nt')] <- 'NT'

# simplify
IUCN_all <- IUCN_all %>% dplyr::select(scientific_name, IUCN_category) %>% 
  rename(., c('scientific_name' = 'valid_name_FishBase'))

# find missing species
missing_IUCN <- sp_name %>% filter(!valid_name_FishBase %in% IUCN_all$scientific_name)
IUCN_full <- full_join(IUCN_all, missing_IUCN)

# assign missing species as DD
IUCN_full$IUCN_category[is.na(IUCN_full$IUCN_category)] <- 'DD'

# check perperties
nrow(IUCN_full); nrow(sp_name)

# bind together fishbase and IUCN information ----

sp_conservation_all <- plyr::join_all(list(sp_name, data.frame(vulnerability), all_lengths, resilience_filtered, IUCN_full)) %>% unique()

sp_conservation_all <- left_join(sp_conservation_all, nuts %>% dplyr::select(valid_name_FishBase, Class_FishBase)) %>% unique()

head(sp_conservation_all)

# save object
saveRDS(sp_conservation_all, file = 'processed-data/fisheries-data/species_conservation.RDS')
