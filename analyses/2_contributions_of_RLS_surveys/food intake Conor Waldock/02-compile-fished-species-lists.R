# script for fished species
#devtools::install_github("ropensci/rfishbase@3.0.1")
library(rfishbase) # set the appropriate fishbase version # use rfishbase v.3 
library(tidyverse)


#### SEA AROUND US: create fished and unfished lists from fishbase -----

# get the list of species considered in our modelling
sdm_species <- gsub('_', ' ', gsub('.RDS', '', list.files('processed-data/fish-range-synthesis/occurrence-records')))

# read in data provided by eva maire 
fished_species <- read_csv('processed-data/fisheries-data/SpList_SAU.csv')

# filter to only fishes
fished_species <- fished_species %>% filter(class == 'marine_finfish')

# read in the nutrient species list
master_list <- read_csv('processed-data/nutrient-content/Spp_NutrientPred_REEF_FUTURESJune2021.csv')

# list of species but with traits
master_traits <- read.csv('processed-data/nutrient-content/Spp_NutrientPred_REEF_FUTURES_March2021.csv')

# get the species, genus and family level matches
species_level_fished  <- master_list$valid_name_FishBase[which(master_list$valid_name_FishBase %in% fished_species$taxon_scientific_name)]
genus_level_fished    <- master_list$valid_name_FishBase[which(master_list$Genus_FishBase %in% fished_species$taxon_scientific_name)]
family_level_fished   <- master_list$valid_name_FishBase[which(master_list$Family_FishBase %in% fished_species$taxon_scientific_name)]

# get speices with maxmimum body size > 20cm
large_fish <- master_traits$valid_name_FishBase[which(master_traits$LMax > 20)]

# create the full list of potentially caught species from SAU
all_fished_SAU <- sort(unique(c(species_level_fished, genus_level_fished, family_level_fished)))

# create list that are both in SAU and large
big_fished_SAU <- sort(all_fished_SAU[all_fished_SAU %in% large_fish])

# combine the all and big fish lists
final_fished_SAU <- sort(unique(c(all_fished_SAU, big_fished_SAU)))

# create a list of only potentially caught species that are ID to species or genus level in SAU
genus_fished_SAU     <- sort(unique(c(species_level_fished, genus_level_fished)))



#### FISHBASE: Get commercial values from fishbase ----

# some debugging caused by a parsing error in readr -> https://github.com/ropensci/rfishbase/issues/221
## # download data from fishbase
## fishbase_importance_19.04 <- species(fields = c('Importance', 'Species'), version = '19.04')
## fishbase_importance_21.04 <- species(fields = c('Importance', 'Species'), version = '21.04')
## 
## table(fishbase_importance_19.04$Importance)
## table(fishbase_importance_21.04$Importance)
## 
## # old version: (for species in dataset)
## # commercial     highly commercial      minor commercial        of no interest of potential interest 
## # 545                    62                   607                   294                     9 
## # subsistence fisheries 
## # 135 
## 
## # version 3.1.9:
## # commercial     highly commercial      minor commercial 
## # 178                    21                   234 
## # of no interest subsistence fisheries 
## # 118                    45 
## 
## 
## # combine with species names
## importance_19.04 <- left_join(data.frame(valid_name_FishBase = master_traits$valid_name_FishBase), fishbase_importance_19.04, by = c('valid_name_FishBase' = 'Species'))
## importance_19.04$unique_name <- paste(importance_19.04$valid_name_FishBase, importance_19.04$Importance)
## importance_21.04 <- left_join(data.frame(valid_name_FishBase = master_traits$valid_name_FishBase), fishbase_importance_21.04, by = c('valid_name_FishBase' = 'Species'))
## importance_21.04$unique_name <- paste(importance_21.04$valid_name_FishBase, importance_21.04$Importance)
## 
## # I want a table that shows the species that differ between the two tables and those specific differences
## importance_both <- left_join(importance_19.04, importance_21.04, by = 'valid_name_FishBase')
## View(importance_both[which(importance_both$unique_name.x != importance_both$unique_name.y),])
## 
## # check where discrepancies exist
## importance_19.04 %>% filter(unique_name %in% importance_21.04$unique_name) %>% unique() %>% nrow()
## importance_21.04 %>% filter(unique_name %in% importance_19.04$unique_name) %>% unique() %>% nrow()
## 
## view(importance_19.04 %>% filter(!unique_name %in% importance_21.04$unique_name))


# read in fishbase importance for all species
fishbase_importance_21.04 <- species(fields = c('Importance', 'Species'), version = '21.04') %>% unique()

# filter to our focal species
fishbase_importance_21.04_sdm <- fishbase_importance_21.04 %>% filter(Species %in% sdm_species)

# table them 
fishbase_importance_21.04_sdm$Importance %>% table

# get how many are fished straight away
fishbase_importance_21.04_sdm %>% filter(Importance %in% c('minor commercial', 'subsistence fisheries', 'commercial', 'highly commercial'))

# reassign those with no information as of no interest
fishbase_importance_21.04_sdm[is.na(fishbase_importance_21.04_sdm$Importance), 'Importance'] <- 'unknown' 

# create a genus column to re-assign fishes that are 'unknown' based on genus level median
fishbase_importance_21.04_sdm$Genus <- sub('\\ .*', '', fishbase_importance_21.04_sdm$Species)

# get the most often status for each genus
genus_fishbase <- fishbase_importance_21.04_sdm %>% 
  group_by(Genus) %>% 
  do(fished = as.numeric(any(as.character(.$Importance) == c('minor commercial', 'subsistence fisheries', 'commercial', 'highly commercial')))) %>% 
  unnest()

# for the unknown species assign based on the genus level if-any above
fishbase_importance_21.04_sdm_unknown <- left_join(fishbase_importance_21.04_sdm %>% filter(Importance == 'unknown'), genus_fishbase)

# create two categories
fishbase_unfished_categories <- fishbase_importance_21.04_sdm %>% filter(Importance %in%  c('of no interest', 'of potential interest'))
fishbase_fished_categories   <- fishbase_importance_21.04_sdm %>% filter(!Importance %in% c('of no interest', 'of potential interest', 'unknown'))

fishbase_unfished_categories$fished = 0
fishbase_fished_categories$fished = 1

# combine all fishbase derived fishing information together
fishbase_fished_all <- rbind(fishbase_importance_21.04_sdm_unknown, fishbase_unfished_categories, fishbase_fished_categories)

# check numbers
table(fishbase_fished_all$fished)
# 0    1 
# 918 1396 

fishbase_fished_all %>% filter(fished == 1, Importance == 'unknown')

# define fished and unfished list
fishbase_unfished <- fishbase_fished_all %>% filter(fished == 0)
fishbase_fished   <- fishbase_fished_all %>% filter(fished == 1)

#### Compare fishbase and SAU fished and unfished lists ----

nrow(fishbase_fished)    # 1396 fished in fishbase
length(genus_fished_SAU) # 995 fished in seaaroundus
length(unique(c(genus_fished_SAU, fishbase_fished$Species))) # 1647 fished in both. 

# SAU species also fished in fishbase
length(genus_fished_SAU[genus_fished_SAU%in%fishbase_fished$Species])
# 774 fished in sau and also in fishbase

# SAU species missing from fishbase
length(genus_fished_SAU[!genus_fished_SAU%in%fishbase_fished$Species])
# 221 species in fished SAU but not in fishbase

# fishbase species also in SAU
length(fishbase_fished$Species[fishbase_fished$Species%in%genus_fished_SAU])
# 774 fished in fishbase and also in SAU. 

# fishbase species missing from SAU
length(fishbase_fished$Species[!fishbase_fished$Species%in%genus_fished_SAU])
# 652 speices fished in fishbase but no found on SAU


#### Combine the SAU and the fishbase categorisations ----

# get the fished union of both: using genus SAU
union_fished_genusSAU     <- union(fishbase_fished$Species, genus_fished_SAU)
union_sdm_fished_genusSAU <- intersect(sdm_species, union_fished_genusSAU)
length(union_sdm_fished_genusSAU)  # 1487 fished species with models

# get the sdm species that are not fished
union_sdm_unfished <- sdm_species[!sdm_species %in% union_fished_genusSAU]
length(union_sdm_unfished)         # 853 species are unfished

# get the sdm species that are large but not fished in fishbase or SAU
union_large_sdm_unfished <- intersect(union_sdm_unfished, large_fish)
length(union_large_sdm_unfished)   # 229 species are unfished and also large (>20cm)

final_fished_object <- list(fished_genus         = union_sdm_fished_genusSAU,
                            unfished_size20      = union_large_sdm_unfished, 
                            unfished_genus       = union_sdm_unfished, 
                            all_species          = sdm_species)

saveRDS(final_fished_object, file = '/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/reef-futures-nutrient/processed-data/fisheries-data/fished_species_lists.RDS')

### create an object to share with expert evaluators to classify unfishable of the unfished fishes ----

          
full_list_to_share <- master_traits %>% 
  dplyr::select(Order, Family, valid_name_FishBase) %>% 
  filter(valid_name_FishBase %in% union_sdm_unfished) %>% 
  unique()

# cleaning
nrow(full_list_to_share) # why doesn't this match 853? 
full_list_to_share[duplicated(full_list_to_share$valid_name_FishBase),]
full_list_to_share %>% filter(valid_name_FishBase == 'Ophthalmolepis lineolata')

# seems to be a gap in the family and order data that is duplicated
full_list_to_share <- full_list_to_share %>% filter(Family != '')

unfished_full_list_to_share <- full_list_to_share[order(full_list_to_share$Order, full_list_to_share$Family, full_list_to_share$valid_name_FishBase),]

writexl::write_xlsx(unfished_full_list_to_share, path = '/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/reef-futures-nutrient/processed-data/super-contributors/unfished_list/unfished_list_to_share.xlsx')
