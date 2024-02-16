################################################################################
##
## 
##
## 2c_biomass_distribution.R
##
## 23/01/2024
##
## Ulysse Flandrin
##
###############################################################################"

## cleaning memory
rm(list=ls())

##------------------- loading datasets-------------------
#Species traits
load(file= here::here("outputs", "RLS_species_traits_inferred.Rdata"))

#RLS observations
load(file = here::here("data/derived_data/2_rls_actino_trop.Rdata"))

#Occurence and biomass matrices
load(file=here::here("data", "derived_data", "2_biom_matrix_sp_survey.Rdata"))
load(file=here::here("data", "derived_data", "2_relative_biom_matrix_sp_survey.Rdata"))

#Functional distances between species (cf 2b_biodiversity_indices)
load(file = here::here("data/derived_data/2_funct_gower_distance_species.Rdata"))

#Phylogenetic tree
#chronogram of ray-finned fishes from Rabosky 2018
# distribution of 100 all-taxon assembled (ATA) time-calibrated trees with 31,526 taxa
trees <- lapply(list.files(
  here::here("data/raw_data/phylogeny/dryad_Rabosky_2018/trees/full_trees/"), full.names = T),
  FUN = function(i){ape::read.tree(i)})

#Occurence and biomass matrices for phylogenetic tree (cf 2b_biodiversity_indices)
save(sp_occ_phylogeny, file = here::here("outputs", "2b_sp_occ_phylogeny_matrix.Rdata"))
save(sp_pbiom_phylogeny, file = here::here("outputs", "2b_sp_pbiom_phylogeny_matrix.Rdata"))

#Biodiversity indices
load(file = here::here("outputs", "2b_biodiv_indices_surveys.Rdata" ))

##------------------- 1) biomass per trophic guild -------------------####

summary(inferred_species_traits$trophic_guild)


## Species with low TL = herbivores_microvores_detritivores diet
species_lowTL <- inferred_species_traits |>
  dplyr::filter(trophic_guild == "Herbivores Microvores Detritivores") |>
  rownames() 
length(species_lowTL) # 232 species

biom_lowTL <- rowSums(surveys_sp_biom[,species_lowTL]) 
summary(biom_lowTL) # from 0 to 466219 , median=11986


## Species with medium TL = all type of invertivorous diet (including planktivores and corallivores)
species_mediumTL <- inferred_species_traits |>
  dplyr::filter(trophic_guild %in% c("corallivore", "sessile invertivores", 
                            "microinvertivore", "macroinvertivore", "crustacivore", 
                            "planktivore") &
                  class != "Elasmobranchii") |>
  rownames() 
length(species_mediumTL) # 1148 species

biom_mediumTL <- rowSums(surveys_sp_biom[,species_mediumTL])
summary(biom_mediumTL) # from 7.8 to 458431.2, median=14624.4



## Species with high TL = piscivores
species_highTL <- inferred_species_traits |>
  dplyr::filter(trophic_guild=="piscivore" & class != "Elasmobranchii") |>
  rownames() 
length(species_highTL) # 190 species

biom_highTL <- rowSums(surveys_sp_biom[,species_highTL])
summary(biom_highTL) # from 0 to 467505.5, median=1145.9


## Total biomass
species_NA <- inferred_species_traits |>
  dplyr::filter(is.na(trophic_guild) & class != "Elasmobranchii") |>
  rownames()

all_actino <- c(species_lowTL, species_mediumTL, species_highTL, species_NA) 
length(all_actino) #1609 -> All actino species have been taken.

total_biomass <- rowSums(surveys_sp_biom[,all_actino])
summary(biom_tot) # from 10.7 to 498916.0, median=35424.5



## Merging biomass
surveys_biomTL <- data.frame(biom_lowTL, biom_mediumTL, biom_highTL, total_biomass) |>
  tibble::rownames_to_column("survey_id")
summary(surveys_biomTL) #OK



## Check the importance of NA 
NA_prop <- rls_actino_trop|> 
  dplyr::filter(species_name %in% species_NA) |> 
  dplyr::group_by(survey_id, abundance_tot_survey, biomass_tot_survey) |> #keep totals at the survey scale
  dplyr::summarise(total = sum(total),
                   biomass = sum(biomass)) |> 
  dplyr::mutate(prop_abund = total / abundance_tot_survey,
                prop_biom = biomass / biomass_tot_survey)|> #proportion of species with NA
  dplyr::filter(prop_abund > 0.2 | prop_biom  > 0.2) 
#SELECT WHEN "NAs" SPECIES REPRESENT MORE THAN 20% OF THE SURVEY, IN ABUNDANCE OR BIOMASS:

na_survey <- unique(NA_prop$survey_id)



## REMOVE NON REPRESENTATIVE ESTIMATIONS:
surveys_biomTL <- surveys_biomTL |> 
  dplyr::mutate(across(-survey_id, ~ifelse(survey_id %in% na_survey, NA, .)))




##------------------- 2) functional entropy in surveys -------------------####

surveys_sp_pbiom <- as.matrix(surveys_sp_pbiom)

# functional entropy / richness on species relative biomass with q=1
funct_entropy <- mFD::alpha.fd.hill(asb_sp_w = surveys_sp_pbiom, 
                                    sp_dist = sp_gower,
                                    q=1, 
                                    tau="mean",
                                    details_returned =FALSE)[,1]

summary(funct_entropy)
hist(funct_entropy)




##------------------- 3) Phylogenetic entropy -------------------####

sp_occ_phylogeny <- as.matrix(sp_occ_phylogeny)

##Extract trees
phylo_100 <- list()
for (i in 1:100) {
  phylo_100[[i]] <- picante::match.phylo.comm(trees[[i]], sp_occ_phylogeny)
}

## phylogenetic entropy: Marcon 2015, from Allen 2009  
# -> very long time to run: run only on 10 trees due to negligible variations
library('entropart')

phylo_entropy_raw <- lapply( phylo_100[1:10], function(x) {
  
  #due to numerical rounding, Rabosky's tree is not ultrametric anymore -> coerce
  # it to be ultrametric:
  # phytools::plotTree(x[[1]][["phy"]])
  tree <- phytools::force.ultrametric(x[["phy"]])
    
  cat("Start computation for one tree... \n")
  
  list <- parallel::mclapply(rownames(sp_pbiom_phylogeny), 
                             mc.cores = parallel::detectCores()-5,
                             function(survey){
                               
    community_biom <- entropart::as.ProbaVector(sp_pbiom_phylogeny[survey,])
    phylo_entropy <- entropart::PhyloEntropy(Ps=community_biom,
                                             q = 1,
                                             Tree = tree, 
                                             Normalize = T) #phylogenetic entropy of order 1 of relative biomass vector.
    phylo_entropy[["Total"]]
  })

  do.call(rbind, list)
  
})

save(phylo_entropy_raw, file = here::here("outputs", "2c_phylogenetic_entropy_surveys_10trees.Rdata"))
# load(file = here::here("outputs", "2c_phylogenetic_entropy_surveys_10trees.Rdata"))


phylo_entropy_surveys_10 <- do.call(cbind, phylo_entropy_raw)
# test variability
sd <- apply(phylo_entropy_surveys_10, 1, sd)
summary(sd) # Median 0.0023459      Mean 0.0053854  3rd Qu. 0.0072062      Max. 0.0538961
#

phylo_entropy_summary <- t(apply(phylo_entropy_surveys_10, 1, summary))
colnames(phylo_entropy_summary) <- paste0("phylo_entropy_", colnames(phylo_entropy_summary))
rownames(phylo_entropy_summary) <- rownames(sp_pbiom_phylogeny)
  



## REMOVE NON REPRESENTATIVE ESTIMATIONS: >20% of biomass or abundance of the 
##  species are not in the phylogenetic tree

#NA are in the same surveys as for Evolutionary distinctiveness (cf 2b_biodiv_indices)
na_survey <- biodiv_indices_surveys$survey_id[
  is.na(biodiv_indices_surveys$ED_mean)]


phylo_entropy <- as.data.frame(phylo_entropy_summary) |> 
  tibble::rownames_to_column(var="survey_id") |> 
  dplyr::select(survey_id, phylo_entropy_Mean) |> 
  dplyr::mutate(phylo_entropy_Mean = as.numeric(ifelse(survey_id %in% na_survey, 
                                                       NA, phylo_entropy_Mean)))

hist(phylo_entropy$phylo_entropy_Mean)






##-------------Merge and Save biomass distribution indices-------------
biomass_indices_surveys <- as.data.frame(funct_entropy) |> 
  tibble::rownames_to_column("survey_id") |> 
  dplyr::full_join(surveys_biomTL) |>
  dplyr::full_join(phylo_entropy) 


## Check indices measures
summary(biomass_indices_surveys)

library(funbiogeo)
fb_plot_species_traits_completeness(dplyr::rename(biomass_indices_surveys, species = survey_id))

pca <- FactoMineR::PCA(biomass_indices_surveys[,-1], scale = T, graph=F, ncp=30)
factoextra::fviz_screeplot(pca, ncp = 20)
factoextra::fviz_pca_var(pca, col.var = "contrib", repel = TRUE)
factoextra::fviz_pca_var(pca, col.var = "contrib", axes=c(3,4), repel = TRUE)

# Check distributions
ggplot(data=tidyr::pivot_longer(biomass_indices_surveys,
                                cols = -survey_id,
                                names_to = "index", values_to = "values"), 
       aes(x=values, group=index, fill=index)) +
  geom_histogram(aes(y = ..density..), bins = 20, color = "grey40", fill ="white") +
  geom_density(aes(fill = index), alpha = 0.2) +
  hrbrthemes::theme_ipsum() +
  facet_wrap(~index, scales = "free") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )
#



save(biomass_indices_surveys, file = here::here("outputs", "2b_biomass_indices_surveys.Rdata" ))
# load(file = here::here("outputs", "2b_biomass_indices_surveys.Rdata" ))
