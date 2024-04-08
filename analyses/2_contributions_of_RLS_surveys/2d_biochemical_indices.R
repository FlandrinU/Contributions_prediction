################################################################################
##
## 
##
## 2d_biochemical_indices.R
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

#Survey metadata
load(file = here::here("data", "raw_data", "environmental_covariates",
                       "all_covariates_benthos_inferred_tropical_surveys.Rdata"))

# Metabolic parameters
parameters <- utils::read.csv(here::here("data", "raw_data", "recycling data",
                                         "species_parameters_metabolic.csv") )
metpar <- utils::read.csv(here::here("data", "raw_data", "recycling data",
                              "metpar_fam_smr.csv") )  |>
  dplyr::rename(family = Family)

# combined data of reef services and renato's extraction of fishbase
kmax <- utils::read.csv( here::here("data", "raw_data", "recycling data",
                                    "kmax_combined.csv") ) |>
  dplyr::mutate(species = gsub(" ", "_", Species), sst = sstmean, linf_m = sizemax) 

# phylogenetic tree
load( here::here("data", "raw_data", "recycling data", "fishtree_glob.RData") )


## Load functions ##
source(here::here("R","evaluation_prediction_model.R"))



##------------------- 1) get all parameters -------------------####
spcombo <- all_covariates_benthos_inferred  |> 
  dplyr::select(survey_id, sst = mean_5year_analysed_sst) |>
  dplyr::right_join(rls_actino_trop, by="survey_id") |>
  dplyr::select(sst, species=rls_species_name, size_class) |>
  dplyr::mutate(sst = round(sst),
                species = gsub(" ", "_", species)) |>
  unique()


nrow(spcombo) # 27 162
length(unique(spcombo$species)) #1609 OK

# Get all parameters that are independent from sst
length(which(parameters$species %in% 
               gsub(" ", "_", rownames(inferred_species_traits)))) #1034 species out of 1609
sp_par <- dplyr::left_join(spcombo, parameters)



# Add metabolic parameters
B0_mean <- mean(metpar$B0)
B0_sd <- sd(metpar$B0)
alpha_mean <- mean(metpar$alpha)
alpha_sd <- sd(metpar$alpha)
theta_mean <- mean(metpar$theta)

sp_par <- sp_par |> 
  dplyr::left_join(metpar)

sp_par[is.na(sp_par$alpha), "alpha"] <- alpha_mean
sp_par[is.na(sp_par$alpha_sd), "alpha_sd"] <- alpha_sd
sp_par[is.na(sp_par$B0), "B0"] <- B0_mean
sp_par[is.na(sp_par$B0_sd), "B0_sd"] <- B0_sd
sp_par[is.na(sp_par$theta), "theta"] <- theta_mean


# temperature adjustment of metabolic constant 
sp_par$B0_adj <- 
  sp_par$B0 * exp(0.59 / 8.62e-5 * (1 / (28 + 273.15) - 1 / (sp_par$sst + 273.15)))

sp_par <- dplyr::rename(sp_par, alpha_m = alpha,
                        f0_m = B0_adj, f0_sd = B0_sd, theta_m = theta) |>
  dplyr::select(-B0)



# k parameter

# get fishtree
sub_fishtree <- fishtree::fishtree_complete_phylogeny(kmax$species, mc.cores = 15)
# just use one tree
set.seed(1)
tree <- sub_fishtree[[sample(1:100, 1)]]
# correlation matrix
A <- ape::vcv(tree, cor = TRUE)

# only use names that exist in tree
kmax <- dplyr::filter(kmax, species %in% colnames(A))


fit_kmax <- brms::brm(
  log(kmax) ~ log(linf_m) + sst + (1|gr(species, cov = A)), 
  data = kmax, data2 = list(A = A),
  family = gaussian(), cores = 4
)

summary(fit_kmax)
brms::bayes_R2(fit_kmax)
# phylogenetic signal
hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
brms::hypothesis(fit_kmax, hyp, class = NULL)


# extrapolate
draws <- tidybayes::spread_draws(fit_kmax, r_species[species,Intercept]) |>
  tidybayes::mean_qi() 

rsp <- draws$r_species
names(rsp) <- draws$species
hist(draws$r_species)

# phylogenetic full tree
# just use one tree
set.seed(2)
treeglob <- fishtree[[sample(1:100, 1)]]

rphy_pred <- picante::phyEstimate(treeglob, rsp) |>
  tibble::rownames_to_column("species") |>
  dplyr::mutate(r_phylo = estimate) |>
  dplyr::select(species, r_phylo)

rphy <- data.frame(
  species = names(rsp),
  r_phylo = rsp
)

rphy <- dplyr::bind_rows(rphy, rphy_pred)

tidybayes::get_variables(fit_kmax)

vars <- fit_kmax |>
  tidybayes::spread_draws(b_Intercept, b_loglinf_m, b_sst)

kpred <- lapply(1:nrow(sp_par), function(i){
  data <-
    dplyr::left_join(sp_par[i,],  rphy)
  
  linf <- unique(log(data$linf_m)) 
  sst <-  unique(data$sst)
  
  kmax <- exp(vars$b_Intercept +
                (linf * vars$b_loglinf_m) +
                (sst * vars$b_sst) +
                (data$r_phylo))
  
  data.frame(
    species = data$species,
    sst = data$sst,
    k_m = mean(kmax),
    k_sd = sd(kmax)
  )
}) |> plyr::ldply()

sp_par <- dplyr::left_join(sp_par, unique(kpred))

# Add values for AE
ae <- utils::read.csv(here::here("data", "raw_data", "recycling data", "ae_dietcat.csv"),
                      sep = ",") |>
  dplyr::mutate(diet_cat = as.character(diet_cat)) |>
  dplyr::mutate(
    an_sd = dplyr::case_when(an_sd > 0.2 ~ 0.2, TRUE ~ an_sd),
    ap_sd = dplyr::case_when(ap_sd > 0.2 ~ 0.2, TRUE ~ ap_sd),
    ac_sd = dplyr::case_when(ac_sd > 0.2 ~ 0.2, TRUE ~ ac_sd)
  ) |>
  dplyr::mutate(diet_cat = as.character(diet_cat))

# mean for missing diet category
ae_5 <- ae |>
  dplyr::select(- diet_cat) |>
  dplyr::summarize_all(mean) |>
  dplyr::mutate(diet_cat = "5")

ae <- dplyr::bind_rows(ae, ae_5) 

sp_par <- dplyr::left_join(
  dplyr::mutate(sp_par, diet_cat = as.character(diet_cat)), ae) |>
  dplyr::mutate(v_m = sst) |>
  unique()

nrow(sp_par)
readr::write_csv(sp_par, file=here::here("data", "derived_data","2d_parameters_sp_sst.csv")  )





##------------------- 2) Run fishflux  -------------------####

sp_par <- utils::read.csv( here::here("data", "derived_data","2d_parameters_sp_sst.csv")  )

data <- sp_par |> tidyr::drop_na()
length(unique(data$species)) #1034 species
#data <- data[1:100,]
cnpflux <- pbmcapply::pbmclapply(1:nrow(data), function(x){
  cat(x, "\n")
  
  dt <- data[x,] 
  par <- dt |> dplyr::select(-species, - sst, - size_class, -family, -diet_cat) |> as.list()
  mod <- fishflux::cnp_model_mcmc(TL = dt$size_class,
                                  param = par, iter = 1000)
  
  
  extr <- fishflux::extract(mod, par = c("F0c", "F0n", "F0p", "Gc", "Gn", "Gp", "Sc", "Sn", "Sp", 
                                         "Ic", "In", "Ip", "Wc", "Wn", "Wp", "Fc", "Fn", "Fp"))
  extr <- cbind(dt[,1:5], extr) 
  lim <- fishflux::limitation(mod, plot = FALSE)
  extr$limitation <- dplyr::first(lim[lim$prop_lim == max(lim$prop_lim), "nutrient"])
  
  return(extr)
}, mc.cores = parallel::detectCores()-5) |> plyr::ldply()


cnpflux <- dplyr::select(cnpflux, - Qc_m, - TL)

# saving
readr::write_csv(cnpflux, here::here("outputs", "2d_cnpflux_sp_size_sst.csv"))





##------------------- 3) Assess flows at the survey scale  -------------------####
# details about species: length-weight, age, nutrient contents, diet cat
species_par <- utils::read.csv(here::here("data", "derived_data","2d_parameters_sp_sst.csv")) |>
  dplyr::select(family, species, lwa_m, lwb_m, k_m, linf_m, sst, Qc_m, Qn_m, Qp_m, diet_cat) 



# fluxes from fishes
data_fishflux <- utils::read.csv(here::here("outputs", "2d_cnpflux_sp_size_sst.csv")) |> unique()
head(data_fishflux)
names(data_fishflux)

# merging datasets and computing fluxes for each species*size_class given abundance and sst
data_surveys_fluxes <- rls_actino_trop |> 
  dplyr::mutate(species = gsub(" ", "_", rls_species_name)) |> 
  dplyr::select(-rls_species_name) |> 
  dplyr::left_join(dplyr::select(all_covariates_benthos_inferred, survey_id, 
                                 sst = mean_5year_analysed_sst )) |>
  dplyr::mutate(sst = round(sst)) |>
  dplyr::left_join(data_fishflux) |>
  dplyr::left_join(unique(species_par))  |>
  dplyr::mutate(fishflux = as.factor(dplyr::case_when(is.na(Gc_median) ~ FALSE, 
                                                      Gc_median >= 0 ~ TRUE) )
  ) |>
  dplyr::mutate(storage_C = Gc_median * total,
                storage_N = Gn_median * total,
                storage_P = Gp_median * total,
                excretion_N = Fn_median * total,
                excretion_P = Fp_median * total,
                egestion_N = Wn_median * total,
                egestion_P = Wp_median * total,
                egestion_C = Wc_median * total,
                biomass_estimated = lwa_m * (size_class^lwb_m) * total, #biomass of individuals with fishflux data only
                abundance_estimated = ifelse(fishflux == "TRUE", total,0) )

head(data_surveys_fluxes)
# summary(data_surveys_fluxes)
nrow(data_surveys_fluxes)


# computing for each flux, matrix survey*species (as for biomass) --

# variables of interest with new names
fluxes_var<-c("storage_C", "storage_N", "storage_P",
              "excretion_N", "excretion_P", 
              "egestion_C", "egestion_N", "egestion_P",
              "biomass_estimated", "abundance_estimated"
)

# list to store matrices
surveys_species_fluxes<-list()

# loop on variables
for (k in fluxes_var ) {
  
  # merging size_classes per species in each survey
  mat_k<-data_surveys_fluxes |> 
    dplyr::select(survey_id, species, !!k ) |>
    dplyr::rename( flux_k = !!k ) |>
    dplyr::group_by(survey_id, species) |>
    dplyr::summarize( sum_flux_k=sum(flux_k, na.rm=T) ) |>
    tidyr::pivot_wider(names_from = species, values_from = sum_flux_k, values_fill = 0) |>
    tibble::column_to_rownames(var="survey_id") |>
    as.matrix()  
  
  # storing
  surveys_species_fluxes[[k]]<-mat_k
  
}# end of k

# lapply(surveys_species_fluxes, dim)

# total of fluxes per survey --
surveys_fluxes <- unlist( sapply(surveys_species_fluxes, function(x) rowSums(x, na.rm = TRUE)) ) |> 
  as.data.frame() |>
  tibble::rownames_to_column("survey_id")
head(surveys_fluxes)
dim(surveys_fluxes)

## computing stoichiometric ratios of excretion and egestion
# for each species*size then median at survey level
mC <- 12
mN <- 14
mP <- 31
median_ratio <- data_surveys_fluxes |>
  tidyr::uncount(total) |>
  dplyr::group_by(survey_id) |>
  dplyr::summarise(egestion_CN = median( (egestion_C/mC) / (egestion_N/mN) , na.rm=T),
                   egestion_CP = median( (egestion_C/mC)/ (egestion_P/mP) , na.rm=T),
                   egestion_NP = median( (egestion_N/mN) / (egestion_P/mP) , na.rm=T),
                   excretion_NP = median( (excretion_N/mN) / (excretion_P/mP) , na.rm=T)
  )

# merging and ratio of quality of excretion and of egestion
surveys_fluxes <- surveys_fluxes |>
  dplyr::left_join(median_ratio) |>
  dplyr::mutate(excrNP_egesNP= excretion_NP / egestion_NP)


## computing recycling as sum of excretion plus egestion
surveys_fluxes <- surveys_fluxes |>
  dplyr::mutate(recycling_C = egestion_C) |>
  dplyr::mutate(recycling_N = excretion_N + egestion_N) |>
  dplyr::mutate(recycling_P = excretion_P + egestion_P)

surveys_species_fluxes[["recycling_N"]] <- surveys_species_fluxes[["excretion_N"]] + surveys_species_fluxes[["egestion_N"]]
surveys_species_fluxes[["recycling_P"]] <- surveys_species_fluxes[["excretion_P"]] + surveys_species_fluxes[["egestion_P"]]


## computing contribution of excretion to recycling of N and P
surveys_fluxes <- surveys_fluxes |>
  dplyr::mutate(pexcr_recycling_N = excretion_N / recycling_N) |>
  dplyr::mutate(pexcr_recycling_P = excretion_P / recycling_P)


## computing ratio between recycling and storing
surveys_fluxes <- surveys_fluxes |>
  dplyr::mutate( recyc_stor_C = recycling_C / storage_C ) |>
  dplyr::mutate( recyc_stor_N = recycling_N / storage_N ) |>
  dplyr::mutate( recyc_stor_P = recycling_P / storage_P )


# summary
summary(surveys_fluxes)
nrow(surveys_fluxes)

### Check the proportion of abundance and biomass evaluated by fishflux ###
tot_biom <- rls_actino_trop |> 
  dplyr::group_by(survey_id, abundance_tot_survey) |> 
  dplyr::summarise(biomass_tot_survey = sum(raw_biomass))

surveys_fluxes_final <- surveys_fluxes |> 
  dplyr::left_join(tot_biom) |> 
  dplyr::mutate(prop_biom_fishflux = biomass_estimated/biomass_tot_survey,
                prop_abund_fishflux = abundance_estimated/abundance_tot_survey)

colnames(surveys_fluxes_final)
var <- c("excretion_N", "excretion_P", 
         "recycling_N" , "recycling_P",
         "prop_biom_fishflux", "prop_abund_fishflux")

# Check distributions
distribution_plot(surveys_fluxes_final, longer = T,
                  cols_plot = var)
#

surveys_fluxes_final <- surveys_fluxes_final |> 
  dplyr::select(survey_id, all_of(var))

na_rows <- which(surveys_fluxes_final$prop_biom_fishflux<0.8 |
    surveys_fluxes_final$prop_abund_fishflux<0.8)

surveys_fluxes_final[na_rows, c("excretion_N", "excretion_P",
                                  "recycling_N", "recycling_P")] <- NA



### saving data ###
save(surveys_species_fluxes, file=here::here("outputs", "2d_surveys_species_fluxes.Rdata") )
save(surveys_fluxes_final, file=here::here( "outputs", "2d_surveys_fluxes.Rdata") )
