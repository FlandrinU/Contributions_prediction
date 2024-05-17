################################################################################
##
##  
##
## 3c_fit_model_hmsc.R
##
## 07/05/2024
##
## Ulysse Flandrin
##
################################################################################

##---------------------------- cleaning memory ---------------------------------
rm(list=ls())

##-----------------------------Loading packages---------------------------------
pkgs <- c("here", "Hmsc", "coda", "ggmcmc")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))
library(ggplot2)

##------------------------------- load data ------------------------------------
load(here::here("data/derived_data/3_all_contributions_to_predict.Rdata"))
load(here::here("data/derived_data/3_all_covariates_to_predict.Rdata"))
response =  observations_final#[1:100,] ####################### reduce data
covariates = covariates_final[rownames(response),]


##----------------------------- Set-up parameters ------------------------------
set.seed(0612)

python = file.path(getwd(),"HMSC_package","hmsc-venv", "bin", "python")  # hmsc-venv for Linux and macOS

#fundamental model fitting parameters determining the MCMC sampling: 
# number of samples to obtain per MCMC chain, thinning and number of chains.
# We also define the regularity of progress printing during MCMC sampling. 
# We set the transient phase being equally long as the sampling phase.

nSamples = 500 #1000 
thin = 200 #100
nChains = 2 #2
verbose = 50 #100
transient = 10000 # 1500 #500 * thin
nb_neighbours = 10
noise_magnitude <- 0.00001 #noise in coordinates



# PATHS
save_init <- here::here("outputs/models/hmsc/init_multi/")
save_out <- here::here("outputs/models/hmsc/out_multi")
localDir <- here::here("outputs/models/hmsc/multivariate")





##------------------------ preping data and model ------------------------------

# Create Y matrix
Y <- response


# Create X matrix
# setting model structure with spatial structure ‘Nearest Neighbour Gaussian Process (NNGP)’

# (1) add noise to coordinates for them to differ sltitly
X <- covariates[rownames(Y), ] |>
  dplyr::rowwise() |> 
  dplyr::mutate(latitude = latitude + runif(1, -noise_magnitude, noise_magnitude),
                longitude = longitude + runif(1, -noise_magnitude, noise_magnitude)) |>
  dplyr::ungroup()
rownames(X) <- rownames(covariates)

# (2) create a matrix with coordinates : xycoords is a matrix with 2 columns "x-coordinate","y-coordinate" and row names with spygen_code
xycoords <- X |>
  dplyr::select(longitude, latitude)
rownames(xycoords) <- rownames(covariates)



# (3) add spatial random effect and study design (rows in Y)
rL.nngp = Hmsc::HmscRandomLevel(sData = xycoords,
                                sMethod = 'NNGP',
                                nNeighbours = nb_neighbours,
                                longlat = F) #Should be True but doesn't work after

rL.nngp = Hmsc::setPriors(rL.nngp, nfMin=1,nfMax=1)

studyDesign <- data.frame(associations = as.factor(rownames(X)),
                          # spatial = as.factor(rownames(X)),
                          # year = as.factor(X$year),
                          country = as.factor(X$country),
                          ecoregion = as.factor(X$ecoregion)
)

# rL = Hmsc::HmscRandomLevel(units = studyDesign$associations)
rL_asso = Hmsc::HmscRandomLevel(units = studyDesign$associations)
rL_year = Hmsc::HmscRandomLevel(units = studyDesign$year)
rL_ecoregion = Hmsc::HmscRandomLevel(units = studyDesign$ecoregion)
rL_country =  Hmsc::HmscRandomLevel(units = studyDesign$country)


fixed_effects <- colnames(dplyr::select(covariates,
                                        -longitude, -latitude, 
                                        -year,
                                        -country,
                                        -ecoregion,
                                        -realm,
                                        
                                        # country level covariates
                                        -control_of_corruption,
                                        -hdi,
                                        -marine_ecosystem_dependency,
                                        -ngo,
                                        -no_violence,
                                        -voice))

formula <- as.formula(paste("~ ", paste(fixed_effects, collapse = "+") ) )

# setting model structure
model_fit = Hmsc::Hmsc(Y = Y,
                       XData = X[, fixed_effects],
                       XFormula = formula,
                       studyDesign = studyDesign,
                       ranLevels = list(associations = rL_asso,
                                        # year = rL_year,
                                        ecoregion = rL_ecoregion,
                                        country = rL_country
                                        # , spatial = rL.nngp
                                        ),
                       distr = "normal") # see help for other distributions

# create object for computation on HMSC-HPC
init_obj = Hmsc::sampleMcmc(model_fit,
                            samples=nSamples,
                            thin=thin,
                            transient=transient, 
                            nChains=nChains,
                            verbose=verbose, 
                            engine="HPC") # HPC : try to use the GPU

# save it locally
file_name <- sprintf(paste0("init_multi_without_spatial&year_",
                            paste(nChains, "chains",
                                  thin, "thin",
                                  nSamples, "samples",
                                  transient, "transient",
                                  sep = "_"),
                            ".Rdata"))
init_file_path = file.path(save_init, file_name)  
saveRDS(jsonify::to_json(init_obj), file=init_file_path)



##----------------------- Construct Python command -----------------------------

# Hmsc-HPC for sequential execution of chains

# List all files in the directory
list.files(save_init)
files <- list.files(save_init, full.names = TRUE)

# list to save all arguments for python
python_cmd_args <- list()

for (file_path in files[grep(file_name, files)] ) { ###########################################################"
  # Extract response_name and j from the file name
  file_name <- basename(file_path)
  
  # Define the output file path
  post_file_path <- file.path(save_out, file_name)
  
  # Construct the Python command
  python_cmd_args[file_path] <- paste("-m hmsc.run_gibbs_sampler",
                                      "--input", shQuote(file_path),
                                      "--output", shQuote(post_file_path),
                                      "--samples", nSamples,
                                      "--transient", transient,
                                      "--thin", thin,
                                      "--verbose", verbose)
  
  # Print the Python command
  # cat(paste(shQuote(python), python_cmd_args), "\n")
}




##---------------- Running Python model (/!\ long time) -----------------------

for (file_path in files[grep(file_name,files)] ) { ##########################################################"
  system2(python, python_cmd_args[file_path])
}


##----------------------------- Study result -----------------------------------

### Importing computed posterior to R
# List all files in the directory

# save_out <- here::here("outputs/models/hmsc/out_multi")
file_name
list.files(save_out)
list_files <- list.files(save_out, full.names = TRUE)

# Import posterior probability
file_rds <- readRDS(file = list_files[5])[[1]] ######################################" Choose the right one
importFromHPC <- jsonify::from_json(file_rds)
postList <- importFromHPC[1:nChains]
cat(sprintf("fitting time %.1f min\n", importFromHPC[[nChains+1]] / 60))

model_fit_mcmc <- Hmsc::importPosteriorFromHPC(model_fit,
                                              postList, 
                                              nSamples, 
                                              thin, 
                                              transient)


model_name <- sprintf(paste0("init_multi_without_spatial_",
                             paste(nChains, "chains",
                                   thin, "thin",
                                   nSamples, "samples",
                                   transient, "transient",
                                   sep = "_"),
                             ".Rdata"))
save(model_fit_mcmc, file=file.path(localDir, model_name))



### Evaluate model
# ## MCMC convergence
# # species niches Beta; residual species associations Omega
# # Histograms of effective sample sizes and potential scale reduction factors (psrf) for Beta and Omega parameters
# 
# ns = 22 # nbr espèces
# mpost = Hmsc::convertToCodaObject(model_fit_mcmc)
# par(mfrow=c(2,2))
# ess.beta = coda::effectiveSize(mpost$Beta)
# psrf.beta = coda::gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
# hist(ess.beta)
# hist(psrf.beta)
# sppairs = matrix(sample(x = 1:ns^2, size = 100))
# tmp = mpost$Omega[[1]]
# for (chain in 1:length(tmp)){
#   tmp[[chain]] = tmp[[chain]][,sppairs]
# }
# ess.omega = coda::effectiveSize(tmp)
# psrf.omega = coda::gelman.diag(tmp, multivariate=FALSE)$psrf
# hist(ess.omega)
# hist(psrf.omega)


#Estimates for each chains
mpost <- Hmsc::convertToCodaObject(model_fit_mcmc)
summary(mpost$Beta)
postBeta <- Hmsc::getPostEstimate(model_fit_mcmc, parName = "Beta")

png("figures/models/hmsc/estimate_significance_(different_from_0).png",
    width = 25, height = 15, units = "cm", res = 300)
par(mar = c(10,10,2,2))
Hmsc::plotBeta(model_fit_mcmc, post = postBeta, param = "Support", supportLevel = 0.95)
dev.off()


#Estimate convergence -> effective size should be around nSamples*nChains, and psrf around 1
# par(mar = c(2,2,2,2))
# plot(mpost$Beta)

summary(coda::effectiveSize(mpost$Beta))
png("figures/models/hmsc/convergence_estimates.png",
    width = 20, height = 20, units = "cm", res = 200)
  par(mfrow=c(2,2), mar = c(4,4,4,4))
  hist(coda::effectiveSize(mpost$Beta), main="ess(beta)")
  abline(v= nSamples*nChains, col = "red4", lty = 2)
  hist(coda::gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
  abline(v= 1, col = "red4", lty = 2)
  hist(coda::effectiveSize(mpost$Omega[[1]]), main="ess(omega)") #also check associations between y
  hist(coda::gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega)")
dev.off()

# Test prediction performance: CV
# partition = Hmsc::createPartition(model_fit_mcmc, nfolds = 3, column = "associations")
# preds = Hmsc::computePredictedValues(hM = model_fit_mcmc,
#                                      partition = partition,
#                                      nParallel = nChains
#                                      )
# MF = Hmsc::evaluateModelFit(hM = model_fit_mcmc, predY = preds)
# MF$R2

#Prediction quality
preds <- Hmsc::computePredictedValues(model_fit_mcmc)
MF <- Hmsc::evaluateModelFit(hM=model_fit_mcmc, predY=preds)

png("figures/models/hmsc/explanatory_power.png",
    width = 10, height = 10, units = "cm", res = 300)
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))
dev.off()


head(model_fit_mcmc$X)



# Diagnostic plot




# Hmsc::computeVariancePartitioning creates bug in Rstudio because it can fill the console of warning.
# to clear it you can run 'cat("\014")' or make Ctrl + L , to continue.
# Otherwise, use the modified following function of computeVariancePartitioning
source("R/HMSC_computeVariancePartitioning.R")


VP <- computeVariancePartitioning(model_fit_mcmc)

VP$R2T$Beta ##############################################################????

# svg("figures/models/hmsc/variance_partitioning.svg")
# Hmsc::plotVariancePartitioning(model_fit_mcmc, VP = VP)
# dev.off()

# Check associations
OmegaCor <- Hmsc::computeAssociations(model_fit_mcmc)
supportLevel <- 0.95

toPlot <- ((OmegaCor[[1]]$support>supportLevel) +
             (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

png("figures/models/hmsc/residual_associations.png",
    width = 25, height = 15, units = "cm", res = 300)
# par(mar = c(10,10,2,2))
corrplot::corrplot(toPlot, method = "color",
                   col=colorRampPalette(c("blue","white","red"))(200),
                   tl.cex=.6, tl.col="black",
                   title=paste("random effect level:", model_fit_mcmc$rLNames[1]), 
                   mar=c(0,0,1,0))
#residual associations among species can be generated by correlated responses to missing covariates,
# or by ecological interactions.
dev.off()


# Variance partitioning ggplot
VP_long <- as.data.frame(VP[["vals"]]) |> 
  tibble::rownames_to_column(var = "Covariate") |> 
  tidyr::pivot_longer(
    cols = - Covariate,
    names_to = "Response",
    values_to = "Value"
  )
VP_long$Covariate <- forcats::fct_relevel(VP_long$Covariate, c(
  unique(VP_long$Covariate)[!grepl("Random", unique(VP_long$Covariate))],
  unique(VP_long$Covariate)[grepl("Random", unique(VP_long$Covariate))]
  ))


library(ggplot2)
ggplot(VP_long, aes(x = Response, y = Value, fill = Covariate)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  # scale_fill_manual(values = heat.colors(length(unique(VP_long$Covariate)), alpha = 1))+
  labs(
    title = "Variance partitioning",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 6)
  )
ggsave(filename = "figures/models/hmsc/variance_partitioning_without_spatial.jpg", width = 15, height = 10)



# Partial graph
Gradient = Hmsc::constructGradient(model_fit_mcmc,
                                   focalVariable = "n_fishing_vessels")
predY = predict(model_fit_mcmc, 
                XData = Gradient$XDataNew, 
                studyDesign = Gradient$studyDesignNew,
                ranLevels = Gradient$rLNew, 
                expected=TRUE)
plotGradient(model_fit_mcmc, Gradient, pred=predY, measure="Y", showData = TRUE,
             index=13)
plotGradient(model_fit_mcmc, Gradient, pred=predY, measure="S", showData = TRUE)


png("figures/models/hmsc/partial_plot.png",
    width = 15, height = 15, units = "cm", res = 300)
plotGradient(model_fit_mcmc, Gradient, pred=predY, measure="Y", showData = TRUE,
             index=13)
dev.off()


###### ggmcmc package #####
# install.packages("ggmcmc")
# library("ggmcmc")
# see help at : http://xavier-fim.net/post/using_ggmcmc/

S <- ggmcmc::ggs(mpost$Beta)
cov <- "n_fishing_vessels"
S <- S[grep(paste(cov, collapse = "|"), S$Parameter),]

# quick look on the distribution of the values and the shape of the posterior distribution.
ggmcmc::ggs_histogram(S)+facet_wrap(~ Parameter, ncol = 5)

# ggs_density() allows comparing the target distribution by chains and whether each
# chain has converged in a similar space.
ggmcmc::ggs_density(S)+facet_wrap(~ Parameter, ncol = 5)

# assess convergence and diagnose chain problems. the expected outcome is to produce 
# “white noise”. + a good tool to assess within-chain convergence
ggmcmc::ggs_traceplot(S)+facet_wrap(~ Parameter, ncol = 5)

# The expected output is a line that quickly approaches the overall mean, in addition
# to the fact that all chains are expected to have the same mean
ggmcmc::ggs_running(S)+facet_wrap(~ Parameter, ncol = 5)

# Ideally, the initial and final parts of the chain have to be sampling in the same 
# target distribution
ggmcmc::ggs_compare_partial(S)+facet_wrap(~ Parameter, ncol = 5)

# The autocorrelation plot expectes a bar at one in the first lag, but no autocorrelation
# beyond it. While autocorrelation is not per se a signal of lack of convergence,
# it may indicate some misbehaviour of several chains or parameters, or indicate that
# a chain needs more time to converge. 
ggmcmc::ggs_autocorrelation(S)+facet_wrap(~ Parameter, ncol = 5)

# diagnose potential problems of convergence due to highly correlated parameters
ggmcmc::ggs_crosscorrelation(S)

# compare the between-chain variation with the within-chain variation. It is expected to be close to 1.
ggmcmc::ggs_Rhat(S) + xlab("R_hat")
ggmcmc::ggs_caterpillar(S)



### ridges plot ###
S <- ggmcmc::ggs(mpost$Beta)

# Extract one or some covariate(s)
unique(S$Parameter)
fixed_effects
# cov <- "coral"
cov <- "effectiveness"
cov <- c("gravtot2", "n_fishing_vessels", "median_5year_degree_heating_week")

df <- S[grep(paste(cov, collapse = "|"), S$Parameter),] |> 
  dplyr::mutate(Parameter = as.character(Parameter),
                resp_name = sapply(strsplit(sapply(strsplit(Parameter, ","),
                                   function(x) x[2]), " "), function(x) x[2]),
                cov = sapply(strsplit(sapply(strsplit(sapply(strsplit(Parameter, ","),
                             function(x) x[1]), "\\["), function(x) x[2]), " "),
                             function(x) x[1]))
# df$cov <- factor(df$cov, levels = c("effectivenessLow", "effectivenessMedium", "effectivenessHigh"))

ggplot(df)+
  aes(y=resp_name, x=value,  fill=resp_name) +
  ggridges::geom_density_ridges(alpha=0.3, bandwidth=0.003) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6)+
  hrbrthemes::theme_ipsum( axis_title_size = 0 ) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 15)
  ) +
  xlab(cov) + ylab("Nature Contributions to People and Nature")+
  facet_wrap(~cov, ncol = 3)

ggsave(filename = paste0("figures/models/hmsc/posterior_distribution_of_estimates_",
                         paste(cov, collapse = "-"), ".jpg"),
       width = 12, height = 8)
