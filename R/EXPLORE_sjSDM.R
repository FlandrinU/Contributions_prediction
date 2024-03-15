remotes::install_gitlab( repo = "vtcfwru/occupancyTuts@1.1.0", 
                         auth_token = Sys.getenv("GITLAB_PAT"),
                         host = "code.usgs.gov", 
                         build_vignettes = FALSE,
                         upgrade = "never")


learnr::available_tutorials(package = "occupancyTuts")
learnr::run_tutorial( name = "intro", package = "occupancyTuts" )
learnr::run_tutorial( name = "study_design", package = "occupancyTuts" )


##-------------sj-SDM-------------
# see vignettes: https://cran.r-project.org/web/packages/sjSDM/vignettes/Dependencies.html
# and : https://cran.r-project.org/web/packages/sjSDM/vignettes/sjSDM_Introduction.html

vignette("Dependencies", package = "sjSDM")

sjSDM::install_sjSDM(method = "gpu") # or "cpu" if you do not have a proper gpu device

#create dataset
com = sjSDM::simulate_SDM(env = 3L, species = 5L, sites = 100L)

#run model
model = sjSDM::sjSDM(Y = com$response, env = com$env_weights, iter = 100L, se=TRUE)

#results
coef(model)
summary(model)
sjSDM::Rsquared(model)

#Option: other type of response: ex. Normal
com = sjSDM::simulate_SDM(env = 3L, species = 5L, sites = 100L,
                          link = "identical", response = "count") 
X = com$env_weights
Y = com$response

model = sjSDM::sjSDM(log(Y+0.01), env = sjSDM::linear(X, ~.), se = TRUE, 
              iter = 50L, family = gaussian("identity"))
summary(model)


#SPATIAL DATA
com = sjSDM::simulate_SDM(env = 3L, species = 5L, sites = 100L, 
                   link = "identical", response = "identical")
X = com$env_weights
Y = com$response
# add spatial residuals (create coordinates and use spatial distance matrix to draw autocorrelated residuals for each species)
XYcoords = matrix(rnorm(200), 100, 2)+2
WW = as.matrix(dist(XYcoords))
spatialResiduals = mvtnorm::rmvnorm( 5L, sigma = exp(-WW))

Ysp = Y + t(spatialResiduals)
Y = ifelse(Ysp < 0, 0, 1) # multivariate probit model

# 3 OPTIONS:
# 1) MORAN'S EIGENVECTOR MAP PREDICTOS
SPeigen = sjSDM::generateSpatialEV(XYcoords)

model = sjSDM::sjSDM(Y, env = sjSDM::linear(X, ~.), 
              spatial = sjSDM::linear(SPeigen, ~0+.), iter = 100L)
summary(model)

#2) Trend surface model - linear
colnames(XYcoords) = c("XX", "YY")
model = sjSDM::sjSDM(Y, 
              env = sjSDM::linear(X, ~.), 
              spatial = sjSDM::linear(XYcoords, ~0+XX+YY+XX:YY+I(XX^2)+I(YY^2)), 
              iter = 100L)
summary(model)

#3)Trend surface model - DNN (Deep Neural Network)
colnames(XYcoords) = c("XX", "YY")
model = sjSDM::sjSDM(Y, 
              env = sjSDM::linear(X, ~.), 
              spatial = sjSDM::DNN(XYcoords, ~0+.), 
              iter = 100L)
summary(model)


## ANALYSE RESULTS
an = anova(model)
print(an)
plot(an)


results = sjSDM::plotInternalStructure(an)
print(results$data$Species)

imp = sjSDM::importance(model)
plot(imp)

plot(model)

weights = sjSDM::getWeights(model) # get layer weights and sigma
sjSDM::setWeights(model, weights)


### Using deep neural networks
com = sjSDM::simulate_SDM(env = 3L, species = 5L, sites = 100L)
X = com$env_weights
Y = com$response

# three fully connected layers with relu as activation function
model = sjSDM::sjSDM(Y = Y, 
              env = sjSDM::DNN(data = X, 
                        formula = ~., 
                        hidden = c(10L, 10L, 10L), 
                        activation = "relu"), 
              iter = 50L, se = TRUE)
summary(model)


## PREDICT NEW DATA
pred = predict(model) # predict on fitted data
pred = predict(model, newdata = X) # predict on new data



## HELP FUNCTION sjSDM
library(sjSDM)
com = simulate_SDM(env = 3L, species = 7L, sites = 100L)

## fit model:
model = sjSDM(Y = com$response, env = com$env_weights, iter = 50L) 
# increase iter for your own data 

coef(model)
summary(model)
getCov(model)

## plot results
species=c("sp1","sp2","sp3","sp4","sp5","sp6","sp7")
group=c("mammal","bird","fish","fish","mammal","amphibian","amphibian")
group = data.frame(species=species,group=group)
plot(model,group=group)

## calculate post-hoc p-values:
p = getSe(model)
summary(p)

## predict with model:
preds = predict(model, newdata = com$env_weights)

## calculate R-squared:
R2 = Rsquared(model)
print(R2)

# Deep neural network
## we can fit also a deep neural network instead of a linear model:
model = sjSDM(Y = com$response,
              env = DNN(com$env_weights, hidden = c(10L, 10L, 10L)),
              iter = 2L) # increase iter for your own data 
summary(model)
getCov(model)
pred = predict(model, newdata = com$env_weights)

## extract weights
weights = getWeights(model)




##-------------loading data and functions-------------
#full data
load(file = here::here("data", "derived_data", "3_all_contributions_to_predict.Rdata"))

#datasets to predict
load( file = here::here("data", "derived_data", "3_datasets_for_predict_CV_80_20.Rdata"))

#covariates
load(file = here::here("data", "derived_data", "3_all_covariates_to_predict.Rdata"))

#load functions
source("R/evaluation_prediction_model.R")



##-------------test sjSDM-------------
# vignette("Dependencies", package = "sjSDM")
# sjSDM::install_sjSDM(method = "gpu") # or "cpu" if you do not have a proper gpu device

X = covariates_final[c(1:500),] |> 
  dplyr::select(-longitude, -latitude) |> 
  as.matrix()

Y = as.matrix(observations_final[c(1:500),])

XYcoords = covariates_final[c(1:500),] |> dplyr::select(XX= longitude, YY =latitude) |> 
  as.matrix()

# 1) MORAN'S EIGENVECTOR MAP PREDICTOS
SPeigen = sjSDM::generateSpatialEV(XYcoords)

model = sjSDM::sjSDM(Y, 
                     env = sjSDM::linear(X, ~.),
                     se = T, 
                     spatial = sjSDM::linear(SPeigen, ~0+.),
                     iter = 50L,
                     family = gaussian("identity"))


model = sjSDM::sjSDM(Y, 
                     env = sjSDM::DNN(X,hidden = c(10L, 10L, 10L)),
                     se = F, 
                     control = sjSDM::sjSDMControl(early_stopping_training = 5,
                                                   scheduler = 10,
                                                   lr_reduce_factor = 0.5),
                     spatial = sjSDM::DNN(SPeigen, hidden = c(5L, 5L), ~0+.),
                     iter = 100L,
                     family = gaussian("identity"))
#results
coef(model)
summary(model)
sjSDM::Rsquared(model) #1

## ANALYSE RESULTS
imp = sjSDM::importance(model)
plot(imp)

plot(model)

new_data <- covariates_final[c(500:600),] |> 
  dplyr::select(-longitude, -latitude) |> 
  as.matrix()

prediction <- predict(model, newdata = new_data, SP = XYcoords)


#### test 
com = simulate_SDM(env = 3L, species = 7L, sites = 100L)
Y = com$response
env = com$env_weights

env = covariates_final[c(1:500),] |> 
  dplyr::select(-longitude, -latitude) |> 
  as.matrix()
XYcoords <- covariates_final[c(1:500),] |>
  dplyr::select(longitude, latitude) |> 
  as.matrix()
SPeigen <- sjSDM::generateSpatialEV(XYcoords)

Y = as.matrix(observations_final[c(1:500),])

model = sjSDM(Y = Y, env = env, iter = 20L, 
              spatial = sjSDM::linear(SPeigen[, 1:20], ~0+.))

newdata = covariates_final[c(500:600),] |> 
  dplyr::select(-longitude, -latitude) |> 
  as.matrix()
SP <- sjSDM::generateSpatialEV(covariates_final[c(500:600),] |> 
                                 dplyr::select(longitude,latitude) |> 
                                 as.matrix())

preds = predict(model, newdata = newdata, SP = SP[, 1:20])

