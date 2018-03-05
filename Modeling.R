rm(list = ls())
require(ggplot2)
require(INLA)
library(MASS)
DATA <- paste0(getwd(), "/Data/")
OUT <- paste0(getwd(), "/Fits/")
verbose = FALSE
##############
### Colors ###
##############
myjblue = rgb(20/256, 50/256, 100/256)
#########################
### Data and metadata ###
#########################
load(paste0(DATA,"Wind_data_BID.Rdata"))
#############
### Dates ###
#############
tmp = paste(paste(data$month, data$day, data$year, sep = "/"), data$hr, sep = " ")
dates = strptime(tmp, "%m/%d/%Y %H", tz = "US/Pacific") #vector of dates (in POSIXlt format)

################
###          ###
### MODELING ###
###          ###
################
#------------#
#- Response -#
#------------#
y.wind = data$Speed
#-------------------#
#- Seasonal effect -#
#-------------------#
# 1 = winter; 2 = spring; 3 = summer; 4 = fall
data$season = rep(NA, nrow(data))
for(i in 1:nrow(data)){
  if(data$month[i] == 1 || data$month[i] == 2) data$season[i] = 1
  if(data$month[i] == 3 || data$month[i] == 4) data$season[i] = 2
  if(data$month[i] == 5 || data$month[i] == 6) data$season[i] = 3
  if(data$month[i] == 7 || data$month[i] == 8) data$season[i] = 4
  if(data$month[i] == 9 || data$month[i] == 10) data$season[i] = 5
  if(data$month[i] == 11 || data$month[i] == 12) data$season[i] = 6
}
idx.seasons = data$season
#-----------------#
#- Hourly effect -#
#-----------------#
hours = dates$hour
idx.hours = inla.group(hours, method = "cut", n = 24)
#-------------------------------------------#
#- Initial precisions for temporal effects -#
#-------------------------------------------#
initialsd.hour = initialsd.season = 0.025
initialprec.hours = 1/initialsd.hour^2
initialprec.seasons = 1/initialsd.season^2
#----------------#
#- INLA formula -#
#----------------#
form = y ~ -1 + intercept + 
  f(id.season, model = 'rw2', cyclic = TRUE, hyper = list(prec = list(initial = log(initialsd.season), fixed = FALSE)), constr = TRUE) +
  f(id.hour, model = 'rw2', cyclic = TRUE, hyper = list(prec = list(initial = log(initialprec.hours), fixed = TRUE)), constr = TRUE)
#-------------#
#- INLA data -#
#-------------#
data.wei = data.frame(y = y.wind, intercept = 1, id.season = idx.seasons, id.hour = idx.hours)
###################################
### Weibull quantile regression ###
###################################
fit.qr <- inla(form, data=data.wei, family ="weibull", 
               control.family = list(variant = 1,control.link = list(quantile = 0.95, model = "quantile")), 
               control.inla = list(int.strategy = "ccd"),
               control.predictor = list(compute = TRUE, link = 1),
               control.compute=list(config = TRUE),
               verbose = verbose)
save(fit.qr, file = paste0(OUT, "fit_qr7.Rdata"))
#--------------------------------#
#- Temporal dependent threshold -#
#--------------------------------#
thr.wei = exp(fit.qr$summary.linear.predictor$mean) 
save(thr.wei, file = paste0(OUT, "thr_wei7.Rdata"))
#- Exceeds above threshold
y.exc = as.numeric(y.wind)
for (j in 1:length(y.exc)){
  thr.j = thr.wei[j]
  y.exc[j] = ifelse(y.exc[j] > thr.j, (y.exc[j] - thr.j), NA)
}
#- Define Bernoulli vector y.pexc with TRUE if an exceedance occurs 
y.pexc = logical(length(y.exc))
y.pexc[!is.na(y.exc)] = TRUE; y.pexc[is.na(y.exc)] = FALSE
# --------- #
# INLA data #
# --------- #
data.ber = data.wei
data.ber$y = as.integer(y.pexc)
data.ber$ntrials = rep(1, length(y.pexc))
#####################
### Bernoulli fit ###
#####################
fit.ber <- inla(form, data = data.ber, family = "binomial", 
                Ntrials = data.ber$ntrials, 
                control.inla = list(int.strategy = "ccd"), 
                control.predictor = list(compute = TRUE, link = 1),
                control.compute=list(config = TRUE),
                verbose = verbose)
save(fit.ber, file = paste0(OUT, "fit_ber7.Rdata"))
# ---------- #
# Filter NAs #
# ---------- #
is.na.exc = is.na(y.exc)
y.exc = y.exc[!is.na.exc]
# --------- #
# INLA data #
# --------- #
data.gpd = data.wei[!is.na.exc, ]
data.gpd$y = y.exc
###############
### GPD fit ###
###############
fit.gp <- inla(form, data = data.gpd, family = "gp",
               control.family = list(control.link = list(quantile = 0.5)),
               control.inla = list(int.strategy = "eb"),
               control.predictor = list(compute = TRUE, link = 1),
               control.compute=list(config = TRUE),
               verbose = verbose)
save(fit.gp, file = paste0(OUT, "fit_gp7.Rdata"))



