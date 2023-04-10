library(tseries)
library(vars)
library(zoo)
library(fGarch)
library(urca)
library(tsDyn)

rate1 <- read.csv('eursek_d.csv')
rate1 <- cbind.data.frame(rate1$Data, rate1$Zamkniecie)
colnames(rate1) <- c('Date', 'Closing')

rate2 <- read.csv('eurisk_d.csv')
rate2 <- cbind.data.frame(rate2$Data, rate2$Zamkniecie)
colnames(rate2) <- c('Date', 'Closing')

rate3 <- read.csv('euregp_d.csv')
rate3 <- cbind.data.frame(rate3$Data, rate3$Zamkniecie)
colnames(rate3) <- c('Date', 'Closing')

dataset <- merge(rate1, rate2, by="Date", incomparables = NA, all.x = T, all.y = T)
dataset <- merge(dataset, rate3, by="Date", incomparables = NA, all.x = T, all.y = T)
colnames(dataset) <- c('Date', 'SEK', 'ISK', 'EGP')

# Interpolation
for (i in 2:4){
  dataset[,i] <- na.approx(dataset[,i], rule=2)
}

# Rates of return
rates <- data.frame()
for(i in 2:nrow(dataset)) {
  rates[i-1,1] <- dataset[i,1]
  rates[i-1,2] <- log(dataset[i,2]/dataset[i-1,2])
  rates[i-1,3] <- log(dataset[i,3]/dataset[i-1,3])
  rates[i-1,4] <- log(dataset[i,4]/dataset[i-1,4])
}
colnames(rates) <- colnames(dataset)

# Stationarity
sapply(rates[2:4], adf.test)

# VAR
VARselect(rates[2:4],5)
model1 <- VAR(rates[2:4], p=1)
summary(model1)

# Autocorelation
Box_result <- as.data.frame(matrix(NA,25,3))
colnames(Box_result) <- c('SEK', 'ISK', 'EGP')
for(i in 1:25){
  Box_result[i,1] <- Box.test(model1$varresult$SEK$residuals, i, type = "Ljung-Box")$p.value
  Box_result[i,2] <- Box.test(model1$varresult$ISK$residuals, i, type = "Ljung-Box")$p.value
  Box_result[i,3] <- Box.test(model1$varresult$EGP$residuals, i, type = "Ljung-Box")$p.value
}

# Arch effect
vars::arch.test(model1)

### ---------- ISK -------------
garch_norm_ISK <- garchFit(~garch(1,1), data = model1$varresult$ISK$residuals, trace = FALSE, cond.dist = 'norm')
garch_snorm_ISK <- garchFit(~garch(1,1), data = model1$varresult$ISK$residuals, trace = FALSE, cond.dist = 'snorm')
garch_sstd_ISK <- garchFit(~garch(1,1), data = model1$varresult$ISK$residuals, trace = FALSE, cond.dist = 'sstd')
summary(garch_norm_ISK)
summary(garch_snorm_ISK)
summary(garch_sstd_ISK)

ks.test(garch_snorm_ISK@residuals/garch_snorm_ISK@sigma.t, 'psnorm', xi= 0.87259)
ks.test(garch_sstd_ISK@residuals/garch_sstd_ISK@sigma.t, 'psstd', nu=2.3931, xi=0.9204)
# only sstd meets the assumptions

forecast_GARCH_ISK <- as.data.frame(fGarch::predict(garch_sstd_ISK, n.ahead = 5, plot=T))
forecast_GARCH_ISK

### ---------- EGP -------------
garch_norm_EGP <- garchFit(~garch(1,1), data = model1$varresult$EGP$residuals, trace = FALSE, cond.dist = 'norm')
garch_snorm_EGP <- garchFit(~garch(1,1), data = model1$varresult$EGP$residuals, trace = FALSE, cond.dist = 'snorm')
garch_sstd_EGP <- garchFit(~garch(1,1), data = model1$varresult$EGP$residuals, trace = FALSE, cond.dist = 'sstd')
summary(garch_norm_EGP)
summary(garch_snorm_EGP)
summary(garch_sstd_EGP)

ks.test(garch_sstd_EGP@residuals/garch_sstd_EGP@sigma.t, 'psstd', nu=6.56203, xi=1.04976)

forecast_GARCH_EGP <- as.data.frame(fGarch::predict(garch_sstd_EGP, n.ahead = 5, plot=T))
forecast_GARCH_EGP
# only sstd meets the assumptions

### ---------- SEK -------------
garch_norm_SEK <- garchFit(~garch(1,1), data = model1$varresult$SEK$residuals, trace = FALSE, cond.dist = 'norm')
garch_snorm_SEK <- garchFit(~garch(1,1), data = model1$varresult$SEK$residuals, trace = FALSE, cond.dist = 'snorm')
garch_sstd_SEK <- garchFit(~garch(1,1), data = model1$varresult$SEK$residuals, trace = FALSE, cond.dist = 'sstd')
summary(garch_norm_SEK)
summary(garch_snorm_SEK)
summary(garch_sstd_SEK)

ks.test(garch_snorm_SEK@residuals/garch_snorm_SEK@sigma.t, 'psnorm', xi=1.10444)

forecast_GARCH_SEK <- as.data.frame(fGarch::predict(garch_snorm_SEK, n.ahead = 5, plot=T))
forecast_GARCH_SEK
# only snorm meets the assumptions

# ------------ FORECAST -------------------
# forecast: rates of return (logarithmic) ---- based on VAR model ----
forecast_RoR <- predict(model1, n.ahead=5, level = 0.95)
forecast_ROR_ISK<-data.frame(forecast_RoR$fcst$ISK)
forecast_ROR_ISK<-forecast_ROR_ISK[,1:3]
forecast_ROR_EGP<-data.frame(forecast_RoR$fcst$EGP)
forecast_ROR_EGP<-forecast_ROR_EGP[,1:3]
forecast_ROR_SEK<-data.frame(forecast_RoR$fcst$SEK)
forecast_ROR_SEK<-forecast_ROR_SEK[,1:3]

# ---- VECM model ----
log_prices <- log(dataset[2:nrow(dataset),2:4])

adf.test(log_prices$SEK)
adf.test(log_prices$ISK)
adf.test(log_prices$EGP)
VARselect(log_prices, 5) # p=2, so VECM(1)

test_johansena <- summary(ca.jo(log_prices, type=c("trace"), K=2, ecdet="none"))
cbind(test_johansena@teststat, test_johansena@cval)

test_johansena_2 <- summary(ca.jo(log_prices, type=c("eigen"), K=2, ecdet="none"))
cbind(test_johansena_2@teststat, test_johansena_2@cval)

model_vecm <- VECM(log_prices, lag=1, r=1, estim = "ML")
summary(model_vecm)

Box_result_VECM <- as.data.frame(matrix(NA,25,3))
colnames(Box_result_VECM) <- c('SEK', 'ISK', 'EGP')
for(i in 1:25){
  Box_result_VECM[i,1] <- Box.test(model_vecm$residuals[,1], i, type = "Ljung-Box")$p.value
  Box_result_VECM[i,2] <- Box.test(model_vecm$residuals[,2], i, type = "Ljung-Box")$p.value
  Box_result_VECM[i,3] <- Box.test(model_vecm$residuals[,3], i, type = "Ljung-Box")$p.value
}

ur.df(model_vecm$residuals[,1])
ur.df(model_vecm$residuals[,2])
ur.df(model_vecm$residuals[,3])
# critical value: -1.95

# forecast: prices (logarithmic) ---- based on VECM model ----
forecast_log_prices <- predict(vec2var(ca.jo(log_prices, ecdet = 'none', type  = 'trace', K = 2, 
                                             spec = 'transitory',  dumvar = NULL)), n.ahead=5, level = 0.95)
forecast_log_prices_ISK <- data.frame(forecast_log_prices$fcst$ISK)
forecast_log_prices_ISK<-forecast_log_prices_ISK[,1:3]
forecast_log_prices_EGP <- data.frame(forecast_log_prices$fcst$EGP)
forecast_log_prices_EGP<-forecast_log_prices_EGP[,1:3]
forecast_log_prices_SEK <- data.frame(forecast_log_prices$fcst$SEK)
forecast_log_prices_SEK<-forecast_log_prices_SEK[,1:3]

# forecast: prices 
sd <- (forecast_log_prices_ISK$fcst - forecast_log_prices_ISK$upper)/1.96
forecast_prices_ISK <- cbind(round(exp(forecast_log_prices_ISK$fcst+0.5*sd^2),4))
forecast_prices_ISK  <- as.data.frame(forecast_prices_ISK)
forecast_prices_ISK[,2:3] <- exp(forecast_log_prices_ISK[,2:3])

sd <- (forecast_log_prices_EGP$fcst - forecast_log_prices_EGP$upper)/1.96
forecast_prices_EGP <- cbind(round(exp(forecast_log_prices_EGP$fcst+0.5*sd^2),4))
forecast_prices_EGP  <- as.data.frame(forecast_prices_EGP)
forecast_prices_EGP[,2:3] <- exp(forecast_log_prices_EGP[,2:3])

sd <- (forecast_log_prices_SEK$fcst - forecast_log_prices_SEK$upper)/1.96
forecast_prices_SEK <- cbind(round(exp(forecast_log_prices_SEK$fcst+0.5*sd^2),4))
forecast_prices_SEK  <- as.data.frame(forecast_prices_SEK)
forecast_prices_SEK[,2:3] <- exp(forecast_log_prices_SEK[,2:3])

# Forecast errors (example for SEK)
april_SEK <- read.csv('eursek_d_april.csv')
ME <- round(mean((april_SEK$Zamkniecie-forecast_prices_SEK[,1])),5)
MAE <- round(mean(abs(april_SEK$Zamkniecie-forecast_prices_SEK[,1])),5)
MSE <- round(mean((april_SEK$Zamkniecie-forecast_prices_SEK[,1])^2),5)
MAPE <- round(mean(abs((april_SEK$Zamkniecie-forecast_prices_SEK[,1])/april_SEK$Zamkniecie)),5)

