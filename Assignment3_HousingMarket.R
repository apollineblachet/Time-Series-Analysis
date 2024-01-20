

library(tseries)
library(forecast)
library(car)
library(aTSA)

data <- read.csv("A3Data.csv")

# Initialize an empty list to store the results
ts_list <- list()

#Convert the data to quarterly Time series
ts <- ts(data, frequency = 4)

# Assign the values to separate variables
time <- seq(1992.5,2024.25,0.25)
denmark <- na.omit(ts[,2])
capital <- na.omit(ts[,3])
sealand <- na.omit(ts[,4])
midjutland <- na.omit(ts[,5])
rural <- na.omit(ts[,6])
interest <- na.omit(ts[,7])
inflation <- na.omit(ts[,8])




# plots
par(mfrow = c(1,1))

# plot the quarterly average sale prices 
plot(time[1:122],denmark, type='l',
     main="Quarterly average sales prices in Denmark",
     xlab='Time [year]' ,ylab="House Price [DKK/m2]")
grid()
abline(h=mean(denmark), lty=2, col="#8B0000")
legend("topleft", legend="Average sales prices", col=1, lty=1, bty="n", cex=0.8)

# plot the interest rate
plot(time[1:124],interest[1:124], type='l',
     main="Interest Rate",
     xlab='Time [year]', ylab="Interest Rate [%]")
grid()
abline(h=mean(interest), lty=2, col="#8B0000")
legend("topright", legend="Interest Rate", col=1, lty=1, bty="n", cex=0.8)

# plot the inflation rate
plot(time[1:124],inflation[1:124], type='l',
     main="Inflation Rate",
     xlab='Time [year]', ylab="Inflation Rate [%]")
grid()
abline(h=mean(inflation), lty=2, col="#8B0000")
legend("topleft", legend="Inflation Rate", col=1, lty=1, bty="n", cex=0.8)




# ACF and PACF of the house prices
par(mfrow = c(1,1))
acf(denmark, main = "Autocorrelation funtion of the raw time series data") 
pacf(denmark, main = "Partial autocorrelation fucntion of the raw time series data")

# Differencing
dif_den <- diff(denmark)
dif_cap <- diff(capital)
dif_seal <- diff(sealand)
dif_midj <- diff(midjutland)
dif_rul <- diff(rural)

# ACF and PACF for the transformed data
par(mfrow = c(1,1))
acf(dif_den, main = "Autocorrelation funtion of the differenced time series data")
pacf(dif_den, main = "Partial autocorrelation fucntion of the differenced time series data")

# cross-correlation function between the four time series of the house prices.
par(mfrow = c(2,3)) # raw series
ccf(capital, sealand, main = "CCF: 'Capital' & 'Sealand' (raw)")
ccf(capital, midjutland, main = "CCF: 'Capital' & 'MidJutland' (raw)")
ccf(capital, rural, main= "CCF: 'Capital' & 'Rural' (raw)")
ccf(sealand, midjutland, main = "CCF: 'Sealand' & 'MidJutland' (raw)")
ccf(sealand, rural, main= "CCF: 'Sealand' & 'Rural' (raw)")
ccf(midjutland, rural, main = "CCF: 'MidJutland' & 'Rural' (raw)")

# Pre-whitening
# Define the ARIMA model for Denmark data
arima_model = arima(denmark, order = c(0, 1, 1), seasonal = list(order = c(2, 0, 0), period = 4))

# Prewhiten the Denmark data
pw_den <- residuals(arima_model)

# Prewhiten the capital data
pw_cap <- residuals(update(arima_model, x = capital))

# Prewhiten the sealand data
pw_seal <- residuals(update(arima_model, x = sealand))

# Prewhiten the midjutland data
pw_midj <- residuals(update(arima_model, x = midjutland))

# Prewhiten the rural data
pw_rul <- residuals(update(arima_model, x = rural))

# Estimate the SCCF for the filtered series
acf(na.pass(cbind(pw_den, pw_cap, pw_seal, pw_midj, pw_rul)), type = "correlation")

# Different plot style
par(mfrow=c(2,2))
ccf(pw_den, pw_cap, lag.max = 20, type = "correlation", main = "CCF: 'Denmark' & 'Capital'")
ccf(pw_den, pw_seal, lag.max = 20, type = "correlation", main = "CCF: 'Denmark' & 'Sealand'")
ccf(pw_den, pw_midj, lag.max = 20, type = "correlation", main = "CCF: 'Denmark' & 'MidJutland'")
ccf(pw_den, pw_rul, lag.max = 20, type = "correlation", main = "CCF: 'Denmark' & 'Rural'")

par(mfrow=c(2,3))
ccf(pw_cap, pw_seal, lag.max = 20, type = "correlation", main = "CCF: 'Capital' & 'Sealand'")
ccf(pw_cap, pw_midj, lag.max = 20, type = "correlation", main = "CCF: 'Capital' & 'MidJutland'")
ccf(pw_cap, pw_rul, lag.max = 20, type = "correlation", main = "CCF: 'Capital' & 'Rural'")
ccf(pw_seal, pw_midj, lag.max = 20, type = "correlation", main = "CCF: 'Sealand' & 'MidJutland'")
ccf(pw_seal, pw_rul, lag.max = 20, type = "correlation", main = "CCF: 'Sealand' & 'Rural'")
ccf(pw_midj, pw_rul, lag.max = 20, type = "correlation", main = "CCF: 'MidJutland' & 'Rural'")




# white noise
tsdiag(m1)
par(mfrow=c(1,1))
pacf(residuals(m1), main='PACF of Residuals')
cpgram(m1$residuals, main='Cumulative Periodogram of Residuals')
grid()
qqPlot(residuals(m1), main='Q-Q plot of Residuals', ylab='Residuals') 





forecast_data <- predict(m1, 6) 
print(forecast_data)
plot(time[1:122],denmark, xlim=c(1992.5,2025), type="l",
     xlab='Time [year]', ylab='House Price [DKK/m2]',
     main='Predictions without external inputs')
matlines(time[123:128], forecast_data$pred + 1.96*cbind(0, -forecast_data$se, forecast_data$se), col=2, lty=c(1,2,2))
grid()
legend("topleft", legend = c("Training data", "Predictions with 95% p.i."), lty=1, col=1:2)




interest.pred <- c(interest[123:124], replicate(4,interest[124]))
inflation.pred <- c(inflation[123:124], replicate(4,inflation[124]))
externalreg.pred <- cbind(cumsum(interest.pred) + externalreg[122,1], cumsum(inflation.pred) + externalreg[122,2])

m2.pred <- predict(m2, n.ahead = 6, newxreg = externalreg.pred)
print(m2.pred)

plot(time[1:122],denmark, xlim=c(1992.5,2025), type="l",
     xlab='Time [year]', ylab='House Price [DKK/m2]',
     main='Predictions with external inputs')
matlines(time[123:128], m2.pred$pred + 1.96*cbind(0, -m2.pred$se, m2.pred$se), col=2, lty=c(1,2,2))
grid()
legend("topleft", legend = c("Training data", "Predictions with 95% p.i."), lty=1, col=1:2)




library(marima)
library(car)

data <- read.csv("A3Data.csv")

# Initialize an empty list to store the results
ts_list <- list()

# Convert the data to quarterly time series
ts <- ts(data, frequency = 4)

# Assign the values to separate variables
time <- seq(1992.5,2024.25,0.25)
denmark <- na.omit(ts[,2])
capital <- na.omit(ts[,3])
sealand <- na.omit(ts[,4])
midjutland <- na.omit(ts[,5])
rural <- na.omit(ts[,6])
interest <- na.omit(ts[,7])
inflation <- na.omit(ts[,8])





par(mfrow=c(1,1))
plot(time[1:122],capital,type='l',
     ylab="House Price [DKK/m2]",xlab="Time [year]",
     main='Quarterly average sales prices for each of the four regions')
lines(time[1:122],sealand,col='steelblue')
lines(time[1:122],midjutland,col='tomato')
lines(time[1:122],rural,col='seagreen')
grid()
legend("topleft",c("Capital region","Sealand","Middle Jutland","Rural areas"),
       col=c("black","steelblue","tomato","seagreen"),lty=1,bty='n')

#Differences
dif_cap <- diff(capital)
dif_seal <- diff(sealand)
dif_midj <- diff(midjutland)
dif_rul <- diff(rural)

plot(time[2:122],dif_cap,type='l',
     ylab="House Price [DKK/m2]",xlab="Time [year]",
     main='Quarterly average sales prices for each of the four regions (differenced)')
lines(time[2:122],dif_seal,col='steelblue')
lines(time[2:122],dif_midj,col='tomato')
lines(time[2:122],dif_rul,col='seagreen')
grid()
legend("topleft",c("Capital region","Sealand","Middle Jutland","Rural areas"),
       col=c("black","steelblue","tomato","seagreen"),lty=1,bty='n')




# Raw data
four_series <- data.frame(ts[1:122,3:6])
par(mfrow=c(2,2))
for(i in 1:4) {
  series <- four_series[,i]
  acf(series, lag.max=25, main=colnames(four_series)[i])
}
par(mfrow=c(2,2))
for(i in 1:4) {
  series <- four_series[,i]
  pacf(series, lag.max=25, main=colnames(four_series)[i])
}
# Differenced data
dif_four_series <- cbind(dif_cap,dif_seal,dif_midj,dif_rul)
par(mfrow=c(2,2))
for(i in 1:4) {
  dif_series <- dif_four_series[,i]
  acf(dif_series, lag.max=25, main=colnames(four_series)[i])
}
par(mfrow=c(2,2))
for(i in 1:4) {
  dif_series <- dif_four_series[,i]
  pacf(dif_series, lag.max=25, main=colnames(four_series)[i])
}





dif_data <- cbind(dif_cap, dif_seal, dif_midj, dif_rul, ts[2:122,7:8]) # don't take the first element because of differencing

# Model 1: just ma process with extrenal reg
ar<-c(0)
ma<-c(1)

Model1 <- define.model(kvar=6, ar=ar, ma=ma, reg.var=c(5,6))

par(mfrow=c(1,1))
M1 <- marima(dif_data, means=1, ar.pattern=Model1$ar.pattern,
             ma.pattern=Model1$ma.pattern, Check=FALSE, Plot="log.det")

short.form(M1$ar.estimates, leading=TRUE) # print estimates
short.form(M1$ma.estimates, leading=TRUE)

acf(t(M1$residuals)[,1:4],lag.max=25) # not good, we see 4-periodicity
pacf(t(M1$residuals)[,1:4],lag.max=25)

# Model 2: ma process with seasonality and with extrenal reg
# (0,0,1)(1,0,0)[4]
ar<-c(4)
ma<-c(1)

Model2 <- define.model(kvar=6, ar=ar, ma=ma, reg.var=c(5,6))
M2 <- marima(dif_data, means=1, ar.pattern=Model2$ar.pattern,
             ma.pattern=Model2$ma.pattern, Check=FALSE, Plot="log.det",penalty=2)

short.form(M2$ar.estimates, leading=TRUE) # print estimates
short.form(M2$ma.estimates, leading=TRUE)

acf(t(M2$residuals)[,1:4],lag.max=25, na.action = na.omit) # better but could be improve especially for the 2 first series
pacf(t(M2$residuals)[,1:4],lag.max=25, na.action = na.omit) # I think it is fine

# Model 3: ma process with seasonality and with extrenal reg
# (0,0,1)(2,0,0)[4]
ar<-c(4,8)
ma<-c(1)

Model3 <- define.model(kvar=6, ar=ar, ma=ma, reg.var=c(5,6))
M3 <- marima(dif_data, means=1, ar.pattern=Model3$ar.pattern,
             ma.pattern=Model3$ma.pattern, Check=FALSE, Plot="log.det", penalty=2)


short.form(M3$ar.estimates, leading=TRUE) # print estimates
short.form(M3$ma.estimates, leading=TRUE)

acf(t(M3$residuals)[,1:4],lag.max=25, na.action = na.omit) 
pacf(t(M3$residuals)[,1:4],lag.max=25, na.action = na.omit)


# Model 4: ma process with seasonality and with independant variable
# (0,0,1)(2,0,0)[4]
ar<-c(4,8)
ma<-c(1)

Model4 <- define.model(kvar=6, ar=ar, ma=ma, indep=c(5,6)) # just change indep param
M4 <- marima(dif_data, means=1, ar.pattern=Model4$ar.pattern,
             ma.pattern=Model4$ma.pattern, Check=FALSE, Plot="log.det",penalty = 2)

## Let's improve model 4

source("/Users/user/Desktop/DTU/1. Courses/4. Spring 23/02417 - Time Series Analysis/Week 9/step.slow.marima_2017.R")
M4 <- marima(dif_data, means=1, ar.pattern=Model4$ar.pattern,
             ma.pattern=Model4$ma.pattern, Check=FALSE, Plot="log.det")
M4_sl <- step.slow(M4, dif_data) #penalty = 2


short.form(M4_sl$ar.estimates, leading=FALSE) # print estimates
short.form(M4_sl$ma.estimates, leading=FALSE)


acf(t(M4_sl$residuals)[,1:4],lag.max=25, na.action = na.omit) 
pacf(t(M4_sl$residuals)[,1:4],lag.max=25, na.action = na.omit) 


# Model 4 - Residual Diagnostics

par(mfrow = c(2,2))
titles <- c("Capital Region", "Sealand", "Middle Jutland", "Rural Areas")
for (i in 1:4){
  plot(M4_sl$residuals[i,], type='l', ylab = "Residual", main = titles[i]) 
} 

# conduct the Ljung-Box test with 12 lags (one year)
for (i in 1:4){
    print(Box.test(M4_sl$residuals[i,], type = "Ljung-Box", lag = 12))
}
# good except maybe for 2 first one

# cpgram
par(mfrow = c(2,2))
for (i in 1:4){
    cpgram(M4_sl$residuals[i,], main = titles[i]) 
}
# good 

# means
for (i in 1:4){
    print(mean(M4_sl$residuals[i,9:121]))
}
# good

# qqplot 
par(mfrow = c(2,2))
for (i in 1:4){
    qqPlot(M4_sl$residuals[i,], ylab = "Residual", main = titles[i]) 
}
# fine

# t test
for (i in 1:4){
    print((t.test(M4_sl$residuals[i,])$p.value) )
}
# good

par(mfrow=c(1,1))
plot(M4_sl$trace,type='l')
# nice convergence





## Predicting 6 quarters ahead 
pred.data <- t(data[2:128,3:8])
# add future value of regressors as in part 1
pred.data[5,124:127] <- rep(data[124,7],4)
pred.data[6,124:127] <- rep(data[124,8],4)
time <- seq(1992.5,2024.25,0.25)

#Differences poly
difference <- matrix(c(1,1,2,1,3,1,4,1), nrow=2)
Y <- define.dif(ts[1:122,3:8], difference=difference)
names(Y)
str(Y)
dif_data <- Y$y.dif
y.lost <- Y$y.lost
dif.poly <- Y$dif.poly
averages <- Y$averages

# make prediction
pred <-  arma.forecast(pred.data, nstart=121, nstep=6, marima=M4_sl, dif.poly = dif.poly, check = TRUE) 

## Time to plot
par(mfrow=c(1,1))
plot(time[1:121], pred.data[1,1:121], xlim=c(1992,2025),
     xlab = 'Time [year]', ylab = 'House Price [DKK/m2]',
     main = 'Quarterly average sales prices: Capital region of Denmark')
lines(time[2:128], pred$forecasts[1,], col=2)
pred.int <- pred$forecasts[1,122:127] + cbind(rep(0, 6), -1, 1)*qnorm(0.975)*sqrt(pred$pred.var[1,1,])
matlines(time[123:128], pred.int, lty=c(1,2,2), col=3, lwd=2 )
grid()
legend("topleft", legend = c("Training data", "1-step predictions", "Forecasts with 95% p.i."), 
       lty = c(0, 1, 1), col = 1:3, pch = c(1, NA, NA))

plot(time[1:121], pred.data[2,1:121], xlim=c(1992,2025),
     xlab = 'Time [year]', ylab = 'House Price [DKK/m2]',
     main = 'Quarterly average sales prices: Sealand')
lines(time[2:128], pred$forecasts[2,], col=2)
pred.int <- pred$forecasts[2,122:127] + cbind(rep(0, 6), -1, 1)*qnorm(0.975)*sqrt(pred$pred.var[2,2,])
matlines(time[123:128], pred.int, lty=c(1,2,2), col=3, lwd=2 )
grid()
legend("topleft", legend = c("Training data", "1-step predictions", "Forecasts with 95% p.i."), 
       lty = c(0, 1, 1), col = 1:3, pch = c(1, NA, NA))

plot(time[1:121], pred.data[3,1:121], xlim=c(1992,2025),
     xlab = 'Time [year]', ylab = 'House Price [DKK/m2]',
     main = 'Quarterly average sales prices: Middle Jutland')
lines(time[2:128], pred$forecasts[3,], col=2)
pred.int <- pred$forecasts[3,122:127] + cbind(rep(0, 6), -1, 1)*qnorm(0.975)*sqrt(pred$pred.var[3,3,])
matlines(time[123:128], pred.int, lty=c(1,2,2), col=3, lwd=2 )
grid()
legend("topleft", legend = c("Training data", "1-step predictions", "Forecasts with 95% p.i."), 
       lty = c(0, 1, 1), col = 1:3, pch = c(1, NA, NA))

plot(time[1:121], pred.data[4,1:121], xlim=c(1992,2025),
     xlab = 'Time [year]', ylab = 'House Price [DKK/m2]',
     main = 'Quarterly average sales prices: Rural areas')
lines(time[2:128], pred$forecasts[4,], col=2)
pred.int <- pred$forecasts[4,122:127] + cbind(rep(0, 6), -1, 1)*qnorm(0.975)*sqrt(pred$pred.var[4,4,])
matlines(time[123:128], pred.int, lty=c(1,2,2), col=3, lwd=2 )
grid()
legend("topleft", legend = c("Training data", "1-step predictions", "Forecasts with 95% p.i."), 
       lty = c(0, 1, 1), col = 1:3, pch = c(1, NA, NA))
