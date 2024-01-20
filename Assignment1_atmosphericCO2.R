

# Load the data
data <- read.table("A1_co2.txt", header = TRUE)

# Split the data into training and testing sets
xtrain <- data$time[1:718]
ytrain <- data$co2[1:718]
xtest <- data$time[719:738]
ytest <- data$co2[719:738]

# Get the number of observations in training set
ntrain <- length(xtrain)

# Plot the training and testing sets
par(mfrow=c(1,1))
plot(xtrain, ytrain, type="l", col="blue",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="Atmospheric CO2 as a function of time", ylim=c(300, max(data$co2)))
lines(xtest, ytest, col="red")
legend("topleft", c("Training set", "Testing set"), col=c("blue", "red"), lty=1, bty='n', lwd=2)

# Boxplot of CO2 concentration
par(mfrow=c(1,1))
boxplot(data$co2, main="CO2 Concentration Boxplot", ylab="Atmospheric CO2 [ppm]", col="lightblue",
        medcol="red", medlwd=2, medpch=19, medcex=1.5, 
        boxcol="lightblue", whiskcol="darkgray")
# Add median value to the plot
text(1, median(data$co2) + 5, round(median(data$co2),2), col="red", pos=4, cex=1.2)



# Check one harmonic period
par(mfrow=c(1,1))
plot(xtrain[1:24], ytrain[1:24], type="o", col="blue",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="Atmospheric CO2 as a function of time", ylim=c(min(data$co2[1:24]), max(data$co2[1:24])))
lines(xtest, ytest, col="red")
legend("topleft", c("Training set"), col=c("blue"), lty=1, bty='n', lwd=2)

# Set the period of the harmonic part
p <- 1 # we've set p to 1 because we can see from the observations that the period of the harmonic part is approximately 1 (1 year)

# Construct the design matrices for the training and testing sets, including Fourier terms to capture periodicity in the data
Xtrain <- cbind(1, xtrain, sin(2*pi*xtrain/p), cos(2*pi*xtrain/p))
Xtest <- cbind(1, xtest, sin(2*pi*xtest/p), cos(2*pi*xtest/p))



# Compute the OLS estimate of the regression coefficients
thetahatOLS <- solve(t(Xtrain) %*% Xtrain) %*% t(Xtrain) %*% ytrain

# Compute the predicted values of CO2 concentration for the training and testing sets
yhatOLS_train <- Xtrain %*% thetahatOLS
yhatOLS_test <- Xtest %*% thetahatOLS

# Compute residual for the training and testing sets
epsOLS_train <- ytrain - yhatOLS_train 
epsOLS_test <- ytest - yhatOLS_test 

# Compute measure of uncertainty for each of the estimates
sigma2OLS <- t(epsOLS_train) %*% epsOLS_train / (length(xtrain) - length(thetahatOLS))
varOLS <- sigma2OLS[1] * solve(t(Xtrain) %*% Xtrain)
stdOLS <- sqrt(diag(varOLS))

# Compute MSE 
MSE_OLS_train <- t(epsOLS_train) %*% epsOLS_train / length(epsOLS_train)
MSE_OLS_test <- t(epsOLS_test) %*% epsOLS_test / length(epsOLS_test)

# Plot the training and testing sets
par(mfrow=c(1,1))
plot(xtrain, ytrain, type="l", col="blue",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="Ordinary least squares model (OLS) data fitting", ylim=c(300, max(data$co2)))
lines(xtest, ytest, col="red")

# Add the OLS fit to the plot
lines(xtrain, yhatOLS_train, col="green")
lines(xtest, yhatOLS_test, col="orange")
legend("topleft", c("Training set", "Testing set", "OLS Training set", "OLS Testing set"), col=c("blue", "red", "green", "orange"), lty=1, bty='n', lwd=2)
abline(v=xtest[1], col="black", lty=2)




# Initial guess of correlation structure
rau <- 0.8
Sigmay <- diag(ntrain) 
for (i in 1:ntrain) { 
    for (j in 1:ntrain) { 
        Sigmay[i,j] <- rau^abs(i-j)
    }
}

for(i in 1:5){ 
     # Estimate parameters using currently assumed correlation structure
     thetahatWLS <- solve(t(Xtrain) %*% solve(Sigmay) %*% Xtrain) %*% (t(Xtrain) %*% solve(Sigmay) %*% ytrain) 

     # Compute residuals for these parameter estimates
     epsWLS <- ytrain - Xtrain %*% thetahatWLS 

     # Select the value for Sigamy which reflects the correlation and variance structure
     rau <- cor(epsWLS[1:ntrain-1], epsWLS[2:ntrain])
     for (i in 1:ntrain) { 
          for (j in 1:ntrain) { 
               Sigmay[i,j] <- rau^abs(i-j)
          }
     }
}

# Compute the predicted values of CO2 concentration for the training and testing sets
yhatWLS_train <- Xtrain %*% thetahatWLS
yhatWLS_test <- Xtest %*% thetahatWLS

# Compute residual for the training and testing sets
epsWLS_train <- ytrain - yhatWLS_train 
epsWLS_test <- ytest - yhatWLS_test 

# Compute measure of uncertainty for each of the estimates
sigma2WLS <- t(epsWLS_train) %*% solve(Sigmay) %*% epsWLS_train / (length(xtrain) - length(thetahatWLS))
varWLS <- sigma2WLS[1] * solve(t(Xtrain) %*% solve(Sigmay) %*% Xtrain)
stdWLS <- sqrt(diag(varWLS))

# Compute MSE 
MSE_WLS_train <- t(epsWLS_train) %*% epsWLS_train / length(epsWLS_train)
MSE_WLS_test <- t(epsWLS_test) %*% epsWLS_test / length(epsWLS_test)

# Plot the training and testing sets
par(mfrow=c(1,1))
plot(xtrain, ytrain, type="l", col="blue",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="Weighted least squares model (WLS) data fitting", ylim=c(300, max(data$co2)))
lines(xtest, ytest, col="red")

# Add the WLS fit to the plot
lines(xtrain, yhatWLS_train, col="Black")
lines(xtest, yhatWLS_test, col="Orange")
legend("topleft", c("Training set", "Testing set",  "WLS Training set", "WLS Testing set"), col=c("blue", "red", "black", "orange"), lty=1, bty='n', lwd=2)
abline(v=xtest[1], col="black", lty=2)






## Local linear trend model for the data.
lambda <- 0.9
f <- function(j) rbind(1, j, sin(2*pi*j/p), cos(2*pi*j/p))
L <- t(matrix(c(1,0,0,0, 1/12,1,0,0, 0,0,cos(2*pi/p/12),sin(2*pi/p/12), 0,0,-sin(2*pi/p/12),cos(2*pi/p/12)),ncol=4))
LInv <- solve(L)

# 4 parameters so at 4 observations are needed to estimate theta 
# however one more is needed to get an estimate of the uncertainty

init <- 5 # Skip estimating for the first 10 observations
## FNinit & hNinit (First observations) using equation (3.100)
F <- matrix(0, nrow=4, ncol=4)
h <- matrix(0,nrow=4,ncol=1)
for (j in 0:(init-1)){
  F <- F + lambda^(j) * f(-j/12) %*% t(f(-j/12))
  h <- h + lambda^(j) * f(-j/12) * ytrain[init-j]
}

## Allocating space
np <- length(h)
theta.all <- matrix(NA,ncol=np, nrow=ntrain)
sigma.all <- rep(NA, ntrain)
sd.err.all <- rep(NA, ntrain)
yhat.all <- rep(NA, ntrain)
err.all <- rep(NA, ntrain)

## Solving at time init
theta.hat <- solve(F, h)
theta.all[init,] <- theta.hat
epsilon <- ytrain[1:init] - cbind(1, (-(init-1):0)/12, sin(2*pi*(-(init-1):0)/12/p), cos(2*pi*(-(init-1):0)/12/p)) %*% theta.hat  

T <- 0
for (j in 0:(init-1)){
     T <- T + lambda^(j)
}

Sigma_inv <- diag(init)
for (j in 0:(init-1)){
     Sigma_inv[j,j] <- lambda^((init-1-j))
}

sigma.all[init] <- sqrt(t(epsilon) %*% Sigma_inv %*% epsilon/(T - np))
sd.err.all[init] <- sigma.all[init] * sqrt(1 + t(f(1/12)) %*% solve(F) %*% f(1/12))
yhat.all[init+1] <- t(f(1/12)) %*% theta.hat
err.all[init+1] <- ytrain[init+1] - yhat.all[init+1]

## Looping over the remaining observations
for (i in (init+1):(ntrain)){
  F <- F + lambda^((i-1)) * f(- (i-1)/12) %*% t(f(- (i-1)/12))
  h <- lambda * LInv %*% h + f(0)*ytrain[i]
  theta.hat <-  solve(F, h)  
  theta.all[i,] <- theta.hat

  ## Adding uncertainty information
  epsilon <- ytrain[1:i] - cbind(1, (-(i-1):0)/12, sin(2*pi*(-(i-1):0)/12/p), cos(2*pi*(-(i-1):0)/12/p)) %*% theta.hat
  T <- 0
  for (j in 0:(i-1)){
     T <- T + lambda^(j)
  }

  Sigma_inv <- diag(i)
  for (j in 0:(i-1)){
     Sigma_inv[j,j] <- lambda^((i-1-j))
  }

  sigma.all[i] <- sqrt(t(epsilon) %*% Sigma_inv %*% epsilon/(T - np))

  ## Estimating s.d. of estimated parameters
  sd.err.all[i] <- sigma.all[i] * sqrt(1 + t(f(1/12)) %*% solve(F) %*% f(1/12))

  yhat.all[i+1] <- t(f(1/12)) %*% theta.hat
  err.all[i+1] <- ytrain[i+1] - yhat.all[i+1]
}

## Predictions on test set
theta_pred <- theta.all[718,] # Get last theta of the training set
sigma_pred <- sigma.all[718] # Get last sigma

ntest <- length(xtest)

# Allocate memory
yhat_test.all <- rep(NA, ntest)
err.test.all <- rep(NA, ntest)
sd.err_test.all <- rep(NA, ntest)

for (j in 1:ntest){
     yhat_test.all[j] <- t(f(j/12)) %*% theta_pred # Make prediction
     sd.err_test.all[j] <- sigma_pred * sqrt(1 + t(f(j/12)) %*% solve(F) %*% f(j/12)) # Calculate standard deviation of the error
     err.test.all[j] <- ytest[j+1] - yhat_test.all[j+1]
}





## Plot Q3.3 ##
par(mfrow=c(1,1))
plot(xtrain[10:length(xtrain)], err.all[10:length(xtrain)], type="l", col="blue",
     xlab="Time [year] ", ylab="One step prediction errors",
     main="One step prediction errors for the local linear trend model (λ=0.9)",
     ylim=c(-2.5, 2.5))
lines(xtest, err.test.all, type="l", col="red")
lines(xtrain[10:length(xtrain)], sigma.all[10:length(xtrain)], col='purple')
lines(xtest, sd.err_test.all, col='orange')
legend("topleft", c("One step prediction errors trainind data", "One step prediction errors test data","Sigma estimate training data","Sigma estimate test data"), col=c("blue", "red","purple", "orange"), lty=1, bty='n', lwd=2)
abline(v=xtest[1], col="red", lty=2)

## Plot Q3.4 ##
t.quan <- qt(p = 0.975, df=rep(ntrain-np, ntrain-9)) # Get t distribution for training
t.quan[1:init] <- NA

t_test.quan <- qt(p = 0.975, df=rep(ntest-np, ntest))
t_test.quan[1:init] <- NA

par(mfrow=c(1,1))
plot(xtrain[10:length(xtrain)], ytrain[10:length(xtrain)], type="l", col="blue",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="One step prediction errors of local linear trend model (λ=0.9)", ylim=c(300, max(data$co2)))
lines(xtest, ytest, col="red")
#lines(xtrain[10:length(xtrain)], yhat.all[10:length(xtrain)],  col="darkgreen")
matlines(xtrain[10:length(xtrain)], yhat.all[10:length(xtrain)] + t.quan * cbind(0,-sd.err.all[10:length(xtrain)],sd.err.all[10:length(xtrain)]),type="l",lty=c(1,2,2),lwd=2, col="darkgreen")
matlines(xtest, yhat_test.all + t_test.quan * cbind(0,-sd.err_test.all,sd.err_test.all),type="l",lty=c(1,2,2),lwd=2, col="purple")
legend("topleft", c("Training set", "Testing set", "Predictions training set -- 95% prediction interval", "Predictions test set -- 95% prediction interval"), col=c("blue", "red", "darkgreen", "purple"), lty=1, bty='n', lwd=2)

## Plot Q3.5 ##
t.quan_test <- qt(p = 0.975, df= rep(ntest-np, ntest)) # Get t distribution for testing

par(mfrow=c(1,1))
plot(xtrain[10:length(xtrain)], ytrain[10:length(xtrain)], col="blue",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="Predictions after 2010", 
     xlim=c(2010, 2020), ylim=c(380, 420))
points(xtest, ytest, col="red")
matlines(xtrain[10:length(xtrain)], yhat.all[10:length(xtrain)] + t.quan * cbind(0,-sd.err.all[10:length(xtrain)],sd.err.all[10:length(xtrain)]),type="l",lty=c(1,2,2),lwd=2, col="darkgreen")
matlines(xtest, yhat_test.all + t.quan_test * cbind(0,-sd.err_test.all, sd.err_test.all),type="l",lty=c(1,2,2),lwd=2, col="purple")
legend("topleft", c("Training set", "Testing set", "Train predictions with 95% prediction interval", "Test predictions with 95% prediction interval"), col=c("blue", "red", "darkgreen", "purple"), lty=1, bty='n', lwd=2)

# Zoom in test data
par(mfrow=c(1,1))
plot(xtest, ytest, col="red",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="Predictions 2018 and onwards",
     ylim=c(404, 418))
matlines(xtest, yhat_test.all + t.quan_test * cbind(0,-sd.err_test.all, sd.err_test.all),type="l",lty=c(1,2,2),lwd=2, col="purple")
legend("topleft", c("Testing set", "Test predictions with 95% prediction interval"), col=c("red", "purple"), lty=1, bty='n', lwd=2)

## Question 1.3.6: Predictions 1, 2, 6, 12 and 20 months ahead & Question 1.3.7
j.all <- c(1,2,6,12,20)
yhat_pred.all <- rep(NA, length(j.all))
dif.pred <- rep(NA, length(j.all))

for (j in 1:length(j.all)){
     yhat_pred.all[j] <- t(f(j.all[j]/12)) %*% theta_pred # Make prediction
}

MSE <- sum((ytest-yhat_test.all)^2)/ntest

## Plot Q3.8 ##
par(mfrow=c(1,1))
plot(xtrain[10:length(xtrain)], ytrain[10:length(xtrain)], type="l", col="blue",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="Estimated mean for each time step", ylim=c(300, max(data$co2)))
lines(xtest, ytest, col="red")
lines(xtrain[10:length(xtrain)], theta.all[10:length(xtrain),1], col='darkgreen')
lines(xtest, rep(theta.hat[1], length(xtest)), col='purple')
legend("topleft", c("Training set", "Testing set", "Estimated mean train set", "Estimated mean test set"), col=c("blue", "red", "darkgreen", "purple"), lty=1, bty='n', lwd=2)






### Find optimal lambda ###

# Extract burning period of 100 months
xtrain <- data$time[523:622]
ytrain <- data$co2[523:622]

xtest <- data$time[719:738]
ytest <- data$co2[719:738]

# Get the number of observations in training set
ntrain <- length(xtrain)

## Construct design matrix ##

# Set the period of the harmonic part
p <- 1 # we've set p to 1 because we can see from the observations that the period of the harmonic part is approxiately 1 (1 year)

# Construct the design matrices for the training and testing sets, including Fourier terms to capture periodicity in the data
Xtrain <- cbind(1, xtrain, sin(2*pi*xtrain/p), cos(2*pi*xtrain/p))
Xtest <- cbind(1, xtest, sin(2*pi*xtest/p), cos(2*pi*xtest/p))


## Local linear trend model for the data.
f <- function(j) rbind(1, j, sin(2*pi*j/p), cos(2*pi*j/p))
L <- t(matrix(c(1,0,0,0, 1/12,1,0,0, 0,0,cos(2*pi/p/12),sin(2*pi/p/12), 0,0,-sin(2*pi/p/12),cos(2*pi/p/12)),ncol=4))
LInv <- solve(L)

lambdas <- seq(from = 0.001, to = 0.9999, by = 0.001)
MSE.all <- rep(NA, length(lambdas))
c <- 1
for (lambda in lambdas){

# 4 parameters so at 4 observations are needed to estimate theta 
# however one more is needed to get an estimate of the uncertainty

init <- 5 # Skip estimating for the first 10 observations
## FNinit & hNinit (First observations) using equation (3.100)
F <- matrix(0, nrow=4, ncol=4)
h <- matrix(0,nrow=4,ncol=1)
for (j in 0:(init-1)){
  F <- F + lambda^(j) * f(-j/12) %*% t(f(-j/12))
  h <- h + lambda^(j) * f(-j/12) * ytrain[init-j]
}

## Allocating space
np <- length(h)
theta.all <- matrix(NA,ncol=np, nrow=ntrain)
sigma.all <- rep(NA, ntrain)
sd.err.all <- rep(NA, ntrain)
yhat.all <- rep(NA, ntrain)
err.all <- rep(NA, ntrain)

## Solving at time init
theta.hat <- solve(F, h)
theta.all[init,] <- theta.hat
epsilon <- ytrain[1:init] - cbind(1, (-(init-1):0)/12, sin(2*pi*(-(init-1):0)/12/p), cos(2*pi*(-(init-1):0)/12/p)) %*% theta.hat  

T <- 0
for (j in 0:(init-1)){
     T <- T + lambda^(j)
}

Sigma_inv <- diag(init)
for (j in 0:(init-1)){
     Sigma_inv[j,j] <- lambda^((init-1-j))
}

sigma.all[init] <- sqrt(t(epsilon) %*% Sigma_inv %*% epsilon/(T - np))
sd.err.all[init] <- sigma.all[init] * sqrt(1 + t(f(1/12)) %*% solve(F) %*% f(1/12))
yhat.all[init+1] <- t(f(1/12)) %*% theta.hat
err.all[init+1] <- ytrain[init+1] - yhat.all[init+1]

## Looping over the remaining observations
for (i in (init+1):(ntrain)){
  F <- F + lambda^((i-1)) * f(- (i-1)/12) %*% t(f(- (i-1)/12))
  h <- lambda * LInv %*% h + f(0)*ytrain[i]
  theta.hat <-  solve(F, h)  
  theta.all[i,] <- theta.hat

  ## Adding uncertainty information
  epsilon <- ytrain[1:i] - cbind(1, (-(i-1):0)/12, sin(2*pi*(-(i-1):0)/12/p), cos(2*pi*(-(i-1):0)/12/p)) %*% theta.hat
  T <- 0
  for (j in 0:(i-1)){
     T <- T + lambda^(j)
  }

  Sigma_inv <- diag(i)
  for (j in 0:(i-1)){
     Sigma_inv[j,j] <- lambda^((i-1-j))
  }

  sigma.all[i] <- sqrt(t(epsilon) %*% Sigma_inv %*% epsilon/(T - np))

  ## Estimating s.d. of estimated parameters
  sd.err.all[i] <- sigma.all[i] * sqrt(1 + t(f(1/12)) %*% solve(F) %*% f(1/12))

  yhat.all[i+1] <- t(f(1/12)) %*% theta.hat
  err.all[i+1] <- ytrain[i+1] - yhat.all[i+1]
  
}

# Compute and store MSE 
MSE <- (t(err.all[(init+1):100]) %*% err.all[(init+1):100])/length(err.all[(init+1):100])
MSE.all[c] <- MSE
c <- c + 1
}

plot(lambdas, MSE.all, col="blue",
     xlab= 'λ value', ylab='MSE',
     main='Optimization of λ')

lambda_opt <- lambdas[which.min(MSE.all)] # 0.969

### Question 4.1 ####

# Split the data into training and testing sets
xtrain <- data$time[1:718]
ytrain <- data$co2[1:718]
xtest <- data$time[719:738]
ytest <- data$co2[719:738]

# Get the number of observations in training set
ntrain <- length(xtrain)

## Construct design matrix ##

# Set the period of the harmonic part
p <- 1 # we've set p to 1 because we can see from the observations that the period of the harmonic part is approxiately 1 (1 year)

# Construct the design matrices for the training and testing sets, including Fourier terms to capture periodicity in the data
Xtrain <- cbind(1, xtrain, sin(2*pi*xtrain/p), cos(2*pi*xtrain/p))
Xtest <- cbind(1, xtest, sin(2*pi*xtest/p), cos(2*pi*xtest/p))

## Local linear trend model for the data.
lambda <- lambda_opt
f <- function(j) rbind(1, j, sin(2*pi*j/p), cos(2*pi*j/p))
L <- t(matrix(c(1,0,0,0, 1/12,1,0,0, 0,0,cos(2*pi/p/12),sin(2*pi/p/12), 0,0,-sin(2*pi/p/12),cos(2*pi/p/12)),ncol=4))
LInv <- solve(L)

# 4 parameters so at 4 observations are needed to estimate theta 
# however one more is needed to get an estimate of the uncertainty.

init <- 5 # Skip estimating for the first 10 observations
## FNinit & hNinit (First observations) using equation (3.100)
F <- matrix(0, nrow=4, ncol=4)
h <- matrix(0,nrow=4,ncol=1)
for (j in 0:(init-1)){
  F <- F + lambda^(j) * f(-j/12) %*% t(f(-j/12))
  h <- h + lambda^(j) * f(-j/12) * ytrain[init-j]
}

## Allocating space
np <- length(h)
theta.all <- matrix(NA,ncol=np, nrow=ntrain)
sigma.all <- rep(NA, ntrain)
sd.err.all <- rep(NA, ntrain)
yhat.all <- rep(NA, ntrain)
err.all <- rep(NA, ntrain)

## Solving at time init
theta.hat <- solve(F, h)
theta.all[init,] <- theta.hat
epsilon <- ytrain[1:init] - cbind(1, (-(init-1):0)/12, sin(2*pi*(-(init-1):0)/12/p), cos(2*pi*(-(init-1):0)/12/p)) %*% theta.hat  

T <- 0
for (j in 0:(init-1)){
     T <- T + lambda^(j)
}

Sigma_inv <- diag(init)
for (j in 0:(init-1)){
     Sigma_inv[j,j] <- lambda^((init-1-j))
}

sigma.all[init] <- sqrt(t(epsilon) %*% Sigma_inv %*% epsilon/(T - np))
sd.err.all[init] <- sigma.all[init] * sqrt(1 + t(f(1/12)) %*% solve(F) %*% f(1/12))
yhat.all[init+1] <- t(f(1/12)) %*% theta.hat
err.all[init+1] <- ytrain[init+1] - yhat.all[init+1]

## Looping over the remaining observations
for (i in (init+1):(ntrain)){
  F <- F + lambda^((i-1)) * f(- (i-1)/12) %*% t(f(- (i-1)/12))
  h <- lambda * LInv %*% h + f(0)*ytrain[i]
  theta.hat <-  solve(F, h)  
  theta.all[i,] <- theta.hat

  ## Adding uncertainty information
  epsilon <- ytrain[1:i] - cbind(1, (-(i-1):0)/12, sin(2*pi*(-(i-1):0)/12/p), cos(2*pi*(-(i-1):0)/12/p)) %*% theta.hat
  T <- 0
  for (j in 0:(i-1)){
     T <- T + lambda^(j)
  }

  Sigma_inv <- diag(i)
  for (j in 0:(i-1)){
     Sigma_inv[j,j] <- lambda^((i-1-j))
  }

  sigma.all[i] <- sqrt(t(epsilon) %*% Sigma_inv %*% epsilon/(T - np))

  ## Estimating s.d. of estimated parameters
  sd.err.all[i] <- sigma.all[i] * sqrt(1 + t(f(1/12)) %*% solve(F) %*% f(1/12))

  yhat.all[i+1] <- t(f(1/12)) %*% theta.hat
  err.all[i+1] <- ytrain[i+1] - yhat.all[i+1]
}


## Predictions on test set
theta_pred <- theta.all[718,] # Get last theta of the training set
sigma_pred <- sigma.all[718] # Get last sigma

ntest <- length(xtest)

# Allocate memory
yhat_test.all <- rep(NA, ntest) 
sd.err_test.all <- rep(NA, ntest)

for (j in 1:ntest){
     yhat_test.all[j] <- t(f(j/12)) %*% theta_pred # Make prediction
     sd.err_test.all[j] <- sigma_pred * sqrt(1 + t(f(j/12)) %*% solve(F) %*% f(j/12)) # Calculate standard deviation of the error
}

## Compute MSE
MSE_test <- sum((ytest-yhat_test.all)^2)/ntest
MSE_train <- sum((ytrain[10:length(xtrain)]-yhat.all[10:length(xtrain)])^2)/ntest

## Plot ##
t.quan <- qt(p = 0.975, df=rep(ntrain-np, ntrain-9) ) # Get t distribution for training
t.quan[1:init] <- NA
t.quan_test <- qt(p = 0.975, df= rep(ntest-np, ntest)) # Get t distribution for testing

# Plot the training and testing sets
plot(xtrain[10:length(xtrain)], ytrain[10:length(xtrain)], col="blue",
     xlab="Time [year]", ylab="Atmospheric CO2 [ppm]",
     main="Predictions after 2010 (λ=0.983)", 
     xlim=c(2010, 2020), ylim=c(380, 420))
points(xtest, ytest, col="red")
matlines(xtrain[10:length(xtrain)], yhat.all[10:length(xtrain)] + t.quan * cbind(0,-sd.err.all[10:length(xtrain)],sd.err.all[10:length(xtrain)]),type="l",lty=c(1,2,2),lwd=2, col="darkgreen")
matlines(xtest, yhat_test.all + t.quan_test * cbind(0,-sd.err_test.all, sd.err_test.all),type="l",lty=c(1,2,2),lwd=2, col="purple")
legend("topleft", c("Training set", "Testing set", "Train predictions with 95% prediction interval", "Test predictions with 95% prediction interval"), col=c("blue", "red", "darkgreen", "purple"), lty=1, bty='n', lwd=2)

