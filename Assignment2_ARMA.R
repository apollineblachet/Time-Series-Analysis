
# Q2.1.3: Simulate 10 realizations of the process with 200 observations each and plot them
n <- 200 # Number of observations
r <- 10 # Number of realization
Y.all <- matrix(NA,ncol=n, nrow=r)
for (j in 1:r){
    Y <- arima.sim(list(ar=0.8,ma=c(0.8,-0.5)),n, sd = 0.4)
    Y.all[j,] <- Y
}
matplot(t(Y.all),lty=1, type = "l", col=rainbow(11),
        xlab = "Time [t]", ylab = "Observation [Xt]",
        main="Stochastic process (10 realizations)")




# Q2.1.4: Estimate the autocorrelation function (ACF) for each realization and plot them
par(mfrow=c(1,1))
ACF.all <- matrix(NA,ncol=24, nrow=r)
colors <- rainbow(r)
for (j in 1:r){
  ACF <- acf(Y.all[j,], xlab="Lag", ylab="ACF", plot=FALSE)
  ACF.all[j,] <- ACF$acf
}
ACF.lags <- ACF$lag
matplot(ACF.lags, t(ACF.all), type="l", col=colors, xlab="Lag", ylab="ACF", main="Autocorrelation function for each realization")
for (j in 1:r){
  text(ACF.lags, ACF.all[j,]+0.05, j, pos=1, col=colors[j])
}
abline(h=c(2/sqrt(n),-2/sqrt(n)), lty=2)




# Q2.1.5: Estimate the partial autocorrelation function (PACF) for each realization and plot them
par(mfrow=c(1,1))
PACF.all <- matrix(NA, ncol=23, nrow=r)
colors <- rainbow(r)
for (j in 1:r){
  PACF <- pacf(Y.all[j,], xlab="Lag", ylab="PACF", plot=FALSE)
  PACF.all[j, 1:min(23, length(PACF$acf))] <- PACF$acf[1:min(23, length(PACF$acf))]
}
PACF.lags <- PACF$lag
matplot(PACF.lags, t(PACF.all), type="b", pch=19, col=colors, xlab="Lag", ylab="PACF", main="Partial autocorrelation function for each realization")
abline(h=c(2/sqrt(n),-2/sqrt(n)), lty=2)
legend("topright", legend=paste("Realization", 1:r), col=colors, pch=19, lty=1, ncol=2)





# Q2.1.6: Calculate the variance of each of the realizations
var.all <- rep(NA, r)
for (j in 1:r){
  var.all[j] <- var(Y.all[j,])
}
(var.all)
(mean(var.all))





# Q2.2.1: Predict the values of Yt for t=2019Q1 and 2019Q2, along with 95% prediction intervals
# Load the data
data <- read.table("A2_sales.txt", header = TRUE)
t <- seq(2014.25,2019,by=0.25)
y <- data$Sales

# Model constants
phi <- c(-1.04, 0.2, -0.86,  0.8944, -0.172)
theta <- -0.42
sigma2 <- 36963
mu <- 2070

# Center data
z <- y - mu

estimate <- function(t, z, eps_t_minus_4) { 
  -phi[1]*z[t-1] - phi[2]*z[t-2] - phi[3]*z[t-4] - phi[4]*z[t-5] - phi[5]*z[t-6] + theta*eps_t_minus_4
} 

# Estimate Yt for t = 2019Q1 using epsilon5 = 0
epsilon5 <- 0
z9_hat <- estimate(9, z, epsilon5)
epsilon9 <- z[9] - z9_hat
z13_hat <- estimate(13, z, epsilon9)
epsilon13 <- z[13] - z13_hat
z17_hat <- estimate(17, z, epsilon13)
epsilon17 <- z[17] - z17_hat
z21_hat <- estimate(21, z, epsilon17)
y21_hat <- z21_hat + mu
(y21_hat)

# Estimate Yt for t = 2019Q2 using epsilon6 = 0
epsilon6 <- 0
z10_hat <- estimate(10, z, epsilon6)
epsilon10 <- z[10] - z10_hat
z14_hat <- estimate(14, z, epsilon10)
epsilon14 <- z[14] - z14_hat
z18_hat <- estimate(18, z, epsilon14)
epsilon18 <- z[18] - z18_hat
z22_hat <- -phi[1]*z21_hat - phi[2]*z[22-2] - phi[3]*z[22-4] - phi[4]*z[22-5] - phi[5]*z[22-6] + theta*epsilon18
y22_hat <- z22_hat + mu
(y22_hat)

# 95% prediction intervals
uncertainity_y21 <- 1.96*sqrt(sigma2)
uncertainity_y22 <- 1.96*sqrt(sigma2)*sqrt(1+phi[1]^2)

# Store data
tpred <- c(2019, 2019.25, 2019.50)
ypred <-  matrix(y[20],ncol=3, nrow=3)
ypred[2,] <- c(y[20], y21_hat, y22_hat)
ypred[1,2] <- ypred[2,2]-uncertainity_y21
ypred[3,2] <- ypred[2,2]+uncertainity_y21
ypred[1,3] <- ypred[2,3]-uncertainity_y22
ypred[3,3] <- ypred[2,3]+uncertainity_y22





# Q2.2.2: Plot the actual and the predicted values
plot(t, y, type = "l", col="blue",
     xlab="Time [year]", ylab="Number of apartment sales [n]",
     xlim=c(2014,2020), ylim=c(1700,3000), cex.lab=1.2)
points(t, y, col="blue")
matlines(tpred, t(ypred), type="l",lty=c(2,1,2),lwd=2, col="#FFA500")
matpoints(tpred[2:3], t(ypred[1:3,2:3]), pch = 1, col="#FFA500")
legend("topright", c("Observation", "Predictions with 95% prediction intervals")
       , col=c("blue", "#FFA500"), lty=1, bty='n', lwd=2)
title("Sales of apartments over time")
grid()





# Simulate 300 observations of each of the four processes 100 times
phi2 <- c(0.52, 0.98)
sigma <- c(0.1, 5)
n <- 300 # Number of observations
r <- 100 # Number of realizations
sim1 <- replicate(r,arima.sim(model = list(ar=c(1.5, -phi2[1]), order=c(2,0,0)), n, sd=sigma[1]))
sim2 <- replicate(r,arima.sim(model = list(ar=c(1.5, -phi2[2]), order=c(2,0,0)), n, sd=sigma[1]))
sim3 <- replicate(r,arima.sim(model = list(ar=c(1.5, -phi2[1]), order=c(2,0,0)), n, sd=sigma[2]))
sim4 <- replicate(r,arima.sim(model = list(ar=c(1.5, -phi2[2]), order=c(2,0,0)), n, sd=sigma[2]))






# Q2.3.2: For each process, make a histogram plot of the estimates of parameter phi2 and indicate the 95% quantiles
nbin <- 20
par(mfrow=c(2,2))

# sim1
param.all1 <- matrix(NA,ncol=2, nrow=r)
for (i in 1:r){
    X <- matrix(NA,ncol=2, nrow=(n-2))
    X[,1] <- -sim1[2:(n-1),i]
    X[,2] <- -sim1[1:(n-2),i]
    Y <- sim1[3:n,i]
    param <- solve(t(X) %*% X) %*% t(X) %*% Y
    param.all1[i,] <- param
}
hist(param.all1[,2], breaks=nbin, main="Simulation 1 (phi2=0.52, sigma=0.1)", xlab="phi2 estimate", col='#BFEFFF')
abline(v=quantile(param.all1[,2], 0.025), lty=2, col='red')
abline(v=quantile(param.all1[,2], 0.975), lty=2, col='red')
quantiles <- quantile(param.all1[,2], c(0.025, 0.5, 0.975))
legend("topleft", legend=c("2.5%", "97.5%"), 
       col=c("red", "red"), lwd=c(1, 2, 1),
       lty=c(1, 1, 1), bty="n", cex=0.8)

# sim2
param.all2 <- matrix(NA,ncol=2, nrow=r)
for (i in 1:r){
    X <- matrix(NA,ncol=2, nrow=(n-2))
    X[,1] <- -sim2[2:(n-1),i]
    X[,2] <- -sim2[1:(n-2),i]
    Y <- sim2[3:n,i]
    param <- solve(t(X) %*% X) %*% t(X) %*% Y
    param.all2[i,] <- param
}
hist(param.all2[,2], breaks=nbin, main="Simulation 2 (phi2=0.98, sigma=0.1)", xlab="phi2 estimate", col='#7EB8DA')
abline(v=quantile(param.all2[,2], 0.025), lty=2, col='red')
abline(v=quantile(param.all2[,2], 0.975), lty=2, col='red')
quantiles <- quantile(param.all2[,2], c(0.025, 0.5, 0.975))
legend("topleft", legend=c("2.5%", "97.5%"), 
       col=c("red", "red"), lwd=c(1, 2, 1),
       lty=c(1, 1, 1), bty="n", cex=0.8)

# sim3
param.all3 <- matrix(NA,ncol=2, nrow=r)
for (i in 1:r){
    X <- matrix(NA,ncol=2, nrow=(n-2))
    X[,1] <- -sim3[2:(n-1),i]
    X[,2] <- -sim3[1:(n-2),i]
    Y <- sim3[3:n,i]
    param <- solve(t(X) %*% X) %*% t(X) %*% Y
    param.all3[i,] <- param
}
hist(param.all3[,2], breaks=nbin, main="Simulation 3 (phi2=0.52, sigma=5)", xlab="phi2 estimate", col='#2E73B9')
abline(v=quantile(param.all3[,2], 0.025), lty=2, col='red')
abline(v=quantile(param.all3[,2], 0.975), lty=2, col='red')
quantiles <- quantile(param.all3[,2], c(0.025, 0.5, 0.975))
legend("topleft", legend=c("2.5%", "97.5%"), 
       col=c("red", "red"), lwd=c(1, 2, 1),
       lty=c(1, 1, 1), bty="n", cex=0.8)

# sim4
param.all4 <- matrix(NA,ncol=2, nrow=r)
for (i in 1:r){
    X <- matrix(NA,ncol=2, nrow=(n-2))
    X[,1] <- -sim4[2:(n-1),i]
    X[,2] <- -sim4[1:(n-2),i]
    Y <- sim4[3:n,i]
    param <- solve(t(X) %*% X) %*% t(X) %*% Y
    param.all4[i,] <- param
}
hist(param.all4[,2], breaks=nbin, main="Simulation 4 (phi2=0.98, sigma=5)", xlab="phi2 estimate", col='#1A3B6C')
abline(v=quantile(param.all4[,2], 0.025), lty=2, col='red')
abline(v=quantile(param.all4[,2], 0.975), lty=2, col='red')
quantiles <- quantile(param.all4[,2], c(0.025, 0.5, 0.975))
legend("topleft", legend=c("2.5%", "97.5%"), 
       col=c("red", "red"), lwd=c(1, 2, 1),
       lty=c(1, 1, 1), bty="n", cex=0.8)






# Q2.3.5: Plot all the estimated pairs of parameters (phi1, phi2) for the four variations
par(mfrow=c(2,2))
# set up colors
colors1 <- c("#8B008B", "#FF8C00")  # Dark Magenta and Dark Orange
colors2 <- c("#00BFFF", "#F0E68C")  # Deep Sky Blue and Khaki
colors3 <- c("#32CD32", "#FF6347")  # Lime Green and Tomato
colors4 <- c("#00FA9A", "#FF1493")  # Medium Spring Green and Deep Pinκ

for (i in 1:4) {
  if (i == 1) {
    matplot(param.all1, lty=1, type="l", main="Simulation 1 (φ2=0.52, sigma=0.1)", 
            xlab="Realization", ylab='Parameters (phi1, phi2)', col=colors1, cex.lab=1.2)
    legend("topright", legend=c("phi1", "phi2"), lty=1, col=colors1, bty='n')
  } else if (i == 2) {
    matplot(param.all2, lty=1, type="l", main="Simulation 2 (phi2=0.98, sigma=0.1)", 
            xlab="Realization", ylab='Parameters (phi1,phi2)', col=colors2, cex.lab=1.2)
    legend("topright", legend=c("phi1", "phi2"), lty=1, col=colors2, bty='n')
  } else if (i == 3) {
    matplot(param.all3, lty=1, type="l", main="Simulation 3 (phi2=0.52, sigma=5)", 
            xlab="Realization", ylab='Parameters (phi1,phi2)', col=colors3, cex.lab=1.2)
    legend("topright", legend=c("phi1", "phi2"), lty=1, col=colors3, bty='n')
  } else {
    matplot(param.all4, lty=1, type="l", main="Simulation 4 (phi2=0.98, sigma=5)", 
            xlab="Realization", ylab='Parameters (phi1,phi2)', col=colors4, cex.lab=1.2)
    legend("topright", legend=c("phi1", "phi2"), lty=1, col=colors4, bty='n')
  }
}
