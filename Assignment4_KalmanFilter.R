
# Load data
data <- read.csv("A4_Kulhuse.csv")
str(data)

# Check for Nan values
(num_missing <- colSums(is.na(data))) # 111 missing values per column
# Convert Datetime from char to datetime object
data$DateTime <- as.POSIXct(data$DateTime, format="%Y-%m-%d %H:%M:%S")
str(data)




## Dissolved Oxygen over Time
par(mfrow=c(1,1))
plot(data$DateTime, data$ODO, col='skyblue',
     xlab='Time [month]', ylab='Dissolved Oxygen [mg/L]',
     main = 'Dissolved Oxygen - Time ')
abline(h=median(data$ODO), col='darkred', lty=2)

## Dissolved Oxygen over Time - Zoom In
plot(data$DateTime, data$ODO, col='skyblue', ylim=c(7.5,11.5),
     xlab='Time [month]', ylab='Dissolved Oxygen [mg/L]',
     main = 'Dissolved Oxygen - Time (Zoom In)')
abline(h=median(data$ODO), col='darkred', lty=2)

## Salinity over Time
par(mfrow=c(1,1))
plot(data$DateTime, data$Sal, col='seagreen',
     xlab='Time [month]', ylab='Water Salinity [g/kg]',
     main = 'Water Salinity - Time ')
abline(h=median(data$Sal), col='darkred', lty=2)

## Salinity over Time - Zoom In
plot(data$DateTime, data$Sal, col='seagreen', ylim=c(16,22),
     xlab='Time [month]', ylab='Water Salinity [g/kg]',
     main = 'Water Salinity - Time (Zoom In) ')
abline(h=median(data$Sal), col='darkred', lty=2)




# random walk model
x0 <- data$Sal[1]
A <- matrix(1)
B <- matrix(0)
C <- matrix(1)
# values of sigma do not have to be specified at this question
Sigma1 <- matrix(0.01)
Sigma2 <- matrix(0.005)




# Kalman filtering using manual implementation
kf1 <-  kalman(data$Sal, A= A, B=B, C=C, Sigma.1=Sigma1, Sigma.2=Sigma2, 
             V0=Sigma1, Xhat0=matrix(x0),n.ahead=1,verbose=TRUE)
str(kf1)

# One-step predictions plot
par(mfrow=c(1,1))
plot(data$Sal, col='seagreen',
     xlab='Time [month]', ylab='Water Salinity [g/kg]',
     main = 'Water Salinity - Time', xaxt="n")
with(kf1, matlines(pred[,1]+sqrt(Sigma.yy.pred[1,1,])%*%cbind(0,-1.96,1.96),
                   type="l", lty=c(1,2,2), col='#FF5733'))
legend('topleft', 
       legend = c('Observations', 'One-step predictions'), 
       col = c('seagreen', '#FF5733'),
       lty = c(1, 1),
       bty = 'n')

# Get indices corresponding to the first day of each month
first_day_indices <- which(format(data$DateTime, "%d") == "01")

# Set x-axis labels to month abbreviations
month_labels <- format(data$DateTime[first_day_indices], "%b")
axis(1, at=first_day_indices, labels=month_labels)

# Standatrdized one-step predictions plot
par(mfrow=c(1,1))
plot((data$Sal[-1]-kf1$pred[2:5000])/sqrt(kf1$Sigma.yy.pred[1,1,2:5000]), type='l', col='#FF5733',
     xlab='Time [month]', ylab='Predicition Error [g/kg]',
     main='Standardized One-Step Prediction Errors',xaxt="n") # ylim=c(-2,2)
axis(1, at=first_day_indices, labels=month_labels)

# Zoom in (One-step predictions)
par(mfrow=c(1,1))
plot(data$Sal[800:950], col='seagreen',
     xlab='Time [day]', ylab='Water Salinity [g/kg]',
     main = 'Water Salinity - Time (Zoom In)', xaxt="n")
with(kf1, matlines(pred[800:950,1]+sqrt(Sigma.yy.pred[1,1,800:950])%*%cbind(0,-1.96,1.96),
                   type="l", lty=c(1,2,2), col='#FF5733'))
legend('bottomright', 
       legend = c('Observations', 'One-step predictions with 95% PI'), 
       col = c('seagreen', '#FF5733'),
       lty = c(1, 1),
       bty = 'n')
axis(1, at=seq(1, 150, by=48), labels=format(data$DateTime[seq(800, 950, by=48)], '%m-%d'))

# Zoom in (Standardized errors)
par(mfrow=c(1,1))
plot((data$Sal[800:950]-kf1$pred[800:950])/sqrt(kf1$Sigma.yy.pred[1,1,800:950]), type='l', col='#FF5733',
     xlab='Time [day]', ylab='Predicition Error [g/kg]',
     main='Standardized One-Step Prediction Errors (Zoom In)',xaxt="n")
axis(1, at=seq(1, 150, by=48), labels=format(data$DateTime[seq(800, 950, by=48)], '%m-%d'))





# Skip outliers
sd(na.omit(data$Sal))*6 

# Kalman filtering using manual implementation
kf2 <-  kalman_outliers2(data$Sal, A= A, B=B, C=C, Sigma.1=Sigma1, Sigma.2=Sigma2, 
             V0=Sigma1, Xhat0=matrix(x0),n.ahead=1,verbose=TRUE) 
str(kf2)

# Zoom in (One-step predictions)
par(mfrow=c(1,1))
plot(data$Sal[800:950], col='seagreen', #ylim=c(15,19),
     xlab='Time [day]', ylab='Water Salinity [g/kg]',
     main = 'Water Salinity - Time (Zoom In)', xaxt="n")
with(kf2, matlines(pred[800:950,1]+sqrt(Sigma.yy.pred[1,1,800:950])%*%cbind(0,-1.96,1.96),
                   type="l", lty=c(1,2,2), col='#FF5733'))
legend('bottomright', 
       legend = c('Observations', 'One-step predictions with 95% PI'), 
       col = c('seagreen', '#FF5733'),
       lty = c(1, 1),
       bty = 'n')
axis(1, at=seq(1, 150, by=48), labels=format(data$DateTime[seq(800, 950, by=48)], '%m-%d'))

kf2$NB.outsider
kf2$Outsider.index

# Calculating the log-likelihood
neps1 <- (data$Sal[-1]-kf1$pred[2:5000,1])^2 / kf1$Sigma.yy.pred[1,1,2:5000]
-0.5 * sum(neps1 + log(kf1$Sigma.yy.pred[1,1,2:5000]),na.rm = TRUE)

neps1 <- (data$Sal[-1]-kf2$pred[2:5000,1])^2 / kf2$Sigma.yy.pred[1,1,2:5000]
-0.5 * sum(neps1 + log(kf2$Sigma.yy.pred[1,1,2:5000]),na.rm = TRUE)





# Find the ML estimates of the two parameters using the first 800 observations.
Ypart <- data$Sal[1:800]
my.obj.part <- function(par){
  Kro <- kalman_outliers2(Ypart, A= A, B=B, C=C, Sigma.1=matrix(par[1]), Sigma.2=matrix(par[2]),
                V0=matrix(par[1]), Xhat0=matrix(x0),n.ahead=1,verbose=TRUE)
  nepso <- (Ypart[-1]-Kro$pred[2:800,1])^2 / Kro$Sigma.yy.pred[1,1,2:800]
  return(0.5 * sum(nepso + log(Kro$Sigma.yy.pred[1,1,2:800]), na.rm = TRUE))
}

(Kmopt.part <- optim(c(0.01,0.005),my.obj.part,method = "L-BFGS-B", lower = c(8e-4,4e-4)))

# Filter the data with the optimal parameters and plot
 Sigma1 <- matrix(Kmopt.part$par[1])
 Sigma2 <- matrix(Kmopt.part$par[2])

# Kalman filtering using manual implementation
kf3 <-  kalman_outliers2(data$Sal, A= A, B=B, C=C, Sigma.1=Sigma1, Sigma.2=Sigma2, 
             V0=Sigma1, Xhat0=matrix(x0),n.ahead=1,verbose=TRUE) 
str(kf3)

# Zoom in (One-step predictions)
par(mfrow=c(1,1))
plot(data$Sal[800:950], col='seagreen', #ylim=c(15,19),
     xlab='Time [day]', ylab='Water Salinity [g/kg]',
     main = 'Water Salinity - Time (Zoom In)', xaxt="n")
with(kf3, matlines(pred[800:950,1]+sqrt(Sigma.yy.pred[1,1,800:950])%*%cbind(0,-1.96,1.96),
                   type="l", lty=c(1,2,2), col='#FF5733'))
legend('bottomright', 
       legend = c('Observations', 'One-step predictions with 95% PI'), 
       col = c('seagreen', '#FF5733'),
       lty = c(1, 1),
       bty = 'n')
axis(1, at=seq(1, 150, by=48), labels=format(data$DateTime[seq(800, 950, by=48)], '%m-%d'))

kf3$NB.outsider
kf3$Outsider.index

# Calculating the log-likelihood
neps1 <- (data$Sal[-1]-kf3$pred[2:5000,1])^2 / kf3$Sigma.yy.pred[1,1,2:5000]
-0.5 * sum(neps1 + log(kf3$Sigma.yy.pred[1,1,2:5000]),na.rm = TRUE)




kalman <- function(Y,A,B=NULL,u=NULL,C,Sigma.1=NULL,Sigma.2=NULL,debug=FALSE,V0=Sigma.1,Xhat0=NULL,n.ahead=1,skip=0,verbose=FALSE){

  ## predictions through data are one-step predictions. n.ahead means
  ## how long we must keep predict after data. These are of course
  ## predictions of different step lengths.

  ## Y has to be columns. 
  if(class(Y)=="numeric"){
    dim.Y <- c(length(Y),1)
    Y <- matrix(Y,ncol=1)
  } ## else {
  dim.Y <- dim(Y)
  ##  }
  
  ## Definition of default variables
  ## A and C must be supplied
  nstates <- dim(A)[1]
  
  ## these default values don't make much sense
  if(is.null(Sigma.1)){
    Sigma.1 <- diag(rep(1,nstates))
  }
  if(is.null(Sigma.2)){
    Sigma.2 <- diag(rep(1,dim.Y[2]))
  }

  if(is.null(B)){
    B <- matrix(rep(0,nstates),ncol=1)
  }
  if(is.null(u)){
    u <- matrix(rep(0,dim.Y[1]+n.ahead),ncol=1)
  }
  if(is.null(V0)){
    V0 <- Sigma.1
  }
  if(is.null(Xhat0)){
    Xhat0 <- matrix(rep(0,nstates),ncol=1)
  }
  
  
  ## i stedet for (10.79)
  X.hat <- Xhat0
  ## (10.80)
  Sigma.xx <- V0
  ## (10.78) (8.78)
  Sigma.yy <- C%*%Sigma.xx%*%t(C)+Sigma.2

  ## for saving reconstruction
  X.rec <- array(dim=c(dim.Y[1]+n.ahead,nstates))
  X.pred <- array(dim=c(dim.Y[1]+n.ahead,nstates))

  ## for saving K, Sigmas.
  if(verbose){
      K.out <- array(dim=c(dim(Sigma.xx%*%t(C)%*%solve(Sigma.yy)),dim.Y[1]))
      Sigma.xx.rec <- array(dim=c(dim(Sigma.xx),dim.Y[1]))
      Sigma.yy.rec <- array(dim=c(dim(Sigma.yy),dim.Y[1]))   
      Sigma.xx.pred <- array(dim=c(dim(Sigma.xx),dim.Y[1]+n.ahead))
      Sigma.yy.pred <- array(dim=c(dim(Sigma.yy),dim.Y[1]+n.ahead))

  }


  for(tt in (skip+1):dim.Y[1]){
    ## (10.75) (8.75)
    K <- Sigma.xx%*%t(C)%*%solve(Sigma.yy)
    
    ## (10.73) (8.73) - reconstruction
    if(!any(is.na(Y[tt,]))){ # At first everything is thrown away if one is missing
    X.hat <- X.hat+K%*%(t(Y[tt,])-C %*% as.matrix(X.hat))
    X.rec[tt,] <- X.hat
    ## (10.74) (8.74)
    Sigma.xx <- Sigma.xx-K%*%C%*%Sigma.xx
    }
    
    if(verbose){
        Sigma.xx.rec[,,tt] <- Sigma.xx
        Sigma.yy.rec[,,tt] <- Sigma.yy
    }
    
    ##(10.76) (8.76) - prediction
    X.hat <- A%*%X.hat + B%*%t(matrix(as.numeric(u[tt,]),nrow=1))
    X.pred[tt+1,] <- X.hat
    
    ##(10.77) (8.77)
    Sigma.xx <- A%*%Sigma.xx%*%t(A)+Sigma.1
    ##(10.78) (8.78)
    Sigma.yy <- C%*%Sigma.xx%*%t(C)+Sigma.2

    if(verbose){
        K.out[,,tt] <- K
#### these are the prediction error variance-covariances
        Sigma.xx.pred[,,tt+1] <- Sigma.xx
        Sigma.yy.pred[,,tt+1] <- Sigma.yy
    }

  }

if(n.ahead>1){
    for(tt in dim.Y[1]+(1:(n.ahead-1))){
      X.hat <- A%*%X.hat + B%*%t(matrix(u[tt,],nrow=1))
      X.pred[tt+1,] <- X.hat
      Sigma.xx <- A%*%Sigma.xx%*%t(A)+Sigma.1
      Sigma.xx.pred[,,tt+1] <- Sigma.xx
      Sigma.yy.pred[,,tt+1] <- C%*%Sigma.xx%*%t(C)+Sigma.2
    }
  }
  if(verbose){
      out <- list(rec=X.rec,pred=X.pred,K=K.out,Sigma.xx.rec=Sigma.xx.rec,Sigma.yy.rec=Sigma.yy.rec,Sigma.xx.pred=Sigma.xx.pred,Sigma.yy.pred=Sigma.yy.pred)
  } else {
      out <- list(rec=X.rec,pred=X.pred)
  }
  return(out)
}




kalman_outliers2 <- function(Y,A,B=NULL,u=NULL,C,Sigma.1=NULL,Sigma.2=NULL,debug=FALSE,V0=Sigma.1,Xhat0=NULL,n.ahead=1,skip=0,verbose=FALSE){

  ## predictions through data are one-step predictions. n.ahead means
  ## how long we must keep predict after data. These are of course
  ## predictions of different step lengths.

  ## Y has to be columns. 
  if(class(Y)=="numeric"){
    dim.Y <- c(length(Y),1)
    Y <- matrix(Y,ncol=1)
  } ## else {
  dim.Y <- dim(Y)
  ##  }
  
  ## Definition of default variables
  ## A and C must be supplied
  nstates <- dim(A)[1]
  
  ## these default values don't make much sense
  if(is.null(Sigma.1)){
    Sigma.1 <- diag(rep(1,nstates))
  }
  if(is.null(Sigma.2)){
    Sigma.2 <- diag(rep(1,dim.Y[2]))
  }

  if(is.null(B)){
    B <- matrix(rep(0,nstates),ncol=1)
  }
  if(is.null(u)){
    u <- matrix(rep(0,dim.Y[1]+n.ahead),ncol=1)
  }
  if(is.null(V0)){
    V0 <- Sigma.1
  }
  if(is.null(Xhat0)){
    Xhat0 <- matrix(rep(0,nstates),ncol=1)
  }
  
  
  ## i stedet for (10.79)
  X.hat <- Xhat0
  ## (10.80)
  Sigma.xx <- V0
  ## (10.78) (8.78)
  Sigma.yy <- C%*%Sigma.xx%*%t(C)+Sigma.2

  ## for saving reconstruction
  X.rec <- array(dim=c(dim.Y[1]+n.ahead,nstates))
  X.pred <- array(dim=c(dim.Y[1]+n.ahead,nstates))

  ## for saving K, Sigmas.
  if(verbose){
      K.out <- array(dim=c(dim(Sigma.xx%*%t(C)%*%solve(Sigma.yy)),dim.Y[1]))
      Sigma.xx.rec <- array(dim=c(dim(Sigma.xx),dim.Y[1]))
      Sigma.yy.rec <- array(dim=c(dim(Sigma.yy),dim.Y[1]))   
      Sigma.xx.pred <- array(dim=c(dim(Sigma.xx),dim.Y[1]+n.ahead))
      Sigma.yy.pred <- array(dim=c(dim(Sigma.yy),dim.Y[1]+n.ahead))

  }

  
  noutsider <- 0
  outsider_index <- NA


  for(tt in (skip+1):dim.Y[1]){
    ## (10.75) (8.75)
    K <- Sigma.xx%*%t(C)%*%solve(Sigma.yy)

    ########################################################
    # compute predictions for observations checking if obs is more than six standard deviations away from the prediction

    ## (10.73) (8.73) - reconstruction
    if(!any(is.na(Y[tt,]))){ # At first everything is thrown away if one is missing
    X.hat.test <- X.hat+K%*%(t(Y[tt,])-C %*% as.matrix(X.hat))
    ## (10.74) (8.74)
    Sigma.xx.test <- Sigma.xx-K%*%C%*%Sigma.xx
    }

    ##(10.76) (8.76) - prediction
    X.hat.test <- A%*%X.hat.test + B%*%t(matrix(as.numeric(u[tt,]),nrow=1))
    ##(10.77) (8.77)
    Sigma.xx.test <- A%*%Sigma.xx.test%*%t(A)+Sigma.1
    ##(10.78) (8.78)
    Sigma.yy.test <- C%*%Sigma.xx.test%*%t(C)+Sigma.2

    ###################################################

    # If obs is not NA AND error is too large
    if (tt<length(Y) && !any(is.na(Y[tt+1,])) && abs(X.hat.test-Y[tt+1]) > 6*sqrt(Sigma.yy.test)) {
        #print(X.hat.test-Y[tt])
        #print(6*sqrt(Sigma.yy.test))
        Y[tt+1] <- NA # put obs to NA
        noutsider <- noutsider+1 # count number of outsider
        outsider_index <- c(outsider_index,tt+1) # store index
    }

    ###################################################

    # run loop as usual
    
    ## (10.73) (8.73) - reconstruction
    if(!any(is.na(Y[tt,]))){ # At first everything is thrown away if one is missing
    X.hat <- X.hat+K%*%(t(Y[tt,])-C %*% as.matrix(X.hat))
    X.rec[tt,] <- X.hat
    ## (10.74) (8.74)
    Sigma.xx <- Sigma.xx-K%*%C%*%Sigma.xx
    }
    
    if(verbose){
        Sigma.xx.rec[,,tt] <- Sigma.xx
        Sigma.yy.rec[,,tt] <- Sigma.yy
    }
    
    ##(10.76) (8.76) - prediction
    X.hat <- A%*%X.hat + B%*%t(matrix(as.numeric(u[tt,]),nrow=1))
    X.pred[tt+1,] <- X.hat
    
    ##(10.77) (8.77)
    Sigma.xx <- A%*%Sigma.xx%*%t(A)+Sigma.1
    ##(10.78) (8.78)
    Sigma.yy <- C%*%Sigma.xx%*%t(C)+Sigma.2

    if(verbose){
        K.out[,,tt] <- K
#### these are the prediction error variance-covariances
        Sigma.xx.pred[,,tt+1] <- Sigma.xx
        Sigma.yy.pred[,,tt+1] <- Sigma.yy
    }

  }

if(n.ahead>1){
    for(tt in dim.Y[1]+(1:(n.ahead-1))){
      X.hat <- A%*%X.hat + B%*%t(matrix(u[tt,],nrow=1))
      X.pred[tt+1,] <- X.hat
      Sigma.xx <- A%*%Sigma.xx%*%t(A)+Sigma.1
      Sigma.xx.pred[,,tt+1] <- Sigma.xx
      Sigma.yy.pred[,,tt+1] <- C%*%Sigma.xx%*%t(C)+Sigma.2
    }
  }
  if(verbose){
      out <- list(rec=X.rec,pred=X.pred,K=K.out,Sigma.xx.rec=Sigma.xx.rec,Sigma.yy.rec=Sigma.yy.rec,Sigma.xx.pred=Sigma.xx.pred,Sigma.yy.pred=Sigma.yy.pred,NB.outsider=noutsider,Outsider.index=outsider_index[2:length(outsider_index)])
  } else {
      out <- list(rec=X.rec,pred=X.pred)
  }
  return(out)
}
