library(kernlab)
library(AtmRay)

############################################################
# 1. Implementing GP Regression
############################################################

# Covariance function
SquaredExpKernel <- function(x1, x2, sigmaF=1, l=1){
  n1 <- length(x1)
  n2 <- length(x2)
  K  <- matrix(NA,n1,n2)
  
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}


# Posterior GP of f
posteriorGP <- function(x, y, xstar, sigmaNoise, k, ...) {
  n <- length(x)
  
  kss <- k(x1=xstar, x2=xstar, ...)
  kxx <- k(x1=x, x2=x, ...)
  kxs <- k(x1=x, x2=xstar, ...)
  
  # Cholesky decomposition of prior marginal variance of noisy training points
  L <- t(chol(kxx + sigmaNoise^2*diag(n)))
  
  # Mean vector of f*  
  alpha <- solve(t(L), solve(L, y))
  meanVec <- t(kxs) %*% alpha
  
  # Covariance matrix of f*  
  v <- solve(L, kxs)
  covMat <- kss - (t(v) %*% v)
  
  return(list(meanVec=meanVec, covMat=covMat))
}


# Plot mean and 95% Probability bands for f*  
plotGPPosterior <- function(xstar, meanFstar, covFstar) {
  plot(xstar, meanFstar, ylim=c(-3, 3), col="red", lwd=2)
  lines(xstar, meanFstar + 1.96*sqrt(diag(covFstar)), col="blue", lwd=2)
  lines(xstar, meanFstar - 1.96*sqrt(diag(covFstar)), col="blue", lwd=2)
}

xgrid <- seq(-1, 1, length=20)
GP <- posteriorGP(x=0.4, y=0.719, 
                  xstar=xgrid, sigmaNoise=0.1, k=SquaredExpKernel,
                  sigmaF=1, l=0.3)
plotGPPosterior(xgrid, GP$meanVec, GP$covMat)


# Update posterior with (x,y)=(-0.6,-0.044)
GP <- posteriorGP(x=c(0.4, -0.6), y=c(0.719, -0.044), 
                  xstar=xgrid, sigmaNoise=0.1, k=SquaredExpKernel,
                  sigmaF=1, l=0.3)
plotGPPosterior(xgrid, GP$meanVec, GP$covMat)


# Compute posterior for five data points with l=0.3
GP <- posteriorGP(x=c(-1.0, -0.6, -0.2, 0.4, 0.8), 
                  y=c(0.768, -0.044, -0.940, 0.719, -0.664), 
                  xstar=xgrid, sigmaNoise=0.1, k=SquaredExpKernel,
                  sigmaF=1, l=0.3)
plotGPPosterior(xgrid, GP$meanVec, GP$covMat)


# Compute posterior for five data points with l=1
GP <- posteriorGP(x=c(-1.0, -0.6, -0.2, 0.4, 0.8), 
                  y=c(0.768, -0.044, -0.940, 0.719, -0.664), 
                  xstar=xgrid, sigmaNoise=0.1, k=SquaredExpKernel,
                  sigmaF=1, l=1)
plotGPPosterior(xgrid, GP$meanVec, GP$covMat)


############################################################
# 2. GP Regression with kernlab
############################################################

# Read data and add variables time and day
SthlmTemp <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/
Code/TempTullinge.csv", header=TRUE, sep=";")
SthlmTemp$time <- 1:nrow(SthlmTemp)
SthlmTemp$day  <- rep(c(1:365), 6)

# Subsample the data- use only every fifth observation
data <- SthlmTemp[seq(1, nrow(SthlmTemp), 5), ]


# Square expontial kernel function 
SEkernel <- function(sigmaf=1, ell=1) 
{
  rval <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    K  <- matrix(NA, n1, n2)
    
    for (i in 1:n2){
      K[, i] <- sigmaf^2*exp(-0.5*( (x1-x2[i])/ell)^2 )
    }
    return(K)
  }
  class(rval) <- "kernel"
  return(rval)
} 
k = SEkernel(sigmaf=1, ell=0.5)


# Evaluate in (x,x')=(1,2)
k(x1=1, x2=2)
# Covariance matrix K(X, X*)
kxs <- kernelMatrix(kernel=k, x=c(1,3,4), y=c(2,3,4))
print(kxs)


# Estimate traning noise using simple quadratic regression on time variable
quadFit <- lm(temp ~ time + I(time^2), data=data)
sigmaNoise <- sd(quadFit$residuals)


# GP regression model
GPfit <- gausspr(x=data$time, y=data$temp, 
                 kernel=SEkernel, kpar=list(sigmaf=20, ell=0.2), 
                 var=sigmaNoise^2)
# Predicting the training data
postMean <- predict(GPfit, data$time) 

# Plot data and mean predictions
plot(data$time, data$temp, lwd=2)
lines(data$time, postMean, col="red", lwd=2)


# 95% prediction interval implementation
postGP <- posteriorGP(x=scale(data$time), y=data$temp, 
                      xstar=scale(data$time), sigmaNoise=sigmaNoise^2, k=SquaredExpKernel, 
                      sigmaF=20, l=0.2)
# Plot 95% prediction bands
lines(data$time, postMean + 1.96*sqrt(diag(postGP$covMat)), col="blue", lwd=2)
lines(data$time, postMean - 1.96*sqrt(diag(postGP$covMat)), col="blue", lwd=2)



# Estimate traning noise using simple quadratic regression on day variable
quadFit <- lm(temp ~ day + I(day^2), data=data)
sigmaNoise <- sd(quadFit$residuals)


# GP regression model
k = SEkernel(sigmaf=20, ell=0.2)
GPfit <- gausspr(x=data$day, y=data$temp, kernel=k, var=sigmaNoise^2)
# Predicting the training data
postMean2 <- predict(GPfit, data$day) 

# Plot data and mean predictions from both GPs
plot(data$time, data$temp, lwd=2)
lines(data$time, postMean, col="red", lwd=2)
lines(data$time, postMean2, col="blue", lwd=2)


# 95% prediction interval implementation
postGP <- posteriorGP(x=scale(data$day), y=data$temp, 
                      xstar=scale(data$day), sigmaNoise=sigmaNoise^2, k=SquaredExpKernel, 
                      sigmaF=20, l=0.2)
# Plot 95% prediction bands
lines(data$time, postMean + 1.96*sqrt(diag(postGP$covMat)), col="green", lwd=2)
lines(data$time, postMean - 1.96*sqrt(diag(postGP$covMat)), col="green", lwd=2)


# Generalized periodic kernel 
PeriodicKern <- function(sigmaf=1, l1=1, l2=1, d) 
{
  rval <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    K  <- matrix(NA, n1, n2)
    
    for (i in 1:n2){
      K[, i] <- sigmaf^2 * exp(-2*(sin(pi*abs(x1-x2[i])/d)/l1)^2) * exp(-0.5*((x1-x2[i])/l2)^2)
    }
    return(K)
  }
  class(rval) <- "kernel"
  return(rval)
} 

# GP regression model
GPfit <- gausspr(x=data$time, y=data$temp, 
                 kernel=PeriodicKern, kpar=list(sigmaf=20, l1=1, l2=10, d=365/sd(data$time)), 
                 var=sigmaNoise^2)
# Predicting the training data
postMean3 <- predict(GPfit, data$time) 

# Plot data and mean predictions from all GPs
plot(data$time, data$temp, lwd=2)
lines(data$time, postMean, col="red", lwd=2)
lines(data$time, postMean2, col="blue", lwd=2)
lines(data$time, postMean3, col="green", lwd=2)

# 95% prediction interval implementation
k = PeriodicKern(sigmaf=20, l1=1, l2=10, d=365/sd(data$time))
postGP <- posteriorGP(x=scale(data$time), y=data$temp, 
                      xstar=scale(data$time), sigmaNoise=sigmaNoise^2, k=k)
# Plot 95% prediction bands
lines(data$time, postMean + 1.96*sqrt(diag(postGP$covMat)), col="orange", lwd=2)
lines(data$time, postMean - 1.96*sqrt(diag(postGP$covMat)), col="orange", lwd=2)


############################################################
# 3. GP Classification with kernlab
############################################################

# Read the data
data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/
GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])

# Train-test split
set.seed(111) 
SelectTraining <- sample(1:dim(data)[1], size=1000, replace=FALSE)
train <- data[SelectTraining, ]
test  <- data[-SelectTraining, ]


# GP classification model using 2 covariates
GPfit     <- gausspr(fraud ~ varWave + skewWave, data=train)
fraudPred <- predict(GPfit, train[, 1:2])

# Confusion matrix and accuracy of the classifier
confMat <- table(fraudPred, train$fraud)
acc     <- sum(diag(confMat)) / sum(confMat)



# Plot contours of the prediction probabilities
x1 <- seq(min(train$varWave), max(train$skewWave), length=100)
x2 <- seq(min(train$skewWave), max(train$skewWave), length=100)

gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind.data.frame(c(gridPoints$x), c(gridPoints$y))
names(gridPoints) <- names(train)[1:2]

fraudProb <- predict(GPfit, gridPoints, type="probabilities")

# Plotting for Prob(fraud = 0)
contour(x1, x2, matrix(fraudProb[, 1], 100, byrow=TRUE), 20, 
        xlab="varWave", ylab="skewWave", main='Prob(fraud = 0)')
points(train[train$fraud==1, 1], train[train$fraud==1, 2], col="blue")
points(train[train$fraud==0, 1], train[train$fraud==0, 2], col="red")

# Plotting for Prob(fraud = 1)
contour(x1, x2, matrix(fraudProb[, 2], 100, byrow=TRUE), 20, 
        xlab="varWave", ylab="skewWave", main='Prob(fraud = 1)')
points(train[train$fraud==1, 1], train[train$fraud==1, 2], col="blue")
points(train[train$fraud==0, 1], train[train$fraud==0, 2], col="red")



# Predict on the test set
fraudPred <- predict(GPfit, test[, 1:2])
# Confusion matrix and accuracy of the classifier
confMat <- table(fraudPred, test$fraud)
acc     <- sum(diag(confMat)) / sum(confMat)



# GP classification model using all 4 covariates
GPfit     <- gausspr(fraud ~ ., data=train)
fraudPred <- predict(GPfit, test[, -5])

# Confusion matrix and accuracy of the classifier
confMat <- table(fraudPred, test$fraud)
acc     <- sum(diag(confMat)) / sum(confMat)
