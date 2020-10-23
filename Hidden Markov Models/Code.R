library("HMM")
library("entropy")
library("gridExtra")
library("ggplot2")

# Q1---------------------------------------------------------
# Build a HMM 
states    <- c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z9", "Z10")
emission  <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")
nstates   <- length(states)

# Transition probabilities
transProb <- matrix(data=0, nrow=nstates, ncol=nstates,
                    dimnames=list(states, states))

for (i in 1:(nstates-1)) {
  # Assign equal probability to current and next state
  transProb[i, c(i, i+1)] <- 0.5    
}
transProb[nstates, c(nstates, 1)] <- 0.5

# Emission probabilities
emissionprob <- matrix(data=0, nrow=nstates, ncol=nstates,
                       dimnames=list(states, emission))

for (i in 1:nstates) {
  # Chain the emissions starting with the emission Xi
  emission_chain <- c(emission[i:nstates], emission[-(i:nstates)])
  
  # Assign equal probability to the current, previous two and next two emissions in the chain
  emissionprob[i, emission_chain[c(9:10, 0:3)]] <- 0.2
}

# Initialize a general discrete time HMM
hmmInit <- initHMM(States=states, 
                   Symbols=emission, 
                   startProbs=c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), 
                   transProbs=transProb, 
                   emissionProbs=emissionprob)

# Q2---------------------------------------------------------
# Simulate the HMM for 100 time steps
T <- 100
simHmm <- simHMM(hmm=hmmInit, length=T)

# Q3---------------------------------------------------------
# Compute forward and backward probabilities
alpha  <- exp(forward(hmmInit, simHmm$observation))
beta   <- exp(backward(hmmInit, simHmm$observation))

# Compute the filtered and smoothed distributions
filtered  <- matrix(data=0, nrow=nstates, ncol=T)  
smoothed  <- matrix(data=0, nrow=nstates, ncol=T)  
for (i in 1:T) {
  filtered[, i] <- alpha[, i]/sum(alpha[, i])
  smoothed[, i] <- alpha[, i]*beta[, i] / sum(alpha[, i]*beta[, i])
}

# Plot the filtered and smoothed distributions for each of the data points
filter_plot <- 
  lapply(X = seq_along(1:T), FUN = function(i) {
    ggplot(as.data.frame(filtered[, i]), aes(1:nstates, filtered[, i])) +
      geom_line() +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),           
            axis.ticks.x=element_line(),
            axis.ticks.y=element_line()) +
      scale_x_continuous(n.breaks=10) + 
      scale_y_continuous(n.breaks=10)   
  })
do.call(grid.arrange, filter_plot)

smooth_plot <- 
  lapply(X = seq_along(1:T), FUN = function(i) {
    ggplot(as.data.frame(smoothed[, i]), aes(1:nstates, smoothed[, i])) +
      geom_line() + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.x=element_line(),
            axis.ticks.y=element_line()) +
      scale_x_continuous(n.breaks=10) + 
      scale_y_continuous(n.breaks=10)
  })
do.call(grid.arrange, smooth_plot)

# Q4---------------------------------------------------------
# Forward-Backward algorithm
FB = function(hmmInit, observations) {
  
  # Compute forward and backward probabilities
  alpha  <- exp(forward(hmmInit, observations))
  beta   <- exp(backward(hmmInit, observations))
  
  # Filtered and Smoothed Distributions
  filtered <- prop.table(x=alpha, margin=2)
  smoothed <- prop.table(x=alpha*beta,  margin=2)
  
  return(list(filtered=filtered, smoothed=smoothed))
} 

# Compare accuracy of different methods
compare_accuracy = function(hmmInit, observations, true_states) {
  
  # Forward-backward algorithm to compute the alpha and beta and the
  # filtering and smoothing distributions
  fb <- FB(hmmInit, observations)
  
  # Most Probable States 
  filter_MPS <- hmmInit$States[apply(X=fb$filtered, MARGIN=2, FUN=which.max)]
  smooth_MPS <- hmmInit$States[apply(X=fb$smoothed, MARGIN=2, FUN=which.max)]
  
  # Most Probable Path
  MPP <- viterbi(hmmInit, observations)
  
  # Confusion matrices
  filter_t <- table(true_states, filter_MPS)
  smooth_t <- table(true_states, smooth_MPS)
  MPP_t    <- table(true_states, MPP) 
  
  # Accuracy
  filter_acc <- sum(diag(filter_t)) / sum(filter_t)
  smooth_acc <- sum(diag(smooth_t)) / sum(smooth_t)
  MPP_acc    <- sum(diag(MPP_t)) / sum(MPP_t)
  
  results <- data.frame(Accuracy=c(filter_acc, smooth_acc, MPP_acc),
                        row.names=c("Filtered", "Smoothed", "MPP"))
  return(results)
}

# Compare the accuracy of most probable states obtained from filtering and smoothing
# distributions with the most probable path 
results <- compare_accuracy(hmmInit, simHmm$observation, simHmm$states)
print(results)

# Q5---------------------------------------------------------
# Simulate the HMM for 150, 200, 250, 300, 350 time steps and compare accuracy 
# for different methods
results <- data.frame(matrix(ncol=5, nrow=3), 
                      row.names=c("Filtered", "Smoothed", "MPP"))
names(results) <- c("T=150", "T=200", "T=250", "T=300", "T=350")

i <- 1
for (T in c(150, 200, 250, 300, 350)) {
  sample  <- simHMM(hmm=hmmInit, length=T)
  results[, i] <- compare_accuracy(hmmInit, sample$observation, sample$states)
  i <- i+1
}
print(results)

# Q6---------------------------------------------------------
# Entropy
fb <- FB(hmmInit, simHmm$observation)
entropyFilter <- apply(X=fb$filtered, MARGIN=2, FUN=entropy.empirical)
graphics::plot(entropyFilter, type="l")

# Q7---------------------------------------------------------
# Prediction
postState <- posterior(hmmInit, simHmm$observation)

z_101 <- postState[, 100] %*% transProb
print(z_101)
