library("bnlearn")
library("gRain")

# Q1---------------------------------------------------------
data("asia")
n = nrow(asia) # number of samples

# Multiple runs of Hill Climbing to learn the BN structure
bn1  <- list()
comp <- c()


for (i in 1:10) {
  bn1[[i]] <- hc(x=asia, score="bde", iss=100, restart=10)
  
  if(i!=1)  
    comp[i] <- all.equal(bn1[[i-1]], bn1[[i]])
}

# Comparison of different BN structures learnt
print(comp)
par(mfcol = c(2,5))
do.call(graphviz.compare, bn1)

# Plot one of the learned BN structure
dev.off()
print(bn1[[1]])
plot(bn1[[1]])

bn_arc <- arcs(bn1[[1]])
print(bn_arc)
bn_vstruct <- vstructs(bn1[[1]])
print(bn_vstruct)

# Q2---------------------------------------------------------
# 80-20 split of data for train-test
set.seed(12345)
index <- sample(1:n, n*0.8)
asia_train <- asia[index, ]
asia_test  <- asia[-index, ]

# Function to fit the parameters of BN and predict S based on the structure of BN
predict_S <- function(bn, evidence_nodes, predict_node, data_train=asia_train, 
                      data_test=asia_test) {
  
  # Fit the parameters of the BN conditional on its structure
  bn_fit <- bn.fit(bn, data=data_train, method="mle")
  #print(bn_fit)
  
  # Convert fit object to grain object
  bn_grain <- as.grain(bn_fit)
  #print(bn_grain)
  #plot(bn_grain)
  
  # Prediction of S
  n_test <- nrow(data_test)  # Number of test samples
  S_predict <- character(n_test) # Vector to hold predictions of S
  
  for (i in 1:n_test) {
    # Set the evidence i.e. the state of all other nodes except S 
    evidence = setEvidence(bn_grain, 
                           nodes=evidence_nodes, 
                           states=sapply(data_test[i, evidence_nodes], toString))
    
    # Conditional distribution of S given all the other nodes in the evidence
    S_cond = querygrain(evidence, nodes=predict_node, type="conditional")
    
    # Classify S as yes or no depending on the conditional posterior of S
    if (S_cond[1] > S_cond[2]) {
      S_predict[i] = names(S_cond)[1] 
    } else {
      S_predict[i] = names(S_cond)[2] 
    }
  }
  
  # Compute the confusion matrix and misclassifcation rate
  conf_mat = table(S_predict, data_test[, predict_node])
  mc_rate = 1 - (sum(diag(conf_mat)) / sum(conf_mat))
  
  return(list(confusion_matrix=conf_mat, mc_rate=mc_rate))
}

# Learn the BN structure from train dataset
bn2 <- hc(x=asia_train)
print(bn2)
plot(bn2)

# Predict S based on the learned BN structure
S_predict = predict_S(bn=bn2,                         # Learned BN
                      evidence_nodes=names(asia)[-2], # All nodes except S
                      predict_node=names(asia)[2])    # Node S
print(S_predict)

# True Asia bayesian network
dag <- model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
print(dag)
plot(dag)

# Predict S based on the true BN structure
S_predict_true = predict_S(bn=dag,                         # True BN
                           evidence_nodes=names(asia)[-2], # All nodes except S
                           predict_node=names(asia)[2])    # Node S
print(S_predict_true)

graphviz.compare(bn2, dag)
graphviz.compare(cpdag(bn2), cpdag(dag))
# Q3---------------------------------------------------------
# Predict S given the markov-blanket of S
S_predict_mb = predict_S(bn=bn2,                           # Learned BN
                         evidence_nodes=mb(bn2, node="S"), # markov blanket of S
                         predict_node=names(asia)[2])      # Node S
print(S_predict_mb)

# Q4---------------------------------------------------------
# Naive-bayes classifier
nb = model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
plot(nb)

S_predict_nb = predict_S(bn=nb,                            # Bayes-classifier
                         evidence_nodes=names(asia)[-2],   # All nodes except S
                         predict_node=names(asia)[2])      # Node S
print(S_predict_nb)

