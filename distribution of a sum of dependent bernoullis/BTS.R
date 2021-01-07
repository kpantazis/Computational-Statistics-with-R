#Compute the distribution of a sum of dependent bernoulli random variables with same probability using Markov property
BTS <- function(n,lambda,p){
  #Create transition probabilities
  alpha <- matrix(0,n+1,n)
  beta <- matrix(0,n+1,n)
  
  
  #n == 1 & restriction on lambda values
  if ((2*p-1)/p <= lambda & lambda <= 1){
  p11 <- lambda
  p01 <- 1 - lambda
  p00 <- ( 1 - 2*p + lambda*p )/ (1-p)
  p10 <- 1-p00
  alpha[1,1] = 1-p
  alpha[2,1] = 0
  beta[1,1] = 0
  beta[2,1] = p
  
  #Transition probabilities for n>=2
  if (n > 1){
  for (l in 2:n){
  alpha[1,l] = (1-p)*p00^(l-1)
  beta[1,l] = 0
  
  for (k in 1:n) {
   if ( k == n) {
     alpha[k+1,k] = 0
     beta[k+1,k] = p*p11^(k-1)
   } else {
     alpha[k+1,l] = p00*alpha[k+1,l-1]+p01*beta[k+1,l-1]
     beta[k+1,l] = p10*alpha[k,l-1]+p11*beta[k,l-1]
   }
  }
  }
}
  #print(alpha[,n])
  #print(beta[,n])
  
  
  #Distribution of the Bernoulli distribution
  pii <- rep(0,n+1)
  pii[1] <- alpha[1,n]
  for (k in 1:n+1){
    pii[k] = alpha[k,n]+beta[k,n]
  }
  return(list(pii))
  }else {
    
    #Gives error when p and lambda values are inappropriate
    library(crayon)
    #print("Error: Choose appropriate p and lambda values")
    cat(red("Error"), ": Choose appropriate", blue("p")," and ",blue("lambda")," values.\n")
  }
}
  