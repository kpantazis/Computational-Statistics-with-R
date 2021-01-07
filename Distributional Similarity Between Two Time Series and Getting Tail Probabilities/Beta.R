#Beta distribution

SP2BBSQK <- function(x1,x2,Increment,BandWidth,T){
  ### m=2,h(x)=(log(x),log(1-x). K means Kernel Density Estimate.
  ### Increment controls the grid at which est. g,g1,G,G1 are evaluated
  ### BandWidth controls smoothness of kernel est. of g,g1,
  ### Reasonable values are Increment=0.05, BandWidth=0.3,0.5
  ### T is an upper probability point: We get 1-G1(T), 1-G(T)
  ### G,g reference; G1,g1 tilted distributions.
  ### We get two tests of equi-distribution: Chi1 and LR (Likelihood Ratio).
  ### We get plots of g,g1,G,G1. The equi-distribution tests match the plots.
  
  
  ###The data: m=2, x2 is the reference sample. 
  n1 <- length(x1)
  n2 <- length(x2)
  rho <- n1/n2
  n <- n1+n2
  level <- 0.05
  t <- c(x1,x2) #Data fusion.
  
  ###MINUS log-likelihood 
  minusloglike <- function(theta) {
    sum(log(1+rho*exp(theta[1] + theta[2]*log(t) + theta[3]*log(1-t))))-
      sum(theta[1]+theta[2]*log(x1)+theta[3]*log(1-x1))}
  ###Maximizing loglikelihood by minimizing MINUS loglikelihood
  min.func <- nlminb( start=c(-0.02,.2,.2),obj = minusloglike)
  
  ###Parameter estimates
  delta <- min.func$par[1]
  epsilon <- min.func$par[2]
  zeta <- min.func$par[3]
  
  ###Reference dist. p=dG and its distortion p1=w*dG
  p <- 1/(n2*(1+rho*exp(delta + epsilon*log(t) + zeta*log(1-t))))
  p1 <- exp(delta + epsilon*log(t) + zeta*log(1-t))*p
  
  #png(file="gEstimate.png")
  ###Plots of est. G,G1,g,g1
  
  par(mfrow=c(2,2), oma=c(0,0,4,0))
  
  ###Estimated ref. cdf G(x):
  G <- function(x){sum(p[t <= x])}
  U <- 1-G(T)
  
  ###To Plot G over the range (min(t),max(t)):
  #Modify Increment as needed. 0.05 reasonable
  x <- seq(min(t),max(t),Increment) 
  
  #Creating a vector out of G for plotting.
  cdf <- length(x)
  for(i in 1:length(x)){
    cdf[i] <- G(x[i])}
  
  ###The distorted cdf of x1:
  G1 <- function(x){sum(p1[t <= x])}
  U1 <- 1-G1(T)
  
  
  cdf1 <- length(x)
  for(i in 1:length(x)){
    cdf1[i] <- G1(x[i])}
  
  plot(x,cdf, type="l")
  #plot(x,cdf, type="l", xlab="x") 
  #NOTE: "xlab" could not be set in high-level plot() function
  lines(x,cdf1, type="l", lty=2)
  title("Estimated G, G1")
  
  
  ###Kernel Density Estimates of g(x) and g1(x).
  ghat <- function(x){
    K <- function(x){(1/sqrt(2*pi))*exp(-x^2/2)}
    sum(p*K(( x-t)/BandWidth)/BandWidth)}
  
  g1hat <- function(x){
    K <- function(x){(1/sqrt(2*pi))*exp(-x^2/2)}
    sum(p1*K(( x-t)/BandWidth)/BandWidth)}
  
  
  ###Plots of ghat and g1hat
  #For plots: Creating vectors from ghat and Distortted g1hat
  x <- seq(min(t),max(t),Increment)
  g <- rep(0,length(x))
  g1 <- rep(0,length(x))
  for(i in 1:length(x)){
    g[i] <- ghat(x[i])
    g1[i] <- g1hat(x[i])} ### Estimated pdf.
  
  #plots
  maxgg1 <- max(c(g,g1))
  plot(x,g, type="l",ylim=c(0,maxgg1))
  #plot(x,g, type="l",ylim=c(0,maxgg1),xlab="x")
  lines(x,g1, type="l",lty=2)
  title("Kernel Est g, g1")
  
  ###Hist of Reference Data vs Estimated Reference pdf
  maxghist <- max(hist(x2,plot=FALSE)$density,g)
  hist(x2,prob=T, xlab="x",ylim=c(0,maxghist),main="Ref Hist & Est g")
  lines(x,g, type="l") 
  
  ###Hist of Distorted Data vs Estimated g1 pdf
  maxg1hist <- max(hist(x1,plot=FALSE)$density,g1)
  hist(x1,prob=T, xlab="x",ylim=c(0,maxg1hist),main="Dist Hist & Est g1")
  lines(x,g1, type="l") 
  
  
  
  #mtext("SP2XXSQ, m=2, h(x)=(x,x^2)", cex=1.2, line=1,side=3,outer=TRUE) 
  mtext("m=2, h(x)=(x,x^2)", cex=1.2, line=1,side=3,outer=TRUE)
  
  #dev.off() #Need dev.off for png(file="gEstimate.png")     
  
  ###The LR Test gamma=0
  logL <- -sum(log(1+rho*exp(delta + epsilon*log(t) + zeta*log(1-t)))) +
    sum(delta + epsilon*log(x1) + zeta*log(1-x1)) 
  
  #MINUS log-likelihood under H_0: gamma=0
  minusloglike0 <- function(theta) {
    sum(log(1+rho*exp(theta[1] + theta[2]*log(t) + 0*log(1-t)))) -
      sum(theta[1]+theta[2]*log(x1)+0*log(1-x1))}
  
  #Maximizing loglikelihood under H_0 by minimizing minusloglike0
  min.func0 <- nlminb( start=c(0,0),obj = minusloglike0)
  
  #Parameter estimates under H_0: gamma=0
  delta0 <- min.func0$par[1]
  epsilon0 <- min.func0$par[2]
  
  #Max log-likelihood under H_0: gamma=0
  logL0 <- -sum(log(1+rho*exp(delta0 + epsilon0*log(t) + 0*log(1-t)))) +
    sum(delta0 + epsilon0*log(x1) + 0*log(1-x1)) 
  
  #Liklihood Ratio Test Statistic for H_0: gamma=0
  LR0 <- -2*(logL0 - logL)
  
  ###The Chi1 Test of Equidistribution (beta,gamma)=(0,0):
  tm1 <- sum(t*p); tm2 <- sum(p*(t)^2); tm3 <- sum(p*(t)^3); tm4 <- sum(p*(t)^4)
  vart <- tm2-(tm1)^2
  covtt2 <- tm3-tm1*tm2
  vart2 <- tm4-(tm2)^2
  VARh <- matrix(c(vart,covtt2,covtt2,vart2),ncol=2)
  A11 <- rho/((1+rho)^2)
  BETA <- c(epsilon,zeta)
  chi1 <- n*A11*BETA%*%VARh%*%BETA
  
  ###The LR Test of Equidistribution (beta,gamma)=(0,0)
  logLR00 <- -n*log(1+rho)
  LR <- -2*(logLR00 - logL)
  
  #LR <- -2*sum(log(1+rho*exp(alpha + beta*(t) + gamma*(t)^2))) +
  #2*sum(alpha + beta*(x1) + gamma*(x1)^2) + 2*n*log(1+rho)
  
  ###Check Sum p, Sum p1
  Sum_p <- sum(p); Sum_p1 <- sum(p1)
  
  ###Check: If integral relation holds
  ALPHA <- -log(sum(exp(epsilon*log(t) + zeta*log(1-t))*p))
  
  list(delta=delta,epsilon=epsilon,zeta=zeta,
       pval_LRzeta = 1-pchisq(LR0,1),
       pval_chi1=1-pchisq(chi1,2), 
       pval_LR=1-pchisq(LR,2),
       Sum_p=Sum_p, Sum_p1=Sum_p1, ALPHA=ALPHA,Upper_Threshold=T,Upper_G1=U1,Upper_G=U)
}

###############################################

#Example with not same distributions
x1 <- rbeta(1000,0.1,0.2)
x2 <- rbeta(1000,.2,.2)
SP2BBSQK(x1,x2,0.05,0.5,1.645)

#Example with almost the same distributions
x1 <- rbeta(1000,.1,.2)
x2 <- rbeta(1000,.1,.2)
SP2BBSQK(x1,x2,0.05,0.5,1.645)


#Note the pval_LR value; in the first case, the pval_LR is 0 rejecting the null hypothesis H_0 (pval_LR<0.05), 
#while, in the second case, the pval_LR is sufficiently large to conclude that the null hypothesis cannot be rejected.
