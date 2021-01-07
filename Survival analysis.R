library(Deriv)
library(crayon)

############################################################
#Theory
#x is age/ F_0(t)/(a,b)
#s_0(t) <- 1-F_0(t)
#F_x(t) <- (F_0(x+t)-F_0(x))/S_0(x)
#S_x(t) <- S_0(x+t)/S_0(x)
#Plot F_x(t),S_x(t),f_x(t),mu_x
#Compute E(T_x),Var(T_x),P(a<T_x<b)=F_x(b)-F_x(a)
#x<- seq(1,10)
#f = deriv(~a^2,'a',func = T)
#attr(f(7),'gradient')
#S_0(x+t)=S_x(t)*S_0(x)
#eg 0<=x<=120 , F_0(x) = 1-(1-x/120)^(1/6) , 0 otherwise
#mu_x = 1/(720-6x)
#####################################################################
#################################################################

#x fixed, t varies
#x life aged
# k corresponds to the exponent in F_0 function. For example, in class we had an example with k=6.
k <- seq(1,10)
x <- 20
t <- seq(1,120)
#a <- 0 
#b <- length(t)-x

F_0 <- function(t,k){
  F_0 <- matrix(0,length(k),length(t))
  for(i in 1:length(k)){
    for(j in 1:length(t)){
      F_0[i,j] <- ifelse(t[j]>0,(1-(1-t[j]/length(t))^(1/k[i])),0)}
  }
  return(F_0)}

#F_0(t,k)

F_x <- function(F_0,t,k,x){
  F_x <- matrix(0,length(k),length(t)-x)
  for(i in 1:length(k)){
    for(j in 1:(length(t)-x)){
      F_x[i,j] <- (F_0(t,k)[i,x+j]-F_0(t,k)[i,x])/(1-F_0(t,k)[i,x])}
  }
  return(F_x)}
#  F_x <- rep(0,length(t)-x)
 # for(i in 1:(length(t)-x)){
#    F_x[i] <- (F_0(t)[x+i]-F_0(t)[x])/(1-F_0(t)[x])  }
  #  return(F_x)}

#F_0(t,k)
#F_x(F_0,t,k,x)

#S_x is survival probability
S_0 <- list()
S_x <- list()
for(i in 1:length(k)){
  S_0[[i]] <- matrix(1,1,length(t))-F_0(t,k)[i,]
  #F_x(F_0,x)[2,]
  S_x[[i]] <- matrix(1,1,(length(t)-x))-F_x(F_0,t,k,x)[i,]}
  #S_x


#Checking s_x(t)=S_0(x+t)/S_0(x)//x=20,t=2,k=1
S_x[[2]][1,2] == S_0[[2]][1,22]/S_0[[2]][1,20]



#Derivative aka density f_x(t)
#x<- seq(1,10)
#f = deriv(~a^2,'a',func = T)
#attr(f(7),'gradient')
#f_0

f_0 <- function(t,k){
  f_0 <- matrix(0,length(k),length(t))
  for(i in 1:length(k)){
    for(j in 1:length(t)){
      f_0[i,j] <- ifelse(t[j]>0,(1/(k[i]*length(t)))*(1-t[j]/length(t))^(-(k[i]-1)/k[i]),0)}
  }
  return(f_0)}

#f_0(t,k)

f_x <- function(f_0,F_0,t,k,x){
  f_x <- matrix(0,length(k),length(t)-x)
  for(i in 1:length(k)){
    for(j in 1:length(t)-x){
      f_x[i,j] <- f_0(t,k)[i,x+j]/(1-F_0(t,k)[i,x])}
  }
  return(f_x)}

#f_x(f_0,F_0,t,k,x)

# force of mortality mu(x) <- -dS_0(x)/S_0(x), mu(x) <- f_0(x)/S_0(x), mu(x+t) <- f_x(t)/S_x(t)
mu_x <- function(f_0,S_0,t,k,x){
  mu_x <- rep(0,length(k))
  for(i in 1:length(k)){
    mu_x[i] <- f_0(t,k)[i,x]/S_0[[i]][1,x]}
  return(mu_x)}

#mu_x(f_0,S_0,t,k,x)

mu_xt <- function(f_0,F_0,f_x,S_x,t,k,x){
  mu_xt <- matrix(0,length(k),(length(t)-x))
  for(i in 1:length(k)){
    for(j in 1:length(t)-x){
      mu_xt[i,j] <- f_x(f_0,F_0,t,k,x)[i,j]/S_x[[i]][1,j]}
  }
  return(mu_xt)}
#mu_xt(f_0,F_0,f_x,S_x,t,k,x)

#Plots( I care about k=2,k=6)
par(oma=c(0,0,2,0)) #<-- To put overall caption
par(mfrow=c(2,2))
#for(i in 1:(length(k)-7)){}
kk <- 6
x1 <- seq(1,length(t)-x)
plot(x1,F_x(F_0,t,k,x)[kk,],xlab="years",ylab = "prob",main ="F_x(t)=P(T_x<t)")
plot(x1,S_x[[kk]],xlab="years",ylab = "prob",main ="S_x(t)=P(T_x>t)-Survival Prob")
#plot(x1,F_x(F_0,x)[2,],xlab="years",ylab = "prob",main ="F_x(t)")
#plot(x1,S_x[[2]],xlab="years",ylab = "prob",main ="S_x(t)=P(T_x>t)-Survival Prob")
#mtext("life aged (x)=2", side=3,line = -1,outer=T,cex=1.5)
plot(x1,f_x(f_0,F_0,t,k,x)[kk,],xlab="years",ylab = "prob",main ="f_x(t)")
plot(x1,mu_xt(f_0,F_0,f_x,S_x,t,k,x)[kk,],xlab="years",ylab = "force",main ="mu_(x+t)")
#plot(x1,S_x[[6]],xlab="years",ylab = "prob",main ="S_x(t)=P(T_x>t)-Survival Prob")
#plot(x1,f_x(f_0,x)[2,],xlab="years",ylab = "prob",main ="f_x(t)")
#plot(x1,mu_xt(f_x,F_x)[2,],xlab="years",ylab = "prob",main ="mu_(x+t)")
#plot(x1,S_x[[2]],xlab="years",ylab = "prob",main ="S_x(t)=P(T_x>t)-Survival Prob")
mtext("life aged (x)=20,k=6", side=3,line = -1,outer=T,cex=1.5)






#Computations!
# P(a < T_x < b)=F_x(b)-F_x(a)
a<- 1
b<- 50
P_ab <- function(F_x,F_0,t,k,x,a,b){
  if(b <= (length(t)-x) && a>0 ){
    Pab <- rep(0,length(k))
    for(i in 1:length(k)){
      Pab[i] <- F_x(F_0,t,k,x)[i,b]-F_x(F_0,t,k,x)[i,a]}}
  else{Pab <- cat(red("Error:"),"Choose a,b appropriately.")}
  return(Pab)}
#P_ab(F_x,F_0,t,k,x,a,b)[6]

#Expectation of T_X = int_{0 to \infty}S_x(t)dt
E_T <- function(S_x,k){
  E_T <- rep(0,length(k))
  for(i in 1:length(k)){
    E_T[i] <- (length(t)-x)*mean(S_x[[i]])}
  return(E_T)}
E_T(S_x,k)[6]

#Variance of T_x = 2*int_{0 to \infty} t*S_x(t)dt
varT <- function(S_x,t,k,x){
  varT <- rep(0,length(k))
  for(i in 1:length(k)){
    tt <- t[1:100]
    varT[i] <- 2*mean(tt*S_x[[i]])}
  return(varT)}
#varT(S_x,t,k,x)[6]
