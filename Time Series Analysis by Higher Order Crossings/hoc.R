set.seed(121)
N <- 1000
NN <- c(1:N)
ts <- c(1:5000)
x1 <- matrix(0,length(ts),length(NN))
x2 <- matrix(0,length(ts),length(NN))
phi <- matrix(0,length(ts),length(NN))
phizc <- rep(0,length(ts))
z <- rep(0,N)
d <- rep(0,N)
for (i in ts){
myts <- ts(rnorm(N))
for (j in NN){
  x1[i,j] <- ifelse(z[j]>=0,1,0)
  z[j] <- 0.8*z[j] + rnorm(1)
  x2[i,j] <- ifelse(z[j]>=0,1,0)
  d[j] <- d[j] + (x2[i,j]-x1[i,j])^2
  phi[i,j] <- cos(3.14159*d[j]/999)
  if (j == N){
    phizc[i] <- phi[i,N]
  }
}
}
Reg <- lm(formula = phizc ~ phi[,1])
summary(Reg)
c <-seq(0,1,1/(5000-1))
length(c)
plot(c,abs(phi[,1]-phizc),main = "Absolute difference between the values of phi and phizc")

#We used the 1st column,but we can use any; we obserbe that as c goes to 1, their difference goes to 0 and this is true 
#for every column.
