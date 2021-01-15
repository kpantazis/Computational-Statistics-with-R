set.seed(12345)


library(igraph)
#library(MCMCpack)
#library(mvtnorm)
library(qualityTools)
#library(ellipse)
#library(matlib)
#library(MASS)


############################################################################################################3
#############################################################################################################
NN <- seq(200, 500, by=10)
#NN <- c(1,2,3,4,5,6)
m <- 3
X <- list()
P <- list()
M <- list()
#discrepancies in one row
graphs <- list()
adj_graphs <- list()
for(i in 1:length(NN)){
  M[[i]] <- matrix(0,m*NN[i],m*NN[i])
}
p <- 0.7
q <- 0.2
m <- 3
K <- 2
#Generate SBMs from the same block prob matrix
bm <- cbind(c(p,q),c(q,p))
for(k in 1:length(NN)){
  for(j in 1:m){
    #rdpg[[3*(k-1)+j]] <- sample_dot_product(t(X[[k]]))
    graphs[[3*(k-1)+j]]<- sample_sbm(NN[k],pref.matrix=bm, block.sizes=c(NN[k]/2,NN[k]/2))}
}

#Get the adjacency matrices
for(k in 1:length(NN)){
  for(j in 1:m)
    adj_graphs[[3*(k-1)+j]] <- as_adjacency_matrix(graphs[[3*(k-1)+j]],type = "upper",sparse = FALSE)}

#Omnibus matrix
#M <- matrix(0,N*m,N*m)
c <- c(0.1,0.4,0.7)
B <- list()
for(k in 1:length(NN)){
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      B[[(3*(k-1)+i+j-2)]] <- (c[i]*adj_graphs[[3*(k-1)+i]]+c[j]*adj_graphs[[3*(k-1)+j]])/(c[i]+c[j])
      #print(i)
    }
  }
  
  for(i in 1:m){
    for(j in 1:m){
      if(k == 1){
        if(i == j){
          M[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- adj_graphs[[3*(k-1)+i]]}
        else if(i<j){
          M[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B[[(3*(k-1)+i+j-2)]]}
        else{M[[k]][(((i-1)*NN[k]+1)):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B[[(3*(k-1)+i+j-2)]]}
      }
      else{
        if(i == j){
          #print(i)
          M[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- adj_graphs[[3*(k-1)+i]]}
        else if(i<j){
          M[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B[[(3*(k-1)+i+j-2)]]}
        else{M[[k]][(((i-1)*NN[k]+1)):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B[[(3*(k-1)+i+j-2)]]}
      }
    }
  }
}

#Embeddings of  M

embed_M <- list()

conv_M <- list()


for(i in 1:length(NN)){
  conv_M[[i]] <- graph_from_adjacency_matrix(M[[i]],mode = "undirected")
  embed_M[[i]] <- embed_adjacency_matrix(conv_M[[i]], 2)
}


#Discrepancy of the rows// Check Normality (First row with n+1 th row)
dis_ncy12 <- matrix(0,length(NN),2)
for(i in 1:length(NN)){
  dis_ncy12[i,] <- sqrt(NN[i])*c(embed_M[[i]]$X[1,]-embed_M[[i]]$X[1+NN[i],])}

#Discrepancy of the rows// Check Normality (First row with 2n+1 th row)
dis_ncy13 <- matrix(0,length(NN),2)
for(i in 1:length(NN)){
  dis_ncy13[i,] <- sqrt(NN[i])*c(embed_M[[i]]$X[1,]-embed_M[[i]]$X[1+2*NN[i],])}

#Discrepancy of the rows// Check Normality (n+1 th row with 2n+1 th row)
dis_ncy23 <- matrix(0,length(NN),2)
for(i in 1:length(NN)){
  dis_ncy23[i,] <- sqrt(NN[i])*c(embed_M[[i]]$X[1+NN[i],]-embed_M[[i]]$X[1+2*NN[i],])}

par(oma=c(0,0,2,0)) #<-- To put overall caption
par(mfrow=c(3,3))
#plot(dis_ncy12)
#qqPlot(omni1)
qqPlot(dis_ncy12)
qqPlot(dis_ncy13)
qqPlot(dis_ncy23)
#qqnorm(dis_ncy23)
#qqline(dis_ncy23)
#ppPlot(dis_ncy23,yaxis = FALSE)
#hist(omni1)
hist(dis_ncy12,breaks = 5)
hist(dis_ncy13,breaks = 5)
hist(dis_ncy23,breaks = 5)
#xxx <- dmvnorm(dis_ncy23 ,mean = c(0,0),sigma = sqrt(covM))
#lele <- function(dis_ncy23){
# return(dmvnorm(dis_ncy23 ,mean = c(0,0),sigma = sqrt(covM)))
#}
#ellipse(lele)
#boxplot(omni1)
boxplot(dis_ncy12)
boxplot(dis_ncy13)
boxplot(dis_ncy23)
mtext("Differences between the rows of the latent positions and their 
      corresponding estimates/SBM(K=2),m=3,c(0.1,0.4,0.7))", 
      side=3,line = -1,outer=T,cex=1.5)
#shapiro.test(omni1)
shapiro.test(dis_ncy12)
shapiro.test(dis_ncy13)
shapiro.test(dis_ncy23)
