set.seed(12345)


library(igraph)
#library(MCMCpack)
#library(mvtnorm)
#library(qualityTools)
library(ellipse)
#library(matlib)
#library(MASS)

#SBM with M1,M2
NN <- seq(200, 500, by=10)
#NN <- c(1,2,3,4,5,6)
m <- 3
P <- list()
M1 <- list()
M2 <- list()
graphs <- list()
adj_graphs <-list()
c1 <- c(.5,.5,.5)
c2 <- c(0.1,0.4,0.7)

p <- 0.7
q <- 0.2
m <- 3
K <- 2


for(i in 1:length(NN)){
  M1[[i]] <- matrix(0,m*NN[i],m*NN[i])
  M2[[i]] <- matrix(0,m*NN[i],m*NN[i])
}

#Generate SBMs from the same block prob matrix
bm <- cbind(c(p,q),c(q,p))
for(k in 1:length(NN)){
  for(j in 1:m)
    
    #rdpg[[3*(k-1)+j]] <- sample_dot_product(t(X[[k]]))
    graphs[[3*(k-1)+j]]<- sample_sbm(NN[k],pref.matrix=bm, block.sizes=c(NN[k]/2,NN[k]/2))}


#Get the adjacency matrices
for(k in 1:length(NN)){
  for(j in 1:m)
    adj_graphs[[3*(k-1)+j]] <- as_adjacency_matrix(graphs[[3*(k-1)+j]],type = "upper",sparse = FALSE)}

#Create prob matrix P
#P1<- cbind(r*matrix(1,N/K,N/K),p*matrix(1,N/K,N/K))
#P2 <- cbind(p*matrix(1,N/K,N/K),r*matrix(1,N/K,N/K))
#P3 <- cbind(q*matrix(1,N/K,N/K),p*matrix(1,N/K,N/K),r*matrix(1,N/K,N/K))

#P <- rbind(P1,P2)
#diag(P) <- 0
#P
#graphs[[70]]

#Omni P--> P_hat
#P_hat <- matrix(0,N*m,N*m)
#for(i in 1:m){
#  for(j in 1:m){
#    P_hat[((i-1)*N+1):(i*N),((j-1)*N+1):(j*N)] <- P}
#}


#Omnibus matrix
#M <- matrix(0,N*m,N*m)
B1 <- list()
B2 <- list()
for(k in 1:length(NN)){
  
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      B1[[(3*(k-1)+i+j-2)]] <- (c1[i]*adj_graphs[[3*(k-1)+i]]+c1[j]*adj_graphs[[3*(k-1)+j]])/(c1[i]+c1[j])
      B2[[(3*(k-1)+i+j-2)]] <- (c2[i]*adj_graphs[[3*(k-1)+i]]+c2[j]*adj_graphs[[3*(k-1)+j]])/(c2[i]+c2[j])
      #print(i)
    }
  }
  
  for(i in 1:m){
    for(j in 1:m){
      if(k == 1){
        if(i == j){
          M1[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- adj_graphs[[3*(k-1)+i]]
          M2[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- adj_graphs[[3*(k-1)+i]]}
        else if(i<j){
          M1[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B1[[(3*(k-1)+i+j-2)]]
          M2[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B2[[(3*(k-1)+i+j-2)]]}
        else{M1[[k]][(((i-1)*NN[k]+1)):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B1[[(3*(k-1)+i+j-2)]]
        M2[[k]][(((i-1)*NN[k]+1)):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B2[[(3*(k-1)+i+j-2)]]}
      }
      else{
        if(i == j){
          #print(i)
          M1[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- adj_graphs[[3*(k-1)+i]]
          M2[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- adj_graphs[[3*(k-1)+i]]}
        else if(i<j){
          M1[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B1[[(3*(k-1)+i+j-2)]]
          M2[[k]][((i-1)*NN[k]+1):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B2[[(3*(k-1)+i+j-2)]]}
        else{M1[[k]][(((i-1)*NN[k]+1)):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B1[[(3*(k-1)+i+j-2)]]
        M2[[k]][(((i-1)*NN[k]+1)):(i*NN[k]),((j-1)*NN[k]+1):(j*NN[k])] <- B2[[(3*(k-1)+i+j-2)]]}
      }
    }
  }
}

#Embeddings of  M

embed_M1 <- list()
embed_M2 <- list()
conv_M1 <- list()
conv_M2 <- list()


for(i in 1:length(NN)){
  conv_M1[[i]] <- graph_from_adjacency_matrix(M1[[i]],mode = "undirected")
  embed_M1[[i]] <- embed_adjacency_matrix(conv_M1[[i]], 2)
  conv_M2[[i]] <- graph_from_adjacency_matrix(M2[[i]],mode = "undirected")
  embed_M2[[i]] <- embed_adjacency_matrix(conv_M2[[i]], 2)
}
#plot(conv_M[[5]])
#Discrepancy of the rows// Check Normality (1st row with n+1 th row)
dis_ncy12c1 <- matrix(0,length(NN),2)
dis_ncy12c2 <- matrix(0,length(NN),2)
for(k in 1:length(NN)){
  dis_ncy12c1[k,] <- sqrt(NN[k])*c(embed_M1[[k]]$X[1,]-embed_M1[[k]]$X[1+NN[k],])
  dis_ncy12c2[k,] <- sqrt(NN[k])*c(embed_M2[[k]]$X[1,]-embed_M2[[k]]$X[1+NN[k],])}

#Discrepancy of the rows// Check Normality (1st row with 2n+1 th row)
dis_ncy13c1 <- matrix(0,length(NN),2)
dis_ncy13c2 <- matrix(0,length(NN),2)
for(k in 1:length(NN)){
  dis_ncy13c1[k,] <- sqrt(NN[k])*c(embed_M1[[k]]$X[1,]-embed_M1[[k]]$X[1+2*NN[k],])
  dis_ncy13c2[k,] <- sqrt(NN[k])*c(embed_M2[[k]]$X[1,]-embed_M2[[k]]$X[1+2*NN[k],])}


#Discrepancy of the rows// Check Normality (n+1 th row with 2n+1 th row)
dis_ncy23c1 <- matrix(0,length(NN),2)
dis_ncy23c2 <- matrix(0,length(NN),2)
for(k in 1:length(NN)){
  dis_ncy23c1[k,] <- sqrt(NN[k])*c(embed_M1[[k]]$X[1+NN[k],]-embed_M1[[k]]$X[1+2*NN[k],])
  dis_ncy23c2[k,] <- sqrt(NN[k])*c(embed_M2[[k]]$X[1+NN[k],]-embed_M2[[k]]$X[1+2*NN[k],])}


#Mean value of the Covariance matrix
#This is for p=0.7,q=0.2
vect1 <- bm*matrix(cbind(.2390,.8017)%*%rbind(.2390,.8017)-(cbind(.2390,.8017)%*%rbind(.2390,.8017))^2,2,2)

#This is for p=0.3,q=0.2
#vect1 <- bm*matrix(cbind(0.3651,.4143)%*%rbind(0.3651,.4143)-(cbind(0.3651,.4143)%*%rbind(0.3651,.4143))^2,2,2)

invSBM <- bm^{-1}


covM1 <- (2/m^2)*invSBM%*%vect1%*%invSBM*(1+(c1[1]/(c1[1]+c1[3]))+(c1[2]/(c1[2]+c1[3]))+(c1[1]^2/(c1[1]+c1[3])^2)+
                                            (c1[2]^2/(c1[2]+c1[3])^2)-c1[1]*c1[2]*(c1[1]*c1[2]+c1[1]*c1[3]+c1[2]*c1[3]+c1[3]^2)/(((c1[1]+c1[3])^2)*(c1[2]+c1[3])^2))
#12
covM12 <- (2/m^2)*invSBM%*%vect1%*%invSBM*(1+(c2[1]/(c2[1]+c2[3]))+(c2[2]/(c2[2]+c2[3]))+(c2[1]^2/(c2[1]+c2[3])^2)+
                                             (c2[2]^2/(c2[2]+c2[3])^2)-c2[1]*c2[2]*(c2[1]*c2[2]+c2[1]*c2[3]+c2[2]*c2[3]+c2[3]^2)/(((c2[1]+c2[3])^2)*(c2[2]+c2[3])^2))

#13
covM13 <- (2/m^2)*invSBM%*%vect1%*%invSBM*(1+(c2[1]/(c2[1]+c2[2]))+(c2[3]/(c2[3]+c2[2]))+(c2[1]^2/(c2[1]+c2[2])^2)+
                                             (c2[3]^2/(c2[3]+c2[2])^2)-c2[1]*c2[3]*(c2[1]*c2[3]+c2[1]*c2[2]+c2[3]*c2[2]+c2[2]^2)/(((c2[1]+c2[2])^2)*(c2[3]+c2[2])^2))

#23
covM23 <- (2/m^2)*invSBM%*%vect1%*%invSBM*(1+(c2[2]/(c2[2]+c2[1]))+(c2[3]/(c2[3]+c2[1]))+(c2[2]^2/(c2[1]+c2[2])^2)+
                                             (c2[3]^2/(c2[3]+c2[1])^2)-c2[2]*c2[3]*(c2[2]*c2[3]+c2[1]*c2[2]+c2[3]*c2[1]+c2[1]^2)/(((c2[1]+c2[2])^2)*(c2[3]+c2[1])^2))


#(2/m^2)*invSBM%*%vect1%*%invSBM#x <- rmvnorm(1000,mean = c(0,0),sigma =covM)
#plot(ellipse(covM),type="l",col = "brown")
#del <- t(t(X[[1]][1000,]))%*%X[[1]][1000,]
#invdel <- del^{-1}
#covMM <- (2/m^2)*invdel%*%del%*%invdel*matrix(X[[1]][1000,]%*%t(t(X[[1]][1000,]))-(X[[1]][1000,]%*%t(t(X[[1]][1000,])))^2,2,2)%*%invdel*(1+(c1[1]/(c1[1]+c1[3]))+(c1[2]/(c1[2]+c1[3]))+(c1[1]^2/
#(c1[1]+c1[3])^2)+(c1[2]^2/(c1[2]+c1[3])^2)-c1[1]*c1[2]*(c1[1]*c1[2]+c1[1]*c1[3]+c1[2]*c1[3]+c1[3]^2)/(((c1[1]+c1[3])^2)*(c1[2]+c1[3])^2))
#covM
par(oma=c(0,0,2,0)) #<-- To put overall caption
par(mfrow=c(1,3))
#covM<- cbind(c(0.52,0.1),c(0.1,0.52))
plot(dis_ncy12c1,col = "blue",main = "(12)blue pts c(1/2,1/2,1/2), red pts c(0.1,0.4,0.7)",xlim=range(dis_ncy12c1[,1]-3,dis_ncy12c2[,1]+3)
     , ylim=range(dis_ncy12c1[,2]-3, dis_ncy12c2[,2]+3))
#plot(dis_ncy12,col = "blue",xlim=range(dis_ncy12[,1]-2,dis_ncy12[,1]+2), ylim=range(dis_ncy12[,2]-10, dis_ncy12[,2]+10))
par(new=TRUE)
points(dis_ncy12c2,col = "red")
points(ellipse(covM1),lwd = 3,type="l",col = "cyan",lty= "dotted")
points(ellipse(covM12),lwd = 3,type="l",col = "brown",lty = "dotted")
legend("topleft",NULL,c("Classical OMNI","weighted OMNI c(0.1,0.4,0.7)"),col=c("cyan","brown"),lty=1:2,bty="n")
mtext("SBM with 2 communities(p=0.7,q=0.2)",cex = 1.5,side = 3, outer = T,line = -1)


#covM<- cbind(c(0.52,0.1),c(0.1,0.52))
plot(dis_ncy13c1,col = "blue",xlim=range(dis_ncy13c1[,1]-3,dis_ncy13c2[,1]+3), ylim=range(dis_ncy13c1[,2]-3, dis_ncy13c2[,2]+3))
#plot(dis_ncy12,col = "blue",xlim=range(dis_ncy12[,1]-2,dis_ncy12[,1]+2), ylim=range(dis_ncy12[,2]-10, dis_ncy12[,2]+10))
par(new=TRUE)
points(dis_ncy13c2,col = "red")
points(ellipse(covM1),lwd = 3,type="l",col = "cyan",lty="dotted")
points(ellipse(covM13),lwd = 3,type="l",col = "brown",lty="dotted")
legend("topleft",NULL,c("Classical OMNI","weighted OMNI c(0.1,0.4,0.7)"),col=c("cyan","brown"),lty=1:2,bty="n")
mtext("n=51,SBM with 2 communities(p=0.7,q=0.2)",cex = 1.5,side = 3, outer = T,line = -1)

#par(oma=c(0,0,2,0)) #<-- To put overall caption
#par(mfrow=c(3,2))
#covM<- cbind(c(0.52,0.1),c(0.1,0.52))
plot(dis_ncy23c1,col = "blue",main = "23-blue pts c(1/2,1/2,1/2), red pts c(0.1,0.4,0.7)",xlim=range(dis_ncy23c1[,1]-3,dis_ncy23c2[,1]+3)
     , ylim=range(dis_ncy23c1[,2]-3, dis_ncy23c2[,2]+3))
#plot(dis_ncy12,col = "blue",xlim=range(dis_ncy12[,1]-2,dis_ncy12[,1]+2), ylim=range(dis_ncy12[,2]-10, dis_ncy12[,2]+10))
par(new=TRUE)
points(dis_ncy23c2,col = "red")
points(ellipse(covM1),lwd = 3,type="l",col = "cyan",lty = "dotted")
points(ellipse(covM23),lwd = 3,type="l",col = "brown",lty = "dotted")
legend("topleft",NULL,c("Classical OMNI","weighted OMNI c(0.1,0.4,0.7)"),col=c("cyan","brown"),lty=1:2,bty="n")
mtext("n=51,SBM with 2 communities(p=0.7,q=0.2)",cex = 1.5,side = 3, outer = T,line = -1)
shapiro.test(dis_ncy23c1)

