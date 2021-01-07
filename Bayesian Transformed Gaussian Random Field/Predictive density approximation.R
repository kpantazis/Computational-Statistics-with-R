Data <- read.table("C:/Users/Pantazis/Desktop/Fall 2019/Comp Statistics(STAT705)/cattxt.cgi")
x <- seq(0, 255)
y <- seq(0, 255)
z <- matrix(0,256,256)

for(i in 1:256) {for(j in 1:256){ z[i,j] <- Data[256*(i-1)+j,3]}}

-------------------------------------------------
  
  
  ## These lines plot the necessary graphs
  
  par(mfrow=c(1,2))
  par(fig=c(0,0.5,0,1))
  
  #persp(x, y, z, main="Matern (10,2)", 
   #     col="blue1", phi=20, theta=50,r=50, d=0.1,expand=0.5, ltheta=30, 
    #    lphi=180, shade=0.15, ticktype="detailed", nticks=1)
  
  
  persp(x, y, z, main="Matern (10,2)", 
        col="red", phi=20, theta=50,r=50, d=0.1,expand=0.5, ltheta=30, 
        lphi=180, shade=0.15, ticktype="detailed", nticks=1)
  

par(fig=c(0.5,1,0.15,0.85),new=TRUE)

contour(x,y,z,)


###With color
#filled.contour(x, y, z, main="Matern(10,2)", color=rainbow)

#Or 
#filled.contour(x, y, z, main="Matern(10,2)", color=terrain.colors)

