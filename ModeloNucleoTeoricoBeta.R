set.seed(12345)
m=49 # each function is observed at m+1 points, including 0 and 1
burnin=50 # the first 50 functions are a burn in period
N=100 # number of functions to simulate
N1=N+burnin
alpha=9/4
# Create 2 matrices, whose entries are 0s.
# Each column represents one function.
X<- matrix(0,50,100)

epsilon<- matrix(0,50,100)
epsilon[,1]<-rnorm(1,sd = 0.02)*sin(pi*(0:m/m))+0.5*rnorm(1,sd = 0.02)*cos(2*pi*(0:m/m))
# the following loop simulates FAR(1).
beta<-matrix(0,ncol = 50,nrow = 50)
t<-seq(0,1,length.out = 50)
for(i in 1:50){
  for (j in 1:50) {
    beta[i,j]<-exp(-(t[i]^2+s[j]^2)/2)
  }
 
}
X[,1]<-epsilon[,1]
for(i in 2:N){
  epsilon[,i]<-rnorm(1)*sin(pi*(0:m/m))+0.5*rnorm(1)*cos(2*pi*(0:m/m))
  X[,i]<-1/2*beta%*%X[,i-1]/50+epsilon[,i]
}
 # Remove the burn in period functions


last=100
N=100
plot.ts(c(X[,(N-last+1):N]),ylim=c(min(X[,(N-last):N])-0.5,0.5
                                   +max(X[,(N-last):N])),axes=F,xlab="",ylab="",lwd=2)
axis(2); axis(1,tick=F,labels=F); abline(h=0)
abline(v=seq(0,last*(m+1),by=m+1), lty=2); box()




phi=9/4*(0:m/m)%*%t(0:m/m)
# True surface evaluated on a discrete bivariate grid
par(mar=c(1.5, 1.5, 3.5, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2
                                                     , 0.2))
par(mfrow=c(2,2))
# 4 panels 2 rows and 2 columns arrangement
axislabelsize=1.5 # controls the size of labels of axes
axisticksize=0.8 # controls the size of ticks of axes
persp((0:m/m),(0:m/m),z=beta,cex.axis = axisticksize,
      cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,
      phi=30,
      ticktype="detailed", main="True",col = "blue")
