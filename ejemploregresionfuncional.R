set.seed(12345)
t = seq(0, 1, length = nt <- 100)
covexp = function(t1, t2) {
  3 * exp(-abs(t1 - t2)/0.5)
}
Sigma = outer(t, t, covexp)
X = rproc2fdata(n <- 200, t, sigma = Sigma)
Y<-matrix(0,nrow = 200,ncol = 100)

beta<-matrix(0,ncol = 100,nrow = 100)
t<-seq(0,1,length.out = 100)
s<-seq(0,1,length.out = 100)
for(i in 1:100){
  for (j in 1:100) {
    beta[j,i]<-exp(-(t[i]^2+s[j]^2)/2)
  }
  
}

epsilon<- matrix(0,200,100)
m<-199
for(i in 1:100){
  epsilon[,i]<-0.5*rnorm(1)*sin(pi*(0:m/m))+0.5*rnorm(1)*cos(2*pi*(0:m/m))
}

Y<-X$data%*%beta/100+epsilon
opt.bsY<-optim.basis(fdata(Y),type.basis = "bspline")
opt.bsX<-optim.basis(X,type.basis = "bspline")



basisx <- create.bspline.basis(c(0, 1),nbasis = 45,norder = 5 )
basisy<-  create.bspline.basis(c(0,1),nbasis = 6,norder = 5 )
datosf<-Data2fd(t(X$data),seq(0,1,length=100),basisobj = basisx)
datosfN<-Data2fd(t(Y),seq(0,1,length=100),basisobj = basisy)

lambdas<-1
lambdat<-1

par(mfrow=c(3,4))
for (i in 1:10) {
  lambdas<-lambdas/10
  lambdat<-lambdat/10
  res.fr = fregre.basis.fr(datosf,datosfN,basisx,basisy,lambda.s = lambdas,lambda.t = lambdat,Lfdobj.s = 1,Lfdobj.t = 1)
  persp((0:44/44),(0:5/5),z=res.fr$beta.estbifd$coefs,cex.axis = axisticksize,
        cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,
        phi=30,
        ticktype="detailed", main="Estimación de Beta(t,s)",col = "red")
  
  
  }

phi=9/4*(0:m/m)%*%t(0:m/m)
# True surface evaluated on a discrete bivariate grid
par(mar=c(1.5, 1.5, 3.5, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2
                                                     , 0.2))
par(mfrow=c(2,2))
# 4 panels 2 rows and 2 columns arrangement
axislabelsize=1.5 # controls the size of labels of axes
axisticksize=0.8 # controls the size of ticks of axes
persp((0:99/99),(0:199/199),z=beta,cex.axis = axisticksize,
      cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,
      phi=30,
      ticktype="detailed", main="Valor teórico de Beta(t,s)",col = "blue")

persp((0:30/30),(0:30/30),z=res.fr$beta.estbifd$coefs,cex.axis = axisticksize,
      cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,
      phi=30,
      ticktype="detailed", main="True",col = "blue")




