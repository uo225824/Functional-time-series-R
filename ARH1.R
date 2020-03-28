library(fda.usc)

pc<-50
t<-50
Tipo=0
serie.f<-Xfs
ARH1<-function(serie.f,t,pc,Tipo){
  serie<-t(matrix(serie.f,nrow=t))
  
  if (Tipo==1) {
    serie<-serie.f
    
  }
  N1<-2
  if (length(serie.f)%%t==0) {
    N1<-1
  }

  fdx<-Data2fd(t(serie[1:(nrow(serie)-N1+1),]),argvals = 0:m/m)
  
  fdXpca=pca.fd(fdx,pc)
  eigenvalues=fdXpca$values; scoresX=fdXpca$scores
  harmonicsX=fdXpca$harmonics
  eigenvectors=eval.fd(evalarg = 0:m/m,harmonicsX)
  eigenvectors<-Data2fd(eigenvectors,argvals = 0:m/m)
  
  l<-1
  for (i in 1:pc) {
    l<-l-fdXpca$varprop[i]
    if (l<0.01) {
      nbasis<-i
      break()
    }
  }
  print(nbasis)
  nbasis<-4
  N<-nrow(serie)-N1
  rho<-0
  rho1<-0
  rho2<-0
  media<-mean.fd(fdx[1:N])
  for (l in 1:nbasis) {
    for (i in 1:(N-1)) {
      for (j in 1:nbasis ) {
        rho[j]<-(1/eigenvalues[j])*inprod(fdx[N]-media,eigenvectors[j])*inprod(fdx[i]-media,eigenvectors[j])*inprod(fdx[i+1]-media,eigenvectors[l])
      }
      rho1[i]<-sum(rho)/2
    }
    rho2[l]<-sum(rho1)/(N-1)
    print(l)
  }
  xn<-0
  xn1<-0
  for (l in 1:nbasis) {
    xn<-rho2[l]*eigenvectors[l]
    xn1<-xn1+xn
  }
  XN1<-xn1+media
  par(mfrow=c(2,2))
  par(mar = rep(2, 4))
  plot(fdx[N+1],col="blue")
  plot(XN1,col="red")

  plot(fdx[N+1],col="blue")
  lines(XN1,col="blue")
  A<-matrix(0,ncol = 2,nrow = length(fdx[N+1]$coefs))
  A[,1]<-XN1$coefs
  A[,2]<-fdx[N+1]$coefs
  myfdata<-fdata(A)
  dis<-metric.lp(myfdata, lp = 2)[2]
  return(XN1)
}

ARH1(serie.f,48,12,2)


serie.1 <- window(x, end=end(x) - c(0,t))
serie.2 <- window(x, start=end(serie.1) + c(0,1))
plot.ts(serie.2,col="blue")
par(new=T)