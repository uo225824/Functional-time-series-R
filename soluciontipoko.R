library(fda.usc)

ARH1<-function(serie,t,pc,N1){
  fdx<-Data2fd(t(serie),argvals = c(1:ttiempo))
  fdXpca=pca.fd(fdx,pc)
  eigenvalues=fdXpca$values; scoresX=fdXpca$scores
  harmonicsX=fdXpca$harmonics
  eigenvectors=eval.fd(evalarg = 1:ttiempo,harmonicsX)
  eigenvectors<-Data2fd(eigenvectors,argvals = c(1:ttiempo))
  
  l<-1
  for (i in 1:pc) {
    l<-l-fdXpca$varprop[i]
    if (l<0.01) {
      nbasis<-i
      break()
    }
  }
  print(nbasis)
  N<-nrow(serie.f)-N1
  rho<-0
  rho1<-0
  rho2<-0
  media<-mean.fd(fdx)
  for (l in 1:nbasis) {
    for (i in 1:(N-1)) {
      for (j in 1:nbasis ) {
        rho[j]<-(1/eigenvalues[j])*inprod(fdx[N]-media,eigenvectors[j])*inprod(fdx[i]-media,eigenvectors[j])*inprod(fdx[i+1]-media,eigenvectors[l])
      }
      rho1[i]<-sum(rho)
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
  plot(fdx[N+1],col="blue")
  par(new=T)
  plot(XN1,col="red")
}

ARH1(serie.f,48,12,2)



