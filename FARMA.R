library(matrixcalc)
t<-12
serie.f<-dserie.funcional
serie.f<-serie.funcional
ARMAH1<-function(serie.f,t,pc,Tipo){
  serie<-t(matrix(serie.f,nrow=t))
  
  if (Tipo==1) {
    serie<-serie.f
    
  }
  N1<-2
  if (length(serie.f)%%t==0) {
    N1<-1
  }
  
  fdx<-Data2fd(t(serie[1:(nrow(serie)-N1+1),]),argvals = c(1:t))
  
  fdXpca=pca.fd(fdx,pc)
  eigenvalues=fdXpca$values; scoresX=fdXpca$scores
  harmonicsX=fdXpca$harmonics
  eigenvectors=eval.fd(evalarg = 1:t,harmonicsX)
  eigenvectors<-Data2fd(eigenvectors,argvals = c(1:t))
  
  l<-1
  for (i in 1:pc) {
    l<-l-fdXpca$varprop[i]
    if (l<0.01) {
      nbasis<-i
      break()
    }
  }
  print(nbasis)
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
  
  print(sqrt(inprod(xn1,xn1)))
  if (sqrt(inprod(xn1,xn1))<1/2) {
    A<-matrix(0,nrow = length(xn1$coefs),ncol = 1)
    for (i in 1:1000) {
      A=hadamard.prod(A^2,xn1$coefs)+xn1$coefs
    }
    basisy<- create.fourier.basis(c(1, t), 51)
    MA<-Data2fd(A,seq(1,12,length=14),basisobj = basisy)
    XN1<-xn1+media+MA
    print("ARMA")
  }
  if (sqrt(inprod(xn1,xn1))>1/2) {
    sc<-sqrt(inprod(xn1,xn1))
    xn1p<-xn1
    xn1p$coefs<-xn1p$coefs/(3*sc[1])
    A<-matrix(0,nrow = length(xn1$coefs),ncol = 1)
    for (i in 1:1000) {
      A=hadamard.prod(A^2,xn1p$coefs)+xn1p$coefs
    }
    basisy<- create.fourier.basis(c(1, t), 51)
    MA<-Data2fd(A,seq(1,t,length=t+2),basisobj = basisy)
    MA$coefs<-MA$coefs*(3*sc[1])
    XN1<-xn1+media+MA
    print("ARMA")}

  
  par(mfrow=c(2,2))
  par(mar = rep(2, 4))
  plot(fdx[N+1],col="blue")
  plot(XN1,col="red")
  plot(MA)
  plot(fdx[N+1],col="blue")
  lines(XN1,col="red")
  lines(xn1+media,col="green")
  dis1<-sqrt(inprod(XN1-fdx[N+1],XN1-fdx[N+1]))
  dis2<-sqrt(inprod(xn1+media-fdx[N+1],xn1+media-fdx[N+1]))
  dis<-data.frame(c(dis1[1],dis2[1]))
  return(dis)
}
