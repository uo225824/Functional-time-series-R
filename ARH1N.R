pc<-50
t<-50
Tipo=0
p<-0
serie.f<-X[,1:100]
serie<-t(serie.f)
ARH1N<-function(serie.f,t,pc,Tipo,p){
  m<-t-1
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
  N<-nrow(serie)-N1
  media<-mean.fd(fdx[1:N])
  
  
  l<-1
  for (i in 1:pc) {
    l<-l-fdXpca$varprop[i]
    if (l<0.01) {
      nbasis<-i
      break()
    }
  }
  print(nbasis)
  #nbasis<-4
  FRMSE<-0
  for (base in 1:nbasis) {
    rho1<-Nucleophi(fdx,eigenvalues,scoresX,harmonicsX,m,base,N,p)
    a<-0
    XN1<-0
    for (o in p:(N-1)) {
      for (n  in 0:(p-1)) {
        XN<-Data2fd(rho1%*%eval.fd(evalarg = 0:m/m,fdx[o-n]-media)/t,argvals = 0:m/m)+media
        XN1<-XN1+XN
        }
      a[o]<-inprod(XN1-fdx[o+1],XN1-fdx[o+1])
    }

    FRMSE[base]<-sqrt(mean(a))
    
    print(base)
  }
  
  base1<-which.min(FRMSE)
  base1<-3
  p<-0
  rho1<-Nucleophi(fdx,eigenvalues,scoresX,harmonicsX,m,base1,N,p)
  XN1<-0
  
  for (n  in 0:(p-1)) {
    XN<-Data2fd(rho1%*%eval.fd(evalarg = 0:m/m,fdx[N-n]-media)/t,argvals = 0:m/m)+media
    XN1<-XN1+XN
  }  
  par(mfrow=c(2,2))
  par(mar = rep(2, 4))
  plot(fdx[N+1],col="blue")
  plot(XN1,col="red")
  
  plot(fdx[N+1],col="blue")
  lines(XN1,col="red")
  return(FRMSE)
}

