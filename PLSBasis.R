

PLSbasis1<-function(Xf,Yf,ttiempo){
  nbasis=31
  basisx <- create.bspline.basis(c(0, 1),nbasis = 31,norder = 5 )
  basisy<-create.bspline.basis(c(0, 1),nbasis = 31,norder = 5 )
  ndata<-ncol(Xf$coefs)
  X0<-matrix(0,nrow = 31,ndata)
  for (i in 1:ndata) {
    X0[,i]<-(Xf[i]$coefs-mean.fd(Xf)$coefs)
  }
  Xf0<-Data2fd(X0,seq(0,1,length=nbasis),basisobj = basisx)
  
  Y0<-matrix(0,nrow = 31,ndata)
  for (i in 1:ndata) {
    Y0[,i]<-(Yf[i]$coefs-mean.fd(Yf)$coefs)
    
  }
  Yf0<-Data2fd(Y0,seq(0,1,length=nbasis),basisobj = basisy)  
  ys<-colMeans(Yf0$coefs)
  
  
  A<-matrix(0,nrow = 31,ndata)
  for (j in 1:ndata) {
    A[,j]<-Xf0[j]$coefs*ys[j]  
  }
  
  Af<-mean.fd(Data2fd(A,seq(0,1,length=nbasis),basisobj = basisy))
  
  K1<-Af$coefs/as.numeric(sqrt(inprod(Af,Af)))
  
  K1f<-Data2fd(K1,seq(0,1,length=nbasis),basisobj = basisy)
  
  eta1<-inprod(K1f,Xf0)
  
  
  P<-matrix(0,nrow = 31,ndata)
  for (j in 1:ndata) {
    P[,j]<-Xf0[j]$coefs*eta1[j]  
  }
  
  P<-mean.fd(Data2fd(P,seq(0,1,length=nbasis),basisobj = basisy))
  
  P1f<-Data2fd(P$coefs/as.numeric(var(as.vector(eta1))),seq(0,1,length=nbasis),basisobj = basisy)
  
  Zeta<-matrix(0,nrow = 31,ndata)
  for (j in 1:ndata) {
    Zeta[,j]<-Yf0[j]$coefs*eta1[j]  
  }
  
  Zeta<-mean.fd(Data2fd(Zeta,seq(0,1,length=nbasis),basisobj = basisy))
  
  Zeta1f<-Data2fd(Zeta$coefs/as.numeric(var(as.vector(eta1))),seq(0,1,length=nbasis),basisobj = basisy)
  
  X1<-matrix(0,nrow = 31,ndata)
  for (j in 1:ndata) {
    Xf01<-eval.fd(Xf0[j],seq(0,1,length=nbasis))
    X1[,j]<-Xf01-P1f[j]$coefs*as.numeric(eta1[j])  
  }
  
  Xf1<-Data2fd(X1,seq(0,1,length=nbasis),basisobj = basisx)
  
  Y1<-matrix(0,nrow = 31,ndata)
  for (j in 1:ndata) {
    Yf01<-eval.fd(Yf0[j],seq(0,1,length=nbasis))
    Y1[,j]<-Yf01-Zeta1f[j]$coefs*as.numeric(eta1[j])  
  }
  
  Yf1<-Data2fd(Y1,seq(0,1,length=nbasis),basisobj = basisx)
  return(list(P1f,Zeta1f,Xf1,Yf1,eta1))
  
}
ttiempo<-50
X<-serie
Xf<-Data2fd(t(X[1:99,]),seq(0,1,length=ttiempo),basisobj = basisx)
Yf<-Data2fd(t(X[2:100,]),seq(0,1,length=ttiempo),basisobj = basisy)
Salida<-PLSbasis1(Xf,Yf,ttiempo = 50)

P1<-matrix(0,ncol = 5,nrow = 31)
P1[,1]<-Salida[[1]]$coefs
for (i in 2:5) {
  Salida<-PLSbasis1(Salida[[3]],Salida[[4]])
  P1[,i]<-Salida[[1]]$coefs
}

