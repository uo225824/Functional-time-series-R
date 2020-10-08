nbasis=60
basisx <- create.bspline.basis(c(0, 1),nbasis = nbasis,norder = 5 )
basisy<-create.bspline.basis(c(0, 1),nbasis = nbasis,norder = 5 )
ndata<-ncol(Xf0$coefs)

Xf0<-Data2fd(t(X$data),seq(0,1,length=101),basisobj = basisx)
plot(Xf0)
ys<-y


A<-matrix(0,nrow = nbasis,ndata)
for (j in 1:ndata) {
  A[,j]<-Xf0[j]$coefs*ys[j]  
}

Af<-mean.fd(Data2fd(A,seq(0,1,length=nbasis),basisobj = basisy))

K1<-Af$coefs/as.numeric(sqrt(inprod(Af,Af)))

K1f<-Data2fd(K1,seq(0,1,length=nbasis),basisobj = basisy)

eta1<-inprod(K1f,Xf0)


P<-matrix(0,nrow = nbasis,ndata)
for (j in 1:ndata) {
  P[,j]<-Xf0[j]$coefs*eta1[j]  
}

P<-mean.fd(Data2fd(P,seq(0,1,length=nbasis),basisobj = basisy))

P1f<-Data2fd(P$coefs/as.numeric(var(as.vector(eta1))),seq(0,1,length=nbasis),basisobj = basisy)

Zeta<-matrix(0,nrow = nbasis,ndata)
for (j in 1:ndata) {
  Zeta[,j]<-Yf0[j]$coefs*eta1[j]  
}

Zeta<-mean.fd(Data2fd(Zeta,seq(0,1,length=nbasis),basisobj = basisy))

Zeta1f<-Data2fd(Zeta$coefs/as.numeric(var(as.vector(eta1))),seq(0,1,length=nbasis),basisobj = basisy)

X1<-matrix(0,nrow = nbasis,ndata)
for (j in 1:ndata) {
  Xf01<-eval.fd(Xf0[j],seq(0,1,length=nbasis))
  X1[,j]<-Xf01-P1f[j]$coefs*as.numeric(eta1[j])  
}

Xf1<-Data2fd(X1,seq(0,1,length=nbasis),basisobj = basisx)

Y1<-matrix(0,nrow = nbasis,ndata)
for (j in 1:ndata) {
  Yf01<-eval.fd(Yf0[j],seq(0,1,length=nbasis))
  Y1[,j]<-Yf01-Zeta1f[j]$coefs*as.numeric(eta1[j])  
}

Yf1<-Data2fd(Y1,seq(0,1,length=nbasis),basisobj = basisx)
return(list(P1f,Zeta1f,Xf1,Yf1,eta1))

