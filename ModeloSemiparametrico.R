pc<-50
tt<-50
serie.f<-X
serie.f2<-X2
Tipo=0
p<-0
SemiPara(serie.f,t,pc,0,0)

SemiPara<-function(serie.f,tt,pc,Tipo,p){
  m<-tt-1
serie<-t(matrix(serie.f,nrow=tt))

if (Tipo==1) {
  serie<-t(serie.f)
  
}
N1<-2
if (length(serie.f)%%tt==0) {
  N1<-1
}
serie<-serie[1:100,]

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
nbasis<-1
p<-0
rho1<-Nucleophi(fdx,eigenvalues,scoresX,harmonicsX,m,nbasis,N,p)

estima<-matrix(0,nrow =nrow(fdx$coefs)-2 ,ncol = nrow(serie))
for (j  in 2:ncol(fdx$coefs)) {
    XN<-Data2fd(rho1%*%eval.fd(evalarg = 0:m/m,fdx[j-1]-media)/tt,argvals = 0:m/m)+media
  estima[,j]<-eval.fd(evalarg = 0:m/m,XN)
  print(j)
} 

for (j  in (ncol(fdx$coefs)+1):170) {
  XN<-Data2fd(rho1%*%eval.fd(evalarg = 0:m/m,XN-media)/tt,argvals = 0:m/m)+media
  estima[,j]<-eval.fd(evalarg = 0:m/m,XN)
  print(j)
} 

Y1<-Data2fd(t(serie[2:170,])-estima[,1:169],argvals = 0:m/m)
X1f<-Data2fd(estima,argvals = 0:m/m)
X1fT<-Data2fd(t(serie),argvals = 0:m/m)
X2<-matrix(0,nrow = nrow(X1fT$coefs),ncol = ncol(X1fT$coefs))
for (j in 1:ncol(X1fT$coefs)) {
  X2[,j]<-NPregre(X1f,X1fT[j],Y1,3)
}
X2f<-Data2fd(X2,argvals = 0:(nrow(X2)-1)/(nrow(X2)-1))
estima2<-eval.fd(evalarg = 0:m/m,X2f)

estimaf<-estima+estima2
return(estimaf)
}

last=108
N=108
plot.ts(c(t(serie[(N-last+1):N,])),ylim=c(min(t(serie[(N-last):N,]))-0.5,0.5
                                          +max(t(serie[(N-last):N,]))),axes=F,xlab="",ylab="",lwd=2)
axis(2); axis(1,tick=F,labels=F); abline(h=0)
abline(v=seq(0,last*(m+1),by=m+1), lty=2); 
box()


lines(c(yest[,(N-last+1):N]),col="red")

plot.ts(as.vector(yest))

lines(c(estima[,(N-last+1):N]),col="blue", lty=3)
lines(c(estimaf[,(N-last+1):N]),col="red", lty=2)
