#Regresion con respuesta funcional
library(fda.usc)
library(tensor)
serie.funcional<-X1
ttiempo<-12
regreF<-function(serie.funcional,ttiempo){
  serie<-t(matrix(serie.funcional,nrow=ttiempo))
  serie
  N1<-2
  if (length(serie.funcional)%%ttiempo==0) {
    N1<-1
  }
  serie.r<-t(serie[1:(nrow(serie)-N1+1),])
  
  x<-matrix(0,ttiempo,ncol = ncol(serie.r)-2)
  y<-matrix(0,ttiempo,ncol = ncol(serie.r)-2)
  for (i in 1:(ncol = ncol(serie.r)-2)) {
    x[,i]<-serie.r[,i]
    y[,i]<-serie.r[,i+1]
  }
  
  fx<-Data2fd(x,seq(1,ttiempo))
  fy<-Data2fd(y,seq(1,ttiempo))
  mediax<-mean.fd(fx,seq(1,ttiempo,length=ttiempo))
  mediay<-mean.fd(fy,seq(1,ttiempo,length=ttiempo))
  FDX<-matrix(0,ttiempo,ncol = ncol(serie.r)-2)
  FDY<-matrix(0,ttiempo,ncol = ncol(serie.r)-2)
  for (i in 1:(ncol = ncol(serie.r)-2)) {
    Fx<-eval.fd(fx[i],seq(1,ttiempo,length=ttiempo))-eval.fd(mediax,seq(1,ttiempo,length=ttiempo))
    Fy<-eval.fd(fy[i],seq(1,ttiempo,length=ttiempo))-eval.fd(mediay,seq(1,ttiempo,length=ttiempo))
    FDX[,i]<-Fx
    FDY[,i]<-Fy
  }
  basisx <- create.bspline.basis(c(1, ttiempo),nbasis = 31,norder = 5 )
  basisy<-  create.bspline.basis(c(1, ttiempo),nbasis = 31,norder = 5 )
  datosf<-Data2fd(FDX,seq(1,ttiempo,length=ttiempo),basisobj = basisx)
  datosfN<-Data2fd(FDY,seq(1,ttiempo,length=ttiempo),basisobj = basisy)
  
  a<-Data2fd(serie[nrow(serie)+1-N1,],argvals = seq(1,ttiempo,length=ttiempo))
  lambdas<-1E-5
  lambdat<-1E-5
  
  CV<-matrix(0,ncol = 6,nrow = 6)
  for (i in 1:6) {
    for (j in 1:6) {
          res.fr = fregre.basis.fr(datosf,datosfN,basisx,basisy,lambda.s = lambdas,lambda.t = lambdat,Lfdobj.s = 3,Lfdobj.t = 3)
          lll<-predict(res.fr,Data2fd(serie[nrow(serie)-N1,],argvals = seq(1,ttiempo,length=ttiempo),basisobj = basisx))
          CV[i,j]<-t((a-lll)$coefs)%*%(a-lll)$coefs/24
    }
    print(i)
  }
  res.fr = fregre.basis.fr(datosf,datosfN,basisx,basisy,lambda.s = 1E-5,lambda.t = 1E-5,Lfdobj.s = 2,Lfdobj.t = 2)
  
  
  plot(a)
  lines(lll)
  
  lines(predict(res.fr,Data2fd(serie[nrow(serie)-N1,],argvals = seq(1,ttiempo,length=ttiempo),basisobj = basisx)),col="red")
  
  plot(fdata.deriv(lll))
  plot(a)
}

plot(predict(res.fr,Data2fd(serie.r[,11],argvals = seq(1,10,length=10),basisobj = basisx)),col="red")


library(plsgenomics)

pls<-pls.regression(X[,1:99],X[,2:100])
AA<-pls$B
AA<-AA[1:4,]
