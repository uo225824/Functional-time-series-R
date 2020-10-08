serie.funcional<-datos
ttiempo<-30

serie<-t(matrix(serie.funcional,nrow=ttiempo))
serie
N1<-2
if (length(serie.funcional)%%ttiempo==0) {
  N1<-1
}
serie.r<-t(serie[1:(nrow(serie)-N1+1),])

x<-matrix(0,ttiempo,ncol = ncol(serie.r)-1)
y<-matrix(0,ttiempo,ncol = ncol(serie.r)-1)
for (i in 1:(ncol = ncol(serie.r)-1)) {
  x[,i]<-serie.r[,i]
  y[,i]<-serie.r[,i+1]
}

#Libreria principal del analisis de datos funcionales.

library(fda.usc)

#Base optima
opt.bsY<-optim.basis(fdata(t(y)),numbasis = seq(4,30,2),type.basis = "bspline")
opt.bsX<-optim.basis(fdata(t(x)),numbasis = seq(4,30,2),type.basis = "bspline")

#Número de elementos de la base óptima 
nbx<-opt.bsX$numbasis.opt
nby<-opt.bsY$numbasis.opt

#Creamos las bases y formamos nuestro set de datos funcionales.  
basisx <- create.bspline.basis(c(0, 1),nbasis = nbx,norder = 5 )
basisy<-  create.bspline.basis(c(0,1),nbasis = nby,norder = 5 )
datosf<-Data2fd((x),seq(0,1,length=30),basisobj = basisx)
datosfN<-Data2fd((y),seq(0,1,length=30),basisobj = basisy)


#Inicializacion de los parametros de penalización
lambdas<-1
lambdat<-1

#Estimación rápida del orden de la penalización
FMSE<-0
for (i in 1:100) {
  lambdas<-lambdas/10
  lambdat<-lambdat/10
  res.fr = fregre.basis.fr(datosf[1:150],datosfN[1:150],basisx,basisy,lambda.s = lambdas,lambda.t = lambdat,Lfdobj.s = 1,Lfdobj.t = 1)
  prf<-predict(res.fr,datosf[151:179],basisobj = basisx)
  pr<-eval.fd(prf,seq(0,1,length=ttiempo))
  FMSE[i]<-mean(colMeans((y[,151:179]-pr)^2))
}

#FMSE mínimo
which.min(FMSE)

#Modelo final

lambdas<-1/10^which.min(FMSE)
lambdat<-1/10^which.min(FMSE)
res.fr = fregre.basis.fr(datosf[1:150],datosfN[1:150],basisx,basisy,lambda.s = lambdas,lambda.t = lambdat,Lfdobj.s = 1,Lfdobj.t = 1)
prf<-predict(res.fr,datosf,basisobj = basisx)
pr<-eval.fd(prf,seq(0,1,length=ttiempo))

#Representación Grafica
plot.ts(as.vector(serie.r[,2:180]),lwd=2,lty=1,main="Ratio Dolar/Euro ",xlab="Tiempo (días)", ylab="Ratio")
lines(as.vector(pr),col="red")
lines(as.vector(pr[,1:150]),col="blue")
legend(1,1.6 , legend=c("Test", "Training"),
       col=c("red", "blue"), lty=2:1, cex=0.8)
