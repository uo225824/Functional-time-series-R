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
y<-0
for (i in 1:(ncol = ncol(serie.r)-1)) {
  x[,i]<-serie.r[,i]
  y[i]<-mean(serie.r[,i+1])
}

#Libreria principal del analisis de datos funcionales.

library(fda.usc)

#Base optima
opt.bsX<-optim.basis(fdata(t(x)),numbasis = seq(4,30,2),type.basis = "bspline")

#Número de elementos de la base óptima 
nbx<-opt.bsX$numbasis.opt

#Creamos las bases y formamos nuestro set de datos funcionales.  
basisx <- create.bspline.basis(nbasis=nbx)
datosf<-fdata(t(x),seq(0,1,length=30))

#Modelo de regresión con respuesta escalar

Y<-y[1:150]
X=datosf[1:150]

dataf=as.data.frame(Y)
ldata=list(df=dataf,X=datosf[1:150])
b.x=list(X=basisx)

res.pls=fregre.pls(X,Y,l=3)
dataf1=as.data.frame(y)
ldatanew=list(df=dataf1,X=datosf[150:180])
prepls<-inprod.fdata(datosf[150:180],res.pls$beta.est)

residuospls<-y[150:180]-prepls
#Representación Grafica
plot(prepls,lwd=2,lty=1,main="Ratio Dolar/Euro ",xlab="Tiempo (días)", ylab="Ratio")
points(y[150:180],col="red")
points(y[1:150],col="blue")
legend(1,1.22 , legend=c("Real", "Estimado"),
       col=c("black", "red"), pch=c(1,1), cex=0.8)



lower<-pre-qnorm(0.95,mean = mean(residuos),sd = sd(residuos))
upper<-pre+qnorm(0.95,mean = mean(residuos),sd = sd(residuos))


#Intervalos por wildbootstrap
library(fANCOVA)

#replicas bootstrap
B<-10000
remuestrapls<-wild.boot(residuospls, nboot = B)

alfa <- 0.01
pto_critpls <- apply(remuestrapls, 1, quantile, probs = c(alfa/2, 1 - alfa/2))
upper_bootpls <- prepls + pto_critpls[2, ]
lower_bootpls <- prepls + pto_critpls[1, ]

require(plotrix)
plotCI(seq(1,31,length.out = 31), prepls, ui=as.numeric(upper_bootpls), li=as.numeric(lower_bootpls),main="Ratio Dolar/Euro",xlab="Mes", ylab="Ratio")
points(y[150:180],pch = 16)
points(prepls,col="red",pch = 16)
legend(1,1.25 , legend=c("Real", "Estimado"),
       col=c("black", "red"), pch=c(16,16), cex=0.8)



for (i in 1 : 15){a[i]<-Box.test(residuospls,lag=i,type="Ljung-Box")$p.value}


plot(a,main="Test de Ljung-Box para los primeros 15 retardos",xlab="Retardo", ylab="P-valor", ylim=c(0,1))
abline(h=0.05,col="red")

plot(residuospls,xlab="",ylab="",main = "Residuos")
abline(h=0)
