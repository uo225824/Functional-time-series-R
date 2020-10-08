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

dataf=as.data.frame(Y)
ldata=list(df=dataf,X=datosf[1:150])
b.x=list(X=basisx)

res.lm=fregre.lm(Y~X,ldata,basis.x=b.x)
dataf1=as.data.frame(y)
ldatanew=list(df=dataf1,X=datosf[150:180])
pre<-predict(res.lm,ldatanew)
prec<-predict(res.lm,ldatanew,interval="prediction")

residuos<-y[150:180]-pre
#Representación Grafica
plot(pre[,1],lwd=2,lty=1,main="Ratio Dolar/Euro ",xlab="Tiempo (días)", ylab="Ratio")
points(y,col="red")
points(y[1:150],col="blue")
legend(1,1.6 , legend=c("Test", "Training"),
       col=c("red", "blue"), lty=2:1, cex=0.8)



lower<-pre-qnorm(0.95,mean = mean(residuos),sd = sd(residuos))
upper<-pre+qnorm(0.95,mean = mean(residuos),sd = sd(residuos))

lines(upper)
lines(lower)

new.lstat = seq(min(1), max(31), length.out=31)
plot(y[150:180])
points( as.vector(pre), col='red', lwd=3) 
points(as.vector(pre[1:150]),col="blue")

polygon(c(rev(new.lstat), new.lstat), c(rev(upper), lower), density=15, col = 'blue', border = NA)
lines(new.lstat, upper, lty = 'dashed',col="red")
lines(new.lstat, lower, lty = 'dashed',col="red")

#Intervalos por wildbootstrap
library(fANCOVA)

#replicas bootstrap
B<-10000
remuestra<-wild.boot(residuos, nboot = B)

alfa <- 0.01
pto_crit <- apply(remuestra, 1, quantile, probs = c(alfa/2, 1 - alfa/2))
upper_boot <- pre + pto_crit[2, ]
lower_boot <- pre + pto_crit[1, ]


new.lstat = seq(min(1), max(31), length.out=31)
plot(y[150:180],main="Ratio Dolar/Euro",xlab="Mes", ylab="Ratio")
points( as.vector(pre), col='red', lwd=3) 
legend(1,1.22 , legend=c("Real", "Estimado"),
       col=c("black", "red"), pch=c(1,1), cex=0.8)


polygon(c(rev(new.lstat), new.lstat), c(rev(upper_boot), lower_boot), density=15, col = 'blue', border = NA)
lines(new.lstat, upper_boot, lty = 'dashed',col="red")
lines(new.lstat, lower_boot, lty = 'dashed',col="red")

require(plotrix)
plotCI(seq(1,31,length.out = 31), pre, ui=as.numeric(upper_boot), li=as.numeric(lower_boot),main="Ratio Dolar/Euro",xlab="Mes", ylab="Ratio")
points(y[150:180],pch = 16)
points(pre,col="red",pch = 16)
legend(1,1.25 , legend=c("Real", "Estimado"),
       col=c("black", "red"), pch=c(16,16), cex=0.8)

Box.test(residuos,lag=1,type="Ljung-Box")

a<-0
for (i in 1 : 15){a[i]<-Box.test(residuos,lag=i,type="Ljung-Box")$p.value}

plot(a,main="Test de Ljung-Box para los primeros 15 retardos",xlab="Retardo", ylab="P-valor", ylim=c(0,1))
abline(h=0.05,col="red")

plot(residuos,xlab="",ylab="",main = "Residuos")
abline(h=0)
