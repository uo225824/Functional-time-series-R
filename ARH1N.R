library(fda.usc)
pc<-30
t<-30
ttraining<-150
serie.f<-na.omit(as.numeric(datos))
serie<-t(serie.f)
serie.fd<-diff(serie.f,lag=30)
#Modelo autorregresivo en espacios de Hilbert de orden 1
#serie.f= set de datos en forma de vector
#t= número de elementos del dato funcional
#pc= número de componentes principales a calcular
#ttraining= número de componentes del set de datos

ARH1N<-function(serie.f,t,pc,ttraining){
  m<-t-1
  
  #Reescribir los datos en forma de matriz
  serie<-t(matrix(serie.fd,nrow=t))
  
  #Eliminamos los elementos repetidos de la matriz
  #que se generan cuando t no es multiplo del número de datos
  N1<-2
  if (length(serie.fd)%%t==0) {
    N1<-1
  }
  
  #Generamos nuestros datos funcionales
  fdx<-Data2fd(t(serie[1:(nrow(serie)-N1+1),]),argvals = 0:m/m)

  #Determinamos las componentes principales
  fdXpca=pca.fd(fdx,pc)
  eigenvalues=fdXpca$values; scoresX=fdXpca$scores
  harmonicsX=fdXpca$harmonics
  N<-nrow(serie)-N1
  media<-mean.fd(fdx[1:N])
  
  #Calculamos el número de componentes principales que
  #explican el 0.99% de la varianza
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
  #Función nucleo
  rho1<-Nucleophi(fdx[1:ttraing],eigenvalues,scoresX,harmonicsX,m,nbasis,N)

  #Guardamos nuestros datos funcionales en forma de matriz
  xn<-matrix(0,nrow = t,ncol = N)
  for (l in 1:N) {
    XN<-Data2fd(rho1%*%eval.fd(evalarg = 0:m/m,fdx[l]-media)/t,argvals = 0:m/m)+media
    xn[,l]<-eval.fd(XN,seq(0,1,length=t))
    
  }
  return(xn)
}

#Representación Grafica
plot.ts(as.vector(serie.r[,2:180]),lwd=2,lty=1,main="Ratio Dolar/Euro ",xlab="Tiempo (días)", ylab="Ratio")
lines(as.vector(xn[,1:179]),col="red")
lines(as.vector(xn[,1:150]),col="blue")
legend(1,1.6 , legend=c("Test", "Training"),
       col=c("red", "blue"), lty=2:1, cex=0.8)


plot.ts(serie.fd,lwd=2,lty=1,main="Ratio Dolar/Euro ",xlab="Tiempo (días)", ylab="Ratio")
lines(as.vector(xn[,1:179]),col="red")
lines(as.vector(xn[,1:150]),col="blue")
legend(1,1.6 , legend=c("Test", "Training"),
       col=c("red", "blue"), lty=2:1, cex=0.8)


plot.ts(diffinv(serie.fd,lag = 30),lwd=2,lty=1,main="Ratio Dolar/Euro ",xlab="Tiempo (días)", ylab="Ratio")
lines(diffinv(as.vector(xn[,1:179]),lag=30),col="red")
lines(as.vector(xn[,1:150]),col="blue")
legend(1,1.6 , legend=c("Test", "Training"),
       col=c("red", "blue"), lty=2:1, cex=0.8)