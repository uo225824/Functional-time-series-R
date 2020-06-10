
serie.funcional<-datos
ttiempo<-30

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


#Base optima
opt.bsY<-optim.basis(fdata(t(y)),numbasis = seq(4,30,2),type.basis = "bspline")
opt.bsX<-optim.basis(fdata(t(x)),numbasis = seq(4,30,2),type.basis = "bspline")

#Número de elementos de la base óptima 
nbx<-opt.bsX$numbasis.opt
nby<-opt.bsY$numbasis.opt

#Creamos las bases y formamos nuestro set de datos funcionales.  
basisx <- create.bspline.basis(c(0, 1),nbasis = nbx,norder = 5 )
basisy<-  create.bspline.basis(c(0,1),nbasis = nby,norder = 5 )
X<-Data2fd((x),seq(0,1,length=30),basisobj = basisx)
Y<-Data2fd((y),seq(0,1,length=30),basisobj = basisy)



#Nucleo Gausiano

W <- function(d, h) {
  k <- (1 / (2 * pi) ^ (1 / 2)) * exp(-((d) / h) ^ 2 / 2) 
  return(k)
}

#Regresion noparametrica funcional

NPregreSo <- function(X,Y1,h) {
  
  d<-0
  K<-0
  n<-ncol(X$coefs)
  for (i in 1:n) {
    d[i]<-mean((Y1$coefs-X[i]$coefs)^2)^1/2
    
  }
  
  for (i in 1:n) {
    K[i]<-W(d[i],h)
  }
  w<-K/sum(K)
  
  
  X1<-matrix(0,ncol = 1,nrow = nrow(X$coefs))
  for (i in 1:(n-1)) {
    X1<-X1+w[i]*X[i+1]$coefs
  }
  
  return(X1)
}


#Determinanción de la ventana por CV

h<-0.01
CVf<-0
for (j in 1:10) {
  
  CV<-0
  h<-h/1.5
  print(h)
  for (i in 1:179) {
    npr<-NPregreSo(X[-i],Y[i],h)
    CV[i]<-(colMeans((Y[i]$coefs-npr)^2))
  }
  CVf[j]<-mean(CV)
  print(j)
}


#Modelo final no parametrico
h<-0.01/(1.5^which.min(CVf))

nprf<-matrix(0,nrow = 30,ncol = 179)
for (i in 1:179) {
  npr<-NPregreSo(X,Y[i],h)
  nprf[,i]<-eval.fd(Data2fd(npr),seq(0,1,length=ttiempo))
}


#Representación Grafica
plot.ts(as.vector(datos[1:(length(datos)-60)]),lwd=2,lty=1,main="Ratio Dolar/Euro ",xlab="Tiempo (días)", ylab="Ratio")
lines(as.vector(nprf),col="red")
lines(as.vector(nprf[,1:150]),col="blue")
legend(1,1.6 , legend=c("Test", "Training"),
       col=c("red", "blue"), lty=2:1, cex=0.8)
