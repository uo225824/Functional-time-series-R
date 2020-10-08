# Función de activación Sigmoide
sigmoid <- function(x) return(1/(1+exp(-x)))
#Derivada de la función sigmoide
d.sigmoid <- function(x) return(exp(-x)/(1+exp(-x))^2)

set.seed(1234)

#Inicialización de los pesos
hl<-5
W1 <- matrix(runif(4*hl[1], 0, 1), nrow = 4)
W2<-W1
#Discretización de la componente temporal
s<-seq(0,1,length.out = 30)
t<-seq(0,1,length.out = 30)

u<-0

#Ratio de aprendizaje
rate<-1E-2
X<-(serie.r)

#Red neuronal
for (l in 1:5000) {
  
#Matriz de pesos  
P<-matrix(0,nrow = 30,ncol = 30)
Pf<-matrix(0,nrow = 30,ncol = 30)
for (k in 1:hl) {
  for (i in 1:30) {
    for (j in 1:30) {
      P[i,j]<-W1[4,k]*sigmoid(W1[1,k]+W1[2,k]*s[i]+W1[3,k]*t[j])
    }
  }
  Pf<-Pf+P
}

#Estimacion
yest<-t(Pf)%*%X[,1:150]/(30)

#Calculo FMRSE
er<-0
for (i in 1:150) {
  er[i]<-t(X[,i+1]-yest[,i])%*%(X[,i+1]-yest[,i])/30
}
u[l]<-mean(er)
if (u[l]==u[which.min(u)]) {
  W3<-W1
}

#Dencenso del gradiente
Pd<-matrix(0,nrow = 30,ncol = 30)
Pdf<-list(Pd)
for (k in 1:hl) {
  for (i in 1:30) {
    for (j in 1:30) {
      Pd[i,j]<-W1[4,k]*d.sigmoid(W1[1,k]+W1[2,k]*s[i]+W1[3,k]*t[j])
    }
  }
  Pdf[[k]]<-Pd
}

STB<-matrix(1,nrow = 30,ncol = 4)
STB[,2]<-s
STB[,3]<-t
grad<-matrix(ncol = hl,nrow = 4)
L<-0
for (k in 1:hl) {
  for (j in 1:4) {
    for (i in 1:150) {
      if (j==1) {
        L[i]<-t((X[,i+1]-yest[,i]))%*%Pdf[[k]]%*%hadamard.prod(STB[,j],X[,i+1])/150^2
      }
      if (j==2) {
        L[i]<-t((X[,i+1]-yest[,i]))%*%Pdf[[k]]%*%hadamard.prod(STB[,j],X[,i+1])/150^2
      }
      if (j==3) {
        L[i]<-t(hadamard.prod((X[,i+1]-yest[,i]),STB[,j]))%*%Pdf[[k]]%*%X[,i+1]/150^2
      }
      if (j==4) {
        L[i]<-t((X[,i+1]-yest[,i]))%*%Pf%*%hadamard.prod(STB[,j],X[,i+1])/150^2}
    if (j==5) {
      L[i]<-mean(X[,i+1]-yest[,i])*mean(X[,i+1])}
    }
    grad[j,k]<-sum(L)
}}

W1<-W1+rate*grad
if (l %% 100==0) {
  plot(u)
}
print(l)

}

uu<-u
#Plot training
plot.ts(as.vector(u[1:5000]),lwd=2,lty=1,main="Entrenamiento de la red neuronal",xlab="Número de iteraciones", ylab="FRMSE")

#Representación Grafica
yest<-t(Pf)%*%X/(30)
plot.ts(as.vector(serie.r[,2:180]),lwd=2,lty=1,main="Ratio Dolar/Euro ",xlab="Tiempo (días)", ylab="Ratio")
lines(as.vector(yest),col="red")
lines(as.vector(yest[,1:150]),col="blue")
legend(1,1.6 , legend=c("Test", "Training"),
       col=c("red", "blue"), lty=2:1, cex=0.8)
