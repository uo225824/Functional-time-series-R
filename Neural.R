


# activation function
# sigmoid
sigmoid <- function(x) return(1/(1+exp(-x)))
d.sigmoid <- function(x) return(exp(-x)/(1+exp(-x))^2)


sigmoid <- function(x) return(tanh(x))
d.sigmoid <- function(x) return(1-tanh(x)^2)

set.seed(123456)
W1 <- matrix(runif(4*hl[1], 0, 0.1), nrow = 4)
W2<-W1
u<-0
s<-seq(0,1,length.out = 50)
t<-seq(0,1,length.out = 50)
hl<-5
rate<-1E-6

for (l in 1:1000) {
  
P<-matrix(0,nrow = 50,ncol = 50)
Pf<-matrix(0,nrow = 50,ncol = 50)

for (k in 1:hl) {
  for (i in 1:50) {
    for (j in 1:50) {
      P[i,j]<-W1[4,k]*sigmoid(W1[1,k]+W1[2,k]*s[i]+W1[3,k]*t[j])
    }
  }
  Pf<-Pf+P

}


yest<-t(Pf)%*%X[,1:99]/(50)
er<-0

for (i in 1:99) {
  er[i]<-t(X[,i+1]-yest[,i])%*%(X[,i+1]-yest[,i])/50
}

u[l]<-mean(er)

if (u[l]==u[which.min(u)]) {
  W3<-W1
}



Pd<-matrix(0,nrow = 50,ncol = 50)
Pdf<-list(Pd)

for (k in 1:hl) {
  for (i in 1:50) {
    for (j in 1:50) {
      Pd[i,j]<-W1[4,k]*d.sigmoid(W1[1,k]+W1[2,k]*s[i]+W1[3,k]*t[j])
    }
  }
  Pdf[[k]]<-Pd
}


STB<-matrix(1,nrow = 50,ncol = 4)
STB[,2]<-s
STB[,3]<-t
grad<-matrix(ncol = hl,nrow = 4)
L<-0
for (k in 1:hl) {
  for (j in 1:4) {
    for (i in 1:99) {
      if (j==1) {
        L[i]<-t((X[,i+1]-yest[,i]))%*%Pdf[[k]]%*%hadamard.prod(STB[,j],X[,i+1])/250
      }
      if (j==2) {
        L[i]<-t((X[,i+1]-yest[,i]))%*%Pdf[[k]]%*%hadamard.prod(STB[,j],X[,i+1])/250
      }
      if (j==3) {
        L[i]<-t(hadamard.prod((X[,i+1]-yest[,i]),STB[,j]))%*%Pdf[[k]]%*%X[,i+1]/250
        
      }
      if (j==4) {
        L[i]<-t((X[,i+1]-yest[,i]))%*%Pf%*%hadamard.prod(STB[,j],X[,i+1])/250}
    if (j==5) {
      L[i]<-mean(X[,i+1]-yest[,i])*mean(X[,i+1])}
    }
    grad[j,k]<-sum(L)
}}


grad

W1<-W1+rate*grad

if (l %% 100==0) {
  plot(u)
  persp((0:49/49),(0:49/49),z=t(Pf),cex.axis = axisticksize,
        cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,
        phi=30,
        ticktype="detailed", main="P=3",col = "blue")
}
print(l)

}


plot(u)
persp((0:49/49),(0:49/49),z=t(Pf),cex.axis = axisticksize,
      cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,
      phi=30,
      ticktype="detailed", main="P=3",col = "blue")

XN<-Data2fd(Pfr%*%eval.fd(evalarg = 0:m/m,fdx[99]-media)/t,argvals = 0:m/m)+media


estima<-matrix(0,nrow =nrow(fdx$coefs)-2 ,ncol = ncol(fdx$coefs))
for (j  in 1:ncol(fdx$coefs)) {
  XN<-Data2fd(Pfr%*%eval.fd(evalarg = 0:m/m,fdx[j]-media)/(m+1),argvals = 0:m/m)+media
  estima[,j]<-eval.fd(evalarg = 0:m/m,XN)
  print(j)
} 
