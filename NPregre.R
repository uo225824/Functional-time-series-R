NPregre(X12f,X12fT[j],Y12,3)
X<-X1f
r2<-X1fT[1]
Y<-Y1
knn<-3
NPregre <- function(X,r2,Y,knn) {
  
  d<-0
  K<-0
  n<-ncol(X$coefs)
  for (i in 1:n) {
    d[i]<-mean((r2$coefs-X[i]$coefs)^2)^1/2

  }
  
  h<-sort(d)[knn]
  
  for (i in 1:n) {
    K[i]<-W(d[i],h)
  }
  w<-K/sum(K)

  
  X1<-matrix(0,ncol = 1,nrow = nrow(X$coefs))
  for (i in 1:(n-1)) {
    X1<-X1+w[i]*Y[i]$coefs
  }
  
  return(X1)
}

W <- function(d, h) {
  k <- (1 / (2 * pi) ^ (1 / 2)) * exp(-((d) / h) ^ 2 / 2) 
  return(k)
}


