library(fds)
data("Electricityconsumption")
serie.f<-as.vector(Electricityconsumption$y)

rho<-0
rho1<-0
P1<-plsteo$rotation$data
eigenvectors<-Data2fd(t(P1),argvals = 0:49/49)
#plot(eigenvectors)
base<-3
  L3<-0
L1<-0
for (k in 1:base) {
  for (i in 1:N) {
    L<-(inprod(fdx[i],eigenvectors[k]))
    L1<-L1+L^2
  }
 L3[k]<-L1/N
 L1<-0
}
  
L3<-sort(L3,decreasing = T)

L3<-L3/20  
base<-3
for (l in 1:base) {
  for (i in 1:(N-1)) {
    for (j in 1:base) {
      rho<-as.numeric((1/L3[j])*(t(eval.fd(evalarg=0:m/m, fdx[i]))%*%(eval.fd(evalarg=0:m/m, eigenvectors[j])))/length(0:m/m)*(t(eval.fd(evalarg=0:m/m, fdx[i+1]))%*%(eval.fd(evalarg=0:m/m, eigenvectors[l])))/length(0:m/m))*eval.fd(evalarg=0:m/m, eigenvectors[j])%*%t(eval.fd(evalarg=0:m/m, eigenvectors[l]))
      rho1<-rho+rho1
    }
  }
  print(l)
}
rho1<-rho1/(N-1)


persp((0:49/49),(0:49/49),z=t(rho1),cex.axis = axisticksize,
      cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,
      phi=30,
      ticktype="detailed", main="P=3",col = "blue")




t(eval.fd(evalarg=0:m/m, fdx[1]))%*%(eval.fd(evalarg=0:m/m, eigenvectors[2]))/50
       
inprod(fdx[1],harmonicsX[2])
scoresX[1,2]


m <- matrix(1:8,ncol=4)
#      [,1] [,2] [,3] [,4]
# [1,]    1    3    5    7
# [2,]    2    4    6    8


rot <- function(x) "[<-"(x, , rev(x))

rot(Pf)

Pfr<-(rot(Pf))
