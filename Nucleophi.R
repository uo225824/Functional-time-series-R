Nucleophi<-function(fdx,eigenvalues,scoresX,harmonicsX,m,base,N,p){
  
  rho<-0
  rho1<-0
  if (p>1) {
    L3<-0
    for (k in 1:base) {
      L2<-0

      for (j in p:N) {
        L1<-0
        for (n in 0:(p-1)) {
          L<-(t(eval.fd(evalarg=0:m/m, fdx[j-n]))%*%(eval.fd(evalarg=0:m/m, harmonicsX[k])))/length(0:m/m)
          L1<-L1+L
        }
        L2<-L2+L1^2
      }
      L3[k]<-L2/N
    }
    for (l in 1:base) {
      for (i in p:(N-1)) {
        for (j in 1:base) {
          for (n in 0:(p-1)) {
            rho<-as.numeric((1/L3[j])*scoresX[i-n,j]*scoresX[i+1-n,l])*eval.fd(evalarg=0:m/m, harmonicsX[j])%*%t(eval.fd(evalarg=0:m/m, harmonicsX[l]))
            rho1<-rho1+rho
            }

        }
      }
      print(l)
    }
    rho1[[1]]<-rho1[[1]]/(N-1)
    rho1[[2]]<-rho1[[2]]/(N-1)
    return(rho1)
    
  }else{
    for (l in 1:base) {
      for (i in 1:(N-1)) {
        for (j in 1:base) {
          rho<-as.numeric((1/eigenvalues[j])*scoresX[i,j]*scoresX[i+1,l])*eval.fd(evalarg=0:m/m, harmonicsX[j])%*%t(eval.fd(evalarg=0:m/m, harmonicsX[l]))
          rho1<-rho+rho1
        }
      }
      print(l)
    }
    rho1<-rho1/(N-1)
    return(rho1)
  }
  
}
base<-3
rho1<-Nucleophi(fdx,eigenvalues,scoresX,harmonicsX,m,base,N,1)
rho2<-Nucleophi(fdx,eigenvalues,scoresX,harmonicsX,m,base,N,1)


persp((0:49/49),(0:49/49),z=rho1,cex.axis = axisticksize,
      cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,
      phi=30,
      ticktype="detailed", main="P=3",col = "red")

L1<-0
for (j in 1:N) {
  L<-(t(eval.fd(evalarg=0:m/m, fdx[j]))%*%(eval.fd(evalarg=0:m/m, harmonicsX[3]))/50)^2
  L1<-L1+L
}
L1/N
