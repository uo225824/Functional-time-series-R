#AR2 Final
serie.f<-x
t<-12
pc<-12
  serie<-t(matrix(serie.funcional,nrow=t))
  
  if (Tipo==1) {
    serie<-serie.f
    
  }
  N1<-2
  if (length(serie.funcional)%%t==0) {
    N1<-1
  }
  
  fdx<-Data2fd(t(serie[1:(nrow(serie)-N1+1),]),argvals = c(1:t))
  
  fdXpca=pca.fd(fdx,pc)
  eigenvalues=fdXpca$values; scoresX=fdXpca$scores
  harmonicsX=fdXpca$harmonics
  eigenvectors=eval.fd(evalarg = 1:t,harmonicsX)
  eigenvectors<-Data2fd(eigenvectors,argvals = c(1:t))
  
  

  l<-1
  for (i in 1:pc) {
    l<-l-fdXpca$varprop[i]
    if (l<0.01) {
      nbasis<-i
      break()
    }
  }
  print(nbasis)
  N<-nrow(serie)-N1
  rho<-0
  rho1<-0
  rho2<-0
  media<-mean.fd(fdx[1:N])

  for (l in 1:nbasis) {
    for (i in 1:(N-1)) {
      for (j in 1:nbasis ) {
        rho[j]<-(1/eigenvalues[j])*inprod(fdx[N]-media,eigenvectors[j])*inprod(fdx[i]-media,eigenvectors[j])*inprod(fdx[i+1]-media,eigenvectors[l])
      }
      rho1[i]<-sum(rho)
    }
    rho2[l]<-sum(rho1)/(N-1)
    print(l)
  }
  xn<-0
  xn1<-0
  for (l in 1:nbasis) {
    xn<-rho2[l]*eigenvectors[l]
    xn1<-xn1+xn
  }
  XN1<-xn1+media
  
  AA<-xn1+fdx[N+1]
  
  for (l in 1:nbasis) {
    for (i in 1:(N-2)) {
      for (j in 1:nbasis ) {
        rho[j]<-(1/eigenvalues[j])*inprod(fdx[N-1]-media,eigenvectors[j])*inprod(fdx[i]-media,eigenvectors[j])*inprod(fdx[i+1]-media,eigenvectors[l])
      }
      rho1[i]<-sum(rho)
    }
    rho2[l]<-sum(rho1)/(N-1)
    print(l)
  }
  xn<-0
  xn1<-0
  for (l in 1:nbasis) {
    xn<-rho2[l]*eigenvectors[l]
    xn1<-xn1+xn
  }
  Num<-inprod(AA,xn1)
  
  lambda<--Num/mean(xn1$coefs)
  
  
  aut<-matrix(0,nrow = (nrow(serie)-N1-1),ncol = pc)
  for (i in 2:(nrow(serie)-N1)) {
    for (j in 1:pc) {
      aut[i-1,j]<-inprod(fdx[i]+lambda*fdx[i-1],eigenvectors[j])*inprod(fdx[i],eigenvectors[j])
      
    }
  }
  autovalores<-colMeans(aut)
  
  la<-1
  for (i in 1:pc) {
    la<-la-autovalores[i]/sum(autovalores)
    if (la<0.01) {
      nbasisa<-i
      break()
    }
  }
  print(nbasisa)
  for (l in 1:nbasisa) {
    for (i in 1:(N-1)) {
      for (j in 1:nbasisa ) {
        rho[j]<-(1/autovalores[j])*inprod(fdx[N]-media,eigenvectors[j])*inprod(fdx[i]-media,eigenvectors[j])*inprod(fdx[i+1]-media,eigenvectors[l])
      }
      rho1[i]<-sum(rho)
    }
    rho2[l]<-sum(rho1)/(N-1)
    print(l)
  }
  xn<-0
  xn1<-0
  for (l in 1:nbasisa) {
    xn<-rho2[l]*eigenvectors[l]
    xn1<-xn1+xn
  }
  XNAR2<-xn1+media
  
  par(mfrow=c(2,2))
  plot(fdx[N+1],col="blue")
  plot(XN1,col="red")
  plot(XNAR2,col="green")
  
  plot(fdx[N+1],col="blue")
  lines(XN1,col="red")
  lines(XNAR2,col="green")
  
  dis1<-sqrt(inprod(XN1-fdx[N+1],XN1-fdx[N+1]))
  dis2<-sqrt(inprod(XNAR2-fdx[N+1],XNAR2-fdx[N+1]))
  dis<-data.frame(c(dis1[1],dis2[1]))

  