#Modelo reescalado


ttiempo<-12
ndatos<-15
MRescala<-function(serie,ttiempo,ndatos){
  Ar<-ARH1(serie,ttiempo,12)
  serie1<-t(matrix(serie,nrow=ttiempo))
  N1<-2
  if (length(serie)%%ttiempo==0) {
    N1<-1
  }
  d<-(nrow(serie1)-ndatos-1-N1)
  df<-0
  for (i in d:(nrow(serie1)-1-N1)) {
    RoYt<-ARH1(serie[1:(ttiempo*i)],ttiempo,12)
    
    e<-Data2fd(serie1[i+1,],argvals = seq(0,ttiempo,length=length(RoYt$coefs)-2))$coef
    A<-matrix(0,nrow = length(e),ncol = 2)
    A[,2]<-RoYt$coefs
    A[,1]<-e
    e[which.min(e)]-RoYt$coefs[which.max(RoYt$coefs)]
    myfdata<-fdata(t(A))
    df[i-d+1]= metric.lp(myfdata, lp = 2)[2]
    print(i)
  }
  
  basisx <- create.fourier.basis(c(0, ttiempo), 5)
  basisy<- create.fourier.basis(c(0, ttiempo), 21)
  Ycon<-df
  Xcon<-Data2fd(t(serie1[d:(nrow(serie1)-1-N1),]),argvals = seq(0,ttiempo,length=ncol(serie1)),basisobj = basisx)
  res.pc = fregre.pc(Xcon, Ycon)
  cor<-predict(res.pc,new.fdataobj =Data2fd(serie1[(nrow(serie1)-N1),],basisobj = basisx))
  lines(fdata(t(cor+ Ar$coefs)))
  plot(fdata(t(cor+ Ar$coefs)))
}


MRescala(serie.funcional,ttiempo = 12,15)
