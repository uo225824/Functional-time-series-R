#Modelo Conjunto
ttiempo<-12
ndatos<-5
serie<-serie.funcional
MConjunto<-function(serie,ttiempo,ndatos){
  Ar<-ARH1(serie,ttiempo,12,0)
  df<-matrix(0,nrow = length(Ar$coefs),ncol=ndatos+1)
  serie1<-t(matrix(serie,nrow=ttiempo))
  N1<-2
  if (length(serie)%%ttiempo==0) {
    N1<-1
  }
  d<-(nrow(serie1)-ndatos-1-N1)
  for (i in d:(nrow(serie1)-1-N1)) {
    RoYt<-ARH1(serie[1:(ttiempo*i)],ttiempo,12,0)
    df[,i-(d-1)]<-Data2fd(serie1[i+1,],argvals = seq(0,ttiempo,length=nrow(df)-2))$coef-RoYt$coefs
    print(i)
  }
  basisx <- create.fourier.basis(c(0, ttiempo), 5)
  basisy<- create.fourier.basis(c(0, ttiempo), 51)
  Ycon<-Data2fd(df,argvals = seq(0,ttiempo,length=nrow(df)),basisobj = basisy)
  Xcon<-Data2fd(t(serie1[d:(nrow(serie1)-1-N1),]),argvals = seq(0,ttiempo,length=ncol(serie1)),basisobj = basisx)
  res.fr = fregre.basis.fr(Xcon,Ycon,basisx,basisy)
  cor<-predict(res.fr,new.fdataobj =Data2fd(t(serie1),basisobj = basisx))
  lines(fdata(t(eval.fd(seq(0,ttiempo,length.out = length(Ar$coefs)),cor)+ Ar$coefs)))
  plot(Data2fd((eval.fd(seq(0,ttiempo,length.out = length(Ar$coefs)),cor[28])+ Ar$coefs),basisobj=basisy))
  return(fdata(t(eval.fd(seq(0,ttiempo,length.out = length(Ar$coefs)),cor[33])+ Ar$coefs)))
}

corf<-MConjunto(serie.funcional,ttiempo,10)
