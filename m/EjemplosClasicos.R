serie.funcional <- ts(scan(file ="series de tiempo/data.csv"), start=c(1972), frequency=12)
datos<-read.csv2(file = "elec/energy_dataset.csv",header = T,sep = ",",dec = ".")
datos2<-read.csv2(file = "elec/weather_features.csv",header = T,sep = ",",dec = ".")
serie.funcional<-(as.numeric(as.vector(datos$price.actual)))
serie.funcional2<-(as.numeric(as.vector(datos2$temp)))
X<-na.omit(serie.funcional)
X2<-na.omit(serie.funcional2)
plot.ts(serie.funcional)

ttiempo<-6



ARH1(serie.funcional,12,12,0)
ARH1A(serie.f,21,21,0)
ARH2A(serie.f,21,21,0)
ARH1N(X,50,3,0,1)
regreF(serie.funcional,ttiempo = 21)

p<-ARMAH1(serie.funcional,12,12,0)

dserie.funcional<-diff(serie.funcional)
plot.ts(dserie.funcional)
ddserie.funcional<-diff(dserie.funcional)
plot.ts(ddserie.funcional)
p<-ARMAH1(dserie.funcional,12,12,0)
p<-ARH2(serie.funcional,12,12,0)

p<-ARH2(dserie.funcional,12,12,0)
p
plot.ts(dserie.funcional)
Mo<-ARH1(dserie.funcional,12,12,0)


regreF(dserie.funcional,ttiempo = 12)

MConjunto(dserie.funcional,12,5)


lines(Data2fd((eval.fd(seq(0,ttiempo,length.out = length(Ar$coefs)),cor[28])+ Ar$coefs),basisobj = basisy))

eval.fd(seq(1,12,length.out = 12),Mo)

prueba<-dserie.funcional
prueba[(length(prueba)-11):length(prueba)]=ts(as.vector(t(eval.fd(seq(1,12,length.out = 12),Mo))))

pr<-diffinv(prueba)
plot(Data2fd(pr[(length(pr)-11):length(pr)]),basisobj = basisy)

prf<-pr+315.71
