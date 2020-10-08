#FMAE
#ARH(1)
FMAEarh1<-mean(colMeans(abs(serie.r[,151:181]-xn[,150:180])))

#Red neuronal
FMAEred<-mean(colMeans(abs(serie.r[,151:181]-yest[,150:180])))

#Parametrico
FMAEpara<-mean(colMeans(abs(serie.r[,151:181]-pr[,150:180])))

#No parametrico
FMAEnopara<-mean(colMeans(abs(serie.r[,151:181]-nprf[,150:180])))




#FMSE

#ARH(1)
FMSEarh1<-mean(colMeans((serie.r[,151:181]-xn[,150:180])^2))

#Red neuronal
FMSEred<-mean(colMeans((serie.r[,151:181]-yest[,150:180])^2))

#Parametrico
FMSEpara<-mean(colMeans((serie.r[,151:181]-pr[,150:180])^2))

#No parametrico
FMSEnopara<-mean(colMeans((serie.r[,151:181]-nprf[,150:180])^2))


#FRMSE


#ARH(1)
FMRSEarh1<-sqrt(mean(colMeans((serie.r[,151:181]-xn[,150:180])^2)))

#Red neuronal
FMRSEred<-sqrt(mean(colMeans((serie.r[,151:181]-yest[,150:180])^2)))

#Parametrico
FMRSEpara<-sqrt(mean(colMeans((serie.r[,151:181]-pr[,150:180])^2)))

#No parametrico
FMRSEnopara<-sqrt(mean(colMeans((serie.r[,151:181]-nprf[,150:180])^2)))
