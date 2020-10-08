library(fds)
data("Electricityconsumption")
serie.f<-Electricityconsumption$y


Fx<-fdata(t(X[,1:99]))
Fx0<-fdata.cen(Fx)

ys<-rowMeans(t(X[,2:100]))

dataf=as.data.frame(ys)

ldata=list("df"=dataf,"x.d2"=Fx)
gls<-fdata2pls(Fx0$Xcen,ys,ncomp = 3)
P1<-gls$rotation$data

P1<-t(P1)

plot((as.numeric(gls$x[1,1])*gls$rotation[1]+as.numeric(gls$x[1,2])*gls$rotation[2]+as.numeric(gls$x[1,3])*gls$rotation[3]))
lines(Fx[1])
plot(Fx[1])




last=99
N=99
plot.ts(c(t(serie[(N-last+1):N,])),ylim=c(min(t(serie[(N-last):N,]))-0.5,0.5
                                   +max(t(serie[(N-last):N,]))),axes=F,xlab="",ylab="",lwd=2)
axis(2); axis(1,tick=F,labels=F); abline(h=0)
abline(v=seq(0,last*(m+1),by=m+1), lty=2); box()



plot.ts(c(t(estima[(N-last+1):N,])),ylim=c(min(t(estima[(N-last):N,]))-0.5,0.5
                                          +max(t(estima[(N-last):N,]))),axes=F,xlab="",ylab="",lwd=2)
axis(2); axis(1,tick=F,labels=F); abline(h=0)
abline(v=seq(0,last*(m+1),by=m+1), lty=2); box()



plot.ts(c((yest[,(N-last+1):N])),ylim=c(min(t(yest[,(N-last):N]))-0.5,0.5
                                           +max(t(yest[,(N-last):N]))),axes=F,xlab="",ylab="",lwd=2)
axis(2); axis(1,tick=F,labels=F); abline(h=0)
abline(v=seq(0,last*(m+1),by=m+1), lty=2); box()

lines(c(X[(N-last+1):N,]),ylim=c(min(X[(N-last):N,])-0.5,0.5
                                   +max(X[(N-last):N,])),axes=F,xlab="",ylab="",lwd=2,col="red")
axis(2); axis(1,tick=F,labels=F); abline(h=0)
abline(v=seq(0,last*(m+1),by=m+1), lty=2); box()











