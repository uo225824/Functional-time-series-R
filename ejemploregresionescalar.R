library(fda.usc)
t = seq(0, 1, length = nt <- 101)
covexp = function(t1, t2) {
  3 * exp(-abs(t1 - t2)/0.5)
}
Sigma = outer(t, t, covexp)
X = rproc2fdata(n <- 200, t, sigma = Sigma)
plot(X,lwd=2,lty=1,main="Procesos gaussianos",xlab="t", ylab="X(t)")
betaf = t + log(t + 0.1)
betaf = fdata(betaf, argvals = t) #Theo. Beta
plot(betaf,lwd=2,lty=1,col="blue",main="Beta",xlab="t", ylab="Beta(t)")
y = drop(inprod.fdata(X, betaf)) + rnorm(n, sd = 0.5) # Simulated response
opt.bs<-optim.basis(X,type.basis = "bspline")


dataf=as.data.frame(y)
ldata=list(df=dataf,abor=X)
b.x=list(abor=create.bspline.basis(nbasis=opt.bs$numbasis.opt))
res.lm=fregre.lm(y~abor,ldata,basis.x=b.x)

plot(betaf,lwd=2,lty=1,col="blue",main="Modelo teórico vs estimado",xlab="t", ylab="Beta(t)")
lines(res.lm$beta.l[["abor"]],lwd=2,lty=2,col="red")


plot(opt.bs$numbasis,opt.bs$gcv,type="l")
