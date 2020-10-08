
library(quantmod)

#############################################################################################################
# Import data from S&P500

getSymbols("SAN.MC",from="2005-01-01",to="2020-03-11")

SAN.prices <- SAN.MC$SAN.MC.Adjusted
XS<-as.vector(SAN.prices)
plot.ts(XS)
X<-log(XS)

dXS<-diff(X)

getSymbols("BBVA.MC",from="2005-01-01",to="2020-03-11")

BBVA.prices <- BBVA.MC$BBVA.MC.Adjusted
XB<-as.vector(BBVA.prices)

dXB<-diff(X)


plot.ts(SAN.prices)
plot.ts(X)

plot.ts(diffinv(as.vector(estima)))
plot.ts(diffinv(as.vector(estimaf))+1.2)

lines(X,col="red")
lines(diffinv(as.vector(estima))+3,col="red")
