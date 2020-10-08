x<-seq(1,18)
y<-c(623,506,822,1259,1544,2000,1630,1890,2443,3431,2822,4946,3646,4517,6584,7937,8578,7817)
modelo<-lm(log(y)~x)
summary(modelo)
plot(modelo)
shapiro.test(modelo$residuals)
exp(predict(modelo,newdata = data.frame(x=19),interval = "prediction",level=0.85))

plot(y)
lines(exp(predict(modelo)))






##############-------Madrid


x<-seq(1,15)
y<-c(8,21,31,56,81,213,355,390,498,628,1756,781,873,1777,2245)
2569
modelo<-lm(log(y)~x)
summary(modelo)
plot(modelo)
shapiro.test(modelo$residuals)
exp(predict(modelo,newdata = data.frame(x=16),interval = "prediction",level=0.9))

plot(y)
lines(exp(predict(modelo)))


###########-------Barcelona



x<-seq(1,11)
y<-c(4,12,18,41,55,82,501,933,1211,1939,2073)
2569
1655
modelo<-lm(log(y)~x)
summary(modelo)
plot(modelo)
shapiro.test(modelo$residuals)
exp(predict(modelo,newdata = data.frame(x=12),interval = "prediction",level=0.95))

plot(y)
lines(exp(predict(modelo)))
