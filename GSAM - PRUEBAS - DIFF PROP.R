
#Library

library(TSA)
library(tseries)
library(sm)
library(dplyr)
library(lubridate)
library(dlnm)
library(fda.usc)
library(utils)
library(PerformanceAnalytics)


#C:/Users/rosal/OneDrive/Escritorio/Master- Tecnicas Estadisticas/tfm/Trabajo TFM/FDA

data_aire<-read.csv2("C:/Users/rsanchez/Desktop/PROSAE/Temp - Admisiones/Data_aire.csv",sep=",")
data_temp<- read.csv2("C:/Users/rsanchez/Desktop/PROSAE/Temp - Admisiones/Data_temp.csv",sep = ",")

# Variables of interest as numeric: 

var<-c(5:9,12:16)

for(i in var) {
  data_temp[,i]<-as.numeric(as.character(data_temp[,i]))}

for(j in 3:9 ){
  data_aire[,j]<-as.numeric(as.character(data_aire[,j]))}

df<-merge(data_temp,data_aire)
df$date<-as.Date(df$date,format="%Y-%m-%d")

#str(df)

######################################################
# Calculating the correlation at lag = 0
# There is a clear negative correlation. 
######################################################

corr_pear<-cor.test(df$prop_hosp,df$Temperatura, method = c("pearson"))$estimate;corr_pear
corr_ken<-cor.test(df$prop_hosp,df$Temperatura, method = c("kendall"))$estimate;corr_ken
corr_spear<-cor.test(df$prop_hosp,df$Temperatura, method = c("spearman"))$estimate;corr_spear


######################################################
#Drawing that RAW data as Time Series:
######################################################

# Omit NA:  Registers of -9999 and -99999.90

ind<-which(df$Humedad==-9999)
df[ind,"Humedad"]<-mean(df[ind-6,"Humedad"])

# Replacing it for the mean

ind<-which(df$O3==-99999.90)
df[ind,"O3"]<-round(mean(df[-ind,"O3"]),0) 

ind<-which(df$PM10==-99999.90)
df[ind,"PM10"]<-round(mean(df[-ind,"PM10"]),0)

ind<-which(df$CO==-99999.90)
df[ind,"CO"]<-round(mean(df[-ind,"CO"]),2)

ind<-which(df$NO2==-99999.90)
df[ind,"NO2"]<-round(mean(df[-ind,"NO2"]),0)

ind<-which(df$NOX==-99999.90)
df[ind,"NOX"]<-round(mean(df[-ind,"NOX"]),0)

ind<-which(df$SO2==-99999.90)
df[ind,"SO2"]<-round(mean(df[-ind,"SO2"]),0)

# Time Series 
interest<-c(4,6,8,9,18,20,21:24)
serie <- ts(cbind(df[,interest]),start=c(2017,1,1), frequency=365)
plot(serie, type="o", lwd=3)



###################################################################
# Arranging the data for Functional Regression with diferent lags
###################################################################

# t represents the register day
# n is the number of days we are doing the lag


# Create variables that incorporate information on time t: 
#...........................................................
# We can incorporate the prediction of meteorological 
# variables as values at time t. Example: Temperature
#
# We want to create a matrix where the rows are X_{t-i} for i=max_lag,...,0
# Example: for max_lag = 7 each row (register day) of the matrix will be the series : 
# X_{t-7},X_{t-6},X_{t-5},X_{t-4},X_{t-3},X_{t-2},X_{t-1},X_{t-0}


lagmatrix_t<-function(n,x){  
  max_lag=n
  lag_matrix<-matrix(NA,nrow =length(x) ,ncol=max_lag+1) 
  for (i in 0:max_lag){lag_matrix[,(i+1)]<-lag(x,max_lag-i)}
  return(lag_matrix[-(1:max_lag),])
}


# Create variables that only have information from the past: 
#...........................................................
# This is beacuse in practice we'll never know the value of a variable at
# time t if the day  hasn't ended yet. Example: Proportion of 
# red cases .
#
# We want to create a matrix where the rows are X_{t-i} for i=max_lag,...,1
# Example: for max_lag = 7 each row (register day) of the matrix will be the series : 
# X_{t-7},X_{t-6},X_{t-5},X_{t-4},X_{t-3},X_{t-2},X_{t-1}


lagmatrix_t1<-function(n,x){  
  max_lag=n
  lag_matrix<-matrix(NA,nrow =length(x) ,ncol=max_lag+1) 
  for (i in 0:max_lag){lag_matrix[,(i+1)]<-lag(x,max_lag-i)}
  return(lag_matrix[-(1:max_lag),-1])
}


# Resume of variables:

# Escalar Variables: 
#...................
# Meteoriological varibales at time t 
# Meteoriological varibales at time t-1
# TRIAJE varibales at time t-1

# Functional Variables: 
#...................
# Meteoriological varibales that start at time t 
# Meteoriological varibales that start at t-1
# TRIAJE varibales that strat at time t-1


# This loop creates all the functional  variables up to a maxlag and 
# acomodates them with the scalar variables in a list: 

list_create<-function(maxlag){
  
names_t<-colnames(df)[c(6:9)]
names_t1<-colnames(df)[c(6:9,12:16)]
  
for(i in 2:maxlag){ # It has to start in 2 because of the t-1 scalar variables
  
  for(j in 1:length(names_t)){
    name<-names_t[j]
    assign(paste0(name,"_t_",i),fdata(lagmatrix_t(i,df[,name]),argvals=0:i)) 
  }
  
  for(j in 1:length(names_t1)){
    name<-names_t1[j]
    assign(paste0(name,"_t1_",i),fdata(lagmatrix_t1(i,df[,name]))) # Question : argvals ?? 
    }
  
  #Scalar variables  (df):
  

  df$diff_prop<-c(NA,df[-1,"num_hosp"]/df[-365,"num_hosp"])
  df_t= df[-(1:i),] 
  df_t1= df[-(c(1:(i-1),nrow(df[-(1:i),]))),-c(1:2,10:11,17)]
  colnames(df_t1)<-c(paste0(colnames(df_t)[-c(1:2,10:11,17)],"_t1"))
  DF=cbind(df_t,df_t1)


  assign(paste0("varlist_",i),list(df = DF,
                                   precipitacion_fd_t=get(paste0(names_t[1],"_t_",i)),
                                   horasdefrio_fd_t=get(paste0(names_t[2],"_t_",i)),
                                   humedad_fd_t=get(paste0(names_t[3],"_t_",i)),
                                   temperatura_fd_t=get(paste0(names_t[4],"_t_",i)),
                                   precipitacion_fd_t1=get(paste0(names_t1[1],"_t1_",i)),
                                   horasdefrio_fd_t1=get(paste0(names_t1[2],"_t1_",i)),
                                   humedad_fd_t1=get(paste0(names_t1[3],"_t1_",i)),
                                   temperatura_fd_t1=get(paste0(names_t1[4],"_t1_",i)),
                                   num_azul_fd_t1=get(paste0(names_t1[5],"_t1_",i)),
                                   num_verde_fd_t1=get(paste0(names_t1[6],"_t1_",i)),
                                   num_amarillo_fd_t1=get(paste0(names_t1[7],"_t1_",i)),
                                   num_naranja_fd_t1=get(paste0(names_t1[8],"_t1_",i)),
                                   num_rojo_fd_t1=get(paste0(names_t1[9],"_t1_",i))),envir = globalenv())

}}


list_create(20) 

varlist_10$df

####################################
# Variable selection for the models
####################################

#Family  = Gaussian

####################################################################################################
#GSAM : --------------------------------------------------------------------------------------------
####################################################################################################

# GSAM has its own funtion to find the best combination of covariables to model , 
# so we are just going to loop for each lag
# and then select the best R^2 

gsam_best<-function(lags){
  
  excluir<-c("X","date","num_admin","num_hosp","time","dow","dates",
             "num_azul","num_verde","num_amarillo","num_naranja","num_rojo",
             "CO","NO","NO2","NOX","O3","PM10","SO2","diff_prop_t1","prop_hosp",
            "CO_t1","NO_t1","NO2_t1","NOX_t1","O3_t1","PM10_t1","SO2_t1")
  
  # Don't know why it doesn`t work for some lags and for Poisson family
  # This loop takes a lot of time to calculate 
  # Aproximate of time: 
  
  works<-lags
  R_gsam<-matrix(0,nrow = 30,ncol=2)
  pb <- txtProgressBar(min = 0, max = length(works), style = 3);k=1
  
  for (i in works){
    #fit<-fregre.gsam.vs(get(paste0("varlist_",i)),y="prop_hosp",exclude = excluir)
    print(i)
    fit<-fregre.gsam.vs(get(paste0("varlist_",i)),y="diff_prop",exclude = excluir)
    R_gsam[i,1]<- as.numeric(summary(fit)$r.sq)
    R_gsam[i,2]<- i
    setTxtProgressBar(pb, k)
    k=k+1
  }
  
  close(pb)
  
  ind<-order(as.numeric(R_gsam[,1]),decreasing = TRUE)
  Best_gsam<<-head(R_gsam[ind,],5)
  
  for(h in 1:5){
    i=Best_gsam[h,2]
    #assign(paste0("best_gsam_",h),fregre.gsam.vs(get(paste0("varlist_",i)),y="prop_hosp",exclude = excluir),envir = globalenv())
    assign(paste0("best_gsam_",h),fregre.gsam.vs(get(paste0("varlist_",i)),y="diff_prop",exclude = excluir),envir = globalenv())
  }
  
}


gsam_best(c(4:10,13:14))# 1:3, 11,12 do not work

#####################
# Top 5 GSAM Models
#####################


gsam_1<-get(paste0("best_gsam_",1))
gsam_2<-get(paste0("best_gsam_",2))
gsam_3<-get(paste0("best_gsam_",3))
gsam_4<-get(paste0("best_gsam_",4))
gsam_5<-get(paste0("best_gsam_",5))

A<-c(gsam_1$formula.ini,
     gsam_2$formula.ini,
     gsam_3$formula.ini,
     gsam_4$formula.ini,
     gsam_5$formula.ini)

res<-cbind(round(Best_gsam,2),as.character(A));
colnames(res)<-c("R.sq","lag","formula")
res

summary(gsam_1)
summary(gsam_2)
summary(gsam_3)
summary(gsam_4)
summary(gsam_5)

windows(25,10)
par(mfrow=c(1,1))
plot.ts(gsam_1$y,col="gray",main="Fits with GSAM",frame=FALSE,xlim=c(0,350))
lines(gsam_5$fitted.values,col="skyblue")
lines(gsam_4$fitted.values,col="skyblue1")
lines(gsam_3$fitted.values,col="skyblue2")
lines(gsam_2$fitted.values,col="skyblue3")
lines(gsam_1$fitted.values,col="skyblue4")
legend("top", legend=round(c(summary(gsam_5)$r.sq, summary(gsam_4)$r.sq,summary(gsam_3)$r.sq,
                             summary(gsam_2)$r.sq,summary(gsam_1)$r.sq),4),
       col=c("skyblue","skyblue1","skyblue2","skyblue3","skyblue4"),
       title = expression(R^2), lty=1,cex=0.5,horiz = TRUE)

plot(gsam_1$basis.x$temperatura_fd_t$basis)
plot(gsam_1)


























#Note:

# Incorporting  air contamination scalar variables at time t-1 does not make any change : 
# So, these variables aren't interesting. 


gsam_best_aire<-function(lags){
  
  excluir<-c("X","date","num_admin","num_hosp","time","dow","dates",
             "num_azul","num_verde","num_amarillo","num_naranja","num_rojo",
             "CO","NO","NO2","NOX","O3","PM10","SO2",            
             "num_admin_t1","num_hosp_t1","prop_hosp_t1")

  
  works<-lags
  R_gsam<-matrix(0,nrow = 30,ncol=2)
  pb <- txtProgressBar(min = 0, max = length(works), style = 2);k=1
  
  for (i in works){
    print(i)
    fit<-fregre.gsam.vs(get(paste0("varlist_",i)),y="prop_hosp",exclude = excluir)
    R_gsam[i,1]<- as.numeric(summary(fit)$r.sq)
    R_gsam[i,2]<- i
    setTxtProgressBar(pb, k)
    k=k+1
  }
  
  close(pb)
  
  ind<-order(as.numeric(R_gsam[,1]),decreasing = TRUE)
  Best_gsam<<-head(R_gsam[ind,],5)
  
  for(h in 1:5){
    i=Best_gsam[h,2]
    assign(paste0("best_gsam_aire",h),fregre.gsam.vs(get(paste0("varlist_",i)),y="prop_hosp",exclude = excluir),envir = globalenv())
  }
  
}


gsam_best_aire(c(4:10,12:15))


##########################################
# Top 5 GSAM Models - AIR CONTAMINATION
##########################################

gsam_aire_1<-get(paste0("best_gsam_aire",1))
gsam_aire_2<-get(paste0("best_gsam_aire",2))
gsam_aire_3<-get(paste0("best_gsam_aire",3))
gsam_aire_4<-get(paste0("best_gsam_aire",4))
gsam_aire_5<-get(paste0("best_gsam_aire",5))

summary(gsam_aire_1)
summary(gsam_aire_2)
summary(gsam_aire_3)
summary(gsam_aire_4)
summary(gsam_aire_5)

plot.ts(df$prop_hosp,col="gray",main="Fits with GSAM",frame=FALSE,xlim=c(0,350))
lines(gsam_aire_5$fitted.values,col="skyblue")
lines(gsam_aire_4$fitted.values,col="skyblue1")
lines(gsam_aire_3$fitted.values,col="skyblue2")
lines(gsam_aire_2$fitted.values,col="skyblue3")
lines(gsam_aire_1$fitted.values,col="skyblue4")
legend("top", legend=round(c(summary(gsam_aire_5)$r.sq, summary(gsam_aire_4)$r.sq,summary(gsam_aire_3)$r.sq,
                             summary(gsam_aire_2)$r.sq,summary(gsam_aire_1)$r.sq),4),
       col=c("skyblue","skyblue1","skyblue2","skyblue3","skyblue4"),
       title = expression(R^2), lty=1,cex=0.5,horiz = TRUE)



getAnywhere(fregre.gsam.vs)

