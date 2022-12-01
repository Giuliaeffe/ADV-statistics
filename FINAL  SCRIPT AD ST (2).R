###library 
library(readxl)
library(regclass)
library(car)
library(mctest)
library(stats)
library(lmtest)
library(sandwich)
library(ggplot2)
library(rvest)
library(gap)
library(tidyr)
library(MASS)
library(tseries)

#data 
datBI <- read.table("datBI.txt",sep = ";")

# TRANSFORM PA0200 IN A DUMMY 
PNA0600a<-as.numeric(datBI$PNA0600a)
pos = c(which(PNA0600a==5))
v1_INT = replace(PNA0600a, pos, 1)
posna = which(is.na(PNA0600a))
v12 = replace(v1_INT,posna, 0)
v2 = replace(v12, -pos, 0)
PA0200<-as.numeric(datBI$PA0200)
pos200 = c(which(PA0200==5))
v1_INT200 = replace(PA0200, pos200, 1)
posna200 = which(is.na(PA0200))
v12200 = replace(v1_INT200,posna, 0)
v2200 = replace(v12200, -pos200, 0)
#same for PA0100
datBI$PA0100
PA0100<-as.numeric(datBI$PA0100)
PA0100 = replace(PA0100, which(PA0100==1), 0)
PA0100 = replace(PA0100, which(PA0100==2), 1)
PA0100 = replace(PA0100, which(PA0100==3), 1)
PA0100 = replace(PA0100, which(PA0100==4), 0)
PA0100 = replace(PA0100, which(PA0100==5), 0)

#### DATA FRAME ####
DATA <- data.frame(datBI$DA1110,datBI$HD1800,datBI$PA0100,v2200,datBI$DI2000,datBI$HD1210,datBI$HB4710,datBI$HB4400,datBI$PE0600)
logDATA <- data.frame(datBI$DA1110),datBI$HD1800,datBI$PA0100,v2200,datBI$DI2000,datBI$HD1210,datBI$HB4710,datBI$HB4400,datBI$PE0600)
### SCATTERPLOT ###
pairs(DATA,lower.panel = panel.smooth, upper.panel = panel.smooth)

####  LET'S TRY TO SPECIFY THE FIRST MODEL ###

LM <- lm(datBI$DA1110~datBI$PE0600+datBI$HD1800+datBI$PA0100+v2200+datBI$DI2000+datBI$HD1210+datBI$HB4710+datBI$HB4400,data = DATA)
summary(LM)

LM0 <- update(LM, .~.-datBI$PE0600 )
summary(LM0)
anova(LM,LM0)

LM1 <- update(LM0, .~.-datBI$HD1800)
summary(LM1)
anova(LM0,LM1)

LM2 <-  update(LM1, .~.-datBI$PA0100)
summary(LM2)
anova(LM1,LM2)


#Confidence intervals 
confint(LM2, level = 0.95)

###########GRAPHICAL ANALYSIS

#Plot of standardized residuals vs each explanatory variable
par(mfrow = c(2,3))
VARIABLIX = model.matrix(LM2)
VARIABLIX = VARIABLIX[,c(2,3,4,5,6,7)]
NAMEVAR = c("2","3","4","5","6","7")
for(i in 1:6){
  plot(VARIABLIX[,i], rstandard(LM2), ylab = "Standardized Residuals", xlab = NAMEVAR[i], main = paste("Standardized residuals vs", NAMEVAR[i], sep = " ") )
  panel.smooth(x = VARIABLIX[,i], y = rstandard(LM2))
}



#Plot Residuals vs fitted values
par(mfrow=c(1,3))
plot(resid(LM2)~fitted(LM2), xlab=expression(paste('estimated ',hat('y'), sep=' ')), ylab='residuals')
abline(h = 0, col = "red")
#Plot Standardized residuals Vs. Fitted values
plot(LM2$fitted.values, rstandard(LM2), ylab = "Standardized residuals", xlab = "Fitted values", main = "Standardized residuals Vs. Fitted values")
abline(h = 0, col = "red")
#Plot studentized  residuals Vs. Fitted values
plot(LM2$fitted.values, rstudent(LM2), ylab = "Studentized residuals", xlab = "Fitted values", main = "Standardized residuals Vs. Fitted values")
abline(h = 0, col = "red")
par(mfrow=c(1,1))




library(car)
qqPlot(LM2, main = "Student-t Q-Q Plot of residuals", col = "red")
par(mfrow=c(1,1))

############NORMALITY

### qq plot and histograms
par(mfrow=c(1,3))
hist(LM2$residuals, main='Frequency distribution of residuals')
qqnorm(resid(LM2), main = "Normal Q-Q Plot of residuals", col = "red")
qqline(resid(LM2), col = "blue", lwd = 2)
library(moments)
skewness(LM2$residuals)
kurtosis(LM2$residuals)
jarque.bera.test(LM2$residuals)
shapiro.test(LM2$residuals)

# Leverage points and outliers

lev<-hat(model.matrix(LM2))
plot(lev,ylab="Leverages",main="Index plot of Leverages") 
lev.t<-2*ncol(model.matrix(LM2))/nrow(model.matrix(LM2)
abline(h=lev.t, col='red')
# units with high leverage:
h.l<-cbind(which(lev > lev.t),lev[c(which(lev > lev.t))])

# Cook's distance
ckd<-cooks.distance(LM2)
plot(ckd)
abline(h=4/length(ckd), col='green')
d.inf<-ckd<= 4/length(ckd)
table(d.inf)

# there are 69 influential observations 

# OTULIER 
library(car)
print(outl<-outlierTest(LM2,cutoff = 0.05, n.max = 69))
outliers<-c(6584,283,4030,2909,7129,85,5884,2687,7452,3962,6037,712)
DATA_senza_out <- DATA[-outliers,]

##### We take away the outliers
LM2_noout <- lm(formula = DATA_senza_out$datBI.DA1110 ~ DATA_senza_out$v2200 + DATA_senza_out$datBI.DI2000 + DATA_senza_out$datBI.HD1210 + 
                  DATA_senza_out$datBI.HB4710 + DATA_senza_out$datBI.HB4400, data = DATA_senza_out)
summary(LM2_noout)

LM3_noout <-  update(LM2_noout, .~.-DATA_senza_out$datBI.HD1210)
summary(LM3_noout)
jarque.bera.test(LM3_noout$residuals)
shapiro.test(LM3_noout$residuals)
plot(LM3_noout$fitted.values,LM3_noout$residuals)

# let's then try to exclude all influential observations:
LM3 = update(LM3_noout, subset = ckd<=4/length(ckd)) 
summary(LM3)
shapiro.test(LM3$residuals)
jarque.test(LM3$residuals)
# Considering the high sample size, we might also appeal to the CLT 

# Let's try first to see, though, if it's possible to overcome the problem 
# by a suitable transformation of the dependent variable:
BOX <- boxCox(LM2) 
boxcox(LM2)
lambda <- BOX$x[which.max(BOX$y)]#log is suggested

DI2000.log<-log(datBI$DI2000)
DI2000.log
DI2000.log.inf = which(DI2000.log== -Inf)
DI2000.log.inf
DI2000.log = replace(DI2000.log,DI2000.log.inf, 0)
DI2000.log

HD1210.log<-log(datBI$HD1210)
HD1210.log
HD1210.log.inf = which(HD1210.log== -Inf)
HD1210.log.inf
HD1210.log = replace(HD1210.log,HD1210.log.inf, 0)
HD1210.log

DATALOG <- data.frame(datBI$DA1110,v2200,DI2000.log,HD1210.log,datBI$HB4710,datBI$HB440,datBI$DI2000,datBI$HD1210)
DATALOG <- na.omit(DATALOG)

LM2.a <- lm(formula =log(DATALOG$datBI.DA1110) ~ DATALOG$v2200 +DATALOG$DI2000.log+ DATALOG$HD1210.log + 
              log(DATALOG$datBI.HB4710) + DATALOG$datBI.HB440, data = DATALOG)

summary(LM2.a)
LM2aa <- update(LM2.a, .~.-HD1210.log)
summary(LM2aa)

jarque.bera.test(LM2aa$residuals)
shapiro.test(LM2aa$residuals)
plot(LM2aa$residuals)
#the problems persist

#so we try to take away the outliers again
outlierTest(LM2aa)
outliers2 <- c(5884,6512,5232,934,4857,1943,5203,3875)
DATALOG2 <- DATALOG[-(outliers2),]

LM2.aa_without_o <- lm(formula =log(DATALOG2$datBI.DA1110)  ~ DATALOG2$v2200 + DATALOG2$DI2000.log+ log(DATALOG2$datBI.HB440 )+ 
                         log(DATALOG2$datBI.HB4710) , data = DATALOG2)
summary(LM2.aa_without_o)
jarque.bera.test(LM2.aa_without_o$residuals)
shapiro.test(LM2.aa_without_o$residuals)
plot(LM2.aa_without_o%)

plot(LM2.aa_without_o$fitted.values,LM2.aa_without_o$residuals)
plot(LM2$fitted.values,LM2$residuals)

par(mfrow = c(2,3))
VARIABLIX2 = model.matrix(LM2.aa_without_o)
VARIABLIX2 = VARIABLIX2[,c(2,3,4,5,6)]
NAMEVAR2 = c("2","3","4","5","6")
for(i in 1:6){
  plot(VARIABLIX2[,i], rstandard(LM2.aa_without_o), ylab = "Standardized Residuals", xlab = NAMEVAR2[i], main = paste("Standardized residuals vs", NAMEVAR2[i], sep = " ") )
  panel.smooth(x = VARIABLIX2[,i], y = rstandard(LM2.aa_without_o))
}


############# Other MIS-SPECIFICATION TESTS

##### RESET <- linearity -> check if some regressors enter in the model in a non linear way 
lmtest::resettest(LM2, power = 2:3, type='fitted')
lmtest::resettest(LM2, power = 2:3, type="regressor")
lmtest::resettest(LM2, power = 2:3, type="princomp")


############### VIF <- irrilevant variable and multicollinearity check 
library(car)
vif(LM3)


############### HETEROSKEDASTICITY check
e2.LM3<-resid(LM2)^2
plot(e2.LM3~fitted(LM2), xlab=expression(paste('fitted values',hat('y'), sep=' ')), ylab=expression(paste('squared residuals',hat('e')^2, sep=' ')))
x.LM3<-as.data.frame(model.matrix(LM2))
par(mfrow=c(2,3))
for (i in 2:6){
  plot(e2.LM3~ x.LM3[,i], ylab=expression(paste('squared residuals',hat('e')^2, sep=' ')), xlab=colnames(x.LM3)[i] )
}
par(mfrow=c(1,1))
# no strong evidence for a particular patterns

## Breush-Pagan test
lmtest::bptest(LM2,studentize=T)       
lmtest::bptest(LM2,studentize=F)   
## white test
library(lmtest)
lmtest::bptest(LM2, varformula=~v2200 + datBI$DI2000 + datBI$HD1210 + 
                 datBI$HB4710 + datBI$HB4400  + I(datBI$DI2000)^2 + I(datBI$HD1210)^2 + 
                 I(datBI$HB4710)^2 + I(datBI$HB4400)^2, data=DATA)



######### Independence (uncorrelated errors) checl
acf(resid(LM2))
pacf(resid(LM2))
#Durbin-Watson test (H0: No residuals correlation)
durbinWatsonTest(LM2)#independence


par(mfrow = c(1,2))
hist(LM2ab$residuals, prob = TRUE, main = "Histogram of resisuals - model LM2ab")
curve(dnorm(x, mean=mean(LM2ab$residuals), sd=sd(LM2ab$residuals)), add=TRUE, col = "blue")

qqnorm(resid(LM2.aa_without_o), main = "Normal QQPlot - model LM2ab", col = "darkgrey")
qqline(resid(LM2.aa_without_o), col = "blue", lwd = 2)


#Durbin-Watson test
durbinWatsonTest(LM2.aa_without_o)

#Test if the model is correct linear specified
resettest(LM2ab, power = 2:3, type = "fitted")
resettest(LM2ab, power = 2:3, type = "regressor")
resettest(LM2ab, power = 2:3, type = "princomp")
#no good results

SQRDRESIDUALS2 = (LM2.aa_without_o$residuals)
plot(SQRDRESIDUALS2) ~ fitted(LM2.aa_without_o), xlab = "Fitted values", ylab = "Squared residuals", main = "Residuals Vs. Fit")
XLM2ab = as.data.frame(model.matrix(LM2.aa_without_o))

par(mfrow = c(1,3))
for (i in 2:4){
  plot(SQRDRESIDUALS2 ~ XLM2ab[,i], ylab = "Squared residuals", xlab = colnames(XLM2ab)[i])
}

#Breush-Pagan test 
bptest(LM2.aa_without_o)

#Check the est parameter:
coeftest(LM2ab, sandwich::vcovHAC(LM2ab, lag = Inf))

#to try the solve the problem of homosk. residuals we try to fit a Weighted least squares regression (preferred to a general least quares, b/c  our residuals are uncorrellated)
wt <- 1 / lm(abs(LM2$residuals) ~ LM2$fitted.values)$fitted.values^2
LM2_WLS <- lm(DATALOG$datBI.DA1110~DATALOG$v2200+DATALOG$datBI.DI2000+DATALOG$datBI.HB4710+DATALOG$datBI.HB440,data = DATALOG,weights = wt)
summary(LM2_WLS)
shapiro.test(LM2_WLS$residuals)#the residuals are stil not normal
jarque.bera.test(LM2_WLS$residuals)
lmtest::bptest(LM2_WLS)# Breusch-Pagan "studentized" test with P-value = 1, they are now homoskedastik
plot(LM2_WLS$residuals)#as we can see from the plot 
## white test, high P-value pass the test
lmtest::bptest(LM2_WLS, varformula=~DATALOG$v2200+DATALOG$datBI.DI2000+DATALOG$datBI.HB4710+DATALOG$datBI.HB440
               + I(DATALOG$datBI.DI2000)^2 + I(DATALOG$datBI.HB4710)^2 + I(DATALOG$datBI.HB440)^2, data=DATALOG)

lmtest::resettest(LM2_WLS, power = 2:3, type='fitted')
lmtest::resettest(LM2_WLS, power = 2:3, type="regressor")
lmtest::resettest(LM2_WLS, power = 2:3, type="princomp")
#still bad resutls

BOXWLS<-boxcox(LM2_WLS)
lambdawls <- BOX$x[which.max(BOX$y)]


LM2_WLS_pesi <- lm((log(DATALOG$datBI.DA1110))~DATALOG$v2200+DATALOG$DI2000.log+log(DATALOG$datBI.HB4710)+log(DATALOG$datBI.HB440),data = DATALOG)
wt2 <- 1 / lm(abs(LM2_WLS_tr$residuals) ~ LM2_WLS_tr$fitted.values)$fitted.values^2
LM2_WLS_log <- lm(log(DATALOG$datBI.DA1110)~DATALOG$v2200+DATALOG$DI2000.log+log(DATALOG$datBI.HB4710)+log(DATALOG$datBI.HB440),data = DATALOG,weights = wt2)
summary(LM2_WLS_log)
shapiro.test(LM2_WLS_log$residuals)
jarque.bera.test(LM2_WLS_log$residuals)

DA1110_log<-(DATALOG$datBI.DA1110)^(2)
DI2000_log<-(DATALOG$datBI.DI2000)^(2)
HB4710_log<-(DATALOG$datBI.HB4710)^(2)
HB440_log<-(DATALOG$datBI.HB440)^(2)
LM2_WLS_tr <- lm(DA1110_log~DATALOG$v2200+(DI2000_log)+(HB4710_log)+(HB440_log),data = DATALOG,weights = wt)
plot(LM2_WLS_tr$residuals)
shapiro.test(LM2_WLS_tr$residuals)
jarque.bera.test(LM2_WLS_tr$residuals)
