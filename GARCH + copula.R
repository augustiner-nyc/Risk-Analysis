setwd("C:/Users/mhube/OneDrive/Desktop")
library(FRAPO)
showEx("C9R1")

GE = read.csv("GE-1.csv")
GERet<-returnseries(GE$Adj.Close)
date<-as.character(GE$Date)
GERetTS<-timeSeries(GERet, charvec = date)
plot(GERetTS)
acf(GERetTS, na.action = na.pass)
pacf(GERetTS, na.action = na.pass) #ARMA(0,0)
acf((GERetTS)^2, na.action=na.pass)

yyy<-garchFit(formula = ~arma(0,0)+garch(1,1), data=na.omit(GERetTS), cond.dist = "std", include.mean = FALSE)
coef(yyy)
summary(yyy)
sigma <-  as.numeric(predict(yyy, n.ahead = 1)[3])
df <- as.numeric(coef(yyy)["shape"])
rand<-rt(100000, df=df)
hist(rand, breaks=100)
quant<-qt(.01, df=df)
quant2<-sort(rand, decreasing = TRUE)[99000]
VAR99<-sigma*quant
predict(yyy, n.ahead=1)[3]



library(QRM)
library(fGarch)
## Losses
data(EuStockMarkets)
head(EuStockMarkets)
str(EuStockMarkets)
#Generate losses of each of the asset
loss<-as.data.frame(na.omit(-1.0*diff(log(EuStockMarkets))*100.0))
head(loss)
acf(loss$FTSE)
pacf(loss$FTSE)#ARMA(0,0) No autocorrelation (kind but a bit)
acf(loss$DAX)
pacf(loss$DAX)#ARMA(0,0) No Autocorrelation at all
acf(loss$CAC)
pacf(loss$CAC)#ARMA(0,0) No Autocorrelation
acf(loss$SMI)
pacf(loss$SMI)#ARMA(0,0) No Correlation either
head(loss)
str(loss)
## GARCH
# Estimate GARCH model for the whole parameters instead of seperate, including shape
# Step 1
gfit<-lapply(loss,garchFit,formula=~arma(0,0)+garch(1,1), cond.dist="std",trace=FALSE)
gfit
mode(gfit)
# one-step-ahead forecasts of the conditional variance (measure of risk) 
# will be used in the computation 
# of the portfolio loss variates, and both of these are needed for the calculation
# ofthe pseudo-uniform variables

gprog<-unlist(lapply(gfit,function(x) predict(x,n.ahead = 1)[3]))
# Estimate degrees-of-freedom parameters for each European market
gshape<-unlist(lapply(gfit, function(x) x@fit$coef["shape"]))
# Can take a look at all paramaters of the GARCH model
gcoef<-unlist(lapply(gfit, function(x) x@fit$coef))
# Estimates conditional standardized residuals are extracted.(h.t - conditional variance)
# Step 2
gresid<-as.matrix(data.frame(lapply(gfit,function(x) x@residuals / sqrt(x@h.t))))
head(gresid)
#QQ plots of the standardized residuals
par(mfrow=c(2,2))
unlist(lapply(gfit, function(x) plot(x, which=13)))
#ACF of the squared residuals
par(mfrow=c(1,1))
unlist(lapply(gfit, function(x) plot(x, which=11)))

## Copula
# pseudo-uniform variables that generates probabilites for each risk 
# (measured as conditional resid)
# Step 3
U <- sapply(1:4, function(y) pt(gresid[, y], df = gshape[y]))
head(U)
hist(U[,4])
# Student's t copula is estimated based on Kendall's rank correlations.  
# Step 4
cop <- fit.tcopula(Udata = U, method = "Kendall")
# 100,000 random losses simulated for each financial instrument 
# Step 5
rcop <- rcopula.t(100000, df = cop$nu, Sigma = cop$P)
head(rcop)
hist(rcop[,3], breaks=100)
# Compute the quantiles for these Monte Carlo draws.
#Step 6
qcop <- sapply(1:4, function(x) qstd(rcop[, x], nu = gshape[x]))
head(qcop)
hist(qcop[,1], breaks = 100)
# creating a matix of 1 period ahead predictions of standard deviations
ht.mat <- matrix(gprog, nrow = 100000, ncol = ncol(loss), byrow = TRUE)
head(ht.mat)
pf <- qcop * ht.mat
head(pf)
## ES 95 percent
weights <- c(0.4, 0.2, 0.2, 0.2)
# The simulated portfolio losses are then determined 
# as the outcome of the matrix-weight vector product
# Step 7
pfall <- (qcop * ht.mat) %*% weights
head(pfall)
tail(pfall)
hist(pfall,breaks = 100)
# Step 8
pfall.es95 <- mean(tail(sort(pfall), 5000))
pfall.es95
pfall.var95 <- min(tail(sort(pfall), 5000))
pfall.var95
pfall.es99 <- mean(tail(sort(pfall), 1000))
pfall.es99
pfall.var99 <- min(tail(sort(pfall), 1000))
pfall.var99