setwd("C:/Users/mhube/OneDrive/Desktop")
library(FRAPO)
#showEx("C11R1")

library(fPortfolio)
library(lattice)
## Loading data and calculating returns
data(SPISECTOR)
Idx <- interpNA(SPISECTOR[, -1], method = "before")
R <- returnseries(Idx, method = "discrete", trim = TRUE)
V <- cov(R)
head(Idx)

## Portfolio Optimizations
GMVw <- Weights(PGMV(R)) #Gloval minimal variance
MDPw <- Weights(PMD(R))  #Most diversified portfolio
MTDw <- Weights(PMTD(R))  #miminum tail dependent
ERCw <- Weights(PERC(V))  #Equal risk contributed
## Combining solutions
W <- cbind(GMVw, MDPw, MTDw, ERCw)
MRC <- apply(W, 2, mrc, Sigma = V) #Marginal Risk Contribution
rownames(MRC) <- colnames(Idx)
colnames(MRC) <- c("GMV", "MDP", "MTD", "ERC")
head(MRC)
## Plot of allocations
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
dotchart(GMVw, xlim = c(0, 40), main = "GMV Allocation", pch = 19)
dotchart(MDPw - GMVw, xlim = c(-20, 20), main = "MDP vs. GMV",
         pch = 19) #Include relative weight comapre to the Global Minimal Variance Portfolio
abline(v = 0, col = "grey")
dotchart(MTDw - GMVw, xlim = c(-20, 20), main = "MTD vs. GMV",
         pch = 19)
abline(v = 0, col = "grey")
dotchart(ERCw - GMVw, xlim = c(-20, 20), main = "ERC vs. GMV",
         pch = 19)
abline(v = 0, col = "grey")
par(oldpar)

## lattice plots of MRC
Sector <- factor(rep(rownames(MRC), 4),
                 levels = sort(rownames(MRC))) #Repeat all Sectors 4 times (for all portfolios)
Port <- factor(rep(colnames(MRC), each = 9),
               levels = colnames(MRC)) #Repeat 9 times
MRCdf <- data.frame(MRC = c(MRC), Port, Sector)
dotplot(Sector ~ MRC | Port, groups = Port, data = MRCdf,
        xlab = "Percentages",
        main = "MR Contributions by Sector per Portfolio",
        col = "black", pch = 19)
dotplot(Port ~ MRC | Sector, groups = Sector, data = MRCdf,
        xlab = "Percentages",
        main = "MR Contributions by Portfolio per Sector",
        col = "black", pch = 19)


showEx("C11R2")
library(PerformanceAnalytics)
## Portfolio Risk Measures and Characteristics
Rdec <- R / 100
Pret <- apply(W, 2, function(x) Rdec %*% x / 100) #Weighted sum of daily returns in decimals
SD <- apply(Pret, 2, sd) * 100
ES95 <- apply(Pret, 2, function(x)
  abs(ES(R = x, method = "modified") * 100)) #Expected shortfall at 95 level
DR <- apply(W, 2, dr, Sigma = V) #Diversification Ratio
CR <- apply(W, 2, cr, Sigma = V) #Concentration Ratio
ExRet<-apply(Pret,2,function(x) mean(x)*100)
## Summarising results
Res <- rbind(SD, ES95, DR, CR, ExRet)
Res

showEx("C11R3")
library(copula)
## S&P 500 (ccan they outperform S&P 500?)
data(INDTRACK6)
length(INDTRACK6[,1]) #291 Observations
## Market and Asset Returns
RM <- returnseries(INDTRACK6[1:260, 1], method = "discrete",
                   trim = TRUE)
RA <- returnseries(INDTRACK6[1:260, -1], method = "discrete",
                   trim = TRUE)
length(INDTRACK6[,1])
# Select stocks that co-move less proportionally than the market in absolute terms
Beta <- apply(RA, 2, function(x) cov(x, RM) / var(RM))
Tau <- apply(RA, 2, function(x) cor(x, RM, method = "kendall"))
## Clayton Copula: Lower Tail Dependence 
ThetaC <- copClayton@iTau(Tau) # copula parameter Theta 
LambdaL <- copClayton@lambdaL(ThetaC) # lower tail dependence coefficients (low lambda)
## Selecting Stocks below median; inverse log-weighted and scaled
IdxBeta <- Beta < median(Beta)
WBeta <- -1 * log(abs(Beta[IdxBeta])) 
WBeta <- WBeta / sum(WBeta) * 100 #low beta portfolio
## TD
IdxTD <- LambdaL < median(LambdaL)
WTD <- -1 * log(LambdaL[IdxTD])
WTD <- WTD / sum(WTD) * 100
Intersection <- sum(names(WTD) %in% names(WBeta)) /
  length(WBeta) * 100
## Out-of-Sample Performance
RMo <- returnseries(INDTRACK6[260:290, 1], method = "discrete",
                    percentage = FALSE) + 1 #If you replace NA with 100, you know what would happen with your 100 dollar investment
RAo <- returnseries(INDTRACK6[260:290, -1], method = "discrete",
                    percentage = FALSE) + 1 #
## Benchmark

RMo[1] <- 100
RMEquity <- cumprod(RMo)
## Low Beta
LBEquity <- RAo[, IdxBeta]
LBEquity[1, ] <- WBeta
LBEquity <- rowSums(apply(LBEquity, 2, cumprod))
## TD 
TDEquity <- RAo[, IdxTD]
TDEquity[1, ] <- WTD
TDEquity <- rowSums(apply(TDEquity, 2, cumprod))
## Collecting results
y <- cbind(RMEquity, LBEquity, TDEquity)

showEx("C11R4")
## Time series plots of equity curves
par(mfrow = c(1, 1))
plot(RMEquity, type = "l", ylim = range(y), ylab = "Equity Index",
     xlab = "Out-of-Sample Periods")
lines(LBEquity, lty = 2)
lines(TDEquity, lty = 3)
legend("topleft",
       legend = c("S&P 500", "Low Beta", "Lower Tail Dep."),
       lty = 1:3)
## Bar plot of relative performance
RelOut <- rbind((LBEquity / RMEquity - 1) * 100,
                (TDEquity / RMEquity - 1) * 100)
RelOut <- RelOut[, -1]
barplot(RelOut, beside = TRUE, ylim = c(-5, 17),
        names.arg = 1:ncol(RelOut),
        legend.text = c("Low Beta", "Lower Tail Dep."),
        args.legend = list(x = "topleft"))
abline(h = 0)
box()
par(mfrow = c(1, 1))

showEx("C11R5")
library(PerformanceAnalytics)
## Key measures (ex post)
RAdec <- RA / 100
RALB <- RAdec[, names(WBeta)]
RATD <- RAdec[, names(WTD)]
LbStd <- StdDev(rowSums(RALB * WBeta / 100)) * 100
TdStd <- StdDev(rowSums(RATD * WTD / 100)) * 100
LbES95 <- abs(ES(R = rowSums(RALB * WBeta / 100),
                 method = "gaussian")) * 100
TdES95 <- abs(ES(R = rowSums(RATD * WTD / 100),
                 method = "gaussian")) * 100
LbDr <- dr(WBeta, Sigma = cov(RALB))
TdDr <- dr(WTD, Sigma = cov(RATD))
LbCr <- cr(WBeta, Sigma = cov(RALB))
TdCr <- cr(WTD, Sigma = cov(RATD))
## Key measure (ex ante)
LbRetO <- returnseries(LBEquity, percent = FALSE, trim = TRUE)
LbRetO <- timeSeries(LbRetO, charvec = 1:30)
TdRetO <- returnseries(TDEquity, percent = FALSE, trim = TRUE)
TdRetO <- timeSeries(TdRetO, charvec = 1:30)
BmRetO <- timeSeries(RMo[-1] - 1, charvec = 1:30)
km <- function(pf, bm, scale){
  ra <- Return.annualized(pf, scale = scale) * 100
  ir <- InformationRatio(pf, bm, scale = scale)
  upr <- UpDownRatios(pf, bm, method = "Capture", side = "Up")
  dnr <- UpDownRatios(pf, bm, method = "Capture", side = "Down")
  res <- c(ra, ir, upr, dnr)
  names(res) <- c("Return", "IR", "UpRatio", "DownRatio")
  return(res)
}
LbKM <- km(LbRetO, BmRetO, scale = 52)
TdKM <- km(TdRetO, BmRetO, scale = 52)