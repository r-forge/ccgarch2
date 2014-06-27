library(ccgarch2)
data(EuStockMarkets, package="datasets")    # Stock indecies at major stock exchanges in EU
EUS <- EuStockMarkets
EUS <- as.matrix(EUS)
head(EUS)
ind <- 1:nrow(EUS)
head(ind)
# Stock indecies data 
EUS <- zoo(EUS, ind)
head(EUS)
# making log-returns
r <- 100*diff(log(EUS))

head(r)

# initial parameter values
ini.dcc <- c(0.1, 0.8)
ini.a <- diag(cov(as.matrix(r)))
ini.A <- diag(rep(0.1, 4))
ini.B <- diag(rep(0.8, 4))






# Constant Conditional Correlation
# exr.ccc <- estimateCCC(ini.a, ini.A, ini.B, data = r, model="extended")
exr.ccc <- estimateCCC(ini.a, ini.A, ini.B, data = r, model="extended")
out.ccc <- summary(exr.ccc)

# Dynamic Conditional Correlation
exr.dcc <- estimateDCC(data = r)
out.dcc <- summary(exr.dcc)

# Corrected Dynamic Conditional Correlation
exr.cdcc <- estimateCDCC(data = r, model = "extended")
out.cdcc <- summary(exr.cdcc)

abc <- residDiag(out.cdcc)

tx1 <- cbind(abc$lb$stat[,1], abc$lb$p.val[,1])
tx2 <- cbind(abc$lb$stat[,2], abc$lb$p.val[,2])
tx3 <- cbind(abc$lb$stat[,3], abc$lb$p.val[,3])

colnames(tx1) <- c("stat", "Pr()")
colnames(tx2) <- c("stat", "Pr()")
colnames(tx3) <- c("stat", "Pr()")


cbind( printCoefmat(tx1, P.values=TRUE, has.Pvalue=TRUE, dig.tst=2, zap.ind=2),
       printCoefmat(tx2, P.values=TRUE, has.Pvalue=TRUE, dig.tst=2, zap.ind=2),
       printCoefmat(tx3, P.values=TRUE, has.Pvalue=TRUE, dig.tst=2, zap.ind=2))
