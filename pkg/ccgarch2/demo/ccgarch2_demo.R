library(ccgarch2)
data(Garch)     # currency data from the Ecdat package
head(Garch)

# date index for zoo
ind <- as.Date(as.character(Garch$date), format = "%y%m%d")

# exchange rate data
exr <- zoo(Garch[, -c(1, 2, 4)], ind)

# making log-returns
r <- 100*diff(log(exr))

# initial parameter values
ini.dcc <- c(0.1, 0.8)
ini.a <- diag(cov(as.matrix(r)))
ini.A <- diag(rep(0.1, 5))
ini.B <- diag(rep(0.8, 5))






# Constant Conditional Correlation
exr.ccc <- estimateCCC(ini.a, ini.A, ini.B, data = r, model="extended")
out.ccc <- summary(exr.ccc)

# Dynamic Conditional Correlation
exr.dcc <- estimateDCC(data = r, model="extended")
out.dcc <- summary(exr.dcc)

# Corrected Dynamic Conditional Correlation
exr.cdcc <- estimateCDCC(data = r[, 1:2], model="extended")
out.cdcc <- summary(exr.cdcc)



abc <- residDiag(out.dcc)

tx1 <- cbind(abc$lb$stat[,1], abc$lb$p.val[,1])
tx2 <- cbind(abc$lb$stat[,2], abc$lb$p.val[,2])
tx3 <- cbind(abc$lb$stat[,3], abc$lb$p.val[,3])

colnames(tx1) <- c("stat", "Pr()")
colnames(tx2) <- c("stat", "Pr()")
colnames(tx3) <- c("stat", "Pr()")


cbind( printCoefmat(tx1, P.values=TRUE, has.Pvalue=TRUE, dig.tst=2, zap.ind=2),
       printCoefmat(tx2, P.values=TRUE, has.Pvalue=TRUE, dig.tst=2, zap.ind=2),
       printCoefmat(tx3, P.values=TRUE, has.Pvalue=TRUE, dig.tst=2, zap.ind=2))
