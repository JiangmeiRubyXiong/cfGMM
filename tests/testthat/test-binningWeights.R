set.seed(1)
phi <- c(0.3, 0.7) # mixing proportion
a <- c(0.5, 8) #shape
b <- c(1/2, 1/3) #scale
n <- 10000 #data size
nbins = 256
ind <- sample(c(1,2) ,size=n, replace = TRUE, prob = phi)
data.gamma <- c(rgamma(sum(ind==1), shape=a[1], rate=b[1]),
                rgamma(sum(ind==2), shape=a[2], rate=b[2]))
## run the model
out <- cfGMM(data.gamma, k=2, diff.conv = 1e-9)
out2 <- cfGMM(data.gamma, weights=rep(1, length(data.gamma)), k=2, diff.conv = 1e-9)
# faster fit on weighted version of data

out[["gamma.pars"]]
out[["lambda"]]

out2[["gamma.pars"]]
out2[["lambda"]]

# Check that binned fit is the same as regular fit
x = do.call(rbind, by(data.gamma, cut(data.gamma, breaks = nbins, include.lowest = TRUE), function(x) c(mean(x), length(x)) ))
out3 = cfGMM(x[,1], k=2, weights=x[,2], diff.conv = 1e-9)
out4 = cfGMM(rep(x[,1], x[,2]), k=2, weights=rep(1, length(data.gamma)), diff.conv = 1e-9)

out3[["gamma.pars"]]
out3[["lambda"]]

out4[["gamma.pars"]]
out4[["lambda"]]

# Check that constrained fit with weights is close to binned constrained fit with weights
constraint = matrix(c(-Inf,2, 24, Inf), byrow=TRUE, nrow=2)
out5 <- cfGMM(data.gamma, weights=rep(1, length(data.gamma)), k=2, constraint=constraint, diff.conv = 1e-9)
x = do.call(rbind, by(data.gamma, cut(data.gamma, breaks = 5000, include.lowest = TRUE), function(x) c(mean(x), length(x)) ))
out6 = cfGMM(rep(x[,1], x[,2]), k=2, constraint=constraint, diff.conv = 1e-9)

out5[["gamma.pars"]]
out5[["lambda"]]
out6[["gamma.pars"]]
out6[["lambda"]]

# Check constrained binned fit is the same as constrained regular fit
out3 = cfGMM(x[,1], k=2, weights=x[,2], constraint=constraint, diff.conv = 1e-9)
out4 = cfGMM(rep(x[,1], x[,2]), k=2, weights=rep(1, length(data.gamma)), constraint=constraint, diff.conv = 1e-9)

out3[["gamma.pars"]]
out3[["lambda"]]

out4[["gamma.pars"]]
out4[["lambda"]]
