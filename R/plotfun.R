#' Plot Function
#'
#' @param cfGMM.list output from function cfGMM
#' @param zero_include logical. Whether to include 0 in the dataset to histogram, if any.
#' @param breaks integer. Number of breaks in the histogram
#' @importFrom graphics hist lines
#' @importFrom stats dgamma
#' @return
#' \item{out}{The histogram.}
#' @export

plotGammaMix <- function(cfGMM.list, zero_include=FALSE, breaks=40){
  x <- cfGMM.list$x
  xs <- (1:100)/100*max(x)
  pars <- cfGMM.list$gamma.pars
  lambda <- cfGMM.list$lambda

  # option to include zero in histogram, but will not show in curve
  if(pars[1,1]==0 & (!zero_include)){
    x=x[x!=0]
    pars[,-1]}

  num.comp <- ncol(pars)
  ys <- matrix(NA, nrow=num.comp, ncol=length(xs))

  for(i in 1:num.comp){
    ys[i, ] <- lambda[i] * dgamma(xs, shape=pars[1,i],scale=pars[2,i])
  }

  # Generate Plot
  list.hist <- hist(x, freq=F, xlim=range(x), breaks=breaks)
  hist(x, freq=F, xlim=range(x), ylim=c(0,max(ys,  max(list.hist$density))), breaks=breaks)
  for(i in 1:num.comp){
    lines(xs, ys[i,],col=(i+2),lty=1,lwd=2)
  }
}
