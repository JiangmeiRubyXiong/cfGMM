#' Plot Function
#'
#' @param cfGMM.list output from function cfGMM
#' @param zero_include logical. Whether to include 0 in the dataset to histogram, if any.
#' @param breaks integer. Number of breaks in the histogram
#' @importFrom graphics hist lines
#' @importFrom stats dgamma
#' @importFrom ggplot2 ggplot xlim ylab ggtitle geom_histogram stat_function stat_function geom_vline
#' @importFrom hrbrthemes theme_ipsum
#' @return
#' \item{out}{The histogram.}
#' @export

plotGammaMix <- function(cfGMM.list, zero_include=FALSE, breaks=40,main="histogram of x", threshold=0){
  x <- cfGMM.list$x
  xs <- (1:100)/100*max(x)
  pars <- cfGMM.list$gamma.pars
  posteriors <- cfGMM.list$posterior
  lambda <- colMeans(posteriors[posteriors[,2]!=0,])

  fun1 <- function(xs){lambda[1] * dgamma(xs, shape=pars[1,1],scale=pars[2,1])}
  fun2 <- function(xs){lambda[2]  * dgamma(xs, shape=pars[1,2],scale=pars[2,2])}
  plot.x <- as.data.frame(x)
  p <- ggplot(plot.x, aes(x=x)) + xlim(min(xs), max(xs))+
    ylab("")+theme_ipsum()+ggtitle(main)+
    geom_histogram(aes(y = stat(density)), bins = breaks, alpha=0.2)+
    stat_function(fun = fun1, n = 101, alpha=0.3, color="forestgreen") +
    stat_function(fun = fun2, n = 101, alpha=0.3, color="red") +
    geom_vline(xintercept=threshold, linetype="dashed", color = "gray")
  print(p)
}
