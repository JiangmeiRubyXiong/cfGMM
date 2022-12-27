#' Plot Function
#'
#' This function makes plot using cfGMM function output.
#' @param cfGMM.list output from function cfGMM
#' @param zero_include logical. Whether to include 0 in the dataset to histogram, if any.
#' @param breaks integer. Number of breaks in the histogram
#' @param main Title of histogram
#' @param threshold Threshold specified for cfgmm fit
#' @importFrom graphics hist lines
#' @importFrom stats dgamma density
#' @importFrom ggplot2 ggplot xlim ylab ggtitle geom_histogram stat_function stat_function geom_abline aes after_stat
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom rlang .data
#' @return
#' \item{out}{The histogram.}
#' @export

plotGammaMix <- function(cfGMM.list, zero_include=FALSE, breaks=40, main="histogram of x", threshold=0){
  x <- cfGMM.list$x
  xs <- (1:100)/100*max(x)
  pars <- cfGMM.list$gamma.pars
  posteriors <- cfGMM.list$posterior
  lambda <- colMeans(posteriors[posteriors[,2]!=0,])

  fun1 <- function(xs){lambda[1] * dgamma(xs, shape=pars[1,1],scale=pars[2,1])}
  fun2 <- function(xs){lambda[2]  * dgamma(xs, shape=pars[1,2],scale=pars[2,2])}
  plot.x <- as.data.frame(x)
  p <- ggplot2::ggplot(plot.x, ggplot2::aes(x=x)) + ggplot2::xlim(min(xs), max(xs))+
    ggplot2::ylab("")+hrbrthemes::theme_ipsum()+ggplot2::ggtitle(main)+
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(stats::density)), bins = breaks, alpha=0.2)+
    ggplot2::stat_function(fun = fun1, n = 101, alpha=0.3, color="forestgreen") +
    ggplot2::stat_function(fun = fun2, n = 101, alpha=0.3, color="red") +
    ggplot2::geom_vline(xintercept=threshold, linetype="dashed", color = "gray")
  print(p)
}
