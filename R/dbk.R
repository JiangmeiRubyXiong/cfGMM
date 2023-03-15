#' Function with the constraint
#'
#' @param bk scale
#' @param mk median
#' @param x data vector
#' @param phi.vec phi
#' @param weights weights for each observation
#' @export

optimize_derivative <- function(bk, mk, x, phi.vec, weights){
  sum(weights * phi.vec*(mk+bk-mk*log(bk)-mk*digamma(mk/bk+1)+mk*log(x)-x))
}
