#' Function with the constraint
#'
#' @param bk scale
#' @param mk median
#' @param x data vector
#' @param phi.vec phi
#' @export

optimize_derivative <- function(bk, mk, x, phi.vec){
  sum(phi.vec*(mk+bk-mk*log(bk)-mk*digamma(mk/bk+1)+mk*log(x)-x))
}
