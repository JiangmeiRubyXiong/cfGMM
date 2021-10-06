#' Function with the constraint
optimize_derivative <- function(bk, mk, data.vec, phi.vec){
  sum(phi.vec*(mk+bk-mk*log(bk)-mk*digamma(mk/bk+1)+mk*log(data.vec)-data.vec))
}
