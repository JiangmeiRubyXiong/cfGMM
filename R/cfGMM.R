
#--------------------------------------------------------------------------
# Main Function -----------------------------------------------------------
#' @x: data vector
#' @k: number of components. Interger that is greater or equal to 2
#' @init.val.a: (optional) shape starting value. Can be either a vector with length k or a constant. If the input is vector and there is NA in the vector, the NAs will be filled with MOM estimator. If the input is a constant, the initial value of a will be a vector of length k. If not specified, the initial value will be a vector of values near MOM estimator. All numerical values input should be greater than 0.
#' @init.val.b: (optional) scale starting value. Same specification as init.val.a
#' @init.val.phi: (optional) Can be either a matrix with nrow=length(xtor) and ncol=k, or a vector with length k. NA value matrix input is not allowed. If the input is vector, the initial value would be a matrix with nrow=length(xtor) and repeated rows of the vector input. If the input is not specified, the initial phi matrix is a nrow=length(xtor) and ncol=k with all elements equal to 1/k.
#' @diff.conv: (optional) default to 1e-5
#' @max.iter: (optional) default to 1e3
#' @max.comp: (optional) whether k indicates maximum number of components or fixed number of component. Logical, TRUE indicates k is the maximum and FALSE indicates k is fixed. Default to be FALSE.
#' @min.lambda: 1e-5
#' @constraint: NULL or matrix of dimension (k,2), each row correspond to alpha_lower, alpha_upper, beta_lower, beta_upper for one component. If one end is not specified, let it be NA.
#' @inner.loop: integer greater or equal to 1. When there isn't constraint, the more inner.loop, the more accurate the result, the longer the algorithm takes. Defaults to be 4.
#' @importFrom mixtools gammamix.init
#' @export


cfGMM <- function(x, k, alpha=NULL, beta=NULL, lambda=NULL, diff.conv=1e-6, max.iter=1e3, maxrestarts = 20, max.comp=FALSE, min.lambda=1e-4, constraint=NULL, inner.loop=4)
{
  n <- length(x)
  if(any(x==0)){stop('Data must be strictly positive.')}
  init.k <- k
  init.constraint <- constraint
  # Getting started with initial values. If not specified, initial value will be replaced with a value near MOM estimator.
  # removes 4th element, which is k
  result4 <- list()
  likelihood4 <- rep(NA,inner.loop)
  for(j in 1:inner.loop){
    param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
    param_current[,3] <- 1/param_current[,3] #convert to shape/rate
    mode <- (param_current[,2]-1)/param_current[,3]
    param_current <- param_current[order(mode),]
    param.init <- param_current
    print("param.init")
    print(param.init)
    m.diff.conv <- Inf
    #loop till convergence.
    iter=0
    nrestarts=0

    m.log_lik_new <-  -Inf

    # loops
    while(m.diff.conv>diff.conv & iter<=max.iter & nrestarts<=maxrestarts){
      # loop parameters
      iter = iter + 1
      param_prev <- param_current
      m.log_lik_old <- m.log_lik_new

      # compute the mixture component probabilities for each observation
      phi_out <- apply(param_current, 1, function(x) dgamma(x,shape=x[2],rate=x[3])*x[1] )
      # divide by their sum (Bayes rule)
      phi_out <- sweep(phi_out, 1, rowSums(phi_out), FUN = '/' )

      # step 3
      # Compute
      param_current[,1] <- colSums(phi_out)
      # this is the denominator?
      param_current[,2] <- param_current[,1] * colSums(phi_out*x*log(x))-
        colSums(phi_out*log(x))*colSums(phi_out*x)
      param_current[,3] <- ((param_current[,1])^2)/param_current[,2]
      param_current[,2] <- (param_current[,1] * colSums(phi_out*x))/param_current[,2]
      param_current[,1] <- param_current[,1]/n

      print(param_current) # prints old parameters
      if(!any(is.infinite(param_current)|is.nan(param_current))){
        if(!is.null(constraint)){
          #TO ADD:check that lower bound is lower than upper bound
          # what if there is only one sided constraint, or a,b doesn't have same constraint
          for (i in 1:k){
            lower <- (constraint[i,1])
            upper <- (constraint[i,2])
            mode <- (param_current[i,2]-1)/param_current[i,3]
            #used the function at the beginning to solve for beta
            if (mode > upper){
              param_current[i,3] <- 1/(uniroot(f=optimize_derivative, interval = c(1e-5,1e6), mk=upper, phi.vec=phi_out[,i], data.vec=x)$root)
              param_current[i,2] <- param_current[i,3]* upper +1
              print(c(param_current[i,3], mode))
            } else if (mode < lower){
              param_current[i,3] <- 1/(uniroot(f=optimize_derivative, interval = c(1e-5,1e6), mk=lower, phi.vec=phi_out[,i], data.vec=x)$root)
              param_current[i,2] <- param_current[i,3]* lower+1
              print(c(param_current[i,3], mode))
            }
          }
        }

        log_lik_new <- apply(param_current, MARGIN = 1,
                             FUN = function(a, x){a[1]*dgamma(x, shape =a[2] , rate=a[3])}, x=x )
        m.log_lik_new <- mean(log(rowSums(log_lik_new)))

        if(is.nan(m.log_lik_new)){
          message('Likelihood NaN, restart')
          k <- init.k
          constraint <- init.constraint
          message('Bad parameter values, restarting with new initial parameters.')
          param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
          param_current[,3] <- 1/param_current[,3] #convert to shape/rate
          mode <- (param_current[,2]-1)/param_current[,3]
          param_current <- param_current[order(mode),]
          iter <- 0
          nrestarts<-nrestarts+1
          m.log_lik_new <-  -Inf
          m.diff.conv <- Inf
        } else {
          m.diff.conv <- abs(m.log_lik_new-m.log_lik_old)#max(abs((param_current-param_prev)/param_current))
          print("mean log likelihood")
          print(m.log_lik_new)

          if(any(param_current[,1]<=min.lambda) & max.comp==TRUE){
            #remove one component
            message('Probability too small, removing one component')
            idx <- which.min(param_current[,1])
            print(idx)
            k=k-1
            param_current <- param_current[-idx,]
            constraint <- constraint[-idx,]
            print(param_current)
          } else if(any(param_current[,1]<=min.lambda)){
            #restart
            k <- init.k
            constraint <- init.constraint
            message('Bad parameter values, restarting with new initial parameters.')
            param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
            param_current[,3] <- 1/param_current[,3] #convert to shape/rate
            mode <- (param_current[,2]-1)/param_current[,3]
            param_current <- param_current[order(mode),]
            iter <- 0
            nrestarts<-nrestarts+1
            m.log_lik_new <-  -Inf
            m.diff.conv <- Inf
          }

        }


      } else {
        #restart
        k <- init.k
        constraint <- init.constraint
        message('Bad parameter values, restarting with new initial parameters.')
        param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
        param_current[,3] <- 1/param_current[,3] #convert to shape/rate
        mode <- (param_current[,2]-1)/param_current[,3]
        param_current <- param_current[order(mode),]
        iter <- 0
        nrestarts<-nrestarts+1
        m.log_lik_new <-  -Inf
        m.diff.conv <- Inf
      }
      if(k==1){
        message('no more components, restart')
        k <- init.k
        constraint <- init.constraint
        param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
        param_current[,3] <- 1/param_current[,3] #convert to shape/rate
        mode <- (param_current[,2]-1)/param_current[,3]
        param_current <- param_current[order(mode),]
        iter <- 0
        nrestarts<-nrestarts+1
        m.log_lik_new <-  -Inf
        m.diff.conv <- Inf
      }


      if(iter==max.iter){
        k <- init.k
        constraint <- init.constraint
        param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
        param_current[,3] <- 1/param_current[,3] #convert to shape/rate
        mode <- (param_current[,2]-1)/param_current[,3]
        param_current <- param_current[order(mode),]
        iter <- 0
        nrestarts<-nrestarts+1
        m.log_lik_new <-  -Inf
        m.diff.conv <- Inf
      }
    }
    convergence <-  m.diff.conv <= diff.conv
    likelihood4[j] <- m.log_lik_new
    result4[[j]] <- list(param_at_conv=param_current, param.init=param.init, z=phi_out, convergence=convergence, nrestarts=nrestarts)
  }
  final.result <- list(result4[[which.max(likelihood4)]],likelihood4)

  return(final.result)
  #return(list(param_at_conv=param_current, params=param.init, z=phi_out, convergence=convergence,nrestarts=nrestarts))
}

```
