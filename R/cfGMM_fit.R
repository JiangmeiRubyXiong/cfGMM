#' Closed-form Gamma Mixture Model estimation.
#'
#' This function returns estimated parameters for Gamma Mixture Model using EM algoritm with cclosed-form estimators.
#' @param x Data vector
#' @param k Number of components. Integer that is greater or equal to 2.
#' @param weights A vector of length equal to the length of x with weights.
#' @param alpha (optional) shape starting value. Can be either a vector with length k or a constant. If the input is vector and there is NA in the vector, the NAs will be filled with MOM estimator. If the input is a constant, the initial value of a will be a vector of length k. If not specified, the initial value will be a vector of values near MOM estimator. All numerical values input should be greater than 0.
#' @param beta (optional) scale starting value. Same specification as init.val.a
#' @param lambda (optional) Can be either a matrix with nrow=length(xtor) and ncol=k, or a vector with length k. NA value matrix input is not allowed. If the input is vector, the initial value would be a matrix with nrow=length(xtor) and repeated rows of the vector input. If the input is not specified, the initial phi matrix is a nrow=length(xtor) and ncol=k with all elements equal to 1/k.
#' @param nbins (optional) To reduce computing time bin x to nbins bins.
#' @param n.rerun (optional) Integer greater or equal to 1. This is the number of times the EM algorithm starting at different values, and the result returned will be from the run with highest mean log-likelihood. When there isn't constraint, the complexity of model with k>2 calls for more restarts for result accuracy. Defaults to be 4.
#' @param diff.conv (optional) default to 1e-5
#' @param max.iter (optional) default to 1e3
#' @param max.restarts (optional) The number of restart within each run. EM algorithm restarts with another initial value when it reaches \code{max.iter} but has not converged. Default to 20.
#' @param max.comp (optional) Whether k indicates maximum number of components or fixed number of component. Logical, \code{TRUE} indicates k is the maximum and \code{FALSE} indicates k is fixed. Default to be \code{FALSE}.
#' @param min.lambda (optional) The minimum proportion assigned to each component. Defaults to 1e-5.
#' @param constraint NULL or matrix of dimension (k,2), each row is c(lower bound, upper bound) for one component. If one end is not specified, let it be NA.
#' @importFrom mixtools gammamix.init
#' @importFrom stats dgamma uniroot
#' @return
#' \item{x}{The raw data.}
#' \item{lambda}{The final mixing proportions.}
#' \item{gamma.pars}{A 2xk matrix where each column provides the component estimates of \code{alpha} and \code{beta}.}
#' \item{loglik}{The final log-likelihood.}
#' \item{posterior}{An nxk matrix of posterior probabilities for observations.}
#' \item{all.loglik}{A vector of each iteration's log-likelihood.}
#' \item{ft}{Name of the function}
#' \item{convergence}{Indicator variable of whether the last run converges}
#' \item{nrestart}{Number of restarts}
#'
#' @examples
#' ### Analyze a 2-component dataset
#'
#' ## data generation
#' set.seed(1)
#' phi <- c(0.3, 0.7) # mixing proportion
#' a <- c(0.5, 8) #shape
#' b <- c(1/2, 1/3) #scale
#' n <- 10000 #data size
#' ind <- sample(c(1,2) ,size=n, replace = TRUE, prob = phi)
#' data.gamma <- c(rgamma(sum(ind==1), shape=a[1], rate=b[1]),
#' rgamma(sum(ind==2), shape=a[2], rate=b[2]))
#' ## run the model
#' out <- cfGMM(data.gamma, k=2)
#' out[["gamma.pars"]]
#' out[["lambda"]]
#' @export


cfGMM <- function(x, k, weights=NULL, alpha=NULL, beta=NULL, lambda=NULL, nbins=NULL, n.rerun=4, diff.conv=1e-6, max.iter=1e3, max.restarts = 20, max.comp=FALSE, min.lambda=1e-4, constraint=NULL)
{
  #try catch with uniroot
  unirootTryCatch <- function(bound){
    tryCatch(
      {result=1/(uniroot(f=optimize_derivative, interval = c(1e-5,1e6), mk=bound, phi.vec=phi_out[,i], x=x, weights=weights)$root)
      return(result)
      },
      error=function(e){
        return(NA)
      })
  }
  if(any(x==0)){stop('Data must be strictly positive.')}
  if(k<1){stop('k must be greater or equal to 1')}
  init.k <- k
  init.constraint <- constraint

  n <- length(x)
  y = x
  weights.original <- weights
  if(!is.null(nbins)){
    weights = do.call(rbind,
                      by(x,
                         cut(x, breaks = nbins, include.lowest = TRUE),
                         function(sub) c(mean(sub), length(sub)) ))
    x = weights[,1]
    weights = weights[,2][order(x)]
  } else if(is.null(weights)){
    sx = table(x)
    if(length(sx)<n){
      weights=as.numeric(sx)
      x = as.numeric(names(sx))
    }
  }
  if(!is.null(weights.original)) {
    weights = (weights.original/sum(weights.original)*n)[order(x)]
  } else {
    weights = rep(1, n)
  }
  x <- sort(x)
  # Getting started with initial values. If not specified, initial value will be replaced with a value near MOM estimator.
  # removes 4th element, which is k
  result4 <- list()
  likelihood4 <- rep(NA,n.rerun)
  all.loglik4 <- list()
  j=1
  for(j in 1:n.rerun){
    param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
    if(is.null(dim(param_current))){param_current <- matrix(param_current, nrow=1)}
    param_current[,3] <- 1/param_current[,3] #convert to shape/rate
    mode <- (param_current[,2]-1)/param_current[,3]
    param_current <- param_current[order(mode),]
    param.init <- param_current
    #print("param.init")
    #print(param.init)
    m.diff.conv <- Inf
    #loop till convergence.
    iter=0
    nrestarts=0

    m.log_lik_new <-  -Inf
    m.loglik.all <- NULL

    # loops
    while(m.diff.conv>diff.conv & iter<=max.iter & nrestarts<=max.restarts){
      # loop parameters
      iter = iter + 1
      param_prev <- param_current
      m.log_lik_old <- m.log_lik_new

      # compute the mixture component probabilities for each observation
      if(is.null(dim(param_current))){param_current <- matrix(param_current, nrow=1)}
      phi_out <- apply(param_current, 1, function(param) dgamma(x,shape=param[2],rate=param[3])*param[1] )
      # divide by their sum (Bayes rule)
      phi_out <- sweep(phi_out, 1, rowSums(phi_out), FUN = '/' )

      # step 3
      # Compute
      param_current[,1] <- colSums(phi_out*weights)
      # this is the denominator?
      param_current[,2] <- param_current[,1] * colSums(phi_out*x*log(x)*weights)-
        colSums(phi_out*log(x)*weights)*colSums(phi_out*x*weights)
      param_current[,3] <- ((param_current[,1])^2)/param_current[,2]
      param_current[,2] <- (param_current[,1] * colSums(phi_out*x*weights))/param_current[,2]
      param_current[,1] <- param_current[,1]/n

      #print(param_current) # prints old parameters
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
              param_current[i,3] <- unirootTryCatch(upper)
              if(is.na(param_current[i,3]))m.log_lik_new <- NaN
              param_current[i,2] <- param_current[i,3]* upper +1
              #print(c(param_current[i,3], mode))
            } else if (mode < lower){
              param_current[i,3] <- unirootTryCatch(lower)
              if(is.na(param_current[i,3]))m.log_lik_new <- NaN
              param_current[i,2] <- param_current[i,3]* lower+1
              #print(c(param_current[i,3], mode))
            }
          }
        }
        if(is.null(dim(param_current))){param_current <- matrix(param_current, nrow=1)}
        log_lik_new <- apply(param_current, MARGIN = 1,
                             FUN = function(a, x){a[1]*dgamma(x, shape =a[2] , rate=a[3])}, x=x )
        m.log_lik_new <- mean(log(rowSums(log_lik_new)))
        m.loglik.all <- c(m.loglik.all, m.log_lik_new)
        if(is.nan(m.log_lik_new)){
          message('Likelihood NaN, restart')
          k <- init.k
          constraint <- init.constraint
          message('Bad parameter values, restarting with new initial parameters.')
          param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
          if(is.null(dim(param_current))){param_current <- matrix(param_current, nrow=1)}
          param_current[,3] <- 1/param_current[,3] #convert to shape/rate
          mode <- (param_current[,2]-1)/param_current[,3]
          param_current <- param_current[order(mode),]
          iter <- 0
          nrestarts<-nrestarts+1
          m.log_lik_new <-  -Inf
          m.diff.conv <- Inf
        } else {
          m.diff.conv <- abs(m.log_lik_new-m.log_lik_old)#max(abs((param_current-param_prev)/param_current))
          #print("mean log likelihood")
          #print(m.log_lik_new)

          if(any(param_current[,1]<=min.lambda) & max.comp==TRUE){
            #remove one component
            message('Probability too small, removing one component')
            idx <- which.min(param_current[,1])
            #print(idx)
            k=k-1
            param_current <- param_current[-idx,]
            constraint <- constraint[-idx,]
            #print(param_current)
          } else if(any(param_current[,1]<=min.lambda)){
            #restart
            k <- init.k
            constraint <- init.constraint
            message('Bad parameter values, restarting with new initial parameters.')
            param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
            if(is.null(dim(param_current))){param_current <- matrix(param_current, nrow=1)}
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
        if(is.null(dim(param_current))){param_current <- matrix(param_current, nrow=1)}
        param_current[,3] <- 1/param_current[,3] #convert to shape/rate
        mode <- (param_current[,2]-1)/param_current[,3]
        param_current <- param_current[order(mode),]
        iter <- 0
        nrestarts<-nrestarts+1
        m.log_lik_new <-  -Inf
        m.diff.conv <- Inf
      }
      if(k<1){
        message('no more components, restart')
        k <- init.k
        constraint <- init.constraint
        param.init <- param_current <- simplify2array(mixtools::gammamix.init(x, lambda, alpha, beta, k = k)[1:3])
        if(is.null(dim(param_current))){param_current <- matrix(param_current, nrow=1)}
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
        if(is.null(dim(param_current))){param_current <- matrix(param_current, nrow=1)}
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
    all.loglik4[[j]] <- m.loglik.all
    result4[[j]] <- list(param_at_conv=param_current, param.init=param.init, z=phi_out, convergence=convergence, nrestarts=nrestarts)
  }
  convs <- c(result4[[1]][["convergence"]],result4[[2]][["convergence"]],result4[[3]][["convergence"]],result4[[4]][["convergence"]])
  if(any(convs)){which.conv <- which(convs)
  final.result <- result4[[which.max(likelihood4[which.conv])]]
  final.loglik.all <- all.loglik4[[which.max(likelihood4)]]
  final.lambda <- final.result[["param_at_conv"]][,1]
  final.pars <- matrix(final.result[["param_at_conv"]][,2:3], nrow=2, byrow=TRUE)
  final.pars[2,] <- 1/final.pars[2,]
  dimnames(final.pars) <- list(c("alpha", "beta"),paste("comp.", 1:ncol(final.pars), sep=""))
  final.lik <- max(likelihood4)
  final.z <- final.result[["z"]]
  weights <- as.integer(round(weights,0))
  final.z <- apply(final.z, 2, function(x){rep(x, weights)} )
  final.z <- final.z[rank(y, ties.method = "random"),]
  final.conv <- final.result[["convergence"]]
  final.restart <- final.result[["nrestarts"]]

  output <- output <- list(x=y, lambda=final.lambda, gamma.pars=final.pars, loglik=final.lik, posterior=final.z, all.loglik=final.loglik.all, convergence = final.conv, nrestart=final.restart,ft="cfGMM")
  }else{output="Did not converge. Consider increase max.iter or check input."}
  return(output)
}
