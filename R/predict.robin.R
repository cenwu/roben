#' make predictions from a robin object
#'
#' make predictions from a robin object
#'
#' @param object robin object.
#' @param X.new a matrix of new values for X at which predictions are to be made.
#' @param E.new a vector of new values for E at which predictions are to be made.
#' @param clin.new a vector or matrix of new values for clin at which predictions are to be made.
#' @param Y.new a vector of the response of new observations. If provided, the prediction mean squared error (PMSE) will be computed based on Y.new.
#' @param ... other predict arguments
#'
#' @details X.new (E.new) must have the same number of columns as X (E) used for fitting the model. If clin was provided when fit the model, clin.new
#' must not be NULL, and vice versa. The predictions are made based on the posterior estimates of coefficients in the robin object.
#' Note that the main effects of environmental exposures E are not subject to selection.
#'
#' @return  an object of class "robin.pred" is returned, which is a list with components:
#' \item{pmse}{predictions mean squared error. pmse is NULL is Y.new=NULL.}
#' \item{y.pred}{predicted values of the new observations.}
#'
#' @rdname predict.robin
#' @method predict robin
#' @seealso \code{\link{robin}}
#'
#' @examples
#' data(GxE_small)
#' test = sample((1:nrow(X)), floor(nrow(X)/5))
#' fit=robin(X[-test,], Y[-test,], E[-test,], clin[-test,], iterations = 5000)
#' predict(fit, X[test,], E[test,], clin[test,], Y[test,])
#'
#' @export
predict.robin=function(object, X.new, E.new, clin.new=NULL, Y.new=NULL, ...){

  intercept = TRUE
  dat = Data.matrix(X.new, Y.new, E.new, clin.new, intercept)
  xx = dat$xx
  y.new = dat$y
  CLC = dat$CLC

  # n = dat$n; s = dat$s
  # env = dat$env
  # size = dat$size
  # G.names = dat$G.names
  # E.names = dat$E.names
  # clin.names = dat$clin.names

  coeff = c(object$coefficient$GE)
  coeff.clc = c(object$coefficient$Int, object$coefficient$clin, object$coefficient$E)

  if(length(coeff)!=ncol(xx)){
    stop(paste("number of columns of X.new dose not match the length of the estimates."))
  }

  if(length(coeff.clc)!=ncol(CLC)){
    stop(paste("incorrect number of clinical covariates (", ncol(CLC)-1, "), supposed to be ", length(coeff.clc)-1, sep = ""))
  }

  y.pred = xx %*% coeff + CLC %*% coeff.clc
  pmse = NULL

  if("RBVS" %in% class(object)){
    pmse = sum(abs(y.new - y.pred))/length(y.new)
  }else{
    pmse = sum((y.new - y.pred)^2)/length(y.new)
  }

  pred = list(pmse=pmse, y.pred=y.pred)
  class(pred) = "robin.pred"
  pred
}


