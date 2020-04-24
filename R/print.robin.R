#' print a robin object
#'
#' Print a summary of a robin object
#'
#' @param x robin object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{robin}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{robin}}
#' @export
print.robin=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(x$call), "\n")
  # cat("\nLambda:\n")
  # print(x$lambda)
  cat("\nCoefficients:\n")
  print(x$coefficient, digits)
  cat("Class:\n")
  print(class(x))
}


#' print a GxESelection object
#'
#' Print a summary of a GxESelection object
#'
#' @param x GxESelection object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{GxESelection}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{GxESelection}}
#' @export
print.GxESelection=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nMethod:\n")
  print(x$method)
  cat("\n")
  print(x$summary)
}


#' print a robin.pred object
#'
#' Print a summary of a robin.pred object
#'
#' @param x robin.pred object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{robin.pred}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{predict.robin}}
#' @export
print.robin.pred=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nPMSE:\n")
  print(x$error, digits)
  cat("\npredicted ", length(x$y.pred), " y (list component y.pred)", sep = "")
}
