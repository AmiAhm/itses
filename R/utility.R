# Contains utility functions used in implementation.
#########################
#    Mean Estimators
#########################
soft.threshold.estimator <- function(y, lambda) {
  v <- abs(y)-lambda
  v <- v*(v>0)
  sign(y)*v
}
hard.threshold.estimator <- function(y, lambda) {
  (abs(y)>lambda)*y
}
#########################
#    Loss functions
#########################
get.l2.loss <- function(theta, theta.hat) {
  sum((theta-theta.hat)^2)
}
loss.w.ht <- function(theta, y, lambda) {
  get.l2.loss(theta, hard.threshold.estimator(y, lambda))
}
loss.w.st <- function(theta, y, lambda) {
  get.l2.loss(theta, soft.threshold.estimator(y, lambda))
}
#########################
#    Risk functions
#########################
#' Risk with soft-threhsold.
#' Assuming unit variance.
#' @noRd
risk.st <- function(theta, lambda) {
  if(!is.finite(lambda)) return(sum(theta^2))
   f <- function(theta, lambda) 1+lambda^2+ (theta^2-lambda^2-1)*(pnorm(lambda-theta)- pnorm(-lambda-theta)) -
    (lambda-theta)*dnorm(lambda+theta)- (lambda+theta)*dnorm(lambda-theta)
    sum(sapply(theta, function(theta) f(theta, lambda)))
}
#' Derivative of risk w/st wrt. threshold.
#' Assuming unit variance.
#' @noRd
dlambda.risk.st  <- function(theta, lambda) {
  if(!is.finite(lambda)) return(0)
  f <- function(theta, lambda) 2*lambda* (1-pnorm(lambda-theta)+ pnorm(-lambda-theta))-
     2*(dnorm(lambda+theta)+dnorm(lambda-theta))
  sum(sapply(theta, function(theta) f(theta, lambda)))
}
#' Double derivative of risk w/st wrt. threshold.
#' Assuming unit variance.
#' @noRd
d2lambda.risk.st <- function(theta, lambda) {
  if(!is.finite(lambda)) return(0)
  f <- function(theta, lambda) 2*(pnorm(-lambda+theta)+pnorm(-lambda-theta))+
     2*theta*(dnorm(lambda+theta)-dnorm(lambda-theta))
  sum(sapply(theta, function(theta) f(theta, lambda)))
}
#' Risk w/ht.
#' Assuming unit variance.
#' @noRd
risk.ht <- function(theta, lambda) {
  if(!is.finite(lambda)) return(sum(theta^2))
  f <- function(theta, lambda) theta^2 *(pnorm(lambda-theta) - pnorm(-lambda-theta))+
    2-pnorm(lambda-theta)- pnorm(lambda+theta)+ (lambda-theta)*
    dnorm(lambda-theta)+(lambda+theta)*dnorm(lambda+theta)
  sum(sapply(theta, function(theta) f(theta, lambda)))
}
#' Derivative of risk w/ht wrt. threshold.
#' Assuming unit variance.
#' @noRd
dlambda.risk.ht  <- function(theta, lambda) {
  if(!is.finite(lambda)) return(0)
  f <- function(theta, lambda) (-lambda^2+2*lambda*theta)*dnorm(lambda-theta)-
    (lambda^2+2*lambda*theta)*dnorm(lambda+theta)
  sum(sapply(theta, function(theta) f(theta, lambda)))
}
#' Double derivative of risk w/st wrt. threshold.
#' Assuming unit variance.
#' @noRd
d2lambda.risk.ht <- function(theta, lambda) {
  if(!is.finite(lambda)) return(0)
  f <- function(theta, lambda) (lambda^3-3*(lambda^2)*theta+2*lambda*theta^2- 2*lambda+2*theta)*dnorm(lambda-theta)+
     (lambda^3+3*(lambda^2)*theta+2*(theta^2)*lambda-2*lambda-2*theta)* dnorm(lambda+theta)
  sum(sapply(theta, function(theta) f(theta, lambda)))
}
