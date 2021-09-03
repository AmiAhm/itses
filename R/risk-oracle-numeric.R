#' Finds optimal threshold to minimize risk given means (numeric).
#' Assuming unit variance.
#' @noRd
get.risk.oracle.threshold <- function(theta,
                                      max_num_iters = 10,
                                      tol = 1e-8,
                                      method = "ST",
                                      debug = FALSE,
                                      min.threshold = 0,
                                      max.threshold = Inf,
                                      also.check = NULL,
                                      k = 25) {
  # Select risk function and its derivatives.
  if(method == "ST") {
    risk.fun <- risk.st
    df <- dlambda.risk.st
    d2f <- d2lambda.risk.st
  }else if(method == "HT") {
    risk.fun <- risk.ht
    df <- dlambda.risk.ht
    d2f <- d2lambda.risk.ht
  }
  if(is.null(max.threshold)) {
    max.threshold <- get.max.threshold(method, theta)
  }
  # Define grid with k points.
  lambdas <- get.lambda_grid(theta, min.threshold =  min.threshold, max.threshold = max.threshold, k = k)
  # due to some convergence issues
  # we add points that we also want to check to the grid. This ensure, without needing convergence,
  # that estimated risk of the oracle is less than the risk of the added grid thresholds.
  lambdas <- sort(unique(get.thresholded.lambda(c(lambdas, also.check),
                                        min.threshold = min.threshold,
                                    max.threshold = max.threshold)
  ))
  # Calculate risks at grid
  risks <- sapply(lambdas, function(lambda) risk.fun(theta, lambda))
  # Find grid minimum
  j <- which.min(risks[is.finite(lambdas)]) # Check non-finite thresholds later
  lambda <- lambdas[is.finite(lambdas)][j]  # Check non-finite thresholds later
  if(debug) print(paste("Starting newton at:", lambda))
  for(i in 1:max_num_iters) {
    lambda.old <- lambda
    d <- df(theta, lambda)
    d2 <- d2f(theta, lambda)
    # Want to avoid strange behaviour, thus offset if 0...
    if(d2 == 0) {
      d2 <- tol
    }
    # Newton step
    lambda <- lambda - d/d2
    # Enfornce boundary
    lambda <- get.thresholded.lambda(lambda, min.threshold, max.threshold)
    lambda
    # Store data
    lambdas <- c(lambdas, lambda)
    risk <- risk.fun(theta, lambda)
    risks <- c(risks, risk)

    diff <- abs(lambda.old - lambda)  # Checks the difference to last threshold

    if(debug) print(paste0("Newton iter: ", i,
                           ", lambda:", lambda, ", diff to last: ",
                           round(diff, digits = 3)
                           )
    )
    # Return threshold if tolerance reached or return if lambda -> Inf.
    if(!is.finite(lambda)|((diff < tol) & (i > 1))) break()
  }
  # Return results
  ord <- order(lambdas)
  risks <- risks[ord]
  lambdas <- lambdas[ord]
  lambda.i <- which.min(risks)
  lambda <- lambdas[lambda.i]
  risk <- risks[lambda.i]
  list(lambda = lambda, risk = risk, lambdas = lambdas, risks = risks)
}
