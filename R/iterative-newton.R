#' Function that takes an iterative numeric step
#' Assuming y of unit variance.
#' @noRd
iter.newton <- function(y,
                          init.lambda = "median",
                          max_num_iters = 10,
                          tol = 1e-8,
                          method = "ST",
                          debug = FALSE,
                          min.threshold = 0,
                         k = 10,
                         max.threshold = Inf) {

  # Select threshold
  lambda <- get.init.lambda(y, init.lambda = init.lambda)

  if(debug) print(paste("Starting iterative with newton", lambda))

  # Find initial means, assume these to be true.
  if(method == "ST") {
    theta <- soft.threshold.estimator(y, lambda)
  }else if(method == "HT") {
    theta <- hard.threshold.estimator(y, lambda)
  }else{
    stop("Invalid method given.")
  }

  # Find oracle threshold given assumed true means.
  get.risk.oracle.threshold(theta,
                           max_num_iters = max_num_iters,
                          tol = tol,
                          method = method,
                          debug = debug,
                            k = k,
                          min.threshold = min.threshold,
                         max.threshold = max.threshold)


}

