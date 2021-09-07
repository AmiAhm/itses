#' Function that does an iterative simulation iteration
#' This function assumes y~N(theta, 1^2)
#' @noRd
iter.sampling <- function(y,
                          init.lambda = "median",
                          sd = 1,
                          b = 20,
                          k = 15,
                          lambda.grid = NULL,
                          max_num_iters = 10,
                          tol = 1e-8,
                          method = "ST",
                          noisetype = NULL,
                          debug = FALSE,
                          min.threshold = 0,
                          max.threshold = Inf
) {

  # Initiate lambda if needed
  lambda <- get.init.lambda(y, init.lambda = init.lambda)
  # Ensure that the threshold is within bounds.
  lambda <- get.thresholded.lambda(lambda, min.threshold, max.threshold)

  if (debug) print("Starting iterative step using sampling.")

  # Initiate grid of threhsolds
  if (is.null(lambda.grid)) {
    lambda.grid <- get.lambda_grid(y,
                                   k = k,
                                   min.threshold = min.threshold,
                                   max.threshold = max.threshold)
  }

  # Add starting threhsold to grid if needed.
  if (!(lambda %in% lambda.grid)) {
    lambda.grid <- c(lambda.grid, lambda)
    # Ensure the thresholds are sorted.
    lambda.grid <- sort(lambda.grid, method = "quick")
  }

  # Fetch needed estimators.
  if(is.character(method)){
    if(method == "ST") {
    theta.t <- soft.threshold.estimator(y, lambda)
    risk.fun <- loss.w.st
    }else if(method == "HT") {
      theta.t <- hard.threshold.estimator(y, lambda)
      risk.fun <- loss.w.ht
    }else{
      stop("Invalid method.")
    }
    } else{
      theta.t <- method(y, lambda)
      risk.fun <- function(theta, y, lambda) get.l2.loss(theta, method(y, lambda))
    }


  # Sample data
  n <- length(y)
  if (debug) print(paste0("Now sampling: ", b * n,"points"))

  if(noisetype == "gaussian"){
    y.star <- rnorm(b * n, mean = theta.t, sd = sd)
    y.star <- matrix(y.star, ncol = b)
  }else{
    #warning("0 centered")
    #noise <- 1+rnorm(b * n, mean = 0, sd = sd)
    #y.star <- matrix(noise, ncol = b)*theta.t
    if(debug)  print("Using custom noise")
    y.star <- noisetype$sample(b, theta.t)
  }

  # Calculate risk over all samples and add at different thresholds
  risks <- sapply(lambda.grid, function(lambda) {
    mean(sapply(1:b, function(i) risk.fun(theta.t, y.star[, i], lambda)))
  })

  # Spline search needs to be on the fintie grid.
  xmax <- max.threshold
  if (!is.finite(max.threshold)) {
    # If non-finite max-threshold set to maximum of absolute of observations.
    xmax <- max(abs(y))
  }

  # Interpolate with first cubic spline
  finite.lambdas <- is.finite(lambda.grid)  # Fit spline only on finite points
  f <- splinefun(lambda.grid[finite.lambdas], risks[finite.lambdas])

  # Start steps fit spline interpolation.
  last.lambda <- -Inf # Initiate first threshold.
  for (i in 1:max_num_iters) {
    if (debug) print(paste("Search", i))

    # Use optimize to find spline minimizer.
    opt <- optimize(f, interval = c(min(min.threshold), xmax))
    next.lambda <- opt$minimum # Extract minimium

    # If no change then stop runs
    if (abs(next.lambda - last.lambda) < tol) {
      if (debug) print("Same as last iteration in spline search.")
      break()
    }
    # Find true sample risk at threshold minimizing risk.
    next.risk <- mean(sapply(1:b,
                             function(i)
                               risk.fun(theta.t, y.star[, i], next.lambda)))
    # Store results of calculation
    risks <- c(risks, next.risk)
    lambda.grid <- c(lambda.grid, next.lambda) # Add to grid of threhsolds

    # Interpolate with new cubic spline, based on all new data.
    finite.lambdas <- is.finite(lambda.grid)
    f <- splinefun(lambda.grid[finite.lambdas], risks[finite.lambdas])

    last.lambda <- next.lambda # Update
  }
  # End spline interpolation, now wrap up and return

  # Find risk minimizing threshold.
  # Only use thresholds within bounds.
  valid.thresholds <- (lambda.grid >= min.threshold) &
    (lambda.grid <= max.threshold)
  j <- which.min(risks[valid.thresholds]) # Find index.
  lambda <- lambda.grid[valid.thresholds][j] # Find threshold itself
  risk <- risks[valid.thresholds][j] # Find connected risk

  # Return results sorted by threhsold.
  ord <- order(lambda.grid)  # Find order of threhsolds.
  lambda.grid <- lambda.grid[ord] # Sort thresholds
  risks <- risks[ord] # Sort connected risks
  list(lambda = lambda,
       risk = risk,
       lambdas = lambda.grid,
       risks = risks,
       y.star = y.star)

}
