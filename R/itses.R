#' Iterative Sparse Estimation
#'
#' Performs iterative sparse estimation of many normal means either using
#' the soft threshold estimator or the hard-threshold estimator.
#'
#'
#' Given \eqn{Y = (Y1, Y2, \dots, Yn)}, such that:
#' \deqn{Yi ~ N(\theta i, \sigma2), 1 \le i \le n.}
#' \eqn{N(\theta, \sigma2)} denotes the normal distribution with
#' means \eqn{\theta} and variance \eqn{\sigma2}. The means are to
#' be estimated. Variance is either known or to be estimated.
#'
#'
#' Means are either estimated using the soft-threhsold (ST) or the
#' hard-threshold (HT) estimator. ITSES minimize expected square error loss
#' under repeated sampling (risk).
#'
#' Soft-threshold: \eqn{\hat{\theta}(x)= sign(x)(abs(x)-\lambda)(abs(x)>\lambda)}.
#'
#'
#' Hard-threshold: \eqn{\hat{\theta}(x)= x(abs(x)>\lambda)}.
#'
#'
#' @param y A `n` vector of observations.
#' @param method Estimator which will be used. Can be either
#' "ST" (Soft-Threshold estimator) or "HT" (Hard-Threshold estimator).
#' Defaults to "ST".
#' @param m The number of iterative steps. Default is `5`.
#' @param init.lambda Initial threshold to start iterations at. Can be:
#' `numeric`; "median"; "visu" (universal threshold); and "half-visu"
#' (universal threshold divded by 2). Defaults to "median".
#' @param max.length The maximum number of observations `n`. If`n` is above
#' internal calculations will downsample `y` to be of `max.length` size.
#' Defaults to `5e6`.
#' @param minimizationmethod   The risk minimization method which will be used.
#' Can be either "numeric" or "sampling". Default is "numeric".
#' @param debug `logical`. Specify wether or not to print debug code from
#' iterations.
#' @param sd Standard deviation of `y` i.e. noise levels. If `NA` will estiamt
#' e using specified MAD estimator. Default is `NA`.
#' @param sparsity Sparsity of `y` if known. `0 < sparsity < 1`, closer to 1
#' means highly sparse. Is used in noise estimation if sparse noise estimation
#' is used. Default is `NA`.
#' @param sparse.mad `logical`. If `TRUE` use sparsity in MAD noise estimation.
#' Default is `TRUE`.
#' @param remove.zero `logical`. Remove observation that are zero from noise
#' estimatiom. Default is `TRUE`.
#' @param tol A numeric that determines the sensitivity of threshold selection.
#' Default is `1e-8`, iteration will stop if change is less than set `tol`.
#' @param h Parameter for deciding eligibility of sparsity measure.
#' Default is `0.4`, only used with SparseMAD estimator.
#' @param ... parameters passed to risk minmization methods.
#' @param b number of samples to take. Default is `10`.
#' @param k number of thresholds to initially evaluate. Default is `10`.
#' @param max_num_iters maximum number of iterations with Newton's method or
#' spline interpolations (depending on method used). Default is `10`.
#' @return List object with iteration results.
#'
#' @export
itses <- function(y,
                  method = "ST",
                  m = 5,
                  init.lambda = "median",
                  noisetype = "gaussian",
                  max.length = 5e6,
                  minimizationmethod   = "numeric",
                  debug = FALSE,
                  sd = NA,
                  sparsity = NA,
                  sparse.mad = TRUE,
                  remove.zero = TRUE,
                  tol = 1e-8,
                  h = 0.4,
                  max.threshold = NULL,
                  min.threshold = NULL,
                  ...
) {
  # Check input
  n <- length(y)
  if(n < 1) {
    stop("Not enough data points")
  }
    if(debug) print(paste("Starting ITSES. Number of points in data: ", n))

  normalizing <- TRUE
  if(noisetype == "gaussian"){
    if(debug) print("Using gaussian noise")
    # Estimate noise if not available.
    if((!is.numeric(sd)) | (sd <= 0)) {
      if(debug) print("Starting sd estimation.")
      y.temp <- y # Temp. variable

      # If wanted, remove 0s from estimation.
      if(remove.zero) {
        if(debug) print("Removing zeros frem sd estimation.")
       y.temp<- y.temp[y.temp != 0]
      }
      if(sparse.mad) {
        if(debug) print("Using SparseMAD estimator.")
        # Estimate noise by building in sparsity assumption.
        if(is.na(sparsity)) {
          # Use 'SparseMAD' to estimate noise.
          if(debug) print("Estimating noise.")
          sd <- sparse.mad.estimator(y.temp, h = h, debug = debug)
        }else{
          if(debug) print("Using given sparsity level.")
          # If sparsity is available find sparsity threshold.
          t <- get.threshold.for.sparsity(y, sparsity)
          # Use suggested 'SparseMAD' with given threshold.
          sd <- thresholded.mad.estimator(y.temp, t)
        }
      }else{
        if(debug) print("Using MAD estimator to estimate noise.")
        sd <- mad.estimator(y.temp)
      }
    }
    if(debug) print(paste0("Sd:", sd))
    # Normalize data for iterative threshold selection.
    z <- y/sd
  }else {
    if(debug) print("Using custom noise.")
    normalizing <- FALSE
    z <- y
  }

  # Set minimum and maximimum bounds for threshold.
  if(is.null(min.threshold))  min.threshold <- 0
  if(is.null(max.threshold)) max.threshold <- get.max.threshold(method, z)

  if(debug) print(paste("Threshold bounds are set to:", min.threshold, max.threshold))
  # Get starting threshold.
  lambda <- get.init.lambda(z, init.lambda = init.lambda)
  lambda <- get.thresholded.lambda(lambda, min.threshold, max.threshold)

  # Set start parameters.
  if(minimizationmethod == "sampling") lambda.grid <- NULL

  # Initiate objects to store results in.
  sds <- sd
  iteration_result <- list()
  lambdas <- lambda # List storing all threholds.

  # Initaitate variable to store last iteration step
  last.lambda <- -Inf

  # Start iterative estimation.
    for(i in 0:(m-1)) {
      if(debug) print(paste("Iteration:", i, "Threshold:", lambda))

      # Down-sample if data exceeds max length.
      z.sample <- z
      if(n  > max.length) {
        if(debug) print("Downsampling")
        z.sample <- sample(z, max.length)
      }

      # Select minimization method.
      if(minimizationmethod   == "sampling") {
        # Find threshold by using sampling to estimate risk.
        iter.results <- iter.sampling(y = z.sample,
                                    init.lambda = lambda,
                                      noisetype = noisetype,
                                    ...,
                                    lambda.grid  = lambda.grid,
                                    debug = debug,
                                    method = method,
                                    tol = tol,
                                      min.threshold = min.threshold,
                                      max.threshold = max.threshold
                                    )
      lambda.grid <- iter.results$lambdas
      }else if(minimizationmethod   == "numeric") {
        # Find threshold by directly minimizing risk.
        if(noisetype != "gaussian") stop("Non supported noisetype")
        if(!is.character(method)) stop("Invalid method")
        iter.results <- iter.newton(y = z.sample,
                                  init.lambda = lambda,
                                  ...,
                                  method = method,
                                  tol = tol,
                                  debug = debug,
                                             min.threshold = min.threshold,
                                             max.threshold = max.threshold )
      }else{
        stop("Not supported minimization method, only 'sampling' and 'numeric' currently supported.")
      }
      # Store results.
      if(normalizing){
        iteration_result[[i+1]] <- cbind(sd*iter.results$lambdas,sd^2*iter.results$risks)
      }else{
        iteration_result[[i+1]] <- cbind(iter.results$lambdas,iter.results$risks)
      }
      lambda <- iter.results$lambda
      lambdas <- c(lambdas, lambda)

      # Check convergence.
      # First equality is to chek if Inf == Inf as Inf - Inf is NA ...
       if(last.lambda == lambda|abs(last.lambda -  lambda)<tol) {
           if(debug) print("Breaking iterations, same threshold as previous iteration.")
           break()
       }
      last.lambda <- lambda
    }


  # Get final mean estimates
  if(normalizing) {
    lambda <- sd*lambda
    lambdas <- sd*lambdas
  }

  if(is.character(method)){
    if(method == "ST") {
    theta <- soft.threshold.estimator(y, lambda)
    }else if(method == "HT") {
      theta <- hard.threshold.estimator(y, lambda)
    }else{
      stop("Invalid method.")
    }
  } else{
    theta <- method(y, lambda)
  }

  # Storing parameters and results.
  result <- list(theta = theta,
                 lambda = lambda,
                 sd = sd,
                 sds = sds,
                 lambdas = lambdas)
  result[["iteration_result"]] <- iteration_result
  result[["y"]] <- y
  result[["method"]] <- method
  result[["m"]] <- m
  result[["minimizationmethod"]] <- minimizationmethod

  # Return results.
  if(debug) print(paste0("Final threshold is: ",  result$lambda))
  result
}

#result <- itses((1+rnorm(10))*rnorm(10), sd = 0.05, method = "ST", max.threshold = Inf,
#                minimizationmethod = "sampling", noisetype = "speckle")
#result$lambda