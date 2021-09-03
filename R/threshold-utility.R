# Contains utility functions specific to threshold selection.
#' Get universal thresholds
#' Assuming y of unit variance.
#' @noRd
get.visu.threshold <- function(y) {
  n <- length(y)
  lambda <- sqrt(2*log(n))
  lambda
}
#' Refine threshold grid
#' @param lambda.grid : initial grid, is assumed to be sorted
#' @param k : number of points to have in grid, count include min and max.
#' @noRd
refine.lambda_grid <- function(lambda.grid, k = NA) {
  # Ensure arguments passed is correct.
  if(!is.numeric(k)) {
    # Do not do anything if k is non numeric, e.g. NA.
    return(lambda.grid)
  }else if(k < 2) {
    stop("Too few points selected for grid. Only supporting k > 2.")
  }
  # Reduce length of grid, select grid values that are evenly apart by what was given.
  if(length(lambda.grid) > k) {
    # Select k points evenly distributed, ensuring that min and max is included.
    qs <- seq(0, 1, length.out = k)
    lambda.grid <- as.vector(quantile(lambda.grid, probs = qs))
  }else if(k - length(lambda.grid) > 0) {
    # Extend the grid, there is a low chance of overlap, this could make number not equal to k..
    lambda.grid.extension <- seq(from = min(lambda.grid),
                                 to = max(lambda.grid[is.finite(lambda.grid)]),
                                 length.out = k - length(lambda.grid) + 2)
    lambda.grid <- c(lambda.grid, lambda.grid.extension)
  }
  # Return grid and ensure sorted and unique.
  lambda.grid <- sort(lambda.grid)
  lambda.grid <- unique(lambda.grid)
  lambda.grid
}
#' Get grid of threhsolds
#' @noRd
get.lambda_grid <- function(y,
                            min.threshold = 0,
                            max.threshold = Inf,
                            k = NA) {
  # Only want unique positve grid values
  y <- unique(abs(y))
  # Define lambda grid, at observations, add max threshold and 0
  lambda.grid <- c(0, y[y!=0], max.threshold)
  # Ensure increasing
  lambda.grid <- sort(lambda.grid)
  # Apply upper limit
  lambda.grid <- lambda.grid[lambda.grid <= max.threshold]
  # Apply lower limit to grid
  lambda.grid <- lambda.grid[min.threshold <=lambda.grid]
  # Ensure all points are unique
  lambda.grid <- unique(lambda.grid)
  # Refine the grid, and make sure it has about k points
  lambda.grid <- refine.lambda_grid(lambda.grid , k = k)
  # Return grid
  lambda.grid
}
#' Select initial threshold, can either be a numeric or a method on the observations..
#' Supports median, visu (universal threshold) and half-visu (universal threshold divded by 2).
#' Unit variance is assumed with universal threhsolds.
#' @noRd
get.init.lambda <- function(y, init.lambda = "median") {
  lambda <- 0
  if(is.numeric(init.lambda)) {
    # If numeric just return the value
    lambda <- init.lambda
  }else{
    visu <- get.visu.threshold(y) # Universal threshold
    # Select threshold by method.
    if(init.lambda == "median") {
      lambda <- median(abs(y))
    }else if(init.lambda == "visu") {
      lambda <- visu
    } else if(init.lambda == "half-visu") {
      lambda <- visu/2
    }else{
      stop("Not supported threshold initialization.")
    }
  }
  lambda
}
#' Gives maximum bound to threshold by method.
#' Assumes unit variance
#' @noRd
get.max.threshold <- function(method, z) {
  if(method == "ST") {
    max.threshold <- get.visu.threshold(z)
  }else if(method == "HT") {
    max.threshold <- Inf #max(abs(z))
  }else{
    stop("Invalid method. Only soft threshold (ST) and hard-threshold (HT) supported.")
  }
  max.threshold
}
#' Projects threshold onto given bounds.
#' @noRd
get.thresholded.lambda <- function(lambda, min.threshold, max.threshold) {
  sapply(lambda, function(lambda) min(max(lambda, min.threshold), max.threshold))
}
