
#' Median abosulte deviaton (MAD) noise estimator.
#' @param y numeric of `n` length to be used in noise estimation.
#' @return Estimated standard deviation
#' @noRd
mad.estimator <- function(y) {
  median(abs(y)/qnorm(0.75))
}

#' Filtered MAD noise estimatation.
#' @param y numeric of `n` length to be used in noise estimation.
#' @param lambda Exclude points with aboslute values greater than `lambda` from estimation.
#' @return Estimated standard deviation
#' @noRd
thresholded.mad.estimator <- function(y, lambda) {
  y <- y[abs(y)<=lambda]
  mad.estimator(y)
}

#' Given sparsity return threshold that gives the same sparsity.
#' @param y a numeric of `n` length.
#' @param s a numeric measure on the sparsity of `y`. `0<s<1`.
#' @return s-th quantile of abs(y).
#' @noRd
get.threshold.for.sparsity <- function(y, s) {
  n <- length(y)
  n.s <- floor(n*s)
  sort(abs(y))[n.s]
}

#' Estiamte sparsity
#'
#' Estiamtes sparsity of a given vector.
#'
#' @noRd
sparsity.estimator <- function(y) {
  (length(y)-sum(abs(y))^2/sum(y^2))/length(y)
}

#' SparseMAD noise estimator.
#'
#' Estiamtes noise using a filtered MAD estimator.
#'
#' @param y Vector of length `n` estiamte noise from.
#' @param h Numeric to decide eligibility of sparsity estimate. Defaults to `0.4`.
#' @param debug `logical`, print logs if `TRUE`. Defaults to `FALSE`.
#' @export
sparse.mad.estimator <- function(y, h=0.4, debug = FALSE) {
  s <- sparsity.estimator(y)
  if(debug) print(paste0("Estiamted sparsity: ", s, ". h at: ", h))
  t <- Inf
  if(s>h) {
    if(debug) print("Thresholding MAD in noise estimation.")
    t <- get.threshold.for.sparsity(y, s)
  }else{
    if(debug) print("Using regular MAD in noise estimation.")
  }
  sd <- thresholded.mad.estimator(y, t)
  sd
}