rm(list = ls())
setwd("~/GitHub/SMNM-Implementation")

library(imager)
library(waveslim)
library(itses)

sd <- 0.05

image_mse <- function(original, estimated, dims = 3){
  v <- NULL
  for(l in 1:dims) {
     vt <- (original[,,,l]-estimated[,,,l])^2
     vt<- vt[1:length(vt)]
    v <- c(v, vt)
    }
  mean(v)
}

image_true_mse <- function(lambdas, original, noisy, estimator, dims = 3){
  image.layers <- denoised.image.list <- lapply(1:dims, function(l){
    img.layer <- noisy[,,,l]
    details <- dwt.2d(img.layer, wf = wf, J = J, boundary = boundary)
    detail.names <- names(details)

     thresholded.details <- lapply(detail.names, function(name){
      print(name)
      res <- list()
      if(grepl("LL", name, fixed = T)){
        res[[name]] <- details[[name]]
      }else{
        detail.vector <- details[[name]][1:length(details[[name]])]
        treated <- estimator(detail.vector, lambdas[[1]][[name]])
        reconstructed.detail <-  matrix(treated, nrow = nrow(details[[name]]))
        res[[name]] <- reconstructed.detail
      }
      return(res)
    })

    thresholded.details <- Reduce(c, thresholded.details)

    for(name in detail.names){
      details[[name]] <- thresholded.details[[name]]
    }
    reconstructed.image.layer <- idwt.2d(details)
    reconstructed.image.layer
}

clip_image_layer <- function (y, min_bound = 0, max_bound = 1){
  clipped_layer <- apply(y, 1, function(x) sapply(x, function(z) min(max(z, min_bound), max_bound)))
  clipped_layer
}

add_gaussian_noise_to_image_layer <- function(y, sd = 1){
  noisy <- apply(y, 1, function(x) rnorm(length(x), mean = x, sd = sd))
  clipped_noisy <- clip_image_layer(noisy)
  clipped_noisy
}

add_speckle_noise_to_image_layer <- function(y, sd = 1){
  noisy <- apply(y, 1, function(x) (1+rnorm(length(x), mean = 0, sd = sd))*x)
  clipped_noisy <- clip_image_layer(noisy)
  clipped_noisy
}

add_poisson_noise_to_image_layer <- function(y){
  noisy <- apply(y, 1, function(x) rpois(length(x), x*255)/255)
  clipped_noisy <- clip_image_layer(noisy)
  clipped_noisy
}


add_mixed_gaussian_poisson_noise_to_image_layer <- function(y, sd = 0.05){
  noisy <- apply(y, 1, function(x) rpois(length(x), x*255)/255 +  rnorm(length(x), mean = 0, sd = sd))
  clipped_noisy <- clip_image_layer(noisy)
  clipped_noisy
}



add_noise_to_image <- function(image, layer_noise, ..., dims = 3, common_layer_noise = F){
  noisy <- image
  if(common_layer_noise){
    stop("Not supported")
  }else{
    for(l in 1:dims) {
      print(l)
      noisy[,,,l]  <- layer_noise(image[,,,l], ...)
    }

  }

  noisy
}

image <- load.image("../SMNM-Implementation/output/figures/peppers256.png")
noisy <- add_noise_to_image(image, add_gaussian_noise_to_image_layer, sd =sd)

par(mfrow = c(1,2))
plot(image)
plot(noisy)



itses_image <- function(image, method = "ST",
                        dims = 3,
                        wf = "haar",
                        boundary = "periodic",
                        sd = NA,
                        ...,
                        J = 6){

  image.result <- image
  denoised.image.list <- lapply(1:dims, function(l){
    print(l)
    img.layer <- image.result[,,,l]
    details <- dwt.2d(img.layer, wf = wf, J = J, boundary = boundary)
    detail.names <- names(details)

    if(is.na(sd)){
    d <- sapply(detail.names, function(name){
      if(grepl("LL", name, fixed = T)){
        return(NULL)
      }else{
        detail.vector <- details[[name]][1:length(details[[name]])]
        return(detail.vector)
      }
    })
    d <- Reduce(c, d)

    sd <- itses:::mad.estimator(d)
    }

    thresholded.details <- lapply(detail.names, function(name){
      print(name)
      res <- list()
      if(grepl("LL", name, fixed = T)){
        res[[name]] <- details[[name]]
      }else{
        detail.vector <- details[[name]][1:length(details[[name]])]
        treated <- itses(detail.vector, method = method, sd = sd, ...)$theta
        reconstructed.detail <-  matrix(treated, nrow = nrow(details[[name]]))
        res[[name]] <- reconstructed.detail
      }
      return(res)
    })

    thresholded.details <- Reduce(c, thresholded.details)

    for(name in detail.names){
      details[[name]] <- thresholded.details[[name]]
    }
    reconstructed.image.layer <- idwt.2d(details)
    reconstructed.image.layer
  })
  for(l in 1:dims)  image.result[,,,l] <- denoised.image.list[[l]]
  image.result

}



itses_image_custom_noise <- function(image,
                                     noise_fun =
                                       function(image_layer)  add_gaussian_noise_to_image_layer(image_layer, sd =sd),
                                     method = "ST",
                                     b = 10,
                                     dims = 3,
                                     m = 10,
                                     wf = "haar",
                                     boundary = "periodic",
                                     min.threshold = NULL,
                                     max.threshold = NULL,
                                      init.lambda = "median",
                                     tol = 1e-8,
                                     debug = TRUE,
                                     ...,
                                     J = 6){
  if(is.character(method)){
      if(method == "ST") {
        method <- function(y, lambda) itses:::soft.threshold.estimator(y, lambda)
      }else if(method == "HT") {
        method <- function(y, lambda) itses:::hard.threshold.estimator(y, lambda)
      }else{
        stop("Invalid method.")
      }
  }
  image.result <- image
  denoised.image.list <- lapply(1:dims, function(l){
    print(l)
    img.layer <- image.result[,,,l]
    details <- dwt.2d(img.layer, wf = wf, J = J, boundary = boundary)
    detail.names <- names(details)

    thresholds <- lapply(detail.names, function(name){
      res <- list()
      if(is.null(min.threshold)){
        detail.min.threshold <- 0
      }
      else {
        detail.min.threshold <- min.threshold
      }
      if(is.null(max.threshold)) {
        detail.max.threshold <- Inf
      }else{
        detail.max.threshold <- max.threshold
      }
      lambda <- itses:::get.init.lambda(details[[name]], init.lambda = init.lambda)
      lambda <-  itses:::get.thresholded.lambda(lambda, min.threshold, max.threshold)

      res[[name]] <- c(lambda, detail.min.threshold, detail.max.threshold)
      res
    })
    thresholds <- Reduce(c, thresholds)

    lambda.grids <- NULL
    lambda.grids <- sapply(detail.names, function (name) lambda.grids[[name]] <- NULL)

    # Initiate objects to store results in.
    iteration_results <- list()

    # Initaitate variable to store last iteration step
    last.thresholds <- NULL
    last.thresholds <- sapply(detail.names, function (name) last.thresholds[[name]] <- -Inf)

    hist.thresholds <- NULL
    hist.thresholds <- sapply(detail.names, function (name) hist.thresholds[[name]] <- list(thresholds[[name]][1]))

    # Start iterative estimation.
    for(i in 0:(m-1)) {
      thresholded.details <- lapply(detail.names, function(name){
            res <- list()
            if(grepl("LL", name, fixed = T)){
              res[[name]] <- details[[name]]
            }else{
              treated <- method(details[[name]][1:length(details[[name]])], thresholds[[name]][1])
              reconstructed.detail <-  matrix(treated, nrow = nrow(details[[name]]))
              res[[name]] <- reconstructed.detail
            }
            res
        })
      thresholded.details <- Reduce(c, thresholded.details)

      thresholded.image.details <- details
      for(name in detail.names){
          thresholded.image.details[[name]] <- thresholded.details[[name]]
        }
      thresholded.image.layer <- idwt.2d(thresholded.image.details)
      thresholded.image.layer <- t(clip_image_layer(thresholded.image.layer))
      image.layer.samples <- lapply(1:b, function(b) noise_fun(image_layer = thresholded.image.layer))
      detail.samples.temp <- lapply(image.layer.samples, function(image.layer){
        dwt.2d(image.layer, wf = wf, J = J, boundary = boundary)
      })

      detail.samples <- NULL
      detail.samples <- sapply(detail.names, function(name){
        detail.samples[[name]] <- sapply(1:b, function(j) {
          ds <- detail.samples.temp[[j]][[name]]
          ds <- ds[1:length(ds)]
      })
      })
      if(debug) print(paste("Iteration:", i))
      iter.resultss <- lapply(detail.names, function(name){
        # Find threshold by using sampling to estimate risk.
        results<- itses:::iter.sampling(y = details[[name]][1:length(details[[name]])],
                                        init.lambda = thresholds[[name]][[1]],
                                        noisetype = noisetype,
                                        samples = detail.samples[[name]],
                                        ...,
                                        lambda.grid  = lambda.grids[[name]],
                                        debug = debug,
                                        b = b,
                                        method = method,
                                        tol = tol,
                                        min.threshold = thresholds[[name]][[2]],
                                        max.threshold =  thresholds[[name]][[3]]
                                    )

        lambda.grids[[name]] <<- results$lambdas
        thresholds[[name]][1] <<- results$lambda
        hist.thresholds[[name]] <<- c(hist.thresholds[[name]], results$lambda)

        results
      })
    }
     thresholded.details <- lapply(detail.names, function(name){
            res <- list()
            if(grepl("LL", name, fixed = T)){
              res[[name]] <- details[[name]]
            }else{
              treated <- method(details[[name]][1:length(details[[name]])], thresholds[[name]][1])
              reconstructed.detail <-  matrix(treated, nrow = nrow(details[[name]]))
              res[[name]] <- reconstructed.detail
            }
            res
        })
      thresholded.details <- Reduce(c, thresholded.details)

      thresholded.image.details <- details
      for(name in detail.names){
          thresholded.image.details[[name]] <- thresholded.details[[name]]
        }
    thresholded.image.layer <- idwt.2d(thresholded.image.details)
    thresholded.image.layer <- t(clip_image_layer(thresholded.image.layer))


    list(thresholded.image.layer = thresholded.image.layer, thresholds = thresholds, hist.thresholds = hist.thresholds)

  })

  for(l in 1:dims)  image.result[,,,l] <- denoised.image.list[[l]]$thresholded.image.layer
  image.result
}

pdf()
sd <- 0.1
par(mfrow = c(2,2))
plot(image, main = "Original")

ms <- round(image_mse(image, noisy), digits = 5)
plot(noisy, main = paste("Noisy, mse: ", ms))


# these looks very close, so seems to work similar
denoised <- itses_image(noisy, sd = sd)
ms <- round(image_mse(image, denoised), digits = 5)
plot(denoised, main = paste("Non-custom, mse: ", ms))


denoised_custom <- itses_image_custom_noise(noisy, b = 1)
ms <- round(image_mse(image, denoised_custom), digits = 5)
plot(denoised_custom, main = paste("custom, mse: ", ms))

## Speckle test
sd = 0.5
noisy_speckle <- add_noise_to_image(image, add_speckle_noise_to_image_layer, sd =sd)

plot(image, main = "Original")

ms <- round(image_mse(image, noisy_speckle), digits = 5)
plot(noisy_speckle, main = paste("Speckle noisy, mse: ", ms))

denoised <- itses_image(noisy_speckle, sd = NA, sparse.mad = F)
ms <- round(image_mse(image, denoised), digits = 5)
plot(denoised, main = paste("Non-custom, mse: ", ms))


denoised_custom <- itses_image_custom_noise(noisy_speckle, b = 1,
                                            noise_fun =  function(image_layer)  add_speckle_noise_to_image_layer(image_layer, sd =sd))
ms <- round(image_mse(image, denoised_custom), digits = 5)
plot(denoised_custom, main = paste("custom, mse: ", ms))

### Poisson

noisy_poisson <- add_noise_to_image(image, add_poisson_noise_to_image_layer)
plot(image, main = "Original")
ms <- round(image_mse(image, noisy_poisson), digits = 5)
plot(noisy_poisson, main = paste("Noisy poisson, mse: ", ms))
denoised <- itses_image(noisy_poisson, sd = NA, sparse.mad = F)
ms <- round(image_mse(image, denoised), digits = 5)
plot(denoised, main = paste("Non-custom, mse: ", ms))
denoised_custom <- itses_image_custom_noise(noisy_poisson, b = 50,
                                            noise_fun =  function(image_layer)  add_poisson_noise_to_image_layer(image_layer)
)
ms <- round(image_mse(image, denoised_custom), digits = 5)
plot(denoised_custom, main = paste("custom, mse: ", ms))


## Mixed poisson

sd = 0.02
noisy_poisson <- add_noise_to_image(image, add_mixed_gaussian_poisson_noise_to_image_layer, sd = sd)
plot(image, main = "Original")
ms <- round(image_mse(image, noisy_poisson), digits = 5)
plot(noisy_poisson, main = paste("Noisy poisson, mse: ", ms))

denoised <- itses_image(noisy_poisson, sd = NA, sparse.mad = F)
ms <- round(image_mse(image, denoised), digits = 5)
plot(denoised, main = paste("Non-custom, mse: ", ms))


denoised_custom <- itses_image_custom_noise(noisy_poisson, b = 1,
                                            noise_fun =  function(image_layer)  add_mixed_gaussian_poisson_noise_to_image_layer(image_layer, sd = sd)
)
ms <- round(image_mse(image, denoised_custom), digits = 5)
plot(denoised_custom, main = paste("custom, mse: ", ms))