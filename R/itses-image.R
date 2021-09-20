rm(list = ls())
setwd("~/GitHub/SMNM-Implementation")

library(imager)
library(waveslim)
library(itses)

sd <- 0.2

clip_image_layer <- function (y, min_bound = 0, max_bound = 1){
  clipped_layer <- apply(y, 1, function(x) sapply(x, function(z) min(max(z, min_bound), max_bound)))
  clipped_layer
}

add_gaussian_noise_to_image_layer <- function(y, sd = 1){
  noisy <- apply(y, 1, function(x) rnorm(length(x), mean = x, sd = sd))
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

image <- load.image("output/figures/stjerten256.png")
noisy <- add_noise_to_image(image, add_gaussian_noise_to_image_layer, sd =sd)

par(mfrow = c(1,2))
plot(image)
plot(noisy)



itses_image <- function(image, method = "ST",
                        dims = 3,
                        wf = "haar",
                        boundary = "periodic",
                        J = 4){

  image.result <- image
  denoised.image.list <- lapply(1:dims, function(l){
    print(l)
    img.layer <- image.result[,,,l]
    details <- dwt.2d(img.layer, wf = wf, J = J, boundary = boundary)
    detail.names <- names(details)
    thresholded.details <- lapply(detail.names, function(name){
      print(name)
      res <- list()
      if("LL" %in% name){
        res[[name]] <- details[[name]]
      }else{
        detail.vector <- details[[name]][1:length(details[[name]])]
        treated <- itses(detail.vector, method = method, sd = sd)$theta
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

denoised <- itses_image(noisy)
plot(denoised)
l <- 1
img.layer <- image[,,l]
details <- dwt.2d(img.layer, "haar", J = 4, boundary = "periodic")
for(name in names(details)){
  if("LL" %in% name){
    next()
  }
  detail <- details[[name]]

}

idwt.2d(details)
img.layer