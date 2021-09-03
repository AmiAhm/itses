# Iterative Sparse Estimation (ITSES) of Many Normal Means

The  main purpose of this package is to perform sparse estimation of many normal means.  Estimation is done using either the hard- or the soft-threhold estimator. 
An iterative approach is used to select a threshold.

## Set up 

Set up is easily done through the `devtools` package. 

Run the following R-script to install `itses`:

```{r}
library(devtools)
devtools::install_github("AmiAhm/itses")
```

Alternatively, download package material or clone the repository, and run the following R-script. Ensure to have the main folder "itses/" in the working directory:

```{r}
library(devtools)
devtools::install("itses")
```


## Example

Running the main method. 

```{r}
library(itses)

# To get data to work with
y <- rnorm(10) 

# To run at default settings (noise level estiamted + soft-threshold + numeric)
itses.result <- itses::itses(y) 

#  Print result
print(paste("Optimal threshold is:", itses.result$lambda))

```

Refer to the documentation and vignettes for in-depth explanation of parameters and alternatives. E.g. by:
```{r}
?itses::itses
```


## Details

The package uses no external dependencies. It was built on R version 4.0.3 (2020-10-10).


## About

This package was made in part of a research project for the MSc in Statistics of Imperial College London.

By Amir Ahmed, supervised by Alastair Young (September 2021)



