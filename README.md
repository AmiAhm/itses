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

Alternatively, download package material or clone the repository, and run the following R-script. Ensure to have the main folder "itses/" in the working directory.

```{r}
library(devtools)
devtools::install("itses")
```


## Example

Running the main method: 

```{r}
library(itses)

# To get data to work with
y <- rnorm(10) 

# To run at default settings (noise level estiamted + soft-threshold + numerical)
itses.result <- itses::itses(y) 

#  Print result
print(paste("Optimal threshold is:", itses.result$lambda))

```

Refer to the documentation and vignettes for in-depth explanation of parameters and alternatives. E.g. by:
```{r}
?itses::itses
```


## Details

The package uses no external dependencies and was built on R version 4.0.3.


## Changelog

### 0.0.0.9001 06-09-2021
Fix rogue linebreak in utility causing faulty HT risk


## About

This package was made as in part of the research project component of the MSc in Statistics of Imperial College London.

Thesis: Sparse Estimation of Many Normal Means
Implementation used in thesis can be found at: https://github.com/AmiAhm/SMNM-Implementation

By Amir Ahmed, supervised by Alastair Young (September 2021)



