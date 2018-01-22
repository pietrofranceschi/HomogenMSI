# HomogenMSI
This package provides a reference implementation of our algorithm for drug homogeneity assessment in ion images extracted from MSI datasets. 

# Installation

* Depends on

```
spatstat, reshape2
```
* Install devtools: in the R/R Studio shell type
```{r}
install.packages("devtools")
```

* Install homogenMSI: In the shell type
```{r}
library(devtools)
install_github("pietrofranceschi/HomogenMSI", dependencies = TRUE) 
```
To additionally build the vignette the previous command should include `build_vignettes = TRUE`

```{r}
library(devtools)
install_github("pietrofranceschi/HomogenMSI", dependencies = TRUE, build_vignettes = TRUE) 
```

* Load the package: In shell type
```{r}
library(HomogenMSI)
```

# Description
For the deails see the package vignette




## Contact

email: Mridula Prasad (<mridula.prasad@fmach.it>)
