This repository accompanies the article: Drug homogeneity index in mass spectrometry imaging and provides a reference implementation of our pipeline for drug homogeneity assessment in a given MSI data at different user define 
parameters. At the moment our pipeline accept single data format (analyze 7.5). 

The pipeline is developed at FEM, San michele, Italy.

For those wishing to use it for more general use, please consider instructions in General usage. 

## Requirments

R version above 3.0
R packages : MALDIquant Foreign, radiomics, msProcess

### Inputs

To process a dataset two things are needed:
1. MSI data location in user computer
2. m/zs for : drug, internal standard (optional), mask image (optional).
In the absence of m/z for mask image, drug image will be used to create tissue mask image. 


## Contact

email: Mridula Prasad 
