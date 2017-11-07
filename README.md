# homogenMSI

This repository accompanies the article: Drug homogeneity index in mass spectrometry imaging and provides a reference implementation of our pipeline for drug homogeneity assessment in a given MSI data at different user-defined parameters. At the moment, our pipeline accepts data from single analytical platform, i.e. analyze 7.5.  

There are two main functions in our pipeline: 

a)	CalculateDHI

b)	DHI

The CalculateDHI function calculates the DHI value from given  MSI data following our pre-processing steps discussed in the paper. If the user simply wants to calculate DHI on their pre-processed data can directly drive it using DHI function. 

Note: while using DHI function, knowledge about a number of pixels in tumor area required and for efficient calculation assigns background gray level equal to the zero. 

## Requirments
R version above 3.0

R packages: MALDIquantForeign, radiomics, msProcess

### Inputs

DHI formula requires following input parameters:

1.	Filename: MSI data folder location in user computer
2.	Binned: m/z bins which can create using CreateBin function
3.	mzs: list of desired mz values: 1 drug, 2 mask, 3 internal standard (optional).In the absence of m/z for mask image, drug image will be used to create tissue mask image.
4.	QuntLevel: Intensity quantization level to be used, ex: 8, 16, 32
5.	Bkg : ‘T’ In our study, we had tissue slice placed on the glass side, hence we had removed the gray level associated with it before DHI calculation. If this is not a case with the user, use ‘F’. Default value is ‘T’

### Example

binned = CreateBin(folderpath)

mzs = c(284.2, 281.2, 289.2)

DHI = CalculateDHI(folderpath, binned, mzs, QuntLevel=8)

### Outputs

![alt text](https://github.com/pietrofranceschi/homogenPy/blob/master/A2780_AVA%20%23326%201.png)
![alt text](https://github.com/pietrofranceschi/homogenPy/blob/master/A2780_PTX%20%23352%202.png)

## Contact

email: Mridula Prasad 
