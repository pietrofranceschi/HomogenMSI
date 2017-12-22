suppressMessages(library("radiomics"))

## function for 3x3 median filtering 

medianFilterR <- function(sampleMat)                              
{
  B = matrix(0, nrow= dim(sampleMat)[1],ncol=dim(sampleMat)[2])
  modifyA = matrix(0,nrow=(dim(sampleMat)[1]+2),ncol=(dim(sampleMat)[2]+2))
  for(x in 1:dim(sampleMat)[1])
  {
    for(y in 1:dim(sampleMat)[2])
    {
      modifyA[x+1,y+1] = sampleMat[x,y]
    }
  }
  for(i in 1:(dim(modifyA)[1]-2))
  {
    for(j in 1:(dim(modifyA)[2]-2))
    {
      window = rep(0,9)
      inc = 1
      for(x in 1:3)
      {
        for(y in 1:3)
        {
          window[inc] = modifyA[i+x-1,j+y-1]
          inc = inc +1
        }
      }
      med = sort(window)
      B[i,j] = med[5]
    }
  }
  return(B)
}


##### calculate DHI value from given MSI data (unnormalized with tumor area)


CalculateDHI <- function(drug_csvfilepath,mask_csvfilepath,QuantLevel=32,Nu=1)
{
Img_drug = as.matrix(read.csv(drug_csvfilepath,header = T));mask =  as.matrix(read.csv(mask_csvfilepath,header = T))
Img_drug = Img_drug[,-1];mask = mask[,-1]

m = QuantLevel/max(Img_drug)
Img_quantized = Img_drug*m
Img_quantized = ceiling(Img_quantized)
Img_quantized = medianFilterR(Img_quantized)
Img_quantized = mask * Img_quantized

szm = glszm(Img_quantized)
szm = szm@.Data
szm = szm[-1,]   ### Removing sz values for the background of image
szm = szm[,-as.numeric(which(colSums(szm) ==0 ))]
szv = as.numeric(colnames(szm))
id = which(szv >= Nu)
DrugHomo = 0
for(j in 1:length(id))
{     
 DrugHomo = DrugHomo + (szv[id[j]]) * sum(szm[,id[j]])
}
DrugHomo = DrugHomo/sum(szm[,id])
TumorArea = table(mask)[[2]]
DrugHomo = DrugHomo/TumorArea

return(DrugHomo)
}
