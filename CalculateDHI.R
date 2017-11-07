
suppressMessages(library("MALDIquant"))
suppressMessages(library("MALDIquantForeign"))
suppressMessages(library("msProcess"))
suppressMessages(library("splus2R"))
suppressMessages(library("data.table"))
suppressMessages(library("radiomics"))

### local minima search in density plot 

minimums <- function(x) which(x - shift(x, 1) < 0  & x - shift(x, 1, type='lead') < 0)   


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


## Create initial bins for peak selection 

CreateBin <- function(filepath,binsize=0.5)
{
analyfile = importAnalyze(filepath)
binlength = seq(range(analyfile[[1]]@mass)[1]-0.5, range(analyfile[[1]]@mass)[2]+0.5, by = binsize) # for fixed bin size
cuts = cut(analyfile[[1]]@mass,binlength)
duration.freq = table(cuts)
duration.freq = cbind(duration.freq)
cutdur = row.names(duration.freq)
binsplit = function(x){as.numeric(unlist(strsplit(gsub("\\(|\\]", "", x),','))[1])}
bin_even = sapply(cutdur,binsplit)
binsplit = function(x){as.numeric(unlist(strsplit(gsub("\\(|\\]", "", x),','))[2])}
bin_odd = sapply(cutdur,binsplit)
binned = cbind(bin_even,bin_odd)
binname = apply(binned,1,mean)
return(binned)
}

##### calculate DHI value from given MSI data (unnormalized with tumor area)

DHI <- function(img,Nu=1,Bkg='T',TumorArea)
{
szm = glszm(img)
szm = szm@.Data
if(Bkg=='T')
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
if(!missing(TumorArea))
   DrugHomo = DrugHomo/TumorArea
return(DrugHomo)
}


CalculateDHI <- function(filepath,binned,mzs,QuantLevel=8, Bkg='T')
{
if(lenght(mzs) < 1) stop('m/z value is missing')

if(length(mzs) ==3){
	mz_drug = mzs[1];mz_mask = mzs[2]; mz_std = mzs[3]
}else if(length(mzs) ==2){
	mz_drug = mzs[1];mz_mask = mzs[2]
}else{
	mz_drug= mzs[1]; mz_mask = mzs[1]}
	
analyfie1 = importAnalyze(filepath)
IntenMatrix = matrix(0,nrow=length(analyfie1),ncol = length(binned[,1])) 

for(i in 1:length(analyfie1))
{
if(max(analyfie1[[i]]@intensity) < 100)   ### To avoid error due to the noisy spectrum 
{IntenMatrix[i,] =0}
 
else
{
z <- msSet(as.matrix(analyfie1[[i]]@intensity),mz = as.vector(analyfie1[[i]]@mass))
z <- msPeak(z, FUN="simple", use.mean=FALSE, snr=median(sort(unique(analyfie1[[i]]@intensity[2:100]))),span=5)
mspeaks = cbind(z$peak.list[[1]]$mass.loc,z$peak.list[[1]]$mass.right,z$peak.list[[1]]$mass.left)
for(j in 1:dim(mspeaks)[1])
{
if(mspeaks[j,1]> (max(mzs)+1))
{
break
}

else
{
id = which(mspeaks[j,1] > binned[,1] & mspeaks[j,1] < binned[,2])
## id length zero if peaks appear on bin boundary, then substract some value from peak to find m/z bin
if(length(id) ==0)
{id = which((mspeaks[j]-0.01) > binned[,1] & (mspeaks[j]-0.01) < binned[,2])}
id = id[[1]]
if(IntenMatrix[i,id] !=0){
IntenMatrix[i,id] = max(IntenMatrix[i,id],z$peak.list[[1]]$intensity[j])}   #### If already peak present then select the one with high intensity value 
else
{IntenMatrix[i,id] = z$peak.list[[1]]$intensity[j]}
}
}
}
} 

IntenMatrix[is.na(IntenMatrix)] = 0


###### Create ion images
x = analyfie1[[length(analyfie1)]]@metaData$imaging$pos[[1]] ;y = analyfie1[[length(analyfie1)]]@metaData$imaging$pos[[2]]

## Find drug and tissue related ion bins 
drug_bin = which(mz_drug > binned[,1] & mz_drug < binned[,2]);
mask_bin = which(mz_mask > binned[,1] & mz_mask < binned[,2]); 

## Create drug and tissue ion images 

Img_drug =  IntenMatrix[,drug_bin]; dim(Img_drug) = c(y,x); Img_mask =  IntenMatrix[,mask_bin]; dim(Img_mask) = c(y,x)
 
if(!missing(mz_std))
{
instd_bin = which(mz_std > binned[,1] & mz_std < binned[,2]);
Img_std =  IntenMatrix[,mz_std]; dim(Img_std) = c(y,x)
Img_mask = (Img_mask/(Img_std+1))
Img_drug = (Img_drug/(Img_std+1))
}
 

##### Masking based on density plot 

mask = medianFilterR(Img_mask)
d = density(mask)
maskvalue = d$x[minimums(d$y)[1]]
mask = mask > maskvalue

par(mfrow=c(1,3))
image(medianFilterR(Img_drug), axes = FALSE);image(medianFilterR(Img_mask), axes = FALSE);image(mask, axes = FALSE)
savePlot(paste(filename,".png",sep=""),type = "png")

m = QuantLevel/max(data_drug)
Img_quantized = Img_drug*m
Img_quantized = ceiling(Img_quantized)
Img_quantized = medianFilterR(Img_quantized)
Img_quantized = mask * Img_quantized
write.table(mask,paste(filename,"_MaskedImg.csv",sep=""),sep=',');write.table(Img_quantized,paste(filename,"_DrugQuantizedImg.csv",sep=""),sep=',')

TumorArea = table(mask)[[2]]
DHIvalue = DHI(Img_quantized,Bkg='T')/TumorArea

return(DHIvalue)
}



  
