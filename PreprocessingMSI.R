#' Pre-process MSI data
#'
#' Derive drug and mask ion images from MSI data.
#' @name PreprocessingMSI
#' @param
#' @param filepath  path of MSI data folder in Analyze 7.5 format,
#' @param filename  filename to save drug and corresponding mask image
#' @param mzs  mass range for drug, tissue, internal standard (optional) ion images
#' @details
#' \enumerate{ The function performs data pre-processing using follwoing steps:
#' \item Create mass bin of size 0.5
#' \item Extract drug, tissue and internal standard peaks for individual spectrum within their respective bins based on local maxima search.
#' \item Create tumor masked image based on tissue ion image.
#' \item Remove spatial noise using median filtering.
#' }
#' @return
#' Returns csv files of drug ion and corresponding mask image for DHI calculation
#' @author Mridula Prasad \email{mridula.prasad@fmach.it}
#' @references \url{https://github.com/pietrofranceschi/HomogenMSI}
#' @examples
#' ## load library
#' library(HomogenMSI)
#'
#' filepath  = 'C:\\folder\\drug' ## MSI data folder path
#'
#' filename ='drug'  ## filename to save csv files
#' mzs = mzs = c(284.2, 281.2, 289.2)  ## mass list corresponds to the drug, tissue and internal standard mass
#'
#' output = PreprocessingMSI(filepath,filename,mzs)
#' @export

suppressMessages(library("MALDIquant"))
suppressMessages(library("MALDIquantForeign"))
suppressMessages(library("msProcess"))
suppressMessages(library("splus2R"))
suppressMessages(library("data.table"))

PreprocessingMSI <- function(filepath,filename,mzs)
{

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


  if(length(mzs) < 1) stop('m/z value is missing')

  if(length(mzs) ==3){
    mz_drug = mzs[1];mz_mask = mzs[2]; mz_std = mzs[3]
  }else if(length(mzs) ==2){
    mz_drug = mzs[1];mz_mask = mzs[2]
  }else{
    mz_drug= mzs[1]; mz_mask = mzs[1]}

  msifile = importAnalyze(filepath) # read file

  # create bins
  binlength = seq(range(msifile[[1]]@mass)[1]-0.5, range(msifile[[1]]@mass)[2]+0.5, by = 0.5) # for fixed bin size
  cuts = cut(msifile[[1]]@mass,binlength)
  duration.freq = table(cuts)
  duration.freq = cbind(duration.freq)
  cutdur = row.names(duration.freq)
  binsplit = function(x){as.numeric(unlist(strsplit(gsub("\\(|\\]", "", x),','))[1])}
  bin_even = sapply(cutdur,binsplit)
  binsplit = function(x){as.numeric(unlist(strsplit(gsub("\\(|\\]", "", x),','))[2])}
  bin_odd = sapply(cutdur,binsplit)
  binned = cbind(bin_even,bin_odd)
  binname = apply(binned,1,mean)


  IntenMatrix = matrix(0,nrow=length(msifile),ncol = length(binned[,1]))

  for(i in 1:length(msifile))
  {
    if(max(msifile[[i]]@intensity) < 100)   ### To avoid error due to the noisy spectrum
    {IntenMatrix[i,] =0}

    else
    {
      z <- msSet(as.matrix(msifile[[i]]@intensity),mz = as.vector(msifile[[i]]@mass))
      z <- msPeak(z, FUN="simple", use.mean=FALSE, snr=median(sort(unique(msifile[[i]]@intensity[2:100]))),span=5)
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
  x = msifile[[length(msifile)]]@metaData$imaging$pos[[1]] ;y = msifile[[length(msifile)]]@metaData$imaging$pos[[2]]

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
  write.csv(mask,paste(filename,"_MaskedImg.csv",sep=""));write.csv(Img_drug,paste(filename,"_DrugImg.csv",sep=""))

  jpeg(file = paste(filename,".jpeg",sep=""))
  par(mfrow=c(1,3))
  image(medianFilterR(Img_drug), axes = FALSE);image(medianFilterR(Img_mask), axes = FALSE);image(mask, axes = FALSE)
  dev.off()
  output = list(Img_drug,mask)
  return(output)
}



