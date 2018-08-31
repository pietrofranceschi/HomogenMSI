#' Perform a permutation test on DHI
#' @name testDHI
#' @description This function perform a permutation testing to assess the significance of the Drug Homogeneity Index.
#' @usage testDHI(drugImg,permutations = 100, Nu=1, QuantLevel=NULL, maskImg=NULL)
#' @param drugImg Input matrix containing the pixel intensities
#' @param permutations Number of permutations used to construct the null distribution. default = 100
#' @param Nu Lowest size-zone value used for DHI calculation. default = 1
#' @param QuantLevel Maximum possible gray-levels in drug ion image. default = 0, i.e., original image
#' @param maskImg Matrix identifying the pixels belonging to the tissue (1 tissue, 0 background). It must have the same size as drugImg
#'
#' @return A list containing the DHI calculated for the original and the permuted images
#' @export
#' @author Mridula Prasad \email{mridula.prasad@fmach.it}
#' @references \url{https://github.com/pietrofranceschi/HomogenMSI}
#'
#' @examples

#' ## load package
#' library(HomogenMSI)
#' data("DHIimages")
#'
#' ## Input drug and mask images
#'
#' drugImg = DHIimages[[20]]
#' maskImg = DHIimages[[1]]
#' maskImg[maskImg !=0] =1
#'
#' # test the significance of the DHI
#' output <- testDHI(drugImg,maskImg=maskImg)
#' 
#' # Plot the results
#' hist(output$DHIPermutedvalues, 
#'      xlim = range(c(output$DHIoriginalvalue,output$DHIPermutedvalues)), 
#'      col = "steelblue", 
#'      xlab = "DHI",
#'      main = "DHI Permutation test")
#' abline(v = output$DHIoriginalvalue, lty = 2)



testDHI <- function(drugImg,
                    permutations=100,
                    Nu=1,
                    QuantLevel=NULL,
                    maskImg=NULL)
{
  ## Check : if there are only zeroes stop and return a warning
  
  if(all(drugImg == 0)){
    stop( "The image contains only zeroes" )
  }
  
  ## If a mask  is present, multiply it with drug image and estimate tumor area, 
  ## otherwise set the tumor area to the image size
  if(!is.null(maskImg))
  {
    if(!identical(dim(drugImg),dim(maskImg))){
      stop("dimensions of drug and mask image matrices are different")
      }
    
    
  
    DHIOriginalvalue <- CalculateDHI(drugImg,
                                     Nu=Nu,
                                     QuantLevel=QuantLevel,
                                     maskImg=maskImg)
    
    ## Extract intensities of pixels within tumor tissue 
    
    TissuePixels <-  which(as.vector(maskImg) == 1)
    TissueIonIntensities <-  as.vector(drugImg)[TissuePixels]
    DHIPermutedvalues <-c()
    
    
    for(i in 1:permutations)
    {   
      PermutatedImg <-  maskImg
      PermutatedImg[PermutatedImg == 1] <-  sample(TissueIonIntensities,length(TissueIonIntensities))
      DHIPermutedvalues[i] <- CalculateDHI(
        PermutatedImg,
        Nu=Nu,
        QuantLevel=QuantLevel,
        maskImg=maskImg)
      print(i)
    }
    
  } else {
    DHIOriginalvalue <- CalculateDHI(drugImg,
                                     Nu=Nu,
                                     QuantLevel=QuantLevel)
    DHIPermutedvalues <- c()
    
    for(i in 1:permutations)
    {   
      PermutatedImg <-  drugImg
      PermutatedImg <-  sample(PermutatedImg,length(PermutatedImg))
      DHIPermutedvalues[i] <- CalculateDHI(PermutatedImg,
                                           Nu=Nu,
                                           QuantLevel=QuantLevel)
      print(i)
    }		
  }
  Output <-  list(DHIoriginalvalue = DHIOriginalvalue, 
                  DHIPermutedvalues = DHIPermutedvalues)
  return(Output)
}

