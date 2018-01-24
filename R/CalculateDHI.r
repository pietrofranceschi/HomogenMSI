#' Calculate DHI value for drug MSI data
#'
#' @name CalculateDHI
#' @description This function calculates the DHI value for the drug image extracted from MSI data.
#' @usage CalculateDHI(drugImg,Nu=1,QuantLevel=NULL,maskImg=NULL)
#' @param  drugImg  Input matrix containing the pixel intensities
#' @param  maskImg  Matrix identifying the pixels belonging to the tissue (1 tissue, 0 background). It must have the same size as drugImg
#' @param  QuantLevel  Maximum possible gray-levels in drug ion image. default = 0, i.e., original image
#' @param Nu  Lowest size-zone value used for DHI calculation. default = 1
#'
#' @details
#' The algorithm is based on the following workflow:
#' \enumerate{
#' \item Derive new quantized drug image at user-defined qunatization value.
#' \item Derive gray-level size-zone matrix (GLSZM) for new drug quantized image (default value =0).
#' \item Derive homogeneous size-zone value from GLSZM at user-defined Nu value (default value =1).
#' \item Normalized overall value with complete tumor area which is obtained from tumor mask file.
#' }
#' @return DHI value for given drug ms file
#' @author Mridula Prasad \email{mridula.prasad@fmach.it}
#' @references \url{https://github.com/pietrofranceschi/HomogenMSI}
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
#' #Calculate DHI with default input parameters
#' print(CalculateDHI(drugImg,maskImg=maskImg))
#'
#' # calculate DHI with user-defined Nu value
#' print(CalculateDHI(drugImg,maskImg=maskImg,QuantLevel=0,Nu=5))
#'
#' @export
#'

suppressMessages(library("spatstat"))


CalculateDHI <- function(drugImg,Nu=1,QuantLevel=NULL,maskImg=NULL)
{
  ## Quantize image to the user-defined number of gray-levels
  if(!missing(QuantLevel))
  {
    m = QuantLevel/max(drugImg)
    drugImg = drugImg*m
    drugImg = round(drugImg,0)

  }

  ## If a mask  is present, multiply it with drug image and estimate tumor area, otherwise set the tumor area to the image size
  if(!missing(maskImg))
  {
    if((dim(drugImg)[1] != dim(maskImg)[1]) | (dim(drugImg)[2] != dim(maskImg)[2]))
    {stop("dimensions of drug and mask image matrices are different")}

    drugImg= maskImg * drugImg
    TumorArea = table(maskImg)[[2]]

  } else {
    TumorArea <-  length(drugImg)
  }

  ## Additional chech: if there are only zeroes stop and return a warning
  if(all(drugImg == 0)){
    stop( "The image contains only zeroes" )
  }
  
  
  # unique number of gray levels in image
  grey_lvls <- unique(c(drugImg))
  ## keep only non NAs
  grey_lvls <- grey_lvls[!is.na(grey_lvls)]

  #convert to data for use with spatstats functions
  ImgMat = spatstat::as.im(drugImg)
  
  szmlist <- do.call(rbind,lapply(grey_lvls, function(l){
    imBinary <- spatstat::levelset(ImgMat, l, compare="==")
    connections <- spatstat::connected(imBinary)
    counts <- table(as.matrix(connections))
    out <- data.frame("grl" = l,
                      "size" = as.matrix(counts))
    return(out)
  }))
  
  rownames(szmlist) <- NULL
  
  ## initialize the glszm
  glszm <- matrix(0,
                 nrow = length(unique(szmlist$grl)), 
                 ncol = length(unique(szmlist$size)))
  rownames(glszm) <- unique(szmlist$grl)
  colnames(glszm) <- sort(unique(szmlist$size))
  
  ## fill ir with the numbers
  for (i in 1:nrow(szmlist)){
    indices <- as.character(szmlist[i,])
    glszm[indices[1],indices[2]] <- glszm[indices[1],indices[2]] + 1
  }
  
  ## Remove zeroes
  glszm <-  glszm[!rownames(glszm) == "0",, drop = FALSE]   
  glszm <-  glszm[,!colSums(glszm) == 0, drop = FALSE]
  
  ## Calculate the DHI
  szv <-  as.numeric(colnames(glszm))
  idin <- szv >= Nu
  
  DHI <- sum(glszm[,idin] %*% szv[idin])/sum(glszm[,idin])
  DHI <- DHI/TumorArea
  
  return(DHI)
}




