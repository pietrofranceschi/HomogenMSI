#' Calculate DHI value for drug MSI data
#'
#' @name CalculateDHI
#' @description This function calculates the DHI value for drug MSI data.
#' @usage CalculateDHI(drug_csvfilepath,mask_csvfilepath,QuantLevel=32,Nu=1)
#' @param
#' @param  drug_csvfilepath  path of drug ion image file (csv format)
#' @param  mask_csvfilepath  path of tissue masked file (csv format)
#' @param  QuantLevel  Maximum possible gray-levels in drug ion image. default = 32
#' @param Nu  Lowest size-zone value used for DHI calculation. default = 1
#'
#' @details
#' The algorithm is based on the following workflow:
#' \enumerate{
#' \item Derive new quantized drug image at user-defined qunatization value.
#' \item Derive gray-level size-zone matrix (GLSZM) for new drug quantized image (default value =32).
#' \item Derive homogeneous size-zone value from GLSZM at user-defined Nu value (default value =1).
#' \item Normalized overall value with complete tumor area which is obtained from tumor mask file.
#' }
#' @return DHI value for given drug ms file
#' @author Mridula Prasad \email{mridula.prasad@fmach.it}
#' @seealso \code{\link{PreprocessingMSI}}
#' @references \url{https://github.com/pietrofranceschi/HomogenMSI}
#' @examples
#' ## load package
#' library(HomogenMSI)
#' data("DHIimages")
#'
#' ## Input drug and mask images
#'
#' drugImg = DHIimages[[2]]
#' maskImg = DHIimages[[1]]
#' maskImg[maskImg !=0] =1
#'
#' #Calculate DHI with default input parameters
#' print(CalculateDHI(drugImg,maskImg))
#'
#' # calculate DHI with user-defined Nu value
#' print(CalculateDHI(drugImg,maskImg,QuantLevel=0,Nu=5))
#'
#' @export
#'

suppressMessages(library("radiomics"))

##### calculate DHI value from given MSI data


CalculateDHI <- function(drugImg,maskImg,QuantLevel=0,Nu=1)
{
 if((dim(drugImg)[1] != dim(maskImg)[1]) | (dim(drugImg)[2] != dim(maskImg)[2]))
   {stop("dimensions of drug and mask image matrices are different")}
  if(QuantLevel !=0)
  {
  m = QuantLevel/max(drugImg)
  drugImg = drugImg*m
  }
  drugImg= maskImg * drugImg

  szm = glszm(drugImg)
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
  TumorArea = table(maskImg)[[2]]
  DrugHomo = DrugHomo/TumorArea

  return(DrugHomo)
}
