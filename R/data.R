#' Set of Synthetic Images for Homogeneity Assessment
#'
#'A list contains 20 images under four different scenarios.
#' @format A data matrix of size 250 x 250
#' There are five basic images under four different scenarios.
#' Each image has a dimension of 250 x250 with the circular object of radius 100 pixels.
#'
#' Images 1 to 5 belong to the scenario 1 where all images have a single gray level with different number of dark zones.
#' Images 2 to 4 are derived from image 1 by gradually adding dark zones/ drug inaccessibility areas.
#' Image 1 does not contains any dark-zones thus gives an accurate delineation of tumor morphology and used to
#' calculate tumor area in drug homogeneity index formula.
#' In scenarios II and III, images (6 to 10 and 11 to 15) are produced from images in the scenario I (1 to 5) by adding
#' 2 and 4 gray levels, respectively. In last scenario IV, contains all images similar to the one in scenario III
#' with added salt-pepper noise.

