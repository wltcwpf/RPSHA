#' The UCERF3 3.2 branch average model for multi-fault rupture sources
#'
#' A large list containing 3.2 branch average model for all multi-fault rupture sources from UCERF3.
#' The data is originally from USGS nshm-cous-2014 github repository.
#'
#' @format A large list with 305,608 rows (number of unique ruptures) and 9 elements
#' (2 elements for the participating fault sub-section indices, see \code{\link{FM32_Secs}} for more information;
#' following by magnitude, occurrence rate,
#' type, depth, dip, rake, and width). The units of depth and width are km, and units of dip and rake are
#' degree. Note the total number of unique ruptures in UCERF 3 by Field et al (2014) is 305,709
#' (101 ruptures more than the data showed here). The data here is re-formatted from
#' USGS nshm-cous-2014 github repository. We may update this list once this inconsistency is resolved.
#'
#' @source USGS github repository: \url{https://github.com/usgs/nshm-cous-2014}
#'
#' UCERF3 Reference: \url{https://doi.org/10.1785/0120130164}
#'
"FM32_Rups"
