#' The point sources filter
#'
#' This function filters the point sources from UCERF3 within the user-defined maximum distance
#' from a site
#'
#' @param slat The latitude of the interest site
#' @param slon The longitude of the interest site
#' @param max_dist The maximum cut-off distance for considered sources. Default is 300 km.
#' @param branch The indicator specifies which branch of UCERF3 is using. 1 for 3.1 branch and 2 for 3.2 branch.
#'
#' @return A filtered point sources dataframe is returned
#'
#' @examples flts_filter(slat = 36.1, slon = -118.1)
#'
#' @export
pts_filter <- function(slat, slon, max_dist = 300, branch = 1) {

  if (branch == 1) {

    FM_Pts <- FM31_Pts

  } else if (branch == 2) {

    FM_Pts <- FM32_Pts

  }

  Tmp_dist <- geosphere::distHaversine( c( slon, slat ), cbind( FM_Pts$Longitude, FM_Pts$Latitude ) ) / 1000

  Idx_Pts <- which( Tmp_dist < max_dist )

  FM_Pts$Rrup <- Tmp_dist

  FM_Pts$Rjb <- Tmp_dist

  FM_Pts_output <- FM_Pts[ Idx_Pts, ]

  return(FM_Pts_output)
}

