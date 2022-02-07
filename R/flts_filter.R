#' The multi-fault rupture sources filter
#'
#' This function filters the multi-fault ruptures from UCERF3 within the user-defined maximum distance
#' from a site
#'
#' @param slat The latitude of the interest site
#' @param slon The longitude of the interest site
#' @param max_dist The maximum cut-off distance for considered sources. Default is 300 km.
#' @param branch The indicator specifies which branch of UCERF3 is using. 1 for 3.1 branch and 2 for 3.2 branch.
#'
#' @return A filtered rupture sources dataframe is returned
#'
#' @examples flts_filter(slat = 36.1, slon = -118.1)
#'
#' @export
flts_filter <- function(slat, slon, max_dist = 300, branch = 1) {

  if (branch == 1) {

    FM_Secs <- FM31_Secs

    FM_Rups <- FM31_Rups

  } else if (branch == 2) {

    FM_Secs <- FM32_Secs

    FM_Rups <- FM32_Rups

  }

  Secs_Rrup <- rep(NA, nrow(FM_Secs))

  Secs_Rjb <- rep(NA, nrow(FM_Secs))

  for (i in 1:nrow(FM_Secs)) {

    Tmp_Secs <- FM_Secs$Traces[i]

    Tmp_Dip <- as.numeric(FM_Secs$Dip[i])

    Tmp_Topd <- as.numeric(FM_Secs$Depth[i])

    Tmp_Width <- as.numeric(FM_Secs$LowerDepth[i]) - as.numeric(FM_Secs$Depth[i])

    Tmp_Traces <- Get_Traces(Tmp_Secs)

    Secs_Rrup[i] <- Get_MinRrupRjb(Tmp_Traces, Tmp_Dip, Tmp_Topd, Tmp_Width, slon, slat)[1]

    Secs_Rjb[i] <- Get_MinRrupRjb(Tmp_Traces, Tmp_Dip, Tmp_Topd, Tmp_Width, slon, slat)[2]
  }

  idx_sel <- which(Secs_Rrup <= max_dist)

  sel_indices <- as.numeric(FM_Secs$Index[idx_sel])

  sel_Rrup <- Secs_Rrup[idx_sel]

  sel_Rjb <- Secs_Rjb[idx_sel]

  Rrup <- unlist(lapply(FM_Rups$Indices, function (x) {
    tmp_idx <- which(sel_indices %in% x)
    ifelse (length(tmp_idx) > 0, min(sel_Rrup[tmp_idx]), NA)
  }))

  Rjb <- unlist(lapply(FM_Rups$Indices, function (x) {
    tmp_idx <- which(sel_indices %in% x)
    ifelse (length(tmp_idx) > 0, min(sel_Rjb[tmp_idx]), NA)
  }))

  idx <- which(!is.na(Rrup))

  FM_Rups_output <- data.frame(Index = FM_Rups$Index[idx],
                               Mag = FM_Rups$Mag[idx],
                               Rate = FM_Rups$Rate[idx],
                               Depth = FM_Rups$Depth[idx],
                               Dip = FM_Rups$Dip[idx],
                               Rake = FM_Rups$Rake[idx],
                               Width = FM_Rups$Width[idx],
                               Rrup = Rrup[idx],
                               Rjb = Rjb[idx],
                               ID = FM_Rups$ID[idx])

  return(FM_Rups_output)
}



# helper functions
# ------------------------------------------------------------------------------------------------------------
Get_Traces <- function( Traces ){

  # Traces: the traces for section

  # return: a data frame with two columns, lon and lat, for the trace

  Tmp <- strsplit( Traces, split = '\n' )[[ 1 ]][ -1 ]

  lon <- as.numeric( unlist( sapply( Tmp, function( x )
    as.numeric( strsplit( x, split = ',' )[[ 1 ]][ 1 ] ) ) ) )

  lat <- as.numeric( unlist( sapply( Tmp, function( x )
    as.numeric( strsplit( x, split = ',' )[[ 1 ]][ 2 ] ) ) ) )

  res <- as.data.frame( cbind( lon, lat ) )

  return( res )
}


Get_MinRrupRjb <- function( Traces, Dip, Topd, Width, Slon, Slat ){

  # this function calculates the minimum Rrup/Rjb from the trace (each of the sub-rectangles) to the site

  Rrup <- rep( NA, nrow( Traces ) - 1 )

  Rjb <- rep( NA, nrow( Traces ) - 1 )

  for( k in 1:( nrow( Traces ) - 1 ) ){

    if( ( Traces$lat[ k ] == Traces$lat[ k + 1 ] ) & ( Traces$lon[ k ] == Traces$lon[ k + 1 ] ) ){

      Tmp_Dist <- geosphere::distVincentyEllipsoid( c( Traces$lon[ k ], Traces$lat[ k ] ), c( Slon, Slat ) )

      Rrup[ k ] <- Tmp_Dist / 1000

      Rjb[ k ] <- Tmp_Dist / 1000

    }else{

      Tmp_Dist <- RrupRjb( fflat1 = Traces$lat[ k ], fflon1 = Traces$lon[ k ],
                           fflat2 = Traces$lat[ k + 1 ], fflon2 = Traces$lon[ k + 1 ],
                           fdip = Dip, topd = Topd, width = Width,
                           slon = Slon, slat = Slat )

      Rrup[ k ] <- Tmp_Dist$Rrup

      Rjb[ k ] <- Tmp_Dist$Rjb
    }
  }

  return( c(min(Rrup), min( Rjb )) )
}




