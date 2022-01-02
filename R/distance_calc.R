#' The source-to-site distances calculator
#'
#' This function calculates rupture distance and Rjb distance for a given rectangle finite fault to the given sites.
#' For a multi-fault rupture source, the final distances should be calculated for each rectangle subsection
#' and then taken the minimum values.
#'
#' @param fflat1 Latitude of the start point of the trace of a rectangle finite fault
#' @param fflon1 Longitude of the start point of the trace of a rectangle finite fault
#' @param fflat2 Latitude of the end point of the trace of a rectangle finite fault
#' @param fflon2 Longitude of the end point of the trace of a rectangle finite fault
#' @param fdip The dip angle of the rectangle finite fault subsection
#' @param fdipDir The dip direction. The default is NA, which assumes trace follows right-hand rule such
#' that hanging wall on the right side of trace
#' @param topd The depth to the top of the rectangle finite fault, in km
#' @param width The width of the rectangle finite fault (perpendicular to the trace), in km
#' @param slon A list of longitude of the interest sites
#' @param slat A list of latitude of the interest sites
#' @return A data frame of distances with number of rows equals to the number of interest sites and two
#' columns
#'
#' @examples RrupRjb(fflat1 = 35.74054, fflon1 = -117.74953, fflat2 = 35.81038,
#' fflon2 = -117.76365, fdip = 50, topd = 0, width = 13, slon = -118.1, slat = 36.1)
#'
#' RrupRjb(fflat1 = 35.74054, fflon1 = -117.74953, fflat2 = 35.81038,
#' fflon2 = -117.76365, fdip = 50, topd = 0, width = 13,
#' slon = c(-118.1, -118.2), slat = c(36.1, 35.9))
#'
#' @references David Eberly (1999). Distance Between Point and Triangle in 3D. Geometric Tools, Redmond WA 98052.
#' \url{http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf}
#'
#' @export
RrupRjb <- function( fflat1, fflon1, fflat2, fflon2, fdip, fdipDir = NA, topd, width, slon, slat ){
  #     This function calculates Rrup and Rjb from the given finite fault model to stations

  botd <- topd + width * sin( fdip / 180 * pi )  # Depth to the bottom of rupture

  rwh <- width * cos( fdip / 180 * pi )  # Horizontal projection of rupture width

  Sxy1 <- CLL2XY_PW( elon = fflon1, elat = fflat1, slon = slon, slat = slat )

  Sxy2 <- CLL2XY_PW( elon = fflon2, elat = fflat2, slon = slon, slat = slat )

  Rrup <- rep( 0, length( slon ) )

  Rjb <- rep( 0, length( slon ) )

  for( i_sta in 1:length( slon ) ){

    pt <- matrix( 0, nrow = 3, ncol = 4 )

    # comment: change X-Y coordinate to centered at site i
    pt[ 1, 1 ] <- -Sxy1$X[ i_sta ]

    pt[ 2, 1 ] <- -Sxy1$Y[ i_sta ]

    pt[ 3, 1 ] <- -topd

    pt[ 1, 2 ] <- -Sxy2$X[ i_sta ]

    pt[ 2, 2 ] <- -Sxy2$Y[ i_sta ]

    pt[ 3, 2 ] <- -topd

    if( is.na( fdipDir ) ){

      dY <- pt[ 2, 2 ] - pt[ 2, 1 ]

      dX <- pt[ 1, 2 ] - pt[ 1, 1 ]

      Tmp_Strike <- atan( dY / dX ) / pi * 180

      if( dX == 0 & dY > 0 ){

        Strike <- 0

      }else if( dX == 0 & dY < 0 ){

        Strike <- 180

      }else if( dX > 0 ){

        Strike <- 90 - Tmp_Strike

      }else if( dX < 0 ){

        Strike <- 270 - Tmp_Strike
      }

      fdipDir <- ifelse( Strike + 90 >= 360, Strike + 90 - 360, Strike + 90 )
    }

    dX <- rwh * sin( fdipDir / 180 * pi )

    dY <- rwh * cos( fdipDir / 180 * pi )

    pt[ 1, 3 ] <- pt[ 1, 2 ] + dX

    pt[ 2, 3 ] <- pt[ 2, 2 ] + dY

    pt[ 3, 3 ] <- -botd

    pt[ 1, 4 ] <- pt[ 1, 1 ] + dX

    pt[ 2, 4 ] <- pt[ 2, 1 ] + dY

    pt[ 3, 4 ] <- -botd

    a <- DEMDT( pt )

    Rrup[ i_sta ] <- a[[ 1 ]]

    pt[ 3, 1 ] <- 0.0

    pt[ 3, 2 ] <- 0.0

    pt[ 3, 3 ] <- 0.0

    pt[ 3, 4 ] <- 0.0

    a <- DEMDT( pt )

    Rjb[ i_sta ] <- a[[ 1 ]]
  }

  Res <- as.data.frame( cbind( Rrup, Rjb ) )

  return( Res )
}


## helper functions
CLL2XY_PW <- function( elon, elat, slon, slat ){
  #     This subroutine uses the haversine formula to calculate distance, bearing, and
  #     convert latitude/longitude to x/y coordinates on the basis of a spherical
  #     earth (ignoring ellipsoidal effects).
  #
  #     LL to XY centered on elon,elat
  #     outputs: a data.frame of X, Y
  #     Reference: http://www.movable-type.co.uk/scripts/latlong.html

  ifelse( length( slon ) == 1, Tmp_Dist <- geosphere::distVincentyEllipsoid( c( elon, elat ), c( slon, slat ) ),
          Tmp_Dist <- geosphere::distVincentyEllipsoid( c( elon, elat ), cbind( slon, slat ) ) )

  ifelse( length( slon ) == 1, Tmp_Bearing <- geosphere::bearing( c( elon, elat ), c( slon, slat ) ),
          Tmp_Bearing <- geosphere::bearing( c( elon, elat ), cbind( slon, slat ) ) )

  X <- Tmp_Dist / 1000 * sin( Tmp_Bearing / 180 * pi )

  Y <- Tmp_Dist / 1000 * cos( Tmp_Bearing / 180 * pi )

  Res <- as.data.frame( cbind( X, Y ) )

  return( Res )
}


DEMDT <- function(P){
  #     input P[3,4]
  #     return a list with P_clst(3) and dclst
  #     The minimum distance to triangle based on "David Eberly, 'Distance Between Point and Triangle in 3D',
  #     Geometric Tools, LLC, (1999)",
  #     Reference: http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

  P_clst <- c()
  inside <- 0
  # Top triangle (P1P2P4) of trapezoid
  B <- c()
  E0 <- c()
  E1 <- c()
  E_clst <- matrix(data <- 0, nrow = 3, ncol = 4)
  B[1] = P[1,1]
  B[2] = P[2,1]
  B[3] = P[3,1]
  E0[1] = P[1,2]-B[1]
  E0[2] = P[2,2]-B[2]
  E0[3] = P[3,2]-B[3]
  E1[1] = P[1,4]-B[1]
  E1[2] = P[2,4]-B[2]
  E1[3] = P[3,4]-B[3]
  a1 = E0[1]*E0[1]+E0[2]*E0[2]+E0[3]*E0[3]
  b1 = E0[1]*E1[1]+E0[2]*E1[2]+E0[3]*E1[3]
  c1 = E1[1]*E1[1]+E1[2]*E1[2]+E1[3]*E1[3]
  d1 = E0[1]*B[1]+E0[2]*B[2]+E0[3]*B[3]
  e11 = E1[1]*B[1]+E1[2]*B[2]+E1[3]*B[3]
  detinv = 1.0/(a1*c1-b1*b1)
  s1 = (b1*e11-c1*d1)*detinv
  t1 = (b1*d1-a1*e11)*detinv
  if(detinv != Inf & s1 >= 0 & s1 <= 1 & t1 >= 0 & t1 <= 1 & s1+t1 <= 1){
    inside <- 1
    P_clst[1] = B[1]+s1*E0[1]+t1*E1[1]
    P_clst[2] = B[2]+s1*E0[2]+t1*E1[2]
    P_clst[3] = B[3]+s1*E0[3]+t1*E1[3]
    dclst = P_clst[1]^2+P_clst[2]^2+P_clst[3]^2
    dclst = sqrt(dclst)
    return(list(dclst, P_clst))
  }
  # Bottom triangle (P3P2P4) of trapezoid
  B[1] = P[1,3]
  B[2] = P[2,3]
  B[3] = P[3,3]
  E0[1] = P[1,2]-B[1]
  E0[2] = P[2,2]-B[2]
  E0[3] = P[3,2]-B[3]
  E1[1] = P[1,4]-B[1]
  E1[2] = P[2,4]-B[2]
  E1[3] = P[3,4]-B[3]
  a2 = E0[1]*E0[1]+E0[2]*E0[2]+E0[3]*E0[3]
  b2 = E0[1]*E1[1]+E0[2]*E1[2]+E0[3]*E1[3]
  c2 = E1[1]*E1[1]+E1[2]*E1[2]+E1[3]*E1[3]
  d2 = E0[1]*B[1]+E0[2]*B[2]+E0[3]*B[3]
  e2 = E1[1]*B[1]+E1[2]*B[2]+E1[3]*B[3]
  detinv = 1.0/(a2*c2-b2*b2)
  s2 = (b2*e2-c2*d2)*detinv
  t2 = (b2*d2-a2*e2)*detinv
  if (detinv != Inf & s2 >= 0 & s2 <= 1 & t2 >= 0 & t2 <= 1 & (s2+t2) <= 1){
    inside = 2
    P_clst[1] = B[1]+s2*E0[1]+t2*E1[1]
    P_clst[2] = B[2]+s2*E0[2]+t2*E1[2]
    P_clst[3] = B[3]+s2*E0[3]+t2*E1[3]
    dclst = P_clst[1]^2+P_clst[2]^2+P_clst[3]^2
    dclst = sqrt(dclst)
    return(list(dclst, P_clst))
  }

  ## Four edges of trapezoid
  # Edge 1: Top edge (P1P2)
  if(-d1 < 0){
    f <- 0
  }else if(-d1 > a1){
    f <- 1
  }else{
    f <- -d1/a1
  }
  E_clst[1,1] = P[1,1] + f * (P[1,2]-P[1,1])
  E_clst[2,1] = P[2,1] + f * (P[2,2]-P[2,1])
  E_clst[3,1] = P[3,1] + f * (P[3,2]-P[3,1])
  # Edge 2: Side edge (P3P2)
  if(-d2 < 0){
    f <- 0
  }else if(-d2 > a2){
    f = 1
  }else{
    f = -d2/a2
  }
  E_clst[1,2] = P[1,3] + f * (P[1,2]-P[1,3])
  E_clst[2,2] = P[2,3] + f * (P[2,2]-P[2,3])
  E_clst[3,2] = P[3,3] + f * (P[3,2]-P[3,3])
  # Edge 3: Bottom edge (P3P4)
  if(-e2 < 0){
    f = 0
  }else if(-e2 > c2){
    f = 1
  }else{
    f = -e2/c2
  }
  E_clst[1,3] = P[1,3] + f * (P[1,4]-P[1,3])
  E_clst[2,3] = P[2,3] + f * (P[2,4]-P[2,3])
  E_clst[3,3] = P[3,3] + f * (P[3,4]-P[3,3])
  # Edge 4: Side edge (P1P4)
  if(-e11 < 0){
    f = 0
  }else if(-e11 > c1){
    f = 1
  }else{
    f = -e11/c1
  }
  E_clst[1,4] = P[1,1] + f * (P[1,4]-P[1,1])
  E_clst[2,4] = P[2,1] + f * (P[2,4]-P[2,1])
  E_clst[3,4] = P[3,1] + f * (P[3,4]-P[3,1])

  ## Find closest point on the edges of Trapezoid
  dis1 = E_clst[1,1]^2+E_clst[2,1]^2+E_clst[3,1]^2
  dis2 = E_clst[1,2]^2+E_clst[2,2]^2+E_clst[3,2]^2
  dis3 = E_clst[1,3]^2+E_clst[2,3]^2+E_clst[3,3]^2
  dis4 = E_clst[1,4]^2+E_clst[2,4]^2+E_clst[3,4]^2
  dis_seq <- c(dis1, dis2, dis3, dis4)
  dclst = min(dis_seq, na.rm = T)
  iclst = which.min(dis_seq)
  dclst = sqrt(dclst)
  P_clst[1] = E_clst[1,iclst]
  P_clst[2] = E_clst[2,iclst]
  P_clst[3] = E_clst[3,iclst]
  return(list(dclst, P_clst))
}




