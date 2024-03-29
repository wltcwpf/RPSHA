% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance_calc.R
\name{RrupRjb}
\alias{RrupRjb}
\title{The source-to-site distances calculator}
\usage{
RrupRjb(
  fflat1,
  fflon1,
  fflat2,
  fflon2,
  fdip,
  fdipDir = NA,
  topd,
  width,
  slon,
  slat
)
}
\arguments{
\item{fflat1}{Latitude of the start point of the trace of a rectangle finite fault}

\item{fflon1}{Longitude of the start point of the trace of a rectangle finite fault}

\item{fflat2}{Latitude of the end point of the trace of a rectangle finite fault}

\item{fflon2}{Longitude of the end point of the trace of a rectangle finite fault}

\item{fdip}{The dip angle of the rectangle finite fault subsection}

\item{fdipDir}{The dip direction. The default is NA, which assumes trace follows right-hand rule such
that hanging wall on the right side of trace}

\item{topd}{The depth to the top of the rectangle finite fault, in km}

\item{width}{The width of the rectangle finite fault (perpendicular to the trace), in km}

\item{slon}{A list of longitude of the interest sites}

\item{slat}{A list of latitude of the interest sites}
}
\value{
A data frame of distances with number of rows equals to the number of interest sites and two
columns
}
\description{
This function calculates rupture distance and Rjb distance for a given rectangle finite fault to the given sites.
For a multi-fault rupture source, the final distances should be calculated for each rectangle subsection
and then taken the minimum values.
}
\examples{
RrupRjb(fflat1 = 35.74054, fflon1 = -117.74953, fflat2 = 35.81038,
fflon2 = -117.76365, fdip = 50, topd = 0, width = 13, slon = -118.1, slat = 36.1)

RrupRjb(fflat1 = 35.74054, fflon1 = -117.74953, fflat2 = 35.81038,
fflon2 = -117.76365, fdip = 50, topd = 0, width = 13,
slon = c(-118.1, -118.2), slat = c(36.1, 35.9))

}
\references{
David Eberly (1999). Distance Between Point and Triangle in 3D. Geometric Tools, Redmond WA 98052.
\url{http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf}
}
