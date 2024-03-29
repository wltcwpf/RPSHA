% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pts_filter.R
\name{pts_filter}
\alias{pts_filter}
\title{The point sources filter}
\usage{
pts_filter(slat, slon, max_dist = 300, branch = 1)
}
\arguments{
\item{slat}{The latitude of the interest site}

\item{slon}{The longitude of the interest site}

\item{max_dist}{The maximum cut-off distance for considered sources. Default is 300 km.}

\item{branch}{The indicator specifies which branch of UCERF3 is using. 1 for 3.1 branch and 2 for 3.2 branch.}
}
\value{
A filtered point sources dataframe is returned
}
\description{
This function filters the point sources from UCERF3 within the user-defined maximum distance
from a site
}
\examples{
flts_filter(slat = 36.1, slon = -118.1)

}
