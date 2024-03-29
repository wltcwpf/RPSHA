% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MR_dist_converter.R
\name{Event_to_MR_Hazd}
\alias{Event_to_MR_Hazd}
\title{A function to convert an event to the target magnitude and/or distance hazard distribution}
\usage{
Event_to_MR_Hazd(
  MR_Hazd,
  Mag,
  AnnualOccurrenceRate,
  Fault_type,
  region = 1,
  Sitelist,
  distribution = "cdf"
)
}
\arguments{
\item{MR_Hazd}{A dataframe with six or seven (if both magnitude and distance hazard distributions are interested)
columns: "IM_type", "Deagg_level", "SiteName", "IM_level", "Hazd", "Magnitude"/"R_distance"
(or "MagRType" and "MagR"). Note the column names should exactly match these variable names.
\code{MR_Hazd} can be directly provided by the function of \code{Deagg_to_MR_Hazd}.}

\item{Mag}{The magnitude of the event}

\item{AnnualOccurrenceRate}{The annual occurrence rate of the event}

\item{Fault_type}{The fault type of the event: 0 for unspecified fault;
1 for strike-slip fault; 2 for normal fault; 3 for reverse fault}

\item{region}{The region of the event. Default is 1 for California}

\item{Sitelist}{A dataframe of sites with four columns: SiteName (should be consistent with
the \code{SiteName} in \code{MR_Hazd}), Vs30, z1, and Rjb. The example info for
SiteName, Vs30, and z1 can be found on this webpage:
\url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/SiteTable.csv}.
The example of Rjb can be found on this webpage:
\url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/EventSet.csv}.}

\item{distribution}{The distribution representation. "cdf" stands for cumulative hazard
and "pdf" stands for discrete bin hazard.}
}
\value{
A dataframe with six or seven (if both magnitude and distance hazard distributions are interested)
columns is returned, the columns inlcude IM_type, Deagg_level,
SiteName, IM_level, Hazd_event, Magnitude/R_distance (or MagRType and MagR)
}
\description{
This function takes a target magnitude and/or distance hazard distribution and an events to
calculate the hazard of the event at each of the corresponding magnitude/distance/
intensity measure type/intensity measure level or return period.
}
