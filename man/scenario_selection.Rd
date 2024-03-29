% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scenario_selection.R
\name{scenario_selection}
\alias{scenario_selection}
\title{A function for scenario (event or ground motion map) selection}
\usage{
scenario_selection(
  Y,
  X,
  min_hazard = 1e-07,
  output_dir = NULL,
  max_rate_multiplier = NULL,
  num_lambda = 1000,
  Weight = 1
)
}
\arguments{
\item{Y}{The target hazard and target hazard magnitude/distance distributions}

\item{X}{The hazard matrix produced by each of the considered scenarios}

\item{min_hazard}{The minimum hazard value in \code{Y} that is considered in regression.
This is used to avoid having extremely large weights. Note the weights is 1/\code{Y}.
The default value for \code{min_hazard} is 1e-7.}

\item{output_dir}{If not NULL, the results will be written into a few CSV files under the specified directory}

\item{max_rate_multiplier}{The maximum value for the multiplier of annual occurrence rate when running LASSO,
a positive numeric number. It is an additional constraint for scenario selection.
The rate multipliers are not straightly constrained: the rate multipliers could be much larger
when the number of selected scenarios is too small; and the rate multipliers will get better
constrained when the number of selected scenarios gets more.}

\item{num_lambda}{The maximum number of trail penalty parameters. A larger \code{num_lambda}
will search for more possible subset events.}

\item{Weight}{The weight vector for each element in \code{Y}. The default is 1 that indicates
using the default weight of scenario selection algorithm (which is 1/\code{Y}). If the other
weights are preferred, then a vector with the same length as \code{Y} for the desired weights
can be specified. Since the weight of 1/\code{Y} has been set in the function, so you'd need
to multiply your desired weights by \code{Y} to compensate the default weights.}
}
\value{
A list of three elements is returned. The first element is
the corresponding adjusted rate multipliers for each event
(or each column in \code{X}), the second element is recovered hazard by the selected subset events,
and the last element is the logarithmic mean squared error.
The events with positive factors are selected events and the events with factors = 0 are
NOT selected.
\item{betas}{A matrix of the adjusted rate multipliers for each scenario for each set of selected
scenarios. The dimension is total number of candidate scenarios by the number of set of selected
scenarios. The scenarios with zero multipliers are unselected}
\item{errors_mat}{A matrix with four columns: the number of selected scenarios following by
three error metrics (mean of squared logarithmic residuals, mean of absolute logarithmic
relative errors, and mean of absolute arithmetic relative errors)}
\item{Y_hat}{The predicted hazard by each set of selected scenarios}
}
\description{
This function takes the target hazard (and target hazard magnitude and/or distance distributions)
and the hazard (and hazard magnitude and/or distance distributions) produced by each of the
considered event scenarios or ground motion map scenarios to efficiently and
optimally select a series of subset scenarios with adjusted annual
occurrence rates by using LASSO regression.
}
