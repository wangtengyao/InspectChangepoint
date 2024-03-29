% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute.threshold.R
\name{compute.threshold}
\alias{compute.threshold}
\title{Computing threshold used in \code{inspect}}
\usage{
compute.threshold(n, p, nrep = 100, show_progress = TRUE)
}
\arguments{
\item{n}{Time length of the observation.}

\item{p}{Dimension of the multivariate time series.}

\item{nrep}{Number of Monte Carlo repetition to be used.}

\item{show_progress}{whether to show the progress of Monte Carlo simulation}
}
\value{
A numeric value indicating the threshold level that should be used based on the Monte Carlo simulation.
}
\description{
The threshold level to be used in \code{inspect} is computed via Monte Carlo simulation of multivariate time series that do not contain any changepoints.
}
\examples{
compute.threshold(n=200, p=50)
}
