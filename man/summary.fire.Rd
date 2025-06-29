% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fire_summary.R
\name{summary.fire_matrix}
\alias{summary.fire_matrix}
\alias{print.summary.fire_matrix}
\alias{summary.fire_tensor}
\alias{print.summary.fire_tensor}
\title{Summary method for FIRe models}
\usage{
\method{summary}{fire_matrix}(object, ...)

\method{print}{summary.fire_matrix}(x, ...)

\method{summary}{fire_tensor}(object, ...)

\method{print}{summary.fire_tensor}(x, ...)
}
\arguments{
\item{object}{A \code{fire_matrix} or \code{fire_tensor} object}

\item{...}{Not used}

\item{x}{A \code{summary.fire_matrix} or \code{summary.fire_tensor} object}
}
\value{
For \code{summary.fire_matrix} and \code{summary.fire_tensor}:
Returns an object of class \code{summary.fire_matrix} or \code{summary.fire_tensor}
containing:
\itemize{
  \item Estimated parameters
  \item Convergence information
  \item Model fit statistics
  \item Kernel specifications
  \item Computation timing
}

The \code{print} methods for these objects display nicely formatted output to the console.
}
\description{
Provides a comprehensive summary of \code{fire_matrix} or \code{fire_tensor} object, including estimated parameters,
convergence information, model fit statistics, and kernel specifications.
}
\examples{
data(Manure)
mod <- fire(X = Manure$absorp[1:10,], Y = Manure$y$DM[1:10],
 dat_T = list(1:700), stop.eps = 2, maxiter = 4)
summary(mod)
}
