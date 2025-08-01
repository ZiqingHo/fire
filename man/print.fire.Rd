% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fire_print.R
\name{print.fire}
\alias{print.fire}
\alias{print.fire_matrix}
\alias{print.fire_tensor}
\title{Print Method for FIRE Models}
\usage{
\method{print}{fire_matrix}(x, ...)

\method{print}{fire_tensor}(x, ...)
}
\arguments{
\item{x}{A model object of class \code{fire_matrix} or \code{fire_tensor}.}

\item{...}{Not used.}
}
\value{
The input object (invisibly) for piping. Prints to console:
\itemize{
  \item Model type (matrix or tensor input)
  \item Convergence status
  \item Data dimensions
  \item Marginal log-likelihood
  \item Estimated hyperparameters
}
}
\description{
Compact display of FIRE models for either matrix or tensor input data.
}
\examples{
# For matrix input
data(Manure)
mod <- fire(X = Manure$absorp[1:10,], Y = Manure$y$DM[1:10],
            dat_T = list(1:700), control = list(stop.eps = 2, maxiter = 4))
print(mod)

# For tensor input
data(Housing)
dat_T <- list(T1 = 1:4, T2 = 1:9)
mod <- fire(X = Housing$X[1:5,,], Y = Housing$y[1:5,2],
            kernels = list(kronecker_delta, kronecker_delta),
            kernels_params = list(NA, NA),
            dat_T = dat_T, control = list(stop.eps = 2, maxiter = 4))
print(mod)
}
\seealso{
\code{\link{fire}}
}
