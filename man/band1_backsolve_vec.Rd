% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{band1_backsolve_vec}
\alias{band1_backsolve_vec}
\title{Backsolve with band 1 upper Cholesky.}
\usage{
band1_backsolve_vec(U, z)
}
\arguments{
\item{U}{An upper triangular square matrix with non-zero entries only on the
main diagonal and the first superdiagonal.}

\item{z}{A vector with as many elements as the number of rows of U.}
}
\value{
A vector.
}
\description{
Backsolve with band 1 upper Cholesky.
}
\keyword{internal}
