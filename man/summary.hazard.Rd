\name{summary.hazard}
\alias{summary.hazard}

\title{
Summary method for objects of class \code{hazard}
}
\description{
Summary method for objects of class hazard.
}
\usage{
\method{summary}{hazard}(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{hazard}.
}
  \item{\dots}{
Further arguments for summary.
}
}
\details{
For an object of class \code{hazard}, this function gives the informations on the type of hazard rates (sojourn time or semi-Markov process), the chosen model, the distribution of the sojourn times, the covariates and the vector of times.
}
\value{
No value returned.
}

\author{
Agnieszka Listwon-Krol
}

\seealso{
\link{hazard}, \link{print.hazard}
}
 \keyword{documentation}
 