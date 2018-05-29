\name{print.semiMarkov}
\alias{print.semiMarkov}
\title{
Print method for object of class \code{semiMarkov}
}
\description{
Print method for objects of class \code{semiMarkov}.
}
\usage{
\method{print}{semiMarkov}(x, CI=TRUE, Wald.test=TRUE, ...)
}

\arguments{
  \item{x}{
An object of class \code{semiMarkov}.
}
  \item{CI}{
  A logical value indicating if the confidence intervals for each parameter should be returned. Default is \code{TRUE}. The confidence level is chosen in \code{semiMarkov}.
  }
  \item{Wald.test}{
    A logical value indicating if the results of the Wald test for each parameter should be returned. Default is \code{TRUE}.
      }
  \item{\dots}{
Further arguments for print or summary.
}
}
%\details{
%The function prints the results of the parametric estimation obtained with function \code{semiMarkov}. It gives the information on the chosen model as well as the results of the estimation.
%including the estimated maximum likelihood, the parameters of the sojourn distribution, transition probability and covariates.
%}
\value{
No value returned.
}

\author{
Agnieszka Listwon-Krol
}

\seealso{
\link{semiMarkov}, \link{summary.semiMarkov}
}
 \keyword{documentation}
