\name{summary.semiMarkov}
\alias{summary.semiMarkov}
\title{
Summary method for objects of class \code{semiMarkov}
}
\description{
Summary method for objects of class \code{semiMarkov}.
}
\usage{
\method{summary}{semiMarkov}(object, all = TRUE, transitions = NULL, ...)
}
\arguments{
  \item{object}{
An object of class \code{semiMarkov}.
}
  \item{all}{
  A logical value indicating if the results should be displayed for all the possible transitions. If set to \code{FALSE}, the transitions to be displayed must be specified using the argument \code{transitions}. Default is \code{TRUE}.
  }
  \item{transitions}{
   A vector of characters specifying the transitions to be displayed when the argument \code{all} is set to \code{FALSE}.
  }
  \item{\dots}{
  Further arguments for summary.
}
}
\value{
A list of data frames giving 
%the possible transitions, transitions affected by covariates, parameters estimation  
%estimated probabilities for the embedded Markov Chain and estimated parameters for each transition. List items for the item Transition are named after the possible transition.
  
 \item{Transition_matrix}{
A matrix containing the informations on the model definition : the possible transitions and the distribution of waiting times for each transition (Exponential, Weibull or Exponentiated Weibull).
}
 \item{param.init}{
Recall the initial values of the parameters. The third column of this object can be used in \code{hazard} function.
}
\item{table.state}{
A table, with starting states as rows and arrival states as columns, which provides the number of observed transitions between two states. This argument can be used to quickly summarize multi-state data.
}
\item{Ncens}{
Number of individuals subjected to censoring.
}


%\item{Transition.probability}{
% Quadratic matrix giving the estimation of the transition probabilities of the Markov chain.
%}

 %\item{table.dist}{
 %A data frame giving the distribution parameters estimation results. It consists of the estimated values, the standard deviations, the confidence intervals, the Wald test statistics and Wald test p-values.
 %}
 %\item{table.coef}{
 %A data frame giving the regression coefficients estimation results. It consists of the estimated values, the standard deviations, the confidence intervals, the Wald test statistics and Wald test p-values.
 % Likelihood ratio statisitc and p-value of Likelihood ratio test.
% }
 \item{table.param}{
 List of data frames (one for each transition). A data frame includes, for each parameter (distribution parameters, the transition probabilities and the regression coefficients), the estimation, the standard deviation, the lower and upper bounds of confidence interval, the Wald test statistic and Wald test p-value (for the distribution parameters and the regression coefficients).
 
}
}
\author{
Agnieszka Listwon-Krol
}
\seealso{
\link{semiMarkov}, \link{print.semiMarkov}
}
 \keyword{documentation}

