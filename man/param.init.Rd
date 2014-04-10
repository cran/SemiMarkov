\name{param.init}
\alias{param.init}
\title{
Defines the initial values of  parameters for a semi-Markov model
}
\description{
Function defining initial values of parameters of the waiting time distributions, probabilities of the Markov chain and optional regression coefficients associated with covariates. The function can either provides the default initial values (the same as those considered in the function \code{semiMarkov}) or can be used to specify particular initial values.}
\usage{
param.init(data = NULL, cov = NULL, states, mtrans, 
           cov_tra = NULL, cens = NULL, dist_init = NULL, 
           proba_init=NULL, coef_init = NULL)
}
\arguments{
  \item{data}{
data frame of the form \code{data.frame(id,state.h,state.j,time)}, where
\itemize{
 \item \code{id}: the individual identification number
 \item \code{state.h}: state left by the process
 \item \code{state.j}: state entered by the process
 \item \code{time}: waiting time in state \code{state.h}
}
The data.frame containts one row per transition (possibly several rows per patient). The data frame \code{data} is not needed if \code{proba_init} is provided.
}
  \item{cov}{
Optional data frame containing the covariates values.
}
  \item{states}{
A numeric vector giving the names of the states (names are values used in \code{state.h}).
}
  \item{mtrans}{
A quadratic matrix of characters describing the possible transitions and the distributions of waiting time. 
The rows represent the left states, and the columns represent the entered states.   
 If an instantaneous transition is not allowed from state \code{h} to state \code{j}, then \code{mtrans} should have \eqn{(h,j)} entry \code{FALSE}, otherwise it should be \code{"E"} (or \code{"Exp"} or \code{"Exponential"}) for Exponential distribution, \code{"W"} (or \code{"Weibull"}) for Weibull distribution or \code{"EW"} (or \code{"EWeibull"} or \code{"Exponentiated Weibull"}) for Exponentiated Weibull distribution. If \code{TRUE} is used instead of the name of the distribution, then a Weibull distribution is considered. By definition of a semi-Markov model, the transitions into the same state are not possible. The diagonal elements of \code{mtrans} must be set to \code{FALSE} otherwise the function will stop. 
}
  \item{cov_tra}{
Optional list of vectors: a vector is associated with covariates included in the model. For a given covariate, the vector contains the transitions \code{"hj"} for which the covariate have an effect (only the transitions not equal to \code{FALSE} in \code{mtrans} are allowed). The effect of covariates can then be considered only on specific transitions. By default, the effects of covariates on all the possible transitions are included in a model.
}
\item{cens}{
   A character giving the code for censored observations in the column \code{state.j} of the data. Default is \code{NULL} which means that the censoring is defined as a transition fron state \eqn{h} to state \eqn{h}.
  }
  \item{dist_init}{
  Optional numeric vector giving the initial values of the distribution parameters. Default is 1 for each distribution parameter. The length of the vector depends on the chosen distribution, number of transitions and states.
  }
 \item{proba_init}{
  Optional numeric vector giving the initial values of the transition probabilities. The sum of the probabilities in the same raw must be equal to 1. According to semi-Markov model, the probability to stay in the same state must be equal to 0. The default values for the transition probabilities are estimated from the data. If \code{data = NULL}, the argument \code{proba_init} is obligatory.
  }
  \item{coef_init}{
  Optional numeric vector giving the initial values of the regression coefficients associated with the covariates. Default is 0 for each regression coefficient meaning no effect of the covariate.
  }
}

\details{
This function returns a data frame containing the initial values of parameters of a semi-Markov model. 
The model parameters are the distribution parameters, the transition probabilities of the Markov chain and the regression coefficients associated with covariates. The number of parameters depends on the chosen model: the distributions of the sojourn times, the number of states and the transitions between states specified with the matrix \code{mtrans}, the number of covariates (\code{cov}) and their effects or not on transitions (\code{cov_tra}).

The default initial values for the distribution parameters are fixed to 1. As the three possible distributions are nested for respective parameters equal to 1 (See details of the \code{semiMarkov} function), the initial distribution corresponds to the exponential distribution with parameter equal to 1 (whatever the chosen distribution). The default initial values for the regression coefficients are fixed to 0 meaning that the covariates have no effect on the hazard rates of the sojourn times. These initial values may be changed using the arguments \code{dist_init} and \code{coef_init}.

By default, the initial probabilities are calculated by simple proportions. The probability associated to the transition from \eqn{h} to \eqn{j} is estimed by the number of observed transitions from state \eqn{h} to state \eqn{j} divided by the total number of transitions from state \code{h} observed in the data. The results are displayed in matrix \code{matrix.P}. The number of parameters for transition probabilities is smaller than the number of possible transitions as the probabilities in the same row sum up to one. Considering this point and that the probability to stay in the same state is zero, the user can change the initial values using the argument \code{proba_init}. 
}

\value{
This function returns an object of class \code{param.init} to be used in functions \link{semiMarkov} and \link{hazard}. An object of class \code{param.init} consists of
 \item{nstates}{
The length of vector \code{states} interpreted as the number of possible states for the process.
}
 \item{table.state}{
A table, with starting states as rows and arrival states as columns, which provides the number of times that a transition between two states is observed. This argument is only returned when \code{data} is provided. It can be used to quickly summarize multi-state data.
}
\item{Ncens}{
Number of individuals subjected to censoring.
}
  \item{matrix.P}{
Quadratic matrix, with starting states as rows and arrival states as columns, giving the default initial values for the transition propabilities of the embedded Markov chain. All diagonal values are zero. The sum of all probabilities of the same row is equal to one.
}
 \item{last}{
The largest duration observed in \code{data} if \code{data} is given.
}
 \item{Transition_matrix}{
A matrix containing the informations on the model definition : the possible transitions o and the distribution of waiting times for each transition (Exponential, Weibull or Exponentiated Weibull).
} 
 \item{dist.init}{
A data frame giving the names of the parameters, transitions associated with and initial values of the distribution parameters.
}
 \item{proba.init}{
 A data frame giving names of the parameters, transitions associated with and initial values of the probabilities of the embedded Markov chain.
}
 \item{coef.init}{
 A data frame giving the names of covariates, transitions associated with and initial values of the regression coefficients.
}
}
\author{
Agnieszka Listwon-Krol
}
\seealso{
\link{hazard}, \link{semiMarkov}
}
\examples{

## Asthma control data
data(asthma)

## Definition of the model:  states, names, possible transtions and waiting time distributions
states_1 <- c("1","2","3")
mtrans_1 <- matrix(FALSE, nrow = 3, ncol = 3)
mtrans_1[1, 2:3] <- c("W","W")
mtrans_1[2, c(1,3)] <- c("EW","EW")
mtrans_1[3, c(1,2)] <- c("W","W")

## Default initial values in a model without covariates
init_1 <- param.init(data = asthma, states = states_1, mtrans = mtrans_1)

## Definition of initial values in a model without covariates
init_2 <- param.init(data = asthma, states = states_1, mtrans = mtrans_1,
          dist_init=c(rep(1.5,6),rep(1.8,6),rep(2,2)),
          proba_init=c(0.2,0.8,0.3,0.7,0.35,0.65))

## Default initial values with a covariate "Sex" influencing transitions " 1->2" and "3->2"
init_3 <- param.init(data = asthma, cov=as.data.frame(asthma$Sex),
          states = states_1, mtrans = mtrans_1, cov_tra=list(c("12","32")))

## Definition of initial values with a covariate "Sex" influencing transitions " 1->2" and "3->2"
init_4 <- param.init(data = asthma, cov=as.data.frame(asthma$Sex),
          states = states_1, mtrans = mtrans_1, cov_tra=list(c("12","32")),
          dist_init=c(rep(1.5,6),rep(1.8,6),rep(2,2)),
          proba_init=c(0.2,0.8,0.3,0.7,0.35,0.65), coef_init=rep(0.3,2))
          
init_5 <- param.init(data = asthma, cov=as.data.frame(asthma$Sex),
          states = states_1, mtrans = mtrans_1, cov_tra=list(c("12","32")),
          coef_init=c(0.2,0.5))          

## Definition of initial values without dataset in an illness-death model
## 1 - healthy, 2 - illness, 3 - death
states_2 <- c("1","2","3")
mtrans_2 <- matrix(FALSE, nrow = 3, ncol = 3)
mtrans_2[1,c(2,3)] <- c("E","W")
mtrans_2[2,c(1,3)] <- c("EW","EW")
init_6<-param.init(states=states_2, mtrans=mtrans_2, proba_init=c(0.7,0.3,0.2,0.8))

}



\keyword{documentation}
