\name{table.state}
\alias{table.state}
\title{
Table giving the numbers of observed transitions
}
\description{
Function returning a table with numbers of transitions between two states observed in the data set. This table can be a used to summarize a multi-state data or to define the matrix \code{mtrans} required in the \code{semiMarkov} function.
}
\usage{
table.state( data, states = NULL, mtrans = NULL, cens = NULL)
}
\arguments{
  \item{data}{
data frame in form \code{data.frame(id,state.h,state.j,time)}, where
\itemize{
 \item \code{id}: the individual identification number
 \item \code{state.h}: state left by the process
 \item \code{state.j}: state entered by the process
 \item \code{time}: waiting time in state \code{state.h}
}
This data.frame containts one row per transition (possibly several rows per patient).
%obligatory format of \code{data} is the long format i.e., with one row per individual and per transition. 
}

  \item{states}{
A numeric vector giving the names of the states (names are values used in \code{state.h}).% Default values are deduced from the data.
}

  \item{mtrans}{
A quadratic matrix of logical values describing the possible transitions. 
The rows represent the left states, and the columns represent the entered states.   
 If an instantaneous transition is not allowed from state \code{h} to state \code{j}, then \code{mtrans} should have \eqn{(h,j)} entry \code{FALSE},
  otherwise it should be \code{TRUE}. Default value is a matrix which allows all the possible transitions between states. %Instead of entry \code{TRUE}, the name of the distribution can be given.
}

\item{cens}{
   A character giving the code for censored observations in the column \code{state.j} of the data. Default is \code{NULL} which means that the censoring is defined as a transition fron state \eqn{i} to state \eqn{i} (by definition of a semi-Markov model, the transitions into the same state are not possible).
  }
}

\value{

\item{table.state}{
A table, with starting states as rows and arrival states as columns, which provides the number of observed transitions between two states. This argument can be used to quickly summarize multi-state data.
}
\item{Ncens}{
Number of individuals subjected to censoring.
}
}

\author{
Agnieszka Listwon
}
\seealso{
\link{param.init}, \link{semiMarkov}
}

   \examples{

## Asthma control data
data(asthma)

# default description
# censoring is implicitly defined as a transition "h->h"
table.state(asthma)
table.state(asthma)$Ncens

# censoring defined as a transition to state "4"
asthma_bis<-asthma
for(i in 1:dim(asthma)[1]){if(asthma[i,2]==asthma[i,3]) asthma_bis[i,3]<-4}
table.state (asthma_bis, cens = 4)

## Definition of the model: states names and possible transtions
states_1 <- c("1","2","3")
mtrans_1 <- matrix(FALSE, nrow = 3, ncol = 3)
mtrans_1[1, 2:3] <- TRUE
mtrans_1[2, c(1,3)] <- c("W","E")
table.state(asthma, states = states_1, mtrans = mtrans_1)

   }

\keyword{documentation}
