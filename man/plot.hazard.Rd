\name{plot.hazard}
\alias{plot.hazard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for objects of class \code{hazard}
}
\description{
Plot method for one or several (maximum 10) objects of class \code{hazard}. Depending on the hazard rate chosen in the function \code{hazard}, the function plots either the hazard rates of sojourn times or the semi-Markov process hazard rate for each considered transition (one plot for each transition).
}
\usage{
\method{plot}{hazard}(x, x2 = NULL, x3 = NULL, x4 = NULL, x5 = NULL, x6 = NULL, x7 = NULL, 
      x8 = NULL, x9 = NULL, x10 = NULL, transitions = NULL, names = NULL, 
      legend = TRUE, legend.pos = NULL, cex = NULL, colors = NULL,
	  xlab = "Time", ylab = "Hazard function", lwd = 3, type = "p", ...)
}
\arguments{
  \item{x}{
Object of class \code{hazard}. At least one \code{hazard} object is needed.
}
  \item{x2}{
Object of class \code{hazard}. Default is \code{NULL}.
}
  \item{x3}{
Object of class \code{hazard}. Default is \code{NULL}.
}
  \item{x4}{
Object of class \code{hazard}. Default is \code{NULL}.
}
  \item{x5}{
Object of class \code{hazard}. Default is \code{NULL}.
}
  \item{x6}{
Object of class \code{hazard}. Default is \code{NULL}.
}
  \item{x7}{
Object of class \code{hazard}. Default is \code{NULL}.
}
  \item{x8}{
Object of class \code{hazard}. Default is \code{NULL}.
}
  \item{x9}{
Object of class \code{hazard}. Default is \code{NULL}.
}
  \item{x10}{
Object of class \code{hazard}. A maximum of ten \code{hazard} objects is possible.  Default is \code{NULL}.
}

  \item{transitions}{
A character vector giving the transitions to be plotted. Default is \code{NULL} which means that all the possible transitions are displayed.
}
  \item{names}{
Names of the hazard rates. Default is \code{NULL} which means that the names used in the \code{semiMarkov} object are applied. 
}
  \item{legend}{
A logical value specifying if a legend should be added. Default is \code{TRUE}.
}
  \item{legend.pos}{
A vector giving the legend position.
}
  \item{cex}{
character expansion factor relative to current \code{par("cex")}.
}
  \item{colors}{
A vector of colours for the hazard rates.
}
  \item{xlab}{
x-axis label. Default is \code{Time}. 
}
  \item{ylab}{
y-axis label. Default is \code{Hazard function}.
}
  \item{lwd}{
Thickness of lines or points.
}
  \item{type}{
Type of graph. Default are points \code{p}.
}
  \item{\dots}{
Further arguments for plot.
}
}
\value{
No value returned.
}

\author{
Agnieszka Listwon-Krol
}

\seealso{
\link{hazard}, \link{semiMarkov}
}
\examples{
\dontrun{


## Asthma control data
data(asthma)

## Definition of the model:  states, names, possible transtions 
# and waiting times distributions
states_1 <- c("1","2","3")
mtrans_1 <- matrix(FALSE, nrow = 3, ncol = 3)
mtrans_1[1, 2:3] <- c("E","E")
mtrans_1[2, c(1,3)] <- c("E","E")
mtrans_1[3, c(1,2)] <- c("W","E")
fit <- semiMarkov(data = asthma, states = states_1, mtrans = mtrans_1)
lambda<-hazard (fit, type = "lambda")

plot(lambda, names = c("lambda"),legend=FALSE)
plot(lambda, transitions = c("13","31"), names = c("lambda"),
legend.pos=c(2,0.09,2,0.4))

## semi-Markov model in each stratum of Severity
fit0 <- semiMarkov(data = asthma[asthma$Severity==0,],
        states = states_1, mtrans = mtrans_1)
fit1 <- semiMarkov(data = asthma[asthma$Severity==1,],
        states = states_1, mtrans = mtrans_1)
lambda0<-hazard (fit0, type = "lambda",s=0,t=5,Length=1000)
lambda1<-hazard (fit1, type = "lambda",s=0,t=5,Length=1000)
plot(lambda0,lambda1, names = c("lambda0", "lambda1"),
legend.pos=c(4,0.18,4,0.8,4,0.2,4,0.09,4,0.7,4,0.21))

## semi-Markov model with covariate "BMI"
fitcov <- semiMarkov(data = asthma, cov = as.data.frame(asthma$BMI),
        states = states_1, mtrans = mtrans_1)
lambda0<-hazard (fitcov, type = "lambda",cov = c(0))
lambda1<-hazard (fitcov, type = "lambda",cov = c(1))
plot(lambda0,lambda1, names = c("lambda0", "lambda1"))

}
}
\keyword{ documentation }

