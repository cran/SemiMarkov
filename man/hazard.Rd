\name{hazard}
\alias{hazard}
%\alias{hazard.semiMarkov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes hazard rates using an object of class \code{semiMarkov} or \code{param.init}
}
\description{
For a given vector of times, the function computes the hazard rates values of an object of class \code{semiMarkov} or \code{param.init} (which provided the hazard rates). Both, values of hazard rate of waiting time of semi-Markov process can be obtained.
}
\usage{
hazard(object, type = "alpha", time = NULL,
      cov = NULL, s = 0, t = "last", Length = 1000)
}

\arguments{
  \item{object}{
Object of class \code{semiMarkov} or \code{param.init}. 
}
  \item{type}{
Type of hazard to be computed: \code{"alpha"} for the hazard rates of waiting times of the embedded Markov chain and
\code{"lambda"} for the hazard rates of the semi-Markov process. Default is \code{"alpha"}. 
}
  \item{time}{
A vector containing the time values for which the hazard rate is computed.
Default value is a vector \code{seq(0, last, length = Length)} where \code{last} is the largest duration observed in the data set and \code{Length} is the length of the vector.
}
  \item{cov}{
  A list with one component for each covariate. Each component gives values of covariates that are to be used for the hazard rates computation. For a time-fixed covariate a single value can be given whereas a whole vector of values is required for time dependent covariates. Default is \code{NULL} which corresponds to time-independent covariates all equal to 0. Note that the same covariates values are used for all transitions.
}
   \item{s}{
Starting value of the time interval \eqn{[s, t]} which is used to compute the hazard rate. This argument is not considered when the vector \code{time} is defined. Default value is 0.
}
  \item{t}{
Ending value of the time interval \eqn{[s, t]} which is used to compute the hazard rate. This argument is not considered when the vector \code{time} is defined. Default value is \code{last} which is the the largest duration observed in the data set.
}
  \item{Length}{
The number of points of the time interval \eqn{[s, t]} for which the hazard rate is computed. These points are equally distributed in the time interval \eqn{[s, t]}. This argument is not considered when the vector \code{time} is defined. Default value is 1000.
}
}
\details{

This function computes the hazard rates of waiting (or sojourn) times and the hazard rates of semi-Markov process defined in the parametric multi-state semi-Markov model described in Listwon and Saint-Pierre (2013). Additional details about the methodology behind the \code{SemiMarkov} package can be found in Limnios and Oprisan (2001), Foucher \emph{et al.} (2006) and Perez-Ocon and Ruiz-Castro (1999).
 
%Let \eqn{(Z_{t}=J_{N(t)}}, \eqn{N(t)=\textrm{sup}\{m\in \mathbb{N}:S_{m}\leq t\}} be an homogeneous semi-Markov process with finite state space \eqn{E}
%where \eqn{\tau_{0}=S_{0}=0,\ldots,\tau_{m}=S_{1}+\cdots+S_{m},\tau_{\infty}=\textrm{sup}\{\tau_{m},m\in \mathbb{N}\}} are the successive  times of entry 
%to states \eqn{J_{0},J_{1},\ldots,J_{m}} (with \eqn{J_{m-1}\neq J_{m}, m=1,2,\ldots}) and \eqn{S_{m}} is the waiting (sojourn) time in state \eqn{J_{m}}.

%Let \eqn{(J_n,S_n)_{n \geq 0}} an homogeneous semi-Markov process with finite state space \eqn{E} where \eqn{S_0=0<S_1<...<S_n} are the successive time of entry in states \eqn{J_0,J_1,...,J_n} (with \eqn{J_{n-1}\neq J_n}). Let \eqn{X_0=0} et \eqn{X_n=S_n-S_{n-1}} the waiting time (gap time) in state \eqn{J_{n-1}}.

The hazard rate of waiting time at time \eqn{t} represents the conditional probability that a transition from state \eqn{h} to state \eqn{j} is observed given that no event occurs until time \eqn{t}. In a parametric framework, the expression of the hazard rates can easily be obtained as the distributions of waiting time belong to a parametric family. The hazard rate values are calculated using the chosen distribution and the given values of the parameters. The effects of both constant and time-varying covariates on the hazard of waiting time can be studied using a proportional intensities model. The effects of covariates can then be interpreted in terms of relative risk.
 
%The hazard rates of waiting time are defined by
%\deqn{\alpha_{hj}(d)=\lim_{\Delta d \rightarrow 0}\frac{P(d<S_{m+1}\leq d+\Delta d | S_{m+1}>d,J_{m}=h,J_{m+1}=j)}{\Delta d}.}
%where \eqn{d} denotes the time and \eqn{h,j \in E}.

The hazard rate of the semi-Markov process at time \eqn{t} represents the conditional probability that a transition into state \eqn{j} is observed given that the subject is in state \eqn{h} and that no event occurs until time \eqn{t}. The hazard rate of the semi-Markov process can be interpreted as the subject's risk of passing from state \eqn{h} to state \eqn{j}. This quantity can be deduced from the transition probabilities of the embedded Markov chain and from the distributions of waiting times.

%The hazard rates of the semi-Markov process represent the individual's risk of passing from state \eqn{h} to state \eqn{j},
%\deqn{\lambda_{hj}(d)=\lim_{\Delta d \rightarrow 0}\frac{P(J_{m+1}=j,d<S_{m+1}\leq d+\Delta d | S_{m+1}>d,J_m=i)}{\Delta d}.}

This function can be used to compute the hazard rates for different values of the covariates or different values of the parameters. These hazard rates can then be plotted using \link{plot.hazard}.

Objects of classes \code{semiMarkov} and \code{param.init} can be used in the function \code{hazard}. These objects contain informations on the model and the values of the parameters for the waiting time distribution, the transition probability of Markov chain and the regression coefficients.

%See the PDF manual ?msm-manual.pdf? in the ?doc? subdirectory for hazard function definitions of all distributions.
%that includes the estimated values of parameters. The covariates can be incorporated to their definitions by the Cox model of proportional baseline risks. Thus, the covariates effect is interpreted in terms of relative risk (RR) by the regression coefficients.
%The estimation of the waiting time is thus given by the relation between the estimated density function \strong{f} and estimated survival function \strong{S}
%\eqn{ \hat{\alpha}_{hj}(d)=\frac{\hat{f}_{hj}}{\hat{S}_{hj}},
%}
% where \strong{d>0} denotes time and \strong{h, j} the states. 

%process hazard rates are the parameters desired to be estimated. They represent the individual's risk of passing from state \strong{h} to state \strong{j}.
% Given the estimated Markov chain transition probabilities \strong{p} this hazard rate is 
 %\eqn{\hat{\lambda}_{hj}(d)=\frac{\hat{p}_{hj}\hat{S}_{hj}(d)\hat{\alpha}_{hj}(d)}{\sum\limits_{j=1}^{s}\hat{p}_{hj}\hat{S}_{hj}(d)}.}  
 
% These hazard rates are difficult in interpretation for the regression coefficients. Thus, it is recommended to trace calculated vectors, using \link{plot.hazard}.


}
\value{
%The waiting time hazard rates and the semi-Markov process hazard rates are calculated similarly as in the function \link{hazard}. The function \code{hazard.init} returns object of class \code{hazard} which can be entered in the function \code{print.hazard}, \code{summary.hazard} and \code{plot.hazard}.

%(one data frame per transition)

Returns an object of class \code{hazard}.
%consisting of a list of data frames. The list items are hazard rate values, vector of time values, vector of covariates values, distribution, type of hazard as well as basic statistics of calculated vector.
  \item{Type}{
  The type of hazard computed by the function \code{hazard}: the hazard of waiting time (\code{alpha}) or the hazard of the semi-Markov process (\code{lambda}).
  }
  \item{vector}{
A data frame containing one vector for each possible transition. A vector contains values of the hazard rate associated to the vector of times.
}
  \item{Time}{
The vector of times used to compute the hazard rate.
}
  \item{Covariates}{
   A list containing the values of the covariates (fixed or time-dependent).
  }
  \item{Summary}{
  A list of data frames (one for each possible transition). A dataframe contains quantiles, means, minimums and maximums of the hazard rate values.
  }
   \item{Transition_matrix}{
A matrix containing informations on the model: the possible transitions and the distribution of waiting times for each transition (Exponential, Weibull or Exponentiated Weibull).
} 
  \item{call}{
  Recall the name of the model.
  }
}

%For the object of class \code{semiMarkov} for each transition the vectors of hazard rates are calculated. They are all of length chosen within the argument \code{length} and the function returns a list of data frames for each objects where the column are the respective transitions.



\references{

Krol, A., Saint-Pierre P. (2015). \emph{SemiMarkov : An R Package for Parametric Estimation in Multi-State Semi-Markov Models}. 66(6), 1-16.

Limnios, N., Oprisan, G. (2001). \emph{Semi-Markov processes and reliability}. Statistics for Industry and Technology. Birkhauser Boston.

Foucher, Y., Mathieu, E., Saint-Pierre, P., Durand, J.F., Daures, J.P. (2006). A semi-Markov
model based on Generalized Weibull distribution with an illustration for HIV disease.
\emph{Biometrical Journal}, 47(6), 825-833.

Perez-Ocon, R., Ruiz-Castro, J. E. (1999). Semi-markov models and applications,
chapter 14, pages 229-238. \emph{Kluwer Academic Publishers}.

%Ouhbi, B. and Limnios, N. (2005). Non-parametric estimation for semi-Markov processes based on k-sample
%paths with application to reliability. \emph{Proceedings of the 11th Symposium on Applied
%Stochastic Models and Data Analysis (ASMDA'05)} J. Janssen and P. Lenca, Eds., pp. 1061-1068,
%Brest, France.

%Hougaard, P. (2000). \emph{Analysis of Multivariate Survival Data}. New York: Springer - Verlag.

%Gill, R.D., Keiding, N., Andersen, P.K., Borgan, O. (1993). \emph{Statistical Models Based on Counting Processes}. New York: Springer-Verlag.
}


\author{
Agnieszka Listwon-Krol
}

\seealso{
\link{plot.hazard}, \link{semiMarkov}, \link{param.init}, \link{summary.hazard}, \link{print.hazard}
}
\examples{

\donttest{

## Asthma control data
data(asthma)

## Definition of the model:  states, names, 
# possible transtions and waiting times distributions
states_1 <- c("1","2","3")
mtrans_1 <- matrix(FALSE, nrow = 3, ncol = 3)
mtrans_1[1, 2:3] <- c("E","E")
mtrans_1[2, c(1,3)] <- c("E","E")
mtrans_1[3, c(1,2)] <- c("W","E")

## semi-Markov model without covariates
fit1 <- semiMarkov(data = asthma, states = states_1, mtrans = mtrans_1)

## Hazard rates of waiting time
alpha1 <- hazard(fit1)
plot(alpha1)

## Hazard rates of the semi-Markov process
lambda1 <- hazard(fit1, type = "lambda")
plot(lambda1)

## Defining a vector of equally distributed times 
alpha2 <- hazard(fit1, s=0, t=3, Length=300)
plot(alpha2)

## Considering times observed in the data set
alpha3 <- hazard(fit1, time=sort(unique(asthma$time)))
plot(alpha3)

## semi-Markov model with a covariate "BMI"
fit2 <- semiMarkov(data = asthma, cov = as.data.frame(asthma$BMI), 
        states = states_1, mtrans = mtrans_1)

## Time fixed covariate
## Covariate equal to 0 and 1 for each transition
alpha4 <- hazard(fit2)
alpha5 <- hazard(fit2, cov=1)
plot(alpha4,alpha5)

## Time dependent covariate 
## Suppose that the covariate value is known for all times values
Time<-sort(unique(asthma$time))             # observed times in ascending order
Cov1<-sort(rbinom(length(Time), 1, 0.3))    # simulation of binary covariate
Cov2<-sort(rexp(length(Time), 5))           # simulation of numeric covariate
alpha6 <- hazard(fit2, time=Time, cov=Cov1)
plot(alpha6)
alpha7 <- hazard(fit2, time=Time, cov=Cov2)
plot(alpha7)

## semi-Markov model with two covariates 
## "BMI" affects transitions "1->3" and "3->1"
## "Sex" affects transition "3->1"
SEX <- as.data.frame(asthma$Sex)
BMI <- as.data.frame(asthma$BMI)
fit3 <- semiMarkov(data = asthma, cov = as.data.frame(cbind(BMI,SEX)),
                   states = states_1, mtrans = mtrans_1,
                   cov_tra = list(c("13","31"),c("31")))
alpha8 <- hazard(fit3, cov=c(0,0))
alpha9 <- hazard(fit3, cov=c(1,1))
plot(alpha8,alpha9)

}
}
 \keyword{documentation}
  
