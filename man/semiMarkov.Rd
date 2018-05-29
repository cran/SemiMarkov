\name{semiMarkov}
\alias{semiMarkov}
\alias{SemiMarkov}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric estimation in multi-state semi-Markov models
}
\description{
This function computes the parametric maximum likelihood estimation in multi-state semi-Markov models in continuous-time. The effect of time varying or fixed covariates can be studied using a proportional intensities model for the hazard of the sojourn time.
%Covariates can be fitted on the base of Cox proportional risks models and they can be either time varying or fixed.
%Fits a continuous-time semi-Markov multi-state model by maximum likelihood using the optimizing function \code{nloptr}. 
}

\usage{
semiMarkov(data, cov = NULL, states, mtrans, cov_tra = NULL, 
          cens = NULL, dist_init=NULL, proba_init = NULL, coef_init = NULL, 
          alpha_ci = 0.95, Wald_par = NULL, eqfun = NULL, 
          ineqLB = NULL,ineqUB = NULL, control = list() )
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
This data.frame containts one row per transition (possibly several rows per patient).
%obligatory format of \code{data} is the long format i.e., with one row per individual and per transition. 
}

  \item{cov}{
Optional data frame containing covariates values.
}

  \item{states}{
A numeric vector giving names of the states (names are values used in \code{state.h}).% Default values are deduced from the data.
}

  \item{mtrans}{
A quadratic matrix of characters describing the possible transitions and the distributions of waiting time. 
The rows represent the left states, and the columns represent the entered states.   
 If an instantaneous transition is not allowed from state \code{h} to state \code{j}, then \code{mtrans} should have \eqn{(h,j)} entry \code{FALSE}, otherwise it should be \code{"E"} (or \code{"Exp"} or \code{"Exponential"}) for Exponential distribution, \code{"W"} (or \code{"Weibull"}) for Weibull distribution or \code{"EW"} (or \code{"EWeibull"} or \code{"Exponentiated Weibull"}) for Exponentiated Weibull distribution. If \code{TRUE} is used instead of the name of the distribution, then a Weibull distribution is considered. By definition of a semi-Markov model, the transitions into the same state are not possible. The diagonal elements of \code{mtrans} must be set to \code{FALSE} otherwise the function will stop. 
}

  \item{cov_tra}{
Optional list of vectors: a vector is associated with covariate included in the model. For a given covariate, the vector contains the transitions \code{"hj"} for which the covariate have an effect (only the transitions specified in \code{mtrans} are allowed). The effect of covariates can then be considered only on specific transitions. By default, the effects of covariates on all the possible transitions are studied.
}

\item{cens}{
   A character giving the code for censored observations in the column \code{state.j} of the data. Default is \code{NULL} which means that the censoring is defined as a transition fron state \eqn{h} to state \eqn{h}.
  }
  
  \item{dist_init}{
  Optional numeric vector giving the initial values of the distribution parameters. Default is 1 for each distribution parameter. The length of the vector depend on the chosen distribution, number of transitions and states.
  }
  
 \item{proba_init}{
  Optional numeric vector giving the initial values of the transition probabilities. The sum of the probabilities in the same row must be equal to 1. According to semi-Markov model, the probability to stay in the same state must be equal to 0. The default values for the transition probabilities are estimated from the data. If \code{data = NULL}, the argument \code{proba_init} is obligatory.
  }
     
  \item{coef_init}{
  Optional numeric vector giving the initial values of the regression coefficients associated with the covariates. Default is 0 for each regression coefficient which means that the covariate has no effect.
  }
    
%  \item{opt.method}{
%Quoted name of the R function used to perform the likelihood optimization. The functions \code{nloptr} from R package \code{nloptr} and \code{constrOptim.nl} from R package \code{alabama} can be considered. The default is the \code{constrOptim.nl} function.
%}
  
%  \item{iter}{
%An optional value specifying the limit on iterations executed by the optimization method. %\code{nloptr}
%}

  \item{alpha_ci}{
Confidence level to be considered for the confidence intervals. The default value is 0.95. 
}

  \item{Wald_par}{
  Optional numeric vector giving the values to be tested (null hypothesis) by the Wald test for each parameter. The Wald statistics are evaluated only for the parameters 
  of distributions and regression coefficients. The length of this vector must then be equal to the number of those parameters. The order of the values must be as in the 
  parameters table given by objects \code{semiMarkov} or \code{param.init} (excluding the parameters associated to the transition probabilities). The default values for 
  the elements of \code{Wald_par} vector are 1 for the distribution parameters and 0 for the regression coefficients.
  %Optional character of theoretical values that estimated parameters are equal to in the null hypothesis in Wald test. As the Wald statistics are
  %computed for the parameters of distributions and regression coefficients the length of the vector must equal to the number of those parameters. 
  %The order of the values is to be as in the table of parameters of the objects \code{semiMarkov} or \code{param.init} (excluding the probabilities 
  %for transitions ).
  }
  
  \item{eqfun}{
     Optional list given equality constraints between parameters. These constraints are passed using the equality constraint function that can be defined in the \code{solnp} optimization function. See below for details.
   %  Optional list given equality constraint function that are to be passed to the \code{solnp} optimization function. See below for details.  
  }
 
  \item{ineqLB}{
     Optional list given values of lower bound for parameters. These values are used in the inequality constraint that can be defined in the \code{solnp} optimization function. See below for details.         
  }
  
  \item{ineqUB}{ 
     Optional list given values of upper bound for parameters. These values are used in the inequality constraint that can be defined in the \code{solnp} optimization function. See below for details.         
     %Optional list given values for the upper bound of the inequality constraint that are to be passed to the \code{solnp} optimization function. See below for details.              
  }
  
  \item{control}{
   The control list of optimization parameters for \code{solnp} optimization function.
}
}

\details{
This function fits parametric multi-state semi-Markov model described in Listwon and Saint-Pierre (2013) to longitudinal data. Additional details about the methodology behind the \code{SemiMarkov} package can be found in Limnios and Oprisan (2001), Foucher \emph{et al.} (2006) and Perez-Ocon and Ruiz-Castro (1999).
  

%Let \eqn{Z_{t}=J_{N(t)}}, \eqn{N(t)=\textrm{sup}\{m\in \mathbb{N}:S_{m}\leq t\}} be an homogeneous semi-Markov process with finite state space \eqn{E} where \eqn{\tau_{0}=S_{0}=0,\ldots,\tau_{m}=S_{1}+\cdots+S_{m},\tau_{\infty}=\textrm{sup}\{\tau_{m},m\in \mathbb{N}\}} are the successive  times of entry to states \eqn{J_{0},J_{1},\ldots,J_{m}} (with \eqn{J_{m-1}\neq J_{m}, m=1,2,\ldots}) and \eqn{S_{m}} is the waiting (sojourn) time in state \eqn{J_{m}}. In a parametric semi-Markov model, the distribution of the waiting time belongs to a parametric family.
%Let \eqn{(J_n,S_n)_{n \geq 0}} an homogeneous semi-Markov process with finite state space \eqn{E} where \eqn{S_0=0<S_1<...<S_n} are the successive time of entry in states \eqn{J_0,J_1,...,J_n} (with \eqn{J_{n-1}\neq J_n}). \eqn{(J_n)_{n \geq 0}} is the embedded Markov chain. Let \eqn{X_0=0} et \eqn{X_n=S_n-S_{n-1}} the sojourn time (gap time) in state \eqn{J_{n-1}}. In a parametric semi-Markov model, the distribution of the waiting time belongs to a parametric family.
 
Consider an homogeneous semi-Markov process with a finite state space. In a parametric framework, distributions of the waiting time belong to parametric families. The distribution of the waiting time can be chosen between the exponential, the Weibull and the exponentiated Weibull distributions. The exponential distribution with scale parameter \eqn{\sigma>0} has a density defined as follows
\deqn{f(x) =(1/\sigma) exp(-x/\sigma).}

The Weibull distribution with scale parameter \eqn{\sigma>0} and shape parameter \eqn{\nu>0} has a density given by (same as one defined in \code{dweibull})
\deqn{g(x) = (\nu/\sigma)(x/\sigma)^{\nu-1} exp(-(x/\sigma)^\nu).}
                                                      
The exponentiated Weibull distribution (or generalized Weibull) with scale parameter \eqn{\sigma>0}, shape parameter \eqn{\nu>0} and family parameter equal to \eqn{\theta>0} has a density given by (same as one defined in function \code{dgweibull} from the R package \code{rmutil})
\deqn{h(x) = \theta(\nu/\sigma)(x/\sigma)^{\nu-1} exp(-(x/\sigma)^\nu) (1-exp(-(x/\sigma)^\nu))^{\theta-1}.}

%has a density given by
%\deqn{f_E(x) = \frac{1}{\sigma} \exp\left(- \frac{x}{\sigma}\right).}
%The Weibull distribution with scale parameter \eqn{\sigma} and shape parameter \eqn{\nu} density has a density given by (same as \code{dweibull})
%\deqn{f_W(x) = \left(\frac{\nu}{\sigma}\right)\left(\frac{x}{\sigma}\right)^{\nu-1} exp\left(- \left(\frac{x}{\sigma}\right)^{\nu}\right).}
%The exponentiated Weibull distribution (also called generalized Weibull) with scale parameter \eqn{\sigma}, shape parameter \eqn{\nu} and family parameter equal to \eqn{\theta} density has a density given by (same definition as in function \code{dgweibull} from the R package \code{rmutil})
%\deqn{f_{EW}(x) = \theta\left(\frac{\nu}{\sigma}\right)\left(\frac{x}{\sigma}\right)^{\nu-1} \exp\left(- \left(\frac{x}{\sigma}\right)^{\nu}\right)
%\left(1-\exp \left(-\frac{x}{\sigma}\right)^{\nu}\right)^{\theta}.}

These three distributions are nested. The exponentiated Weibull density with \eqn{\theta=1} gives a Weibull distribution and the Weibull density with \eqn{\nu=1} gives the exponential density. 

Note that the effects of both constant and time-varying covariates on the hazards of sojourn time can be studied using a proportional intensities model. The effects of covariates can then be interpreted in terms of relative risk.

The model parameters are the distribution parameters, the transition probabilities of the Markov chain and the regression coefficients associated with covariates. The number of parameters depends on the chosen model: the distributions of the sojourn times, the number of states and the transitions between states specified with the matrix \code{mtrans}, the number of covariates (\code{cov}) and their effects or not on transitions (\code{cov_tra}).
  %%

The default initial values for the distribution parameters are fixed to 1. As the three possible distributions are nested for parameters equal to 1 (See details of the \code{semiMarkov} function), the initial distribution corresponds to an exponential with parameter equal to 1 (whatever the chosen distribution). The default initial values for the regression coefficients are fixed to 0 meaning that the covariates have no effect on the hazard rates of the sojourn times. These initial values may be changed using the arguments \code{dist_init} and \code{coef_init}.

By default, the initial probabilities are calculated by simple proportions. The probability associated to the transition from \eqn{h} to \eqn{j} is estimed by the number of observed transitions from state \eqn{h} to state \eqn{j} divided by the total number of transitions from state \code{h} observed in the data. The results are displayed in matrix \code{matrix.P}. The number of parameters for transition probabilities is smaller than the number of possible transitions as the probabilities in the same row sum up to one. Considering this point and that the probability to stay in the same state is zero, the user can changed the initial values using the argument \code{proba_init}.

The Yinyu Ye optimization solver to nonlinear problem is applied to maximize the log-likelihood using the function \code{solnp} created by A. Ghalanos and S. Theussl. In order to modify the optimization parameters refer to the package \code{Rsolnp} documentation.

Some optimization difficulties can arrise when there is not enough information in the data to estimate each transition rate. The users can change the optimization parameters and the initial values. It may be appropriate to reduce the number of states in the model, the number of allowed transitions, or the number of covariate effects, to ensure convergence.

Some additionals constraints can be introduced using \code{eqfun}, \code{ineqLB} and \code{ineqUB}. These constraints on distribution parameters, transition probabilities 
and regression coefficients can be defined using lists of vectors. 
The argument \code{eqfun} gives the possibility to add constraints of type \eqn{par1 = a*par2} (\eqn{a} is a constant). This equality constraint must be expressed with a vector 
of 3 elements where the first element is the identifier of the parameters type (\code{"dist"} for distribution parameters, \code{"proba"} for the transition probabilities and 
\code{"coef"} for the regression coefficients), the second and the third elements are the index of \eqn{par1} and \eqn{par2}, respectively. The index values of 
distribution parameters, transition probabilities and regression coefficients can be found in the table provided by an object \code{semiMarkov}. The last element of the vector 
corresponds to the constant \eqn{a}. The arguments \code{ineqLB} and \code{ineqUB} allow to add constraints of type \eqn{par \geq a} and \eqn{par \leq a}, respectively. 
These arguments are lists of vectors of length 3 where the first element is the type of the parameter (\code{"dist"}, \code{"proba"} or \code{"coef"}), the second element 
is the index of parameter \eqn{par} and the last one is the constant \eqn{a}.
If a chosen constraint corresponds to a transition probability, it should be considered that the last probabilities in a row of the transition matrix are not estimated but 
obtained directly since the sum of transition probabilities in the same row is equal to 1. Thus, no additional constraints related to these parameters are permitted. Moreover, 
note that the argument \code{eqfun} does not allowed to define relationships between parameters of different types (for instance, a transition probability can not be equal to 
a regression coefficient). The optional constraints on parameters should be used prudently since they may induce problems in the convergence of the optimization method. In particular, 
the Wald statistic and the standard deviation may not be computed for some parameters due to negative values in the hessian matrix. Note that the default constraints induce by the model 
definition are treated in priority.

%It is possible to adjust the constraints to the needs of user using \code{eqfun}, \code{ineqLB} and \code{ineqUB} that should be lists of 
%vectors.  The constraints are possible for the distribution parameters, probabilities and regression coefficients.
%If a chosen constraint corresponds to a probability parameter, it should be considered that the last probabilities in a line from the transition matrix are not estimated but
%  calculated directly after the estimation to maintain the assumption of the sum of probabilities in one line equals to 1. Thus, no additional 
%  constraints related to these parameters are permitted.
%The option \code{eqfun} gives the possibility to add constraints of type \eqn{par1 = a*par2} (\eqn{a} - constant) and that should be expressed with a vector of 3 elements
%where the first element is the identifier of the parameters type (\code{"dist"} for distribution parameters, \code{"proba"} for the probability and \code{"coef"} for the regression
%coefficients), the second and the third elements are the positions of \eqn{par1} and \eqn{par2}, respectively (the information on the positions of distribution parametersand regression coefficients
%can be found in the respective tables of the object \code{param.init} and in the table of the object \code{semiMarkov} and the positions of probabilities can be found in the table of the object \code{semiMarkov}). The last element of the vector corresponds to the constant \eqn{a}. 
%The options \code{ineqLB} and \code{ineqUB} create the possibility to add the constraints of type \eqn{par \geq a} and    \eqn{par \leq a}, respectively.
%These are lists of vectors of length 3 where the first element is the identifier of the parameter  type (\code{"dist"}, \code{"proba"} or \code{"coef"}), the second ele;ent
%is the position of the \eqn{par} and the last one is the constant \eqn{a}. 
%These options
%should be used prudently as they may impose problems for the convergence or reliability of a model. When using the additional constraints the Wald statistic may not be computed for some parameters.
%What is more, the default constraints for the parameters are treated
%with the priority. Finally, in \code{eqfun} no relationships between parameters of different types are allowed. 

%the optimization method,

%A METTRE DANS L'ARTICLE
 %- When an exponential distribution is chosen, the semi-Markov model corresponds to the usual Markov model.
 %- The effect of covariates are indirectly reflected on the hazard functions of the semi-Markov process but this effect can not be interpreted in terms of relative risk (only graphically). 
 %- The parameters of the distribution can then be tested to evaluate the relevance of the chosen distribution (for instance, using an exponentiated weibull distribution, one can test \eqn{H_0:\theta=1}).


%For full details about the methodology behind the semiMarkov package, refer to the PDF manual ?semiMarkov-manual.pdf? in the ?doc? subdirectory of the package. This includes a tutorial in the typical use of \code{semiMarkov}.

%\code{semiMarkov} was created to fit semi-Markov models for which transitions between the states from the finite state space can occur at any time. It is possible that the model is the subject to independent right-censoring.

%To encompass covariates to the model, the Cox risk proportionality is assumed. The additional assumption is that covariates act on the waiting time distributions. Indirectly, their effects are reflected on the hazard functions of the semi-Markov process but they are only interpretable graphically (P. Saint-Pierre \emph{et al.}, 2005).

%If the arguments \code{par_init} and \code{cov_init} are not given, then the likelihood is calculated in terms of initial transition intensity matrix where the distribution parameters are set to 1 and covariates to 0. Matrix of initial values may be displayed using function \link{param.init}. When Exponential distribution is chosen, then the classical Markov model is considered. The initial parameters may be chosen by the utilisator, to interpret their influence on the model it is recommended to test them with functions \link{param.init} and \link{hazard.init}.

%There must be enough information in the data on each state to estimate each transition rate, otherwise the likelihood will be flat and the maximum will not be found.

%The option to choose between the Exponential, Weibull and Exponentiated distributions was designed for their advantage of being nested. Thus, the Likelihood Ratio Statistics
%can be used to evaluate the relevance of a larger number of parameters. In the results \code{sigma} is the scale parameter that occurs in all three distributions, 
%\code{nu} is the shape parameter found in Weibull in Exponentiated Weibull distributions and \code{theta} - the shape parameter only for the Exponentiated Weibull. 

}
\value{
 \item{call}{
The original call to \code{semiMarkov}.
}
 \item{minus2loglik}{
Minus twice the maximized log-likelihood.
}
 \item{solution}{
Etimations of the distribution parameters, probabilities of the embedded Markov chain and regression coefficients (if any considered) for each transition specified in the matrix \code{mtrans}. This is a data frame with three columns: the label of the parameter, the transition associated with and the estimated value.
}
 \item{opt.message}{
The message giving the information on the optimization result returned by the \code{constrOptim.nl} function. %or the \code{nloptr}
}
 \item{opt.iter}{
Number of outer iterations of the optimization method.
}
 \item{nstates}{
The length of vector \code{states} interpreted as the number of possible states for the process.
}
 \item{table.state}{
A table, with starting states as rows and arrival states as columns, which provides the number of observed transitions between two states. This argument can be used to quickly summarize multi-state data.
}
\item{Ncens}{
Number of individuals subjected to censoring.
}
 \item{Transition_matrix}{
A matrix containing the informations on the model definition : the possible transitions and the distribution of waiting times for each transition (Exponential, Weibull or Exponentiated Weibull).
} 
 \item{param.init}{
Recall the initial values of the parameters. The third column of this object can be used in \code{hazard} function.
}
 \item{table.dist}{
Statistics for the estimations of distribution parameters of waiting time distributions. For the exponential distribution one data frame for the parameter \code{sigma} is returned, for the Weibull distribution two data frames for \code{sigma} and \code{nu} are returned, and for the Exponentiated Weibull distribution three data frames for \code{sigma}, \code{nu} and \code{theta} are returned. The columns of each data frame are the possible transitions, the estimations, the standard deviations, the lower and upper bounds of confidence intervals, the Wald test null hypothesis, the Wald test statistics and the p-values of the Wald test when testing hypothesis \code{sigma=1}, \code{nu=1} or \code{theta=1}.
}
 \item{table.proba}{ 
 A data frame giving the estimations of the transition probabilities of the Markov chain and their standard deviations. By definition, the probability associated to the last possible transition of each row of the matrix \code{mtrans} is equal to \eqn{1-pr}, where \eqn{pr} is the sum of all other probabilities from the row.
}  
% To simplify the optimization for the last possible transition of each row of the matrix \code{mtrans},
 \item{table.coef}{
If some covariates are included in the model it returns a data frame with the statistics for the estimated values of the regression coefficients. The columns of the data frame are the transitions associated with the coefficients, the estimations, the standard deviations, the lower and upper bounds of confidence intervals, the Wald test null hypothesis, the Wald test statistics and the p-values of the Wald test when testing hypothesis \code{coef=0}.
}
 \item{table.param}{
Data frame with the statistics for all the model parameters, that is \code{table.dist}, \code{table.proba} and \code{table.coef} in a single data frame.
 }
%Printing a \code{semiMarkov} object by typing the object's name at the command line implicitly invokes \code{print.semiMarkov}. This formats and prints the important information in the model fit. 
%This includes the estimated values, standard deviations, confidence intervals and p-value of the Wald test for all the parameters and the minus twice the maximised log-likelihood. 
}
\note{
Printing a \code{semiMarkov} object by typing the object's name at the command line implicitly invokes \code{print.semiMarkov}.
% By definition of a semi-Markov model, the transitions into a same state are not allowed. If transitions into the same state are detected in the data, the function will stop. The diagonal elements of \code{mtrans} must be set to \code{FALSE}.
}

\references{
  Krol, A., Saint-Pierre P. (2015). \emph{SemiMarkov : An R Package for Parametric Estimation in Multi-State Semi-Markov Models}. 66(6), 1-16.

Ghalanos, A. and Theussl S. (2012). Rsolnp: General Non-linear Optimization Using Augmented Lagrange Multiplier Method. R package version 1.14.

Limnios, N., Oprisan, G. (2001). \emph{Semi-Markov processes and reliability}. Statistics for Industry and Technology. Birkhauser Boston.

Foucher, Y., Mathieu, E., Saint-Pierre, P., Durand, J.F., Daures, J.P. (2006). A semi-Markov
model based on Generalized Weibull distribution with an illustration for HIV disease.
\emph{Biometrical Journal}, 47(6), 825-833.

Perez-Ocon, R., Ruiz-Castro, J. E. (1999). Semi-markov models and applications,
chapter 14, pages 229-238. \emph{Kluwer Academic Publishers}.
}
\author{
Agnieszka Listwon-Krol, Philippe Saint-Pierre
}

\seealso{
\link{param.init}, \link{hazard}, \link{summary.semiMarkov}, \link{print.semiMarkov}
}
\examples{

\donttest{

## Asthma control data
data(asthma)

## Definition of the model:  states, names, possible transtions and waiting time 
## distributions
states_1 <- c("1","2","3")
mtrans_1 <- matrix(FALSE, nrow = 3, ncol = 3)
mtrans_1[1, 2:3] <- c("E","E")
mtrans_1[2, c(1,3)] <- c("E","E")
mtrans_1[3, c(1,2)] <- c("W","E")

## semi-Markov model without covariates
fit1 <- semiMarkov(data = asthma, states = states_1, mtrans = mtrans_1)

## semi-Markov model with one covariate 
## "BMI" affects all transitions
BMI <- as.data.frame(asthma$BMI)
fit2 <- semiMarkov(data = asthma, cov = BMI, states = states_1, mtrans = mtrans_1)
## semi-Markov model with one covariate 
## "BMI" affects the transitions "1->3" and "3->1"
fit3 <- semiMarkov(data = asthma, cov = BMI, states = states_1, mtrans = mtrans_1,
                   cov_tra = list(c("13","31")))

## semi-Markov model with two covariates 
## "BMI" affects the transitions "1->3" and "3->1"
## "Sex" affects the transition "3->1"
SEX <- as.data.frame(asthma$Sex)
fit4 <- semiMarkov(data = asthma, cov = as.data.frame(cbind(BMI,SEX)),
                   states = states_1, mtrans = mtrans_1,
                   cov_tra = list(c("13","31"),c("31")))
                   
## semi-Markov model using specific initial values      
## same model as "fit1" but using different initial values
## "fit5" and "fit6" are equivalent

init <- param.init(data = asthma, states = states_1, mtrans = mtrans_1,
        dist_init=c(rep(1.5,6),c(1.8)), proba_init=c(0.2,0.8,0.3,0.7,0.35,0.65))
fit5 <- semiMarkov(data = asthma, states = states_1, mtrans = mtrans_1,
                   dist_init=init$dist.init[,3], proba_init=init$proba.init[,3])
fit6 <- semiMarkov(data = asthma, states = states_1, mtrans = mtrans_1,
                   dist_init=c(rep(1.5,6),c(1.8)), 
                   proba_init=c(0.2,0.8,0.3,0.7,0.35,0.65))                   


## The Wald test null hypothesis is modified
## Wald statistics when testing nullity of distribution parameters
## and regression coefficients equal to -1
fit7 <- semiMarkov(data = asthma, cov = BMI, states = states_1, mtrans = mtrans_1,
                   Wald_par = c(rep(0,7),rep(-1,6)))
            
## semi-Markov model with additional constraints 
## distribution parameters sigma for transition "1->3" = sigma for transition "2->1" 
fit8 <- semiMarkov(data = asthma, cov = BMI, states = states_1, mtrans = mtrans_1,
                   eqfun = list(c("dist",2,3,1)))

## semi-Markov model with additional constraints 
## regression coefficients beta for transition "1->2" = beta for transition "2->1" 
##                                                    = beta for transition "2->3"  
fit9 <- semiMarkov(data = asthma, cov = BMI, states = states_1, mtrans = mtrans_1,
                   eqfun = list(c("coef",1,3,1),c("coef",1,4,1)))

## semi-Markov model with additional constraints 
## regression coeficient beta for transition "1->2" belongs to [-0.2,0.2]
## and regression coeficient beta for transition "2->3" belongs to [-0.05,0.05]
fit10 <- semiMarkov(data = asthma, cov = BMI, states = states_1, mtrans = mtrans_1,
                    ineqLB = list(c("coef",1,-0.2),c("coef",4,-0.05)),
                    ineqUB = list(c("coef",1,0.2),c("coef",4,0.05)))

}

}
 \keyword{documentation}
