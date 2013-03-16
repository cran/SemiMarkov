\name{asthma}
\alias{asthma}
\docType{ data }
\title{Asthma control data}
\description{

Data from a follow-up study of severe asthmatic patients. At each visit, covariates are recorded and asthma was evaluated using the concept of control scores. 
Such scores reflect a global judgement of the disease gravity based on official criteria. Three levels are considered (optimal, suboptimal and unacceptable control) 
and can be used to define the subject's state at each visit. The aim is to investigate the evolution of asthma control and to evaluate the effect of covariates.
The data contains an extraction of 371 patients with at least two visits. The table is presented in long format with one row for each observed transition between two states. 
The rows corresponding to the same subject are ordered chronologically. The last sojourn time is right-censored by the end of the study and represent the time until censoring. A censored transition is defined as a transition to the same state \code{h->h}. 
}

\usage{data(asthma)}
\format{
A data frame containing 876 rows. 
%Rows are grouped by patient number and ordered by the identification number. 
Each row represents a patient examination and contains several covariates.
\describe{
     \item{\code{id}}{Patient identification number}
     \item{\code{state.h}}{Starting state (1 for optimal, 2 for suboptimal and 3 for unacceptable control state)}
     \item{\code{state.j}}{Arrival state (1 for optimal, 2 for suboptimal and 3 for unacceptable control state)}
     \item{\code{time}}{Waiting (sojourn) time in state \code{state.h}}
     \item{\code{Severity}}{Disease severity (1=severe, 0=mild-moderate asthma)}
     \item{\code{BMI}}{Body Mass Index (1=BMI>=25, 0=otherwise)}
     \item{\code{Sex}}{Sex (1=men, 0=women)}
%\item{\code{CorticoInhDaily}}{the covariate for patient's daily dose of inhaled corticosteroids (1=dose>500\eqn{\mu g}, 0=otherwise)}
%\item{\code{CorticoOralDaily}}{the covariate for patient's daily dose of oral corticosteroids (1=treatment with oral corticostereoids, 0=otherwise)}
%\item{\code{CorticoOralCumul}}{the covariate for patient's cumulated dose of oral corticosteroids (1=dose>2g, 0=otherwise)}
}
This presentation of the data implies that, for a given patient, the visited states are the sequence of \code{state.h} 
and the follow-up time is the the cumulated sum of \code{time}.
}

  \source{
ARIA (Association pour la Recherche en Intelligence Artificielle), France.
}
\references{
Saint-Pierre P., Combescure C., Daures J.P., Godard P. (2003). The analysis of
asthma control under a Markov assumption with use of covariates. \emph{Statistics in Medicine},
22(24):3755-70.
}
\examples{
data(asthma)
head(asthma)
}
\keyword{ datasets }

