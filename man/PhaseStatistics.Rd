% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prog.R
\name{PhaseStatistics}
\alias{PhaseStatistics}
\title{Summary statistics of a phase}
\usage{
PhaseStatistics(PhaseMin_chain, PhaseMax_chain, level = 0.95)
}
\arguments{
\item{PhaseMin_chain}{numeric vector containing the output of the MCMC algorithm for the minimum of the dates included in the phase. The MCMC samples should be in calendar year (BC/AD).}
\item{PhaseMax_chain}{numeric vector containing the output of the MCMC algorithm for the maximum of the dates included in the phase. The MCMC samples should be in calendar year (BC/AD).}
\item{level}{probability corresponding to the level of confidence used for the credible interval and the highest density region}
}
\value{
  Returns a list of values corresponding to the summary statistics:
  \item{1}{Statistics of the minimum of the dates included in the phase}
  \item{2}{Statistics of the maximum of the dates included in the phase}
  \item{3}{Statistics of the duration of the dates included in the phase}
  The results are given in calendar year (in format BC/AD).
}
\description{
Estimation of several summary statistics of the minimum, the maximum and the duration of the dates included in the phase. 
}
\details{
The summary statistics are those given by MarginalStatistics function. The time range is given by PhaseTimeRange function.
The duration is computed as follow duration = maximum - minimum at each iteration of the MCMC output.
}
\author{
  Anne Philippe <Anne.Philippe@univ-nantes.fr> and

  Marie-Anne Vibet <Marie-Anne.Vibet@univ-nantes.fr>
}
\examples{
  data(Phases); attach(Phases)
  PhaseStatistics(Phase.1.alpha, Phase.1.beta, 0.95)
  PhaseStatistics(Phase.2.alpha, Phase.2.beta, 0.95)
}

\keyword{ mean }
\keyword{ maximum a posteriori }
\keyword{ credible interval }
\keyword{ highest posterior density }
\keyword{ phase time range }
\keyword{ summary statistics}
\keyword{ individual phase}