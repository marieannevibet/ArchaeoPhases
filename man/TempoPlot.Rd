\name{TempoPlot}
\alias{TempoPlot}
\title{Plot of the occurence of events}
\usage{
TempoPlot(data, position, level = 0.95 , count = TRUE, Gauss = FALSE, 
title = "Tempo plot", x.label="Calendar Year", y.label="Cumulative events", 
line.types=c("solid", "12", "11", "28", "28"), plot.wi = 7, plot.ht = 7, 
base.font = 11, colors=TRUE, out.file=NULL)
}
\arguments{
   \item{data}{dataframe containing the output of the MCMC algorithm. The MCMC samples should be in calendar year (BC/AD). }
  \item{position}{numeric vector containing the position of the column corresponding to the MCMC chains of interest}
  \item{level}{probability corresponding to the level of confidence used for the credible interval}
  \item{count}{if TRUE the counting process is given as a number, otherwise it is a probability}
  \item{Gauss}{if TRUE, the Gaussian approximation of the CI is used}
  \item{title}{title of the graph}
  \item{x.label}{label of the x-axis}
  \item{y.label}{label of the y-axis}
  \item{line.types}{type of the lines drawn of the graph}
  \item{plot.wi}{width of the graph}
  \item{plot.ht}{height of the graph}
  \item{base.font}{font of the graph}
  \item{colors}{if TRUE, the graph is drawn with colors, otherwise it is drawn in black and white}
  \item{out.file}{the name of the graph + extension that will be saved if chosen. Null by default.}
}
\value{
It calculates the cumulative frequency of specified events by calculating how many events took place before each date in a specified range of dates. The result is given in calendar year (in format BC/AD). 
}
\description{
A statistical graphic designed for the archaeological study of rhythms of the long term that embodies a theory of archaeological evidence for the occurrence of events.
}
\details{
The tempo plot is one way to measure change over time: it estimates the cumulative occurrence of archaeological events in a Bayesian calibration.
The tempo plot yields a graphic where the slope of the plot directly reflects the pace of change: a period of rapid change yields a steep slope and a period of slow change yields a gentle slope. When there is no change, the plot is horizontal. When change is instantaneous, the plot is vertical.
}
\author{
  Anne Philippe <Anne.Philippe@univ-nantes.fr>, Thomas S. Dye  <TSD@tsdye.com> and

  Marie-Anne Vibet <Marie-Anne.Vibet@univ-nantes.fr>
}
\references{
Dye, T.S. (2016) Long-term rhythms in the development of Hawaiian social stratification. Journal of Archaeological Science, 71, 1--9.

}
\examples{
  data(Events); 
  TempoPlot(Events[1:1000,], c(2:5))
  TempoPlot(Events[1:1000,], c(2:5), count = TRUE)
}
\keyword{ tempo plot}
\keyword{ credible interval }
\keyword{ mean }

