\name{MultiDatesPlot}
\alias{MultiDatesPlot}
\title{ Plot of the endpoints of credible intervals or HPD intervals of a series of dates}
\usage{
MultiDatesPlot(data, position, level = 0.95,  intervals = "CI", 
  title = "Plot of intervals", labelXaxis = "Calendar Year", 
  exportFile = NULL, exportFormat = "PNG")
}
\arguments{
  \item{data}{dataframe containing the output of the MCMC algorithm. The MCMC samples should be in calendar year (BC/AD). }
  \item{position}{numeric vector containing the position of the column corresponding to the MCMC chains of interest}
  \item{level}{probability corresponding to the level of confidence used to estimate the credible interval}
  \item{intervals}{"CI" corresponds to the credible intervals, "HPD" to the highest density regions}
  \item{title}{title of the graph}
  \item{labelXaxis}{label of the horizontal axis}
  \item{exportFile}{the name of the file to be saved. If NULL then no graph is saved}
  \item{exportFormat}{the format of the export file : SVG or PNG.}
}
\description{
Draws a plot of segments corresponding to the endpoints of the intervals (CI or HPD) of each selected date. The result is given in calendar year (in format BC/AD).
}
\author{
  Anne Philippe <Anne.Philippe@univ-nantes.fr> and

  Marie-Anne Vibet <Marie-Anne.Vibet@univ-nantes.fr>
}
\examples{
  data(Events)
  MultiDatesPlot(Events, c(2,4,3), level = 0.95, intervals ="CI", title = "Plot of CI intervals")
  MultiDatesPlot(Events, c(2,4,3), level = 0.95, intervals ="HPD", title = "Plot of HPD intervals")
  MultiDatesPlot(Events, c(2,4,3), level = 0.95, intervals ="HPD", exportFile = "MultiDatesPlot", 
  exportFormat = "PNG")
}
\keyword{ credible interval }
\keyword{ highest posterior density regions }
