\name{TempoActivityPlot}
\alias{TempoActivityPlot}
\title{Plot of the activity of events}
\usage{
TempoActivityPlot(data, position, plot.result = NULL, level = 0.95,
                       count = TRUE, Gauss = FALSE, title = "Tempo plot",
                       subtitle = NULL, caption = "ArcheoPhases",
                       legend.title = "Legend",
                       legend.labels = c("Bayes estimate",
                                         "Credible interval, low",
                                         "Credible interval, high",
                                         "Gaussian approx., high",
                                         "Gaussian approx., low"),
                       x.label = "Calendar year",
                       y.label = "Cumulative events",
                       line.types = c("solid", "12", "11", "28", "28"),
                       width = 7, height = 7, units = "in",
                       x.min = NULL, x.max = NULL, colors = TRUE,
                       file = NULL, x.scale = "calendar",
                       elapsed.origin.position = NULL, appliShiny=FALSE, print.data.result = TRUE)
}
\arguments{
  \item{data}{dataframe containing the output of the MCMC algorithm. The MCMC samples should be in calendar year (BC/AD). }
  \item{position}{numeric vector containing the position of the column corresponding to the MCMC chains of interest}
  \item{plot.result}{a list containing the data to plot, typically the result of a previous run of TempoPlot()}
  \item{level}{probability corresponding to the level of confidence used for the credible interval}
  \item{count}{if TRUE the counting process is given as a number, otherwise it is a probability}
  \item{Gauss}{if TRUE, the Gaussian approximation of the CI is used}
  \item{title}{title of the graph}
  \item{subtitle}{subtitle of the graph}
  \item{caption}{caption of the graph}
  \item{legend}{the title of the legend}
  \item{legend.labels}{a vector of strings to label legend entries}  
  \item{x.label}{label of the x-axis}
  \item{y.label}{label of the y-axis}
  \item{line.types}{type of the lines drawn of the graph in the order of legend.labels}
  \item{width}{width of the plot in units}
  \item{height}{height of the plot in units}
  \item{units}{units used to specify width and height.  One of "in", "cm", or "mm".  Default = "in".}
  \item{x.min}{minimum value for x axis}
  \item{y.max}{maximum value for x axis}
  \item{colors}{if TRUE, the graph is drawn with colors, otherwise it is drawn in black and white}
  \item{file}{the name of the graph (+ extension) that will be saved if chosen. Null by default.}
  \item{x.scale}{one of "calendar", "bp", or "elapsed"}
  \item{elapsed.origin.position}{if x.scale is "elapsed", the position of the column corresponding to the occurrence from which elapsed time is calculated}
  \item{appliShiny}{whether the plot is drawn with the shiny app or not}
  \item{print.data.result}{If TRUE, the list containing the data to plot will be given. Default = TRUE.}
}
\value{
It is the derivative of the TempoPlot bayesian estimate. 
It may also return a list containing the data to plot (if print.data.result = TRUE). The result is given in calendar year (in format BC/AD).
}
\description{
A statistical graphic designed for the archaeological study of rhythms of the long term that embodies a theory of archaeological evidence for the occurrence of events.
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
  TempoActivityPlot(Events[1:1000,], c(2:5))
  TempoActivityPlot(Events[1:1000,], c(2:5), count = TRUE)
}
\keyword{ Tempo Activity plot}
\keyword{ mean }

