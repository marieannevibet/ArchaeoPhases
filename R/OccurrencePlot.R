
###################################
###       Occurrence plot ## NEW 2017/08      ###


# Occurrence plot 
# A statistical graphic designed for the archaeological study of rhythms of the long term that embodies a theory of archaeological evidence for the occurrence of events
#' @param data dataframe containing the output of the MCMC algorithm
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence
#' @param intervals one of "CI" for credible intervals or "HPD" for highest posterior density intervals
#' @param title title of the graph
#' @param subtitle subtitle of the graph
#' @param caption caption of the graph
#' @param labelXaxis label of the x-axis
#' @param labelYaxis label of the y-axis
#' @param language a string indicating a language recognized by the toOrdinal library
#' @param occurrence a string to append to each y axis tic label
#' @param width width size in units
#' @param height height size in units
#' @param units a string recognized by the ggsave function, one of "in", "cm", "mm"
#' @param x.min minimum x axis value
#' @param x.max maximum x axis value
#' @param x.scale one of "calendar" for calendar years, "BP" for years before present, or "elapsed" for time elapsed from a specified origin
#' @param elapsed.origin.position position of the column to use as the origin for elapsed time calculations
#' @param dumbbell.size size of the plot symbol
#' @param dot.guide switch for a horizontal guide from the y axis
#' @param dot.guide.size size of the dot guide
#' @param y.grid switch for horizontal grid lines
#' @param count, if TRUE the counting process is a number, otherwise it is a probability
#' @param file the name of the graph (+ extension) that will be saved if chosen. Null by default.
#' @param appliShiny whether the plot is drawn with the shiny app or not
#' @param print.data.result If TRUE, the list containing the data to plot will be given
#' @return a plot
#' @export

OccurrencePlot <- function(data, position, level = 0.95, intervals = "CI",
                           title = "Occurrence plot",
                           subtitle = NULL,
                           caption = "ArchaeoPhases",
                           labelXaxis = "Calendar year",
                           labelYaxis = NULL,
                           language = "English", occurrence = "occurrence",
                           height = 7, width = 7, units = "in",
                           x.min = NULL, x.max = NULL, x.scale = "calendar",
                           elapsed.origin.position = NULL,
                           dumbbell.size = 3, dot.guide = FALSE,
                           dot.guide.size = 0.25, y.grid = FALSE,
                           count = TRUE, file = NULL, appliShiny=FALSE, print.data.result = TRUE)
{
  
  library(ggalt)
  library(toOrdinal)
  sort.rows <- function(x) {
    if (is.numeric(as.matrix(x))) {
      res <- as.data.frame(t(apply(x,1,sort)))
      colnames(res) <- 1:ncol(res)
      res
    }
    else {
      stop("Cannot sort non-numeric data frame")
    }
  }
  group.dates <- function(data, position) {
    L = length(position)
    res = matrix(ncol = L, nrow=nrow(data))
    for (i in 1:L) {
      res[,i] = data[,position[i]]
    }
    res
  }
  if (x.scale == "elapsed") {
    if (is.null(elapsed.origin.position)) {
      stop("Elapsed origin not specified")
    }
    else {
      data <- data - data[,elapsed.origin.position]
    }
  }
  groupOfDates <- group.dates(data, position)
  groupOfEvents <- sort.rows(groupOfDates)
  
  if (intervals == "CI") {
    Bornes = MultiCredibleInterval(groupOfEvents, 1:ncol(groupOfEvents), level = level)
    Ordered.df <- as.data.frame(Bornes)
    Ordered.df$y.labs <- paste(sapply(as.integer(rownames(Bornes)), toOrdinal, language = language), occurrence, sep= " ")
    
  }
  else if (intervals == "HPD") {
    Bornes = MultiHPD(groupOfEvents, 1:ncol(groupOfEvents), level = level)
    Ordered.df <- as.data.frame(Bornes)
    Ordered.df$y.labs <- paste(sapply(as.integer(rownames(Bornes)), toOrdinal, language = language), occurrence, sep= " ")
    
    # In the case of (multiple) two intervals
    x = (dim(Bornes)[2] - 1) /2 
    if ( x > 1 ){
      data = NULL
      for(j in 1:x){
        Bornesj = subset(Bornes,select=-c(2*j,2*j+1))
        data = rbind(Bornesj, data)
      }
      data = data[is.na(data[,2])==FALSE,]
      Ordered.df <- as.data.frame(data)
      Ordered.df$y.labs <- paste(sapply(as.integer(rownames(data)), toOrdinal, language = language), occurrence, sep= " ")
    }
    
  }
  
  if (x.scale == "BP") {
    Ordered.df[,2] <- 1950-Ordered.df[,2]
    Ordered.df[,3] <- 1950-Ordered.df[,3]
  }
  

  h <- ggplot(data = Ordered.df,
              aes(y=factor(y.labs, levels=unique(y.labs), ordered=TRUE),
                  x=Ordered.df[,2],
                  xend=Ordered.df[,3])) +
    geom_dumbbell(size = dumbbell.size, dot_guide = dot.guide,
                  dot_guide_size = dot.guide.size) +
    labs(x = labelXaxis, y = labelYaxis, title = title,
         subtitle = subtitle, caption = caption)
  if (!y.grid) {
    h <- h + theme(panel.grid.major.y=element_blank())
  }
  if (!is.null(x.min) & !is.null(x.max)) {
    h <- h + xlim(x.min, x.max)
  }
  if (!is.null(file)) {
    ggsave(filename = file, plot = h, height = height,
           width = width, units = units)
  }
  if(appliShiny == FALSE) {
    dev.new(height = height, width = width)
  }
  print(h)
  
  ## If the result is desired
  if (print.data.result == TRUE){
    Bornes
  }
}