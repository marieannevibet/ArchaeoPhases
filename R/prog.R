#####################################################
#        Run ArchaeoPhases shiny apps               #
#####################################################

#' Run ArchaeoPhases shiny apps
#' 
#' 
#' @export app_ArchaeoPhases
app_ArchaeoPhases <- function() {
  
  app <- shiny::runApp(system.file(paste0("shiny/"), package = "ArchaeoPhases"), launch.browser = TRUE)
  
}



#####################################################
#             Importing a CSV file                  #
#####################################################

#' Importing a CSV file 
#'
#' Importing a CSV file containing the output of the MCMC algorithm from any software
#'
#' @details 
#' @param file the name of the CSV file containing the output of the MCMC algorithm 
#' @param dec the character used in the file for decimal points for the use of read.csv()
#' @param sep the field separator character for the use of read.csv()
#' @param comment.char a character vector of length one containing a single character or an empty string for the use of read.csv()
#' @param header a character vector of length one containing a single character or an empty string for the use of read.csv()
#' @param iterationColumn the column number containing the iteration number. If specified, the column will be withdrawn. Default = NULL. 
#' @param referenceYear the year of reference of the date format 
#' @param rowToWithdraw the number of the row to be withdrawn. Default = NULL.
#' @return A data frame (data.frame) containing a representation of the data in the file.
#' @export
#'
ImportCSV <- function(file, dec='.', sep=',', comment.char = '#', header = TRUE, iterationColumn = NULL, referenceYear = NULL, rowToWithdraw=NULL){
  
  # importing the CSV file
  data = read.csv(file, dec = dec, sep=sep, comment.char = comment.char, header = header)
  
  # Withdrawing the iterations column
  if (is.null(iterationColumn)){
    return(data)
  } else {
    data = data[,-iterationColumn]
  }
  
  # Withdrawing a row
  if (is.null(rowToWithdraw)){
    return(data)
  } else {
    data = data[-rowToWithdraw, ]
  }
  
  # Conversion of the MCMC samples in date format cal BP or other to BC/AD
    if (is.null(referenceYear)){
      return(data)
    } else {
    data2 = data
    L = length(data)
    conv <- function(value, T0){
      T0 - value
    }
    for (i in 1:L){
      if( is.numeric(data[,i]) == TRUE){
        data2[,i] = sapply(data[,i], conv, referenceYear)
      }
    }
    return(data2)
  }

  
}




#####################################################
#        Create a mcmc.list for CODA users         #
#####################################################

#' Create a mcmc.list for CODA users
#' 
#' 
#' @export coda.mcmc

coda.mcmc <- function(data, numberChains = 1, iterationColumn = NULL){
  
  # Withdrawing the iteration column
  if (!is.null(iterationColumn)){
    data = data[,-iterationColumn]
    }
  
  dim =dim(data)
  L = dim[1]/numberChains 
  
  # select only numeric columns
  vect = NULL
  for(i in 1:dim[2]){
    if(is.numeric(data[,i])==TRUE) { vect = c(vect,i)}
  }
  data2 = data[,vect]
    
  obj <- list(NA); 
  for (i in 1:numberChains){
    obj[[i]] = mcmc(data2[ (L*(i-1)+1):(L*i),], start=1, end=L)
  }
  
  mcmcList = mcmc.list(obj)
  return(mcmcList)
}


#####################################################
#         Constructing the Phases min max            #
#####################################################

#' Constructing the minimum and the maximum for each phase  
#'
#' Constructing a dataframe containing the output of the MCMC algorithm corresponding to the minimum and the maximum of each group of events
#'
#' @details 
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of all dates included in the phase of interest
#' @param name name of the current group of dates or phase
#' @param add the name of the dataframe in which the current minimum and maximum should be added. Null by default. 
#' @param exportFile the name of the final file that will be saved if chosen. Null by default. 

#' @return A dataframe containing the minimum and the maximum of the group of dates included in the phase of interest. These values may be added to an already existing file "addFile" if given. 
#' @export
#'
CreateMinMaxGroup <- function(data, position, name ="Phase", add=NULL, exportFile=NULL){
  
  # importing the CSV file
  dataTemp = data[position]
  Min = apply(dataTemp, 1, min)
  Max = apply(dataTemp, 1, max)
  
  name.Min = paste(name,".Min", sep="")
  name.Max = paste(name,".Max", sep="")
  
  MinMaxCurrentPhase = cbind(Min,Max)
  colnames(MinMaxCurrentPhase) <- c(name.Min,name.Max)
  
  if (is.null(add)){
    MinMaxPhase = MinMaxCurrentPhase
  } else {
    MinMaxPhase = cbind(add, MinMaxCurrentPhase)
  }
  
  if (is.null(exportFile)){
    
  } else {
    write.csv(MinMaxPhase, exportFile, row.names=FALSE)
  }
  
  return(as.data.frame(MinMaxPhase))
  
}


#####################################################
#     Anteriority / posteriority Probability        #
#####################################################

#' Bayesian test for anteriority / posteriority
#'
#' Tests that a_chain < b_chain
#'
#' @param a_chain : numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param b_chain : numeric vector containing the output of the MCMC algorithm for the parameter b
#' @return The baysesian probability that a < b knowing the data
#' @export
MarginalProba <- function(a_chain,b_chain){

  mean(ifelse(a_chain < b_chain, 1, 0))   ## bayesion test : a < b

}


#####################################################
#     Estimation of Credible interval        #
#####################################################

#' Bayesian credible interval
#'
#' Estimation of the shorest credible interval of the output of the MCMC algorithm for the parameter a
#'
#' @details A 100*level % credible interval is an interval that keeps N*(1-level) elements of the sample outside the interval
#' The 100*level % credible interval is the shortest of all those intervals.
#' @param a_chain numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param level probability corresponding to the level of confidence used for the credible interval and the highest density region
#' @return The endpoints of the shortest credible interval
#' @export
#'
CredibleInterval <- function(a_chain, level=0.95){

  sorted_sample <- sort(a_chain)     # ordering the sample
  N = length(a_chain)                # calculation of the sample size
  OutSample = N * (1-level)          # calculation of the number of data to be outside the interval

  I =  cbind(sorted_sample[1:(OutSample+1)] , sorted_sample[(N-OutSample):N])    #   combinasion of all credible intervals

  l = I[,2]-I[,1]   # length of intervals
  i <- which.min(l) # look for the shortest interval

  c(level = level, CredibleIntervalInf=I[i,1],CredibleIntervalSup=I[i,2])   # returns the level and the endpoints

}




#####################################################
#            Marginal  Statistics                   #
#####################################################
#' Summary statistics
#'
#' Estimation of all usual statistics
#'
#' @param a_chain numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param level probability corresponding to the level of confidence used for the credible interval and the highest density region
#' @param title title of the summary statistics
#' @return A list of values corresponding to all the following statistics
#' @export
MarginalStatistics <- function(a_chain, level=0.95){

  # Position
  mean = round(mean(a_chain), 0)
  hdr = hdr(a_chain, prob = c(level * 100))
  map = round(hdr$mode, 0)
  quantiles = round(quantile(a_chain, c(0.25,0.5,0.75)), 0)

  # Dispersion
  sd = round(sd(a_chain), 0)            # standard deviation using the 'sd' function
  CI = c(round(CredibleInterval(a_chain, level)[2], 0), round(CredibleInterval(a_chain, level)[3], 0))           # Credible Interval using the function 'CredibleInterval' from the package 'Rchronomodel'
  HPDR = round(hdr$hdr, 0)              # Highest posterior density function region using the function 'hdr' from the package 'hdrcde'

  # Resulted
  res = c(mean, map, sd, quantiles[1], quantiles[2], quantiles[3], level, CI[1], CI[2], HPDR)
  Mat = matrix(nrow=length(res), ncol=1)
  Mat[,1] = res
  
  nom=c()
  for( k in (1: (length(HPDR)/2) ) ) {
    nom=c(nom,paste("HPDRInf",k))
    nom=c(nom,paste("HPRDSup",k))
  }
  names1 = c("mean", "MAP", "sd", "Q1", "median", "Q2", "level", "CredibleInterval Inf", "CredibleInterval Sup")
  rownames(Mat) = c(names1, nom)
  return(Mat)
}




#####################################################
#          Marginal posterior Density               #
#####################################################
#' Marginal posterior Density
#'
#' Plots the density of a_chain + statistics (mean, CI, HPDR)
#'
#' @param a_chain numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param level probability corresponding to the level of confidence used for the credible interval and the highest density region
#' @param title label of the title
#' @param colors if TRUE  -> use of colors in the graph
#' @param GridLength length of the grid used to estimate the density
#' @param exportFile the name of the file to be saved. If NULL then no graph is saved. 
#' @param exportFormat the format of the export file : PNG or SVG.
#' @return a plot with the density of a_chain + +CI + mean + HDR
#' @export
MarginalPlot <- function(a_chain, level=0.95, title="Characteristics of a date", colors=TRUE, exportFile = NULL, exportFormat = "PNG", GridLength=1024){

  maxValuex <- min(max(density(a_chain, n=GridLength)$x), 2016)
  minValuex <- min(density(a_chain, n=GridLength)$x)
  
  # Computing min and max values
  x = 10^c(0:10)
  if(minValuex!=0){  
    c =0
    for(i in 1:length(x)) { if( abs(minValuex/x[i])>1) {c=c+1}}
    if(c>3){ minValuex = floor(minValuex/x[c-1])*x[c-1]} else {minValuex = floor(minValuex/x[c])*x[c]}
  }
  if(maxValuex!=0){  
    d=0
    for(i in 1:length(x)) { if( abs(maxValuex/x[i])>1) {d=d+1}}
    if(d>3){ maxValuex = ceiling(maxValuex/x[d-1])*x[d-1]} else {maxValuex = ceiling(maxValuex/x[d])*x[d]}
  }
  
  # x-axis
  middleValuex <- minValuex + ( maxValuex - minValuex ) / 2
  P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
  P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
  
  # y-axis
  step <- max(density(a_chain, n=GridLength)$y) /50   # used to draw CI and mean above the curve
  maxValuey <- max(density(a_chain, n=GridLength)$y)
  middleValuey <- maxValuey /2

  if(!is.null(exportFile)) {
    
    if(exportFormat == "PNG") {
      png( filename = paste(exportFile,"png", sep =".") )
    } 
    if(exportFormat == "SVG") {
      svg( filename = paste(exportFile,"svg", sep =".") )
    }
    
  } 
    if (colors==T){
      par(mfrow=c(1,1))
      plot(density(a_chain, n=GridLength), main = title, xlab = "Calendar Year", axes = F, xlim=c(minValuex,maxValuex), ylim=c(0,max(density(a_chain, n=GridLength)$y) + step))
      # abscissa axis
      axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) )
      # ordinate axis
      axis(2, at=c(0, middleValuey , maxValuey),labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )
      
      segments(CredibleInterval(a_chain, level)[2], 0, CredibleInterval(a_chain, level)[3], 0, lwd=6, col = 4)
      points(mean(a_chain), 0 , lwd=6, col = 2)
      
      # legend
      legend(P3Valuex, maxValuey, c("Density", "Credible Interval", "Mean"), lty=c(1, 1, 0), bty="n", pch=c(NA, NA,1), col = c("black","blue","red"), lwd=c(1,6,6), x.intersp=0.5, cex=0.9)
      
    } else {
      par(mfrow=c(1,1))
      plot(density(a_chain, n=GridLength), main = title, xlab = "Calendar Year", axes = F, xlim=c(minValuex,maxValuex), ylim=c(0,max(density(a_chain, n=GridLength)$y) + step))
      # abscissa axis
      axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) )
      # ordinate axis
      axis(2, at=c(0, middleValuey , maxValuey),labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )
      
      segments(CredibleInterval(a_chain, level)[2], 0, CredibleInterval(a_chain, level)[3], 0, lwd=6, lty=1)
      points(mean(a_chain), 0 , lwd=6, pch=1)
      
      # legend
      legend(P3Valuex, maxValuey, c("Density", "Credible Interval", "Mean"), lty=c(1,1,0), pch=c(NA,NA,1), bty="n", lwd=c(1,6,6), x.intersp=0.5, cex=0.9)
    }
  if(!is.null(exportFile)) {

    dev.off()
  } 

}

#####################################################
#          Hiatus between two dates                 #
#####################################################
#'  Test of the hiatus between two dates
#' Finds if it exists a gap between two dates that is the longest interval that satisfies : P(a_chain < IntervalInf < IntervalSup < b_chain | M) = level
#'
#' @param a_chain : numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param b_chain : numeric vector containing the output of the MCMC algorithm for the parameter b
#' @param level probability corresponding to the level of confidence
#' @return The endpoints of the longest gap
#' @export
DatesHiatus <- function(a_chain, b_chain, level=0.95){

  if(length(a_chain) != length(b_chain)) {stop('Error : the parameters do not have the same length')} # test for the length of both chains
   
       gamma = mean((a_chain<b_chain))
       if (gamma < level) {print("No hiatus at this level")
       					   return(c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')) } else # 
       {
    
      interval <- function(epsilon, P1End, P2Beginning, level)
      {
        q1 = quantile(P1End ,probs = 1-epsilon) ;
        indz = (P1End < q1)
        q2 = quantile(P2Beginning[indz],probs= (1-level-epsilon)/(1-epsilon))
        c(q1,q2)
      }
      hia = Vectorize(interval,"epsilon")

      indz = which(a_chain<b_chain)
      epsilon = seq(0,1-level,gamma)
      p = hia(epsilon, a_chain[indz], b_chain[indz], level/gamma)
      rownames(p)<- c("HiatusIntervalInf", "HiatusIntervalSup")

      D<- p[2,]-p[1,]
      DD = D[D>0]

      if (length(DD) > 0){
        I = which(D==max(DD))
        interval2 = round( p[,I], 0)
        if (p[2,I] != p[1,I]) {
          c(level=level, interval2[1], interval2[2])
        } else {
          c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')
        }#end if (p[2,I] != p[1,I])

      } else {
        c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')
        }#end if (length(DD) > 0)

      } # end if( sum(ifelse(PhaseBeginning < PhaseEnd, 1, 0) == length(PhaseBeginning) ) ) {  # test for Beginning < End


}




#####################################################
#                 Phase Time Range                 #
#####################################################

#' Phase Time Range
#'
#' Computes the shortest interval that satisfies : P(PhaseMin_chain =< IntervalInf < IntervalSup =< PhaseMax_chain | M) = level
#'
#' @param PhaseMin_chain : numeric vector containing the output of the MCMC algorithm for the minimum of the events included in the phase
#' @param PhaseMax_chain : numeric vector containing the output of the MCMC algorithm for the maximum of the events included in the phase
#' @param level probability corresponding to the desired level of confidence
#' @return The endpoints of the shortest time range associated with the desired level
#' @export
PhaseTimeRange <- function(PhaseMin_chain, PhaseMax_chain, level=0.95){

  if(length(PhaseMax_chain) != length(PhaseMin_chain)) { print('Error : the parameters do not have the same length')}   # test the length of both chains
  else{

    if( sum(ifelse(PhaseMin_chain <= PhaseMax_chain, 1, 0)) != length(PhaseMin_chain) )  {  # test for Beginning < End
      print('Error : PhaseMin_chain should be older than PhaseMax_chain')
    } else {

      periode <- function(epsilon, PMin, PMax, level){
        q1 = quantile(PMin, probs = epsilon)    # Computes the 'level'th quantile of the minimum of the events included in the phase
        indz = (PMin > q1)
        q2 = quantile(PMax[indz], probs= (level/(1-epsilon)))
        c(q1,q2)
      }   # end periode <- function(epsilon, PMin, PMax, level){
      per = Vectorize(periode,"epsilon")

      epsilon = seq(0,1-level,.001)       # sequence of values used to compute
      p = per(epsilon,PhaseMin_chain, PhaseMax_chain, level)
      rownames(p)<- c("TimeRangeInf", "TimeRangeSup")

      D<- p[2,]-p[1,]     # computes the length of all intervals
      I = which.min(D)    # finds the shortest interval
      range = round(p[,I], 0)
      c(level=level, range[1], range[2]) # returns the endpoints of the shortest interval

    }
    # end if( sum(ifelse(PhaseMin_chain < PhaseMin_chain, 1, 0) == length(PhaseMin_chain) ) ) {  # test for Beginning < End

  }# end if(length(PhaseMin_chain) != length(PhaseMin_chain)) {

} # end PhaseTimeRange <- function(PhaseMin_chain, PhaseMax_chain, level, plot = F){





#####################################################
#             Statistics  for Phases               #
#####################################################

#' Summary statistics for phases
#'
#' Estimation of all usual statistics of the beginning and the end of a phase and the duration of the phase
#'
#' @param PhaseMin_chain numeric vector containing the output of the MCMC algorithm for the minimum of the dates included in the phase
#' @param PhaseMax_chain numeric vector containing the output of the MCMC algorithm for the maximum of the dates included in the phase
#' @param level probability corresponding to the level of confidence used for the credible interval and the highest density region
#' @return A matrix of values corresponding to all the summary statistics
#' @export
PhaseStatistics <- function(PhaseMin_chain, PhaseMax_chain, level=0.95){

  #Statistics according to PhaseMin_chain
  MinStat = MarginalStatistics(PhaseMin_chain, level)

  #Statistics according to PhaseMax_chain parameter
  MaxStat = MarginalStatistics(PhaseMax_chain, level)

  #Statistics according to the duration 
  DurationStat = MarginalStatistics(PhaseMax_chain-PhaseMin_chain, level)

  # Resulted List
  
  if (length(MinStat) > length(MaxStat)) {
    
    NbDiff = length(MinStat) - length(MaxStat)
    Add = rep(NA,NbDiff) 
    MaxStat = c(MaxStat, Add)

  }else if (length(MinStat) < length(MaxStat)) {
    
    NbDiff = length(MaxStat) - length(MinStat)
    Add = rep(NA,NbDiff) 
    MinStat = c(MinStat, Add)
    
  }
  Mat1 = cbind(MinStat, MaxStat)
  
  if (dim(Mat1)[1] > length(DurationStat)) {
    
    NbDiff = dim(Mat1)[1] - length(DurationStat)
    Add = rep(NA,NbDiff) 
    DurationStat = c(DurationStat, Add)
    
  }else if (dim(Mat1)[1] < length(DurationStat))  {
    
    NbDiff = length(DurationStat) - dim(Mat1)[1]
    Add = cbind( rep(NA,NbDiff), rep(NA,NbDiff)) 
    Mat1 = rbind(Mat1, Add)
    
  }
  
  Mat = cbind(Mat1, DurationStat) 
  colnames(Mat) = c("Minimum", "Maximum", "Duration")
  return(Mat)

}


#####################################################
#           Phase marginal density plot             #
#####################################################

#' Phase marginal density plot
#'
#' Plot of the density of the minimum and the maximum of the events included in the phase and summary statistics (mean, CI)
#'
#' @param PhaseMin_chain numeric vector containing the output of the MCMC algorithm for the minimum of the events included in the phase
#' @param PhaseMax_chain numeric vector containing the output of the MCMC algorithm for the maximum of the events included in the phase
#' @param level probability corresponding to the level of confidence used for the credible interval and the time range
#' @param title The Title of the graph
#' @param colors if TRUE  -> use of colors in the graph
#' @param GridLength length of the grid used to estimate the density
#' @param exportFile the name of the file to be saved. If NULL then no graph is saved. 
#' @param exportFormat the format of the export file : PNG or SVG.
#' @return A plot with the density of PhaseMin_chain + PhaseMax_chain + additionnal summary statitsics
#' @export

PhasePlot <- function(PhaseMin_chain, PhaseMax_chain, level=0.95, title = "Characterisation of a group of dates", colors = TRUE, exportFile = NULL, exportFormat = "PNG", GridLength=1024){

  if(length(PhaseMax_chain) != length(PhaseMin_chain)) { print('Error : the parameters do not have the same length')}   # test the length of both chains
  else{

  if( sum(ifelse(PhaseMin_chain <= PhaseMax_chain, 1, 0)) == length(PhaseMin_chain) ) {

    minValuex <- min( density(PhaseMin_chain, n=GridLength)$x)
    maxValuex <- min(max(density(PhaseMax_chain, n=GridLength)$x), 2016)
    
    x = 10^c(0:10)
    if(minValuex!=0){  
      c =0
      for(i in 1:length(x)) { if( abs(minValuex/x[i])>1) {c=c+1}}
      if(c>3){ minValuex = floor(minValuex/x[c-1])*x[c-1]} else {minValuex = floor(minValuex/x[c])*x[c]}
    }
    if(maxValuex!=0){  
      d=0
      for(i in 1:length(x)) { if( abs(maxValuex/x[i])>1) {d=d+1}}
      if(d>3){ maxValuex = ceiling(maxValuex/x[d-1])*x[d-1]} else {maxValuex = ceiling(maxValuex/x[d])*x[d]}
    }
    # x-axis
    middleValuex <- ( maxValuex + minValuex) / 2
    P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
    P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4

    # y-axis
    maxValuey <- max ( max(density(PhaseMin_chain, n=GridLength)$y) , max(density(PhaseMax_chain, n=GridLength)$y))
    middleValuey <- maxValuey /2
    step <- maxValuey  /20

    # Options for export
    if(!is.null(exportFile)) {
      if(exportFormat == "PNG") {
        png( filename = paste(exportFile,"png", sep =".") )
      } 
      if(exportFormat == "SVG") {
        svg( filename = paste(exportFile,"svg", sep =".") )
      }
    } 
    if (colors==TRUE){
    # first graph
    par(las=1, mfrow=c(1,1), cex.axis=0.8)
    plot(density(PhaseMax_chain, n=GridLength), main = title, xlab="Calendar Year", axes = F, ylim=c(0,maxValuey+step), xlim=c(minValuex, maxValuex), lty =1, lwd=2, col="steelblue4")
    lines(density(PhaseMin_chain, n=GridLength), lty =1, lwd=2, col ="steelblue1")

    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor(minValuex), floor(P1Valuex), floor(middleValuex), floor(P3Valuex), floor(maxValuex)))
    # ordinate axis
    axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )

    # segment representing the CredibleInterval of the max of the phase
    CIEnd = CredibleInterval(PhaseMax_chain, level)
    segments(CIEnd[2], 0, CIEnd[3], 0, lty = 1, lwd=6, col = "steelblue4")
    # point in red representing the mean
    points(mean(PhaseMax_chain), 0, lwd=6, col = "steelblue4")
    # segment representing the Time range of the phase in green
    PTR = PhaseTimeRange(PhaseMin_chain, PhaseMax_chain, level)
    segments(PTR[2], maxValuey+step, PTR[3], maxValuey+step, lwd=6, col = "violetred4")

    # segment representing the CredibleInterval in blue
    CIBeginning = CredibleInterval(PhaseMin_chain, level)
    segments(CIBeginning[2], step, CIBeginning[3], step, lty= 1, lwd=6, col = "steelblue1")
    # point in red representing the mean
    points(mean(PhaseMin_chain), step , lwd=6, col = "steelblue1")

    # legend
    legend(P3Valuex, maxValuey, c("Density of the Minimum", "Density of the Maximum", "with Credible Interval" ," and Mean (o)", " Phase Time Range"), lty=c(1,1,0,0,1), bty="n",col = c("steelblue1","steelblue4","black","black","violetred4"), lwd=c(2,2,6,6,6), x.intersp=0.5, cex=0.9)

    } else {

      # first graph
      par(las=1, mfrow=c(1,1), cex.axis=0.8)
      plot(density(PhaseMax_chain, n=GridLength), main = title, xlab="Calendar Year", axes = F, ylim=c(0,maxValuey+step), xlim=c(minValuex, maxValuex), lty =2, lwd=2)
      lines(density(PhaseMin_chain, n=GridLength), lty =3, lwd=2)

      # abscissa axis
      axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor(minValuex), floor(P1Valuex), floor(middleValuex), floor(P3Valuex), floor(maxValuex)))
      # ordinate axis
      axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )

      # segment representing the CredibleInterval in blue
      CIEnd = CredibleInterval(PhaseMax_chain, level)
      segments(CIEnd[2], 0, CIEnd[3], 0, lty = 2, lwd=6, col = 1)
      # point in red representing the mean
      points(mean(PhaseMax_chain), 0, lwd=6, col = 1)
      # segment representing the Time range of the phase in green
      PTR = PhaseTimeRange(PhaseMin_chain, PhaseMax_chain, level)
      segments(PTR[2], maxValuey+step, PTR[3], maxValuey+step, lwd=6, col = 1)

      # segment representing the CredibleInterval in blue
      CIBeginning = CredibleInterval(PhaseMin_chain, level)
      segments(CIBeginning[2], step, CIBeginning[3], step, lty= 3, lwd=6, col = 1)
      # point in red representing the mean
      points(mean(PhaseMin_chain), step , lwd=6, col = 1)

      # legend
      legend(P3Valuex, maxValuey, c("Density of the Minimum", "Density of the Maximum", "with Credible Interval", "and Mean (o)", " Phase Time Range"), lty=c(3,2,0,0,1), bty="n",col = c(1,1,1,1,1), lwd=c(2,2,6,6,6), x.intersp=0.5, cex=0.9)

    }
    
    if(!is.null(exportFile)) {
    dev.off()
    } 

  } else {
    print('Error : PhaseMin_chain should be older than PhaseMax_chain')
  }

  } # end if(length(PhaseMax_chain) != length(PhaseMin_chain))

}


#####################################################
#     Phase duration marginal density plot          #
#####################################################

#' Phase duration marginal density plot 
#'
#' Plot of the density of the time elapsed between the minimum and the maximum of the events included in a phase and summary statistics (mean, CI)
#'
#' @param PhaseMin_chain numeric vector containing the output of the MCMC algorithm for the minimum of the events included in the phase
#' @param PhaseMax_chain numeric vector containing the output of the MCMC algorithm for the maximum of the events included in the phase
#' @param level probability corresponding to the level of confidence used for the credible interval and the time range
#' @param title The Title of the graph
#' @param colors if TRUE  -> use of colors in the graph
#' @param exportFile the name of the file to be saved. If NULL then no graph is saved. 
#' @param exportFormat the format of the export file : PNG or SVG.
#' @param GridLength length of the grid used to estimate the density
#' @return A plot with the density of the duration of the phase + additionnal summary statitsics
#' @export
PhaseDurationPlot <- function(PhaseMin_chain, PhaseMax_chain, level=0.95, title = "Duration of a group of dates", colors = TRUE, exportFile = NULL, exportFormat = "PNG", GridLength=1024){
  
  if(length(PhaseMax_chain) != length(PhaseMin_chain)) { print('Error : the parameters do not have the same length')}   # test the length of both chains
  else{
    if( sum(ifelse(PhaseMin_chain <= PhaseMax_chain, 1, 0)) == length(PhaseMin_chain) ) { 
    
      if(   sum(ifelse(PhaseMin_chain == PhaseMax_chain, 1, 0)) == length(PhaseMin_chain)    ) { print('Error : no duration')}
      else{
      
      a_chain = PhaseMax_chain -PhaseMin_chain 
    
      step <- max(density(a_chain, n=GridLength)$y) /50   # used to draw CI and mean above the curve
      maxValuex <- max(density(a_chain, n=GridLength)$x)
      minValuex <- min(density(a_chain, n=GridLength)$x)
      middleValuex <- minValuex + ( maxValuex - minValuex ) / 2
      P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
      P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
      maxValuey <- max(density(a_chain, n=GridLength)$y)
      middleValuey <- maxValuey /2
      
      # Options for export
      if(!is.null(exportFile)) {
        
        if(exportFormat == "PNG") {
          png( filename = paste(exportFile,"png", sep =".") )
        } 
        if(exportFormat == "SVG") {
          svg( filename = paste(exportFile,"svg", sep =".") )
        }
        
      } 
      if (colors==T){
        par(mfrow=c(1,1))
        plot(density(a_chain, n=GridLength), main = title, xlab = "Calendar Year", axes = F, ylim=c(0,max(density(a_chain, n=GridLength)$y) + step))
        # abscissa axis
        axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor(minValuex), floor(P1Valuex), floor(middleValuex), floor(P3Valuex), floor(maxValuex )))
        # ordinate axis
        axis(2, at=c(0, middleValuey , maxValuey),labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )
        
        segments(CredibleInterval(a_chain, level=level)[2], 0, CredibleInterval(a_chain, level=level)[3], 0, lwd=6, col = 4)
        points(mean(a_chain), 0 , lwd=6, col = 2)
        
        # legend
        legend(P3Valuex, maxValuey, c("Density", "Credible Interval", "Mean"), lty=c(1, 1, 0), bty="n", pch=c(NA, NA,1), col = c("black","blue","red"), lwd=c(1,6,6), x.intersp=0.5, cex=0.9)
        
      }else {
        par(mfrow=c(1,1))
        plot(density(a_chain, n=GridLength), main = title, xlab = "Calendar Year", axes = F, ylim=c(0,max(density(a_chain, n=GridLength)$y) + step))
        # abscissa axis
        axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor(minValuex), floor(P1Valuex), floor(middleValuex), floor(P3Valuex), floor(maxValuex )))
        # ordinate axis
        axis(2, at=c(0, middleValuey , maxValuey),labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )
        
        segments(CredibleInterval(a_chain, level=level)[2], 0, CredibleInterval(a_chain, level=level)[3], 0, lwd=6, lty=1)
        points(mean(a_chain), 0 , lwd=6, pch=1)
        
        # legend
        legend(P3Valuex, maxValuey, c("Density", "Credible Interval", "Mean"), lty=c(1,1,0), pch=c(NA,NA,1), bty="n", lwd=c(1,6,6), x.intersp=0.5, cex=0.9)
      }
      
      # options for export
      if(!is.null(exportFile)) {
        dev.off()
      } 
      
      }
      
    } else {
      print('Error : PhaseMin_chain should be older than PhaseMax_chain')
    }
    
  #} # end if(length(PhaseMax_chain) != length(PhaseMin_chain))
  } # if( sum(ifelse(PhaseMin_chain == PhaseMax_chain, 1, 0)) == length(PhaseMin_chain) )
}


#####################################################
#          Hiatus between two phases             #
#####################################################
#'  Gap/Hiatus between two successive phases (for phases in temporal order constraint)
#'
#' Finds if it exists a gap between two phases that is the longest interval that satisfies : P(Phase1Max_chain < IntervalInf < IntervalSup < Phase2Min_chain | M) = level
#'
#' @param Phase1Max_chain numeric vector containing the output of the MCMC algorithm for the maximum of the events included in the oldest phase
#' @param Phase2Min_chain numeric vector containing the output of the MCMC algorithm for the minimum of the events included in the following phase
#' @param level probability corresponding to the level of confidence
#' @return The endpoints of the longest gap
#' @export
PhasesGap <- function(Phase1Max_chain, Phase2Min_chain, level=0.95){

  if(length(Phase1Max_chain) != length(Phase2Min_chain)) { stop('Error : the parameters do not have the same length')} # test for the length of both chains
    else{

      if( sum(ifelse(Phase1Max_chain <=Phase2Min_chain, 1, 0)) != length(Phase2Min_chain) )  {  # test for Phase1Max_chain < Phase2Min_chain
        stop('Error : Phase1Max_chain should be older than Phase2Min_chain')
      } else {

      interval <- function(epsilon, P1Max, P2Min, level)
      {
        q1 = quantile(P1Max ,probs = 1-epsilon) ;
        indz = (P1Max < q1)
        q2 = quantile(P2Min[indz],probs= (1-level-epsilon)/(1-epsilon))
        c(q1,q2)
      }
      hia = Vectorize(interval,"epsilon")

      epsilon = seq(0,1-level,.001)
      p = hia(epsilon, Phase1Max_chain, Phase2Min_chain, level)
      rownames(p)<- c("HiatusIntervalInf", "HiatusIntervalSup")

      D<- p[2,]-p[1,]
      DD = D[D>0]

      if (length(DD) > 0){
        I = which(D==max(DD))
        interval2 = round( p[,I], 0)
        if (p[2,I] != p[1,I]) {
          c(level=level, interval2[1], interval2[2])
        } else {
          c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')
        }#end if (p[2,I] != p[1,I])

      } else {
        c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')
        }#end if (length(DD) > 0)

      } # end if( sum(ifelse(Phase2Min_chain < Phase1Max_chain, 1, 0) == length(Phase2Min_chain) ) ) {  # test for Phase1Max_chain < Phase2Min_chain

    } # if(length(Phase1Max_chain) != length(Phase2Min_chain)) {

}





#####################################################
#                Phases Transition                  #
#####################################################

#'  Transition range between two successive phases (for phases in temporal order constraint)
#'
#' Finds if it exists the shortest interval that satisfies : P(TransitionRangeInf < Phase1Max_chain < Phase2Min_chain < TransitionRangeSup  | M) = level
#'
#' @param Phase1Max_chain numeric vector containing the output of the MCMC algorithm for the maximum of the events included in the oldest phase
#' @param Phase2Min_chain numeric vector containing the output of the MCMC algorithm for the minimum of the events included in the following phase
#' @param level probability corresponding to the level of confidence
#' @return the endpoints of the transition interval
#' @export
PhasesTransition <- function(Phase1Max_chain, Phase2Min_chain, level=0.95){

  result = as.matrix( PhaseTimeRange(Phase1Max_chain, Phase2Min_chain, level=level))
  rownames(result)<- c(level=level, "TransitionRangeInf", "TransitionRangeSup")
  result <- t(result)

  return(result[1,])
}




#####################################################
#             Succession Plot                  #
#####################################################

#' Density Plots of two successive phases (for phases in temporal order constraint)
#'
#' Plot of the densities of two successive phases + statistics (mean, CI, HPDR)
#'
#' @param Phase1Min_chain numeric vector containing the output of the MCMC algorithm for the minimum of the events included in the oldest phase
#' @param Phase1Max_chain numeric vector containing the output of the MCMC algorithm for the maximum of the events included in the oldest phase
#' @param Phase2Min_chain numeric vector containing the output of the MCMC algorithm for the minimum of the events included in the youngest phase
#' @param Phase2Max_chain numeric vector containing the output of the MCMC algorithm for the maximum of the events included in the youngest phase
#' @param level probability corresponding to the level of confidence
#' @param title title of the graph
#' @param exportFile the name of the file to be saved. If NULL then no graph is saved. 
#' @param exportFormat the format of the export file : PNG or SVG.
#' @param GridLength length of the grid used to estimate the density
#' @return a plot of all densities + CI + mean + HDR
#' @export

SuccessionPlot <- function(Phase1Min_chain, Phase1Max_chain, Phase2Min_chain, Phase2Max_chain, level=0.95,  title = "Characterisation of a succession of groups", 
                           exportFile = NULL, exportFormat = "PNG", GridLength=1024){


  if(length(Phase1Max_chain) != length(Phase2Min_chain)) { stop('Error : the parameters do not have the same length')} # test for the length of both chains
  else{

  if( sum(ifelse(Phase1Min_chain <= Phase1Max_chain, 1, 0)) != length(Phase1Min_chain) ||  sum(ifelse(Phase2Min_chain <= Phase2Max_chain, 1, 0)) != length(Phase1Min_chain) || sum(ifelse( Phase1Max_chain <= Phase2Min_chain, 1, 0)) != length(Phase1Min_chain) ) {
    # test for PhaseMin_chain < PhaseMax_chain and Phase1 < Phase2
    stop('Error : PhaseMin_chain should be older than PhaseMax_chain')
  } else {

    minValuex <- min(density(Phase1Min_chain, n=GridLength)$x, density(Phase2Min_chain, n=GridLength)$x)
    maxValuex <- min( max(density(Phase1Max_chain, n=GridLength)$x, density(Phase2Max_chain, n=GridLength)$x), 2016)
    
    x = 10^c(0:10)
    if(minValuex!=0){  
      c =0
      for(i in 1:length(x)) { if( abs(minValuex/x[i])>1) {c=c+1}}
      if(c>3){ minValuex = floor(minValuex/x[c-1])*x[c-1]} else {minValuex = floor(minValuex/x[c])*x[c]}
    }
    if(maxValuex!=0){  
      d=0
      for(i in 1:length(x)) { if( abs(maxValuex/x[i])>1) {d=d+1}}
      if(d>3){ maxValuex = ceiling(maxValuex/x[d-1])*x[d-1]} else {maxValuex = ceiling(maxValuex/x[d])*x[d]}
    }
    
    # x-axis
    middleValuex <- ( maxValuex + minValuex) / 2
    P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
    P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
    
    # y-axis
    maxValuey <- max ( max(density(Phase1Min_chain, n=GridLength)$y) , max(density(Phase1Max_chain, n=GridLength)$y), max(density(Phase2Min_chain, n=GridLength)$y) , max(density(Phase2Max_chain, n=GridLength)$y))
    middleValuey <- maxValuey /2
    minValuey <- min ( min(density(Phase1Min_chain, n=GridLength)$y) , min(density(Phase1Max_chain, n=GridLength)$y), min(density(Phase2Min_chain, n=GridLength)$y) , min(density(Phase2Max_chain, n=GridLength)$y))

    haut = seq(minValuey,maxValuey,length.out=5)
    middleA <- maxValuey+ (haut[1] + haut[2]) / 2

    # Options for export
    if(!is.null(exportFile)) {
      
      if(exportFormat == "PNG") {
        png( filename = paste(exportFile,"png", sep =".") )
      } 
      if(exportFormat == "SVG") {
        svg( filename = paste(exportFile,"svg", sep =".") )
      }
      
    } 
    
    par(las=1, mfrow=c(1,1), cex.axis=0.8)
    plot(density(Phase1Max_chain, n=GridLength), main = title, ylab="Density", xlab = "Calendar Year", ylim=c(0,maxValuey+maxValuey), xlim=c(minValuex, maxValuex), bty='n',lty =1, lwd=2, axes=F, col = "steelblue")
    lines(density(Phase1Min_chain, n=GridLength), lty =1, lwd=2, col = "steelblue")

    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor( minValuex), floor( P1Valuex), floor( middleValuex), floor( P3Valuex), floor( maxValuex)))
    # ordinate axis
    axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )

    ## Phase2
    lines(density(Phase2Min_chain, n=GridLength), lty =1, lwd=2, col ="violet")
    lines(density(Phase2Max_chain, n=GridLength), lty =1, lwd=2, col ="violet")

    ## Phase Time Range
    PTR1 = PhaseTimeRange(Phase1Min_chain, Phase1Max_chain, level=level)
    PTR2 = PhaseTimeRange(Phase2Min_chain, Phase2Max_chain, level=level)
    segments(PTR1[2],maxValuey+haut[2],PTR1[3],maxValuey+haut[2],lwd=6,col="steelblue")
    segments(PTR2[2],maxValuey+haut[3],PTR2[3],maxValuey+haut[3],lwd=6, col ="violet")
    text(minValuex, middleA,"Time range",srt =90)

    ## Phase Transition
    PTrans = PhasesTransition(Phase1Max_chain, Phase2Min_chain, level=level)
    segments(PTrans[2],maxValuey+haut[4],PTrans[3],maxValuey+haut[4],lwd=6, col = "steelblue")
    segments(PTrans[2],maxValuey+haut[4],PTrans[3],maxValuey+haut[4],lwd=6, col = "violet", lty=4)

    PGap = PhasesGap(Phase1Max_chain, Phase2Min_chain, level=level)
    if (PGap[2] == "NA" || PGap[3] == "NA") {
      points( (PTrans[3]+PTrans[2])/2, maxValuey+haut[5], lwd=2, col = "steelblue", pch=4)
    } else {
      segments(PGap[2],maxValuey+haut[5],PGap[3],maxValuey+haut[5],lwd=6, col = "steelblue")
      segments(PGap[2],maxValuey+haut[5],PGap[3],maxValuey+haut[5],lwd=6, col = "violet", lty=4)
    }

    text(minValuex, maxValuey+haut[4],"Transition",srt =90)
    text(minValuex, maxValuey+haut[5],"Gap",srt =90)

    # options for export
    if(!is.null(exportFile)) {
      dev.off()
    } 
    
  }

  }

}







#####################################################
#     Estimation of Credible interval        #
#####################################################

#' Bayesian credible interval for a series of MCMC chains
#'
#' Estimation of the shorest credible interval of the output of the MCMC algorithm for the parameter a
#'
#' @details A 100*level % credible interval is an interval that keeps N*(1-level) elements of the sample outside the interval
#' The 100*level % credible interval is the shortest of all those intervals.
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence used for the credible interval
#' @return The endpoints of the shortest credible interval
#' @export
#'
MultiCredibleInterval <- function(data, position, level=0.95){

  # number of chains
  L = length(position)

  # matrix of results for each pair of phases
  result = matrix(nrow=L, ncol=3)

  colnames(result) <- c("Level","CredibleIntervalInf", "CredibleIntervalSup")
  
  # names
  rownames(result) <- names(data)[position]

  for (i in 1:L) {

    sorted_sample <- sort(data[,position[i]])     # ordering the sample
    N = length(sorted_sample)                     # calculation of the sample size of the chain
    OutSample = N * (1-level)          # calculation of the number of data to be outside the interval

    I =  cbind(sorted_sample[1:(OutSample+1)] , sorted_sample[(N-OutSample):N])    #   combinasion of all credible intervals

    l = I[,2]-I[,1]   # length of intervals
    j <- which.min(l) # look for the shortest interval

    result[i,] =   c(level, round(I[j,1], digits = 0) , round(I[j,2],digits = 0) )   # returns the level and the endpoints

  }
  return(result)
}





#####################################################
#                   MultiHPD                        #
#####################################################

#' Bayesian HPD regions for a series of MCMC chains
#'
#' Estimation of the HPD region of the output of the MCMC algorithm for the parameter a
#'
#' @details Highest posterior density function region using the function 'hdr' from the package 'hdrcde'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence used for the credible interval
#' @return The endpoints of the shortest credible interval
#' @export
#'
MultiHPD <- function(data, position, level=0.95){
  
  # matrix of results for each pair of phases
    hdr = hdr(data[,position[1]], prob = c(level * 100))$hdr
    HPDR = round(hdr, digits = 0)  
    result = matrix( c(level, HPDR), nrow=1)
    dim = dim(result)[2]
    
    if(length(position) > 1){
      
      for (i in 2:length(position)) {
        
        hdr = hdr(data[,position[i]], prob = c(level * 100))$hdr
        HPDR = round(hdr, digits = 0) 
        res = c(level, HPDR)
        
        if (length(res) > dim) {
          NbCol = length(res) - dim
          AddColum = rep(NA, i-1) 
          Ajout = NULL 
          for (j in 1:NbCol){
            Ajout = cbind(Ajout, AddColum)
          }
          resultTemp = cbind(result, Ajout)
          result =  rbind(resultTemp, res)
          
        }else if (length(res) < dim) {
          NbCol = dim - length(res) 
          Add = rep(NA, NbCol) 
          Ajout = c(res,Add)

          result = rbind(result, Ajout)
          
        }else{
          result =  rbind(result, res)  
        }
        dim = dim(result)[2] 
      }

    }

    nom=c()
    for( k in (1:((dim-1)/2)) ) {
      nom=c(nom,paste("HPDRInf",k))
      nom=c(nom,paste("HPRDSup",k))
    }
    colnames(result) <- c("Level", nom)
    rownames(result) <- names(data)[position]
    
  return(result) # returns a matrix with the level and the endpoints
}

###############################################
#             MultiDatesPlot                  #
###############################################

#' Plot of credible intervals or HPD regions of a series of dates
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence
#' @param title title of the graph
#' @param exportFile the name of the file to be saved. If NULL then no graph is saved. 
#' @param exportFormat the format of the export file : PNG or SVG.
#' @return a plot of the endpoints of the credible intervals of a series of dates
#' @export

MultiDatesPlot <- function(data, position, level=0.95, intervals = "CI", title = "Plot of intervals", labelXaxis = "Calendar Year", exportFile = NULL, exportFormat = "PNG"){
  
  if(intervals =="CI"){
  Bornes = MultiCredibleInterval(data, position, level=level) 
  }else if(intervals =="HPD") {
  Bornes = MultiHPD(data, position, level=level) 
  }
  Ordered = Bornes[order(Bornes[,2]),]
  nbCol = dim(Ordered)[2]
  
  NbDates = length(position)
  DatesNames <- rownames(Ordered)
  
  minValuex <- floor( min(Ordered[,2], na.rm = TRUE)/100) * 100
  maxValuex <- ceiling( max(Ordered[,-c(1,2)],na.rm = TRUE)/100) * 100
  middleValuex <- ( maxValuex + minValuex) / 2
  P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
  P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
  seq = seq(minValuex, maxValuex)
  
  
  ## Graph
  par(las=1, mfrow=c(1,1), cex.axis=0.8, mar=c(5,6,4,2))
  plot(0, main = title, ylab="", xlab = labelXaxis, ylim=c(0, NbDates), xlim=c(minValuex, maxValuex), type="n", axes=F)
  grid()
  # abscissa axis
  axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) ) 
  # ordinate axis
  axis(2, at=1:NbDates, labels =DatesNames, las =2)
  
  ## Phase Time Range
  for (i in 1:NbDates ) { for (j in seq(2,(nbCol-1), by = 2)) { segments(Ordered[i,j], i, Ordered[i,j+1], i, lwd=6) } }
 
  # Options for export
  if(!is.null(exportFile)) {
    
    if(exportFormat == "PNG") {
      png( filename = paste(exportFile,"png", sep =".") )
    } 
    if(exportFormat == "SVG") {
      svg( filename = paste(exportFile,"svg", sep =".") )
    }
    
    par(las=1, mfrow=c(1,1), cex.axis=0.8, mar=c(5,6,4,2))
    plot(0, main = title, ylab="", xlab = labelXaxis, ylim=c(0, NbDates), xlim=c(minValuex, maxValuex), type="n", axes=F)
    grid()
    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) ) 
    # ordinate axis
    axis(2, at=1:NbDates, labels =DatesNames, las =2)
    
    ## Phase Time Range
    for (i in 1:NbDates ) { for (j in seq(2,(nbCol-1), by = 2)) { segments(Ordered[i,j], i, Ordered[i,j+1], i, lwd=6) } }
    
    dev.off()
  } 
  
}



#####################################################
#         Multiple Phase Time Range                 #
#####################################################

#' Phase Time Range for multiple groups
#'
#' Computes the shortest interval that satisfies : P(PhaseMin < IntervalInf < IntervalSup < PhaseMax | M) = level for each phase
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_minimum numeric vector containing the column number corresponding to the minimum of the events included in each phase
#' @param position_maximum numeric vector containing the column number corresponding to the maximum of the phases set in the same order as in position_minimum
#' @param level probability corresponding to the desired level of confidence
#' @return The endpoints of the shortest time range associated with the desired level
#' @export


MultiPhaseTimeRange <- function(data, position_minimum, position_maximum=position_minimum+1, level=0.95){
  
  if (length(position_minimum)!= length(position_maximum)) {
    print('Error : the position vectors do not have the same length')
    } else {

  # number of phases
  L = length(position_minimum)
  
  # names
  names_beginning <- names(data)[position_minimum]
  names_end <- names(data)[position_maximum]
  
  # matrix of results
  result = matrix(nrow=L, ncol=3)
  colnames(result)<- c("Level","TimeRangeInf", "TimeRangeSup")
  
  phasenames <- vector(length = L)
  for (i in 1:L) { phasenames[i] = paste(names_beginning[i], names_end[i] ) }
  rownames(result)<- phasenames
  
  for (i in 1:L){
    result[i,] = PhaseTimeRange( data[,position_minimum[i]], data[,position_maximum[i]], level=level)
  }
  
  return(result)
  
  
  }
}

#####################################################
#       Hiatus between a succession of  groups      #
#####################################################

#'  Gap/Hiatus between a succession of groups (for groups in temporal order constraint)
#'
#' Finds if it exists a gap between two groups that is the longest interval that satisfies : P(Phase1Max < IntervalInf < IntervalSup < Phase2Min | M) = level
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_minimum numeric vector containing the column number corresponding to the minimum of the events included in each group
#' @param position_maximum numeric vector containing the column number corresponding to the end of the phases set in the same order as in position_minimum
#' @param level probability corresponding to the level of confidence
#' @return The endpoints of the longest gap
#' @export

MultiPhasesGap <- function(data, position_minimum, position_maximum = position_minimum+1, level=0.95){
  
  if (length(position_minimum)!= length(position_maximum)) {
    print('Error : the position vectors do not have the same length')
  } else {
    
  # number of phases
  L = length(position_minimum)
  
  #names
  names_Min <- names(data)[position_minimum]
  names_Max <- names(data)[position_maximum]
  
  # matrix of results
  result = matrix(nrow=L-1, ncol=3)
  colnames(result)<- c("Level","HiatusIntervalInf", "HiatusIntervalSup")
  
  phasenames <- vector(length = (L-1))
  for (i in 1:L-1) { phasenames[i] = paste(names_Max[i], "&", names_Min[i+1]) }
  rownames(result)<- phasenames
  
  for (i in 1:(L-1)){
    result[i,] = PhasesGap( data[,position_maximum[i]], data[,position_minimum[i+1]], level=level)
  }
  
  return(result)
  
  }
}



#####################################################
#         Multiple Phases Transition               #
#####################################################

#'  Transition range for a succession of groups (for groups in temporal order constraint)
#'
#' Finds if it exists the shortest interval that satisfies : P(TransitionRangeInf < Phase1Max < Phase2Min < TransitionRangeSup  | M) = level
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_minimum numeric vector containing the column number corresponding to the minimum of the events included in each group
#' @param position_maximum numeric vector containing the column number corresponding to the end of the groups set in the same order as in position_minimum
#' @param level probability corresponding to the level of confidence
#' @return the endpoints of the transition interval
#' @export


MultiPhasesTransition <- function(data, position_minimum, position_maximum = position_minimum+1, level=0.95){

  if (length(position_minimum)!= length(position_maximum)) {
    print('Error : the position vectors do not have the same length')
  } else {
    
  # number of phases
  L = length(position_minimum)
  
  #names
  names_min <- names(data)[position_minimum]
  names_max <- names(data)[position_maximum]

  # matrix of results
  result = matrix(nrow=L-1, ncol=3)
  colnames(result)<- c(level, "TransitionRangeInf", "TransitionRangeSup")
  
  phasenames <- vector(length = L-1)
  for (i in 1:L-1) { phasenames[i] = paste(names_max[i], "&", names_min[i+1]) }
  rownames(result)<- phasenames

  for (i in 1:(L-1)){
    result[i,] = PhaseTimeRange( data[,position_maximum[i]], data[,position_minimum[i+1]], level=level)
  }

  return(result)

  }
}





#####################################################
#        Multiple Phases Density Plots              #
#####################################################

#' Successive Phases Density Plots (for phases in temporal order constraint)
#'
#' Plot of the densities of several successive groups + statistics (mean, CI, HPDR)
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_minimum numeric vector containing the column number corresponding to the minimum of the events included in each group
#' @param position_maximum numeric vector containing the column number corresponding to the end of the groups set in the same order as in position_minimum
#' @param level probability corresponding to the level of confidence
#' @param title title of the graph
#' @param colors vector of colors corresponding to each group of dates
#' @param exportFile the name of the file to be saved. If NULL then no graph is saved. 
#' @param exportFormat the format of the export file : PNG or SVG.
#' @return a plot of all densities + CI + mean + HDR
#' @export

MultiSuccessionPlot <- function(data, position_minimum, position_maximum = position_minimum+1, level=0.95, title = "Characterisation of a succession of groups", 
                                colors = NULL, exportFile = NULL, exportFormat = "PNG"){
  
  if (length(position_minimum)!= length(position_maximum)) {
    print('Error : the position vectors do not have the same length')
  } else {
    
    # construction of the new dataset
    L = length(position_minimum)
    
    GridLength=1024
    phase = matrix(ncol = L*2, nrow=nrow(data))
    densityX = matrix(ncol = L*2, nrow=GridLength)
    densityY = matrix(ncol = L*2, nrow=GridLength)
    
    for (i in 1:L) {
      phase[,2*i-1] = data[,position_minimum[i]]
      phase[,2*i] = data[,position_maximum[i]]
      
      densityX[,2*i-1] = density(data[,position_minimum[i]], n = GridLength)$x
      densityX[,2*i] = density(data[,position_maximum[i]], n = GridLength)$x
      
      densityY[,2*i-1] = density(data[,position_minimum[i]], n = GridLength)$y
      densityY[,2*i] = density(data[,position_maximum[i]], n = GridLength)$y
    }
    
    minValuex <- min (apply(densityX,2,min))
    maxValuex <- min(max( apply(densityX,2,max) ), 2016)
    
    # rounding up x and y values
    x = 10^c(0:10)
    if(minValuex!=0){  
      c =0
      for(i in 1:length(x)) { if( abs(minValuex/x[i])>1) {c=c+1}}
      if(c>3){ minValuex = floor(minValuex/x[c-1])*x[c-1]} else {minValuex = floor(minValuex/x[c])*x[c]}
    }
    if(maxValuex!=0){  
      d=0
      for(i in 1:length(x)) { if( abs(maxValuex/x[i])>1) {d=d+1}}
      if(d>3){ maxValuex = ceiling(maxValuex/x[d-1])*x[d-1]} else {maxValuex = ceiling(maxValuex/x[d])*x[d]}
    }
    # x-axis values
    middleValuex <- ( maxValuex + minValuex) / 2
    P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
    P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
    # y-axis values
    maxValuey <- max( apply(densityY,2,max))
    middleValuey <- maxValuey /2
    minValuey <- min( apply(densityY,2,min))
    # 
    haut = seq(minValuey,maxValuey,length.out=(2*L+2) )
    
    # Options for colors
    if (is.null(colors)) {
      pal = rainbow(L)
    } else {
      pal = colors
    }
    
    ## Graph
    par(las=1, mfrow=c(1,1), cex.axis=0.8)
    plot(density(phase[,1], n=GridLength), main = title, ylab="Density", xlab = "Calendar Year", ylim=c(0,2.5*maxValuey), xlim=c(minValuex, maxValuex), bty='n',lty =1, lwd=2, axes=F, col = pal[1])
    lines(density(phase[,2], n=GridLength), lty =1, lwd=2, col = pal[1])
    
    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor( minValuex), floor( P1Valuex), floor( middleValuex), floor( P3Valuex), floor( maxValuex)))
    # ordinate axis
    axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 5), round(maxValuey, 5)) )
    
    
    ## Following phases
    for(i in 2:L) {
      lines(density(phase[,2*i-1]), col=pal[i], lwd=2, lty=1)
      lines(density(phase[,2*i]), col=pal[i], lwd=2, lty=1)
    }
    
    ## Phase Time Range
    MPTR = MultiPhaseTimeRange(data=data, position_minimum=position_minimum, position_maximum=position_maximum, level=level)
    for (i in 1:(L) ) { segments(MPTR[i,2], maxValuey + haut[i+1], MPTR[i,3], maxValuey + haut[i+1], lwd=6,col=pal[i]) }
    text(minValuex, maxValuey + haut[trunc((L+1)/2)] , "Time range",srt =90)
    
    ## Phase Transition / Gap
    PTrans = MultiPhasesTransition(data=data, position_minimum=position_minimum, position_maximum=position_maximum, level=level)
    PGap = MultiPhasesGap(data=data, position_minimum=position_minimum, position_maximum=position_maximum, level=level)
    
    # Separateur 
    segments(minValuex, maxValuey+haut[L+2], maxValuex, maxValuey+haut[L+2], lwd=.2)
    text(minValuex, maxValuey+haut[2*L+1], "Transition ",srt =90)
    
    segments(minValuex, 2*maxValuey+haut[1], maxValuex, 2*maxValuey+haut[1], lwd=.2)
    text(minValuex, 2*maxValuey+haut[trunc((L+2)/2)]," Gap",srt =90)
    

    for (i in 1:(L-1) ) {
      segments(PTrans[i,2], maxValuey+haut[L+i+2],PTrans[i,3], maxValuey+haut[L+i+2],lwd=6, col = pal[i])
      segments(PTrans[i,2], maxValuey+haut[L+i+2],PTrans[i,3], maxValuey+haut[L+i+2],lwd=6, col = pal[i+1], lty=4)
      
      if (PGap[i,2] == "NA" || PGap[i,3] == "NA") {
        points( (PTrans[i,3]+PTrans[i,2])/2, 2*maxValuey+haut[i+1], lwd=2, col = pal[i], pch=4)
        
      } else {
        segments(as.numeric(PGap[i,2]), 2*maxValuey+haut[i+1], as.numeric(PGap[i,3]), 2*maxValuey+haut[i+1], lwd=6, col = pal[i])
        segments(as.numeric(PGap[i,2]), 2*maxValuey+haut[i+1], as.numeric(PGap[i,3]), 2*maxValuey+haut[i+1], lwd=6, col = pal[i+1], lty=4)
      }
      
    }
    
    # Options for export
    if(!is.null(exportFile)) {
      
      if(exportFormat == "PNG") {
        png( filename = paste(exportFile,"png", sep =".") )
      } 
      if(exportFormat == "SVG") {
        svg( filename = paste(exportFile,"svg", sep =".") )
      }
      
      par(las=1, mfrow=c(1,1), cex.axis=0.8)
      plot(density(phase[,1], n=GridLength), main = title, ylab="Density", xlab = "Calendar Year", ylim=c(0,2.5*maxValuey), xlim=c(minValuex, maxValuex), bty='n',lty =1, lwd=2, axes=F, col = pal[1])
      lines(density(phase[,2], n=GridLength), lty =1, lwd=2, col = pal[1])
      
      # abscissa axis
      axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor( minValuex), floor( P1Valuex), floor( middleValuex), floor( P3Valuex), floor( maxValuex)))
      # ordinate axis
      axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 5), round(maxValuey, 5)) )
      
      
      ## Following phases
      for(i in 2:L) {
        lines(density(phase[,2*i-1]), col=pal[i], lwd=2, lty=1)
        lines(density(phase[,2*i]), col=pal[i], lwd=2, lty=1)
      }
      
      ## Phase Time Range
      MPTR = MultiPhaseTimeRange(data=data, position_minimum=position_minimum, position_maximum=position_maximum, level=level)
      for (i in 1:(L) ) { segments(MPTR[i,2], maxValuey + haut[i+1], MPTR[i,3], maxValuey + haut[i+1], lwd=6,col=pal[i]) }
      text(minValuex, maxValuey + haut[trunc((L+1)/2)] , "Time range",srt =90)
      
      ## Phase Transition / Gap
      PTrans = MultiPhasesTransition(data=data, position_minimum=position_minimum, position_maximum=position_maximum, level=level)
      PGap = MultiPhasesGap(data=data, position_minimum=position_minimum, position_maximum=position_maximum, level=level)
      
      # Separateur 
      segments(minValuex, maxValuey+haut[L+2], maxValuex, maxValuey+haut[L+2], lwd=.2)
      text(minValuex, maxValuey+haut[2*L+1], "Transition ",srt =90)
      
      segments(minValuex, 2*maxValuey+haut[1], maxValuex, 2*maxValuey+haut[1], lwd=.2)
      text(minValuex, 2*maxValuey+haut[trunc((L+2)/2)]," Gap",srt =90)
      
      
      for (i in 1:(L-1) ) {
        segments(PTrans[i,2], maxValuey+haut[L+i+2],PTrans[i,3], maxValuey+haut[L+i+2],lwd=6, col = pal[i])
        segments(PTrans[i,2], maxValuey+haut[L+i+2],PTrans[i,3], maxValuey+haut[L+i+2],lwd=6, col = pal[i+1], lty=4)
        
        if (PGap[i,2] == "NA" || PGap[i,3] == "NA") {
          points( (PTrans[i,3]+PTrans[i,2])/2, 2*maxValuey+haut[i+1], lwd=2, col = pal[i], pch=4)
          
        } else {
          segments(as.numeric(PGap[i,2]), 2*maxValuey+haut[i+1], as.numeric(PGap[i,3]), 2*maxValuey+haut[i+1], lwd=6, col = pal[i])
          segments(as.numeric(PGap[i,2]), 2*maxValuey+haut[i+1], as.numeric(PGap[i,3]), 2*maxValuey+haut[i+1], lwd=6, col = pal[i+1], lty=4)
        }
        
      }
      
      dev.off()
    }
    
  }
}




#####################################################
#           Multiple Phases     Plots            #
#####################################################

#' Several Phases Density Plots
#'
#' Plot of the densities of several groups + statistics (mean, CI, HPDR)
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_minimum numeric vector containing the column number corresponding to the minimum of the events included in each group
#' @param position_maximum numeric vector containing the column number corresponding to the end of the groups set in the same order as in position_minimum
#' @param level probability corresponding to the level of confidence
#' @param title title of the graph
#' @param colors vector of colors corresponding to each group of dates
#' @param exportFile the name of the file to be saved. If NULL then no graph is saved. 
#' @param exportFormat the format of the export file : PNG or SVG.
#' @return a plot of all densities + CI + mean + HDR
#' @export

MultiPhasePlot <- function(data, position_minimum, position_maximum = position_minimum+1, level=0.95, title = "Characterisation of several groups",
                           colors = NULL, exportFile = NULL, exportFormat = "PNG"){

  if (length(position_minimum)!= length(position_maximum)) {
    print('Error : the position vectors do not have the same length')
  } else {
    
  # construction of the new dataset
  L = length(position_minimum) 
  
  phase = matrix(ncol = L*2, nrow=nrow(data))
  GridLength = 1024
  densityX = matrix(ncol = L*2, nrow=GridLength)
  densityY = matrix(ncol = L*2, nrow=GridLength)

  for (i in 1:L) {
    phase[,2*i-1] = data[,position_minimum[i]]
    phase[,2*i] = data[,position_maximum[i]]
    
    densityX[,2*i-1] = density(data[,position_minimum[i]], n=1024)$x
    densityX[,2*i] = density(data[,position_maximum[i]], n=1024)$x
    
    densityY[,2*i-1] = density(data[,position_minimum[i]], n=1024)$y
    densityY[,2*i] = density(data[,position_maximum[i]], n=1024)$y
  }

  minValuex <- min (apply(densityX,2,min))
  maxValuex <- min(max( apply(densityX,2,max) ), 2016)
  
  # rounding up x and y values
  x = 10^c(0:10)
  if(minValuex!=0){  
    c =0
    for(i in 1:length(x)) { if( abs(minValuex/x[i])>1) {c=c+1}}
    if(c>3){ minValuex = floor(minValuex/x[c-1])*x[c-1]} else {minValuex = floor(minValuex/x[c])*x[c]}
  }
  if(maxValuex!=0){  
    d=0
    for(i in 1:length(x)) { if( abs(maxValuex/x[i])>1) {d=d+1}}
    if(d>3){ maxValuex = ceiling(maxValuex/x[d-1])*x[d-1]} else {maxValuex = ceiling(maxValuex/x[d])*x[d]}
  }
  
  # x-axis
  middleValuex <- ( maxValuex + minValuex) / 2
  P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
  P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
  # y-axis
  maxValuey <- max( apply(densityY,2,max))
  middleValuey <- maxValuey /2
  minValuey <- min( apply(densityY,2,min))

  haut = seq(minValuey,middleValuey,length.out=(L+1) )

  # Options for colors
  if (is.null(colors)) {
    pal = rainbow(L)
  } else {
    pal = colors
  }
  

  # Graph
  par(las=1, mfrow=c(1,1), cex.axis=0.8)
  plot(density(phase[,1], n=GridLength), main = title, ylab="Density", xlab = "Calendar Year", ylim=c(0,maxValuey+middleValuey), xlim=c(minValuex, maxValuex), bty='n',lty =1, lwd=2, axes=F, col = pal[1])
  lines(density(phase[,2], n=GridLength), lty =1, lwd=2, col = pal[1])

  # abscissa axis
  axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor( minValuex), floor( P1Valuex), floor( middleValuex), floor( P3Valuex), floor( maxValuex)))
  # ordinate axis
  axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 5), round(maxValuey, 5)) )

  ## Following phases
  for(i in 2:L) {
    lines(density(phase[,2*i-1], n = GridLength), col=pal[i], lwd=2, lty=1)
    lines(density(phase[,2*i], n = GridLength), col=pal[i], lwd=2, lty=1)
  }

  # ## Phase Time Range
   MPTR = MultiPhaseTimeRange(data=data, position_minimum=position_minimum, position_maximum=position_maximum, level=level)
   for (i in 1:L ) { segments(MPTR[i,2], maxValuey+haut[i+1],MPTR[i,3], maxValuey+haut[i+1],lwd=6,col=pal[i]) }
   text(minValuex, maxValuey+haut[2], "Time range", srt =90)
  
   # Options for export
   if(!is.null(exportFile)) {
     
     if(exportFormat == "PNG") {
       png( filename = paste(exportFile,"png", sep =".") )
     } 
     if(exportFormat == "SVG") {
       svg( filename = paste(exportFile,"svg", sep =".") )
     }
     
     # Graph
     par(las=1, mfrow=c(1,1), cex.axis=0.8)
     plot(density(phase[,1], n=GridLength), main = title, ylab="Density", xlab = "Calendar Year", ylim=c(0,maxValuey+middleValuey), xlim=c(minValuex, maxValuex), bty='n',lty =1, lwd=2, axes=F, col = pal[1])
     lines(density(phase[,2], n=GridLength), lty =1, lwd=2, col = pal[1])
     
     # abscissa axis
     axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor( minValuex), floor( P1Valuex), floor( middleValuex), floor( P3Valuex), floor( maxValuex)))
     # ordinate axis
     axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 5), round(maxValuey, 5)) )
     
     ## Following phases
     for(i in 2:L) {
       lines(density(phase[,2*i-1], n = GridLength), col=pal[i], lwd=2, lty=1)
       lines(density(phase[,2*i], n = GridLength), col=pal[i], lwd=2, lty=1)
     }
     
     # ## Phase Time Range
     MPTR = MultiPhaseTimeRange(data=data, position_minimum=position_minimum, position_maximum=position_maximum, level=level)
     for (i in 1:L ) { segments(MPTR[i,2], maxValuey+haut[i+1],MPTR[i,3], maxValuey+haut[i+1],lwd=6,col=pal[i]) }
     text(minValuex, maxValuey+haut[2], "Time range", srt =90)

    dev.off()
  }   # end option export
  
  }
}




  ####################################
###   Tempo plot   NEW Version 2017/03   ###


#library(ggplot2)
#library(ggthemes)
#library(scales)

# The tempo plot introduced by T. S. Dye 
# A statistical graphic designed for the archaeological study of rhythms of the long term that embodies a theory of archaeological evidence for the occurrence of events
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence
#' @param count, if TRUE the counting process is a number, otherwise it is a probability
#' @param Gauss if TRUE, the Gaussian approximation of the CI is used
#' @param title title of the graph
#' @param x.label label of the x-axis
#' @param y.label label of the y-axis
#' @param line.types type of the lines drawn of the graph
#' @param plot.wi width size
#' @param plot.ht heigth size
#' @param base.font font of the text
#' @param colors if TRUE, the graph is drawn with colors, otherwise it is drawn in black and white
#' @param out.file the name of the graph (+ extension) that will be saved if chosen. Null by default. 
#' @return a plot 
#' @export
TempoPlot <- function(data, position, level = 0.95 , count = TRUE, Gauss = FALSE, title = "Tempo plot", x.label="Calendar Year", y.label="Cumulative events",
                           line.types=c("solid", "12", "11", "28", "28"), plot.wi = 7, plot.ht = 7, base.font = 11,
                           colors=TRUE, out.file=NULL) {

  # Construction of a new dataset containing the columns corresponding to the phases of interest
  L = length(position)
  
  groupOfDates = matrix(ncol = L, nrow=nrow(data))
  for (i in 1:L) {
    groupOfDates[,i] = data[,position[i]]
  }
  
  # Horizontal axis : min / max
  min = min(apply(groupOfDates,2, min))
  max = max(apply(groupOfDates,2, max))
  
  ## rounding the minimum value down and the maximum value up.
  x = 10^c(0:10)
  if(min!=0){  
    c =0
    for(i in 1:length(x)) { if( abs(min/x[i])>1) {c=c+1}}
    if(c>3){ min = floor(min/x[c-1])*x[c-1]} else {min = floor(min/x[c])*x[c]}
  }
  if(max!=0){     
    d=0
    for(i in 1:length(x)) { if( abs(max/x[i])>1) {d=d+1}}
    if(d>3){ max = ceiling(max/x[d-1])*x[d-1]} else {max = ceiling(max/x[d])*x[d]}
  }
  
  
  t = seq( min, max, length.out = 50*ncol(groupOfDates))
  
  # Calculation of the Bayes estimate and its credible interval
  f= function(x){
    g=ecdf(x)   
    y=g(t) 
    if (count)  y = y * ncol(groupOfDates)
    y
  } 
  F = t( apply( groupOfDates,1,f ) )
  moy = apply(F,2,mean)
  ec = apply(F,2,sd)
  qu = cbind( apply(F,2, quantile, probs  = (1-level)/2 , type = 8) ,  apply(F,2, quantile, probs  = 1-((1-level)/2) , type = 8)  )
  quG = cbind( moy+qnorm(1-(1-level)/2)*ec, moy-qnorm(1-(1-level)/2)*ec )
  
  result = list(t=t, moy=moy, qu=qu, quG=quG)
  
  if (Gauss) {
    result.mat <- cbind(result$moy, result$qu, result$quG)
    colnames(result.mat) <-  c("Bayes estimate", "Credible interval, low", "Credible interval, high", "Gaussian approx., high",
                             "Gaussian approx., low")
  }
  else {
    result.mat <- cbind(result$moy, result$qu)
    colnames(result.mat) <-  c("Bayes estimate", "Credible interval, low", "Credible interval, high")
  }
  
  plot.result <- as.data.frame.table(result.mat)
  colnames(plot.result) <- c("Var1", "Legend", "Count")
  plot.result$Year <- result$t
  
  if (colors)  {
    h <- ggplot(plot.result, aes(x=plot.result$Year, y=plot.result$Count, group=plot.result$Legend, colour=plot.result$Legend))
  }
  else {
    h <- ggplot(plot.result, aes(x=plot.result$Year, y=plot.result$Count, group=plot.result$Legend))
  }
  
  old.theme <- theme_set(theme_bw(base_size=base.font))
  h <- h + theme(legend.title=element_blank())
  h <- h + geom_line(aes(linetype=plot.result$Legend))
  h <- h + scale_linetype_manual(values=line.types)
  
  if (!is.null(x.label)) {
    h <- h + xlab(x.label)
  }
  if (!is.null(y.label)) {
    h <- h + ylab(y.label)
  }
  h <- h + ggtitle(title)
  
  if (!is.null(out.file)) {
    ggsave(filename=out.file, plot = h, height = plot.ht, width = plot.wi)
  }
  
  #old.par <- par(no.readonly = T)
  #dev.new(width = plot.wi, height = plot.ht)
  #par(new)
  print(h)
  #par(old.par)
  #theme_set(old.theme)
}


##############################################
###   Tempo Activity plot   NEW 2016/09   ###

# A statistical graphic designed for the archaeological study of rhythms of the long term that embodies a theory of archaeological evidence for the occurrence of events
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence
#' @param count, if TRUE the counting process is a number, otherwise it is a probability
#' @param title title of the graph
#' @return a plot 
#' @export
TempoActivityPlot <- function(data, position, level = 0.95, count = TRUE, title = "Activity plot") {
  
  # Construction of a new dataset containing the columns corresponding to the phases of interest
  L = length(position)
  
  groupOfDates = matrix(ncol = L, nrow=nrow(data))
  for (i in 1:L) {
    groupOfDates[,i] = data[,position[i]]
  }
  
  min = min(apply(groupOfDates,2, min))
  max = max(apply(groupOfDates,2, max))
  
  x = 10^c(0:10)
  c =0
  for(i in 1:length(x)) { if( abs(min/x[i])>1) {c=c+1}}
  if(c>3){ min = floor(min/x[c-1])*x[c-1]} else {min = floor(min/x[c])*x[c]}
  d=0
  for(i in 1:length(x)) { if( abs(max/x[i])>1) {d=d+1}}
  if(d>3){ max = ceiling(max/x[d-1])*x[d-1]} else {max = ceiling(max/x[d])*x[d]}
  
  t = seq( min, max, length.out = 50*ncol(groupOfDates))
  
  f= function(x){
    g=ecdf(x)   
    y=g(t) 
    if (count)  y = y * ncol(groupOfDates)
    y
  } 
  F = t( apply( groupOfDates,1,f ) )
  moy = apply(F,2,mean)
  ec = apply(F,2,sd)
  
  plot(x<-t[-1], y<-diff(moy)/diff(t), type = "l", xlab="Time", ylab="Incidence", main =title)  
  grid()
}

