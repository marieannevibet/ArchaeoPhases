
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
          plot(density(a_chain, n=GridLength), main = title, xlab = "Calendar Year", axes = FALSE, ylim=c(0,max(density(a_chain, n=GridLength)$y) + step))
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
          plot(density(a_chain, n=GridLength), main = title, xlab = "Calendar Year", axes = FALSE, ylim=c(0,max(density(a_chain, n=GridLength)$y) + step))
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