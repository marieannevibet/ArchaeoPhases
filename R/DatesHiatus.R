
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
