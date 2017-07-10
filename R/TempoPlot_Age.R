
TempoPlot_Age <- function(data, position, level = 0.95 , count = TRUE, Gauss = FALSE, title = "Tempo plot", x.label="Calibrated year BP", y.label="Cumulative events",
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
  
  #Modified to obtain the graphic in Year BP 
  
  result = list(t=1950-t, moy=(moy), qu=qu, quG=quG)
  
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

