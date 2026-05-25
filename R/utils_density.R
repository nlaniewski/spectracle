peak.bounds <- function(x, adjust = 1, height.threshold = 0.01, plot = F){
  d <- stats::density(x, adjust = adjust)
  pv <- diff(sign(diff(d$y)))
  peaks.i <- which(pv == -2) + 1
  peak.major.i <- peaks.i[which.max(d$y[peaks.i])]
  ##
  threshold.value <- d$y[peak.major.i] * height.threshold
  ##
  bounds.i <- sapply(c('left','right'), function(side){
    i <- peak.major.i
    while(d$y[i] > threshold.value){
      i <- get(if(side=='left'){"-"}else if(side=='right'){"+"})(i,1)
      if(pv[i] != 2){
        i
      }else{
        break
      }
    }
    return(i)
  })
  ##
  if(plot){
    plot(d)
    graphics::abline(v = d$x[peak.major.i], col = "red")
    graphics::abline(v = d$x[bounds.i], col = "red", lty = "dashed")
  }
  ##
  return(d$x[bounds.i])
}
