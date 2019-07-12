RPM <- function(counts) { 
  RPM <-apply(counts,2,function(x) x*1e6/sum(x)) 
  return(RPM)
  }