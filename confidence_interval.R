conf.interval <- function(p,n,z){
  diff <- (sqrt((p*(1-p))/n)*z) #find difference
  plus <- p + diff
  minus <- p - diff
  ci <- paste("[",minus,", ", plus, "]", sep = "")
  return(ci)
}
conf.interval(0.71,21,1.465)
