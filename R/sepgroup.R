#' Title
#'
#' @param x
#' @param threshold
#'
#' @returns generate group according to threshold
#' @export
#'
#' @examples
Sepgroup <- function(x,threshold){
  if(length(x)==1){return(x)}
  group = c()
  x0 = x
  i = 1
  while(i <= length(x)){
    diff = x[i]-x
    loc = which(abs(diff) <= threshold,x)
    set = which(x0 %in% unique(x[loc]))
    group = append(group,list(set))
    x = x[-loc]
    if(i==0){break}
  }
  return(group)
}
