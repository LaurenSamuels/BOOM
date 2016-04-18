GetPValue <- function(TE, SE){
    zval <- abs(TE / SE)
    2 * pnorm(zval, lower.tail= FALSE)
}
