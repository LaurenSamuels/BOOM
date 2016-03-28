GetConfInt <- function(est, SE, conf.level) {
    alpha <- 1 - conf.level
    zval <- qnorm(1 - alpha/2)
    c(est - zval * SE, est + zval * SE)
}
