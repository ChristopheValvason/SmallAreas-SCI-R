#'Bonferroni Correction for Simultaneous Confidence Interval
#'
#'BonferroniSCI is a function to compute a simultaneous confidence interval (SCI) using the Bonferroni correction.
#'@param theta a vector of estimates of size Dx1.
#'@param variances a vector of variances for theta of size Dx1.
#'@param alpha a numeric value, 0 < alpha < 1, indicating the significance level.
#'@param n a numeric value for the overall sample size.
#'@param sidak a boolean value that indicates if the Sidak correction must be computed instead of the Bonferroni correction.
#'@param t_dist a boolean value that indicates if the t-distribution is used instead of the normal distribution for computing the quantiles.
#'@param domains an optional vectors of names for the domains of size Dx1; must be in the same order as theta.
#'@export
BonferroniSCI <- function(theta, variances, alpha, n, sidak = FALSE, t_dist = FALSE, domains = NULL){
  m <- length(theta)
  if(sidak == FALSE){
    alpha_m <- alpha/(2*m)
    if(t_dist == FALSE){c_alpha <- qnorm(1-alpha_m, mean = 0, sd = 1)}
    else{c_alpha  <- qt(1-alpha_m, df = n - m)}
  
    L <- theta - (c_alpha * sqrt(variances))
    U <- theta + (c_alpha * sqrt(variances))
    if(is.null(domains)){ return(data.frame(Lower = L, Upper= U))}
    else{return(data.frame(Domains = domains, Lower = L, Upper= U))}
  }
  else{
    alpha_m <- 1 - (1-alpha/2)^(1/m)
    if(t_dist == FALSE){c_alpha <- qnorm(1-alpha_m, mean = 0, sd = 1)}
    else{c_alpha  <- qt(1-alpha_m, df = n - m)}
    
    L <- theta - (c_alpha * sqrt(variances))
    U <- theta + (c_alpha * sqrt(variances))
    if(is.null(domains)){ return(data.frame(Lower = L, Upper= U))}
    else{return(data.frame(Domains = domains, Lower = L, Upper= U))}
  }
}