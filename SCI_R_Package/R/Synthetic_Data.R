#'Generation of synthethic data from a nested regression model
#'
#'NERMData is a function to generate synthetic data from a nested error regression model.
#'@param N is the size of the population
#'@param dom is a N * 1 vector of domains that partioned the population
#'@param N_x is the number of auxilliary variables
#'@param Dist_x is the distributions of the auxilliary variables in the form r + distribution
#'@param param_x is a list of the parameters for each distribution in Dist_x
#'@param Beta is the N_x * 1 vector of beta's for the model
#'@param s_u is the standard error for the random effects
#'@param s_e is the standard error for the error terms
#'@param intercept is a boolean that indicates if the model has an intercept. By default it is set to TRUE.
#'@export
NERMData <- function(N, dom, N_x, Dist_x, param_x, Beta, s_u, s_e, intercept = TRUE){
  X <- matrix(0, nrow = N, ncol = N_x)
  for(i in 1:N_x){
    X[,i] <- Dist_x[[i]](N, param[[i]])
  }
  
  if(isTRUE(intercept)){
    X <- cbind(rep(1, N), X)
  }
  
  D <- length(levels(factor(dom)))
  u <- rnorm(D, mean = 0, sd = s_u)
  e <- rnorm(N, mean = 0, sd = s_e)
  
  y <- X%*%Beta + u[dom] + e
  
  return(cbind(y, dom))
}