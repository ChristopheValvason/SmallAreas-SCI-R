#' Gross Bootstrap for finite population and simple random sampling 
#' 
#' bootGross is used to obtain a bootstrap replicate of a sample according to the method developped in Gross(1981).
#' @param y a vector or a matrix of sampled data
#' @param PI a matrix of first and second order inclusion probabilities
#' @param N the size of the population
#' @return a list with a bootstrap sample and the matrix of inclusion probabilities
#' @export
bootGross <- function(y, PI, N){
  #sample size
  if(is.vector(y)){
    n <- length(y)
  }
  else{
    n <- dim(y)[1]
  }
  #sampling rate
  f <- n/N
  #inclusion probabilities
  pik <- diag(PI)
  pikl <- PI[1,2]
  
  #reconstruction of a population
  if(is.vector(y)){
    U_star <- rep(y, each = round(1/f))
    N_star <- length(U_star)
  }
  else{
    n_col <- dim(y)[2]
    U_star <- matrix(rep(y, each = round(1/f)), n*round(1/f), n_col)
    N_star <- dim(U_star)[1]
  }
  PI_star <- diag(rep(pik, each = round(1/f)))
  PI_star[upper.tri(PI_star)] <- pikl
  PI_star[lower.tri(PI_star)] <- pikl
  
  #bootstrap sample
  s_star <- sampling::srswor(n, N_star)
  if(is.vector(y)){
    y_star_s <- U_star[s_star == 1]
  }
  else{
    y_star_s <- U_star[s_star == 1, ]
  }
  PI_star_s <- PI_star[s_star == 1, s_star == 1]
  
  return(list(sample = y_star_s, PI = PI_star_s))
}

#'Booth and al. (1994) bootstrap for finite population and simple random sampling
#'
#'bootBooth is used to obtain a bootstrap replicate of a sample according to the method developped in Booth and al. (1994).
#' @return a list with a bootstrap sample and the matrix of inclusion probabilities
bootBooth <- function(y, PI, N){
  #sample size
  if(is.vector(y)){
    n <- length(y)
  }
  else{
    n <- dim(y)[1]
  }
  #sampling rate
  f <- n/N
  #inclusion probabilities
  pik <- diag(PI)
  pikl <- PI[1,2]
  
  p <- round(1/f - 0.5)
  r <- N - (n*p)
  r_s <- sample(1:n, size = r, replace = FALSE)
  
  if(is.vector(y)){
    U_star <- c(rep(y, each = p), y[r_s])
    N_Star <- length(U_star)
  }
  else{
    n_col <- dim(y)[2]
    U_star <- rbind(matrix(rep(y, each = round(1/f)), n*round(1/f), n_col), y[r_s,])
    N_star <- dim(U_star)[1]
  }
  PI_star <- diag(c(rep(pik, each = p), pik[r_s]))
  PI_star[upper.tri(PI_star)] <- pikl
  PI_star[lower.tri(PI_star)] <- pikl
  
  #bootstrap sample
  s_star <- sampling::srswor(n, N_star)
  if(is.vector(y)){
    y_star_s <- U_star[s_star == 1]
  }
  else{
    y_star_s <- U_star[s_star == 1, ]
  }
  PI_star_s <- PI_star[s_star == 1, s_star == 1]
  
  return(list(sample = y_star_s, PI = PI_star_s))
}