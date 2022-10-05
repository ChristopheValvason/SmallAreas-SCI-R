#' Max-Type bootstrap 
#' 

MaxTypeBoot <- function(y, PI, N, R, estimator, theta, method, mean){
  Smax <- numeric(R)
  for(b in 1:R){
    if(method == "Gross"){
      B <- bootGross(y, PI, N)
    }
    else{
      B <- bootBooth(y, PI, N)
    }
    if(estimator == "HorvitzThompson"){
      theta_hat <- HTestimator(y = B$sample[,1], weight = diag(B$PI), domains = B$sample[,2] , mean = mean)
      Nd <- theta_hat[,3]
      theta_hat_var <- HTVar(y = B$sample[,1], Weights = B$PI, domains = B$sample[,2], mean = mean, Nd = Nd)
      S0 <- (theta_hat[,2] - theta)/ theta_hat_var[,2]
    }
    Smax[b] <- max(abs(S0))
  }
  return(Smax)
}