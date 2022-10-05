#' Horvitz-Thompson Estimator
#' 
#' HTestimator is used to compute Horvitz-Thompson estimators of the total or of the mean for domains.
#' @param y a vector of size nx1 (n is the sample size) of variables of interest.
#' @param weight a vector of size nx1 of sampling weights.
#' @param domains a vector of size nx1 of domains. If set to NULL, the function returns an estimator for the population.
#' @param mean a boolean variable that indicates if the function has to compute the mean instead of the total.
#' @param Nd a vector of size Dx1 (D is the number of domains) of domains population size for computing the means. If set to null, the function compute the Hajek estimators of the mean. 
#' @return A data frame with 3 columns : Domains, Estimates and Sample Sizes.
#' @examples
#' y <- rnorm(100)
#' dom <- sample(1:5, size = 100, replace = TRUE)
#' Nd <- c(30,40,60,30,20)
#' weight <- sum(Nd) / 100
#' t_d <- HTestimator(y = y, weight = weight, domains = dom, mean = TRUE, Nd = Nd)
#' @export
HTestimator <- function(y, weight, domains = NULL, mean = FALSE, Nd = NULL){
  
  if(mean == FALSE){
    if(is.null(domains)){
      t <- sum(y/weight)
      return(t)
    }
    else{
      nd <- aggregate(y, by = list(domains), length)
      t_d <- aggregate(y*weight, by = list(domains), sum)
      t_d <- cbind(t_d, nd$x)
      colnames(t_d) <- c("Domains", "Totals", "Sample Size")
      return(t_d)
    }
  }
  else{
    if(is.null(domains)){
      if(is.null(Nd)){
        Nd <- sum(weight)
        tbar <- sum(y*weight) / Nd
        return(tbar)
      }
      else{
        tbar <- sum(y*weight) / Nd
        return(tbar)
      }
    }
    else{
      if(is.null(Nd)){
        Nd <- aggregate(weight, by = list(domains), sum)
        nd <- aggregate(y, by = list(domains), length)
        tbar_d <- aggregate(y*weight, by = list(domains), sum)
        tbar_d$x <- tbar_d$x / Nd$x
        tbar_d <- cbind(tbar_d, nd$x)
        colnames(tbar_d) <- c("Domains", "Means", "Sample Size")
        return(tbar_d)
      }
      else{
        nd <- aggregate(y, by = list(domains), length)
        tbar_d <- aggregate(y*weight, by = list(domains), sum)
        tbar_d$x <- tbar_d$x / Nd
        tbar_d <- cbind(tbar_d, nd$x)
        colnames(tbar_d) <- c("Domains", "Means", "Sample Size")
        return(tbar_d)
      }
    }
  }
}

#' Horvitz-Thompson Estimator of the Variance
#' 
#' HTvar is used to compute Horvitz-Thompson estimators of the variance for the total or for the mean for domains.
#' @param y a vector of size nx1 (n is the sample size) of variables of interest.
#' @param Weights a matrix of size nxn of sampling weights.
#' @param domains a vector of size nx1 of domains. If set to NULL, the function returns an estimator for the population.
#' @param mean a boolean variable that indicates if the function has to compute the mean instead of the total.
#' @param Nd a vector of size Dx1 (D is the number of domains) of domains population size for computing the means.
#' @return A data frame with 3 columns : Domains, Variances and Sample Sizes.
#' @examples
#' n <- 100
#' y <- rnorm(n)
#' dom <- sample(1:5, size = n, replace = TRUE)
#' Nd <- c(30,40,60,30,20)
#' N <- sum(Nd)
#' pik <- n/N
#' pikl <- (n *(n-1)) / (N *(N-1))
#' Weights <- matrix(0, n, n)
#' diag(Weights) <- pik
#' Weights[upper.tri(Weights)] = PI[lower.tri(Weights)]  <- pikl
#' t_d_var <- HTvar(y = y, Weights = Weights, domains = dom, mean = TRUE, Nd = Nd)
#' @export
HTVar <- function(y, Weights, domains = NULL, mean = FALSE, Nd, method = 1){
  #1. No domains
  if(is.null(domains)){
    #1.1 Total
    if(mean == FALSE){
      if(length(y) == 1){
        return((y/Pi)^2 * (1 - Pi))
      }
      else{
        return(varHT(y, Weights, method = method))
      }
    }
    #1.2 Mean
    else{
      if(length(y) == 1){
        return((y/Pi)^2 * (1 - Pi) / (Nd^2) )
      }
      else{
        return(varHT(y, Weights, method = method) / (Nd^2))
      }
    }
  }
  #2. Domains
  else{
    index_domains <- unique(domains)
    D <- length(index_domains)
    nd <- aggregate(y, by = list(domains), length)
    #2.1 Total
    if(mean == FALSE){
      t_d_var <- numeric(D)
      for(d in 1:D){
        t_d_var[d] <- varHT(y[domains == index_domains[d]], 
                         Weights[domains == index_domains[d],domains == index_domains[d]],
                         method = method)
      }
      t_d_var <- data.frame(Domains = index_domains, 
                            Variances = as.numeric(t_d_var), 
                            Samples_Sizes = nd$x)
      return(t_d_var[order(index_domains),])
    }
    #2.2 Means
    else{
      t_d_var <- numeric(D)
      for(d in 1:D){
        t_d_var[d] <- varHT(y[domains == index_domains[d]], 
                            Weights[domains == index_domains[d],domains == index_domains[d]],
                            method = method) / (Nd[d]^2)
      }
      t_d_var <- data.frame(Domains = index_domains, 
                            Variances = as.numeric(t_d_var), 
                            Samples_Sizes = nd$x)
      return(t_d_var[order(index_domains),])
    }
  }
}