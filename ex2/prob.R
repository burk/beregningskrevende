library(spam);
library(MASS);

etaProb <- function(eta, y, kv, ku, u, e) {
  eta = c(eta);
  y = c(y);
  n = length(eta);
  
  return(-0.5 * kv * t(eta) %*% eta + kv * t(eta) %*% u + t(eta) %*% y - t(exp(eta)) %*% e);
}

xProb <- function(x, y, kv, ku, e) {
  x = c(x);
  y = c(y);
  n = length(x) / 2;
  
  Q <- blockMatrix(kv, ku, R);
  
  return(-0.5 * t(x) %*% Q %*% x + t(x) %*% c(rep(0,n), y) - exp(x) %*% c(rep(0,n), e));
}

xQ <- function(x, x0, y, kv, ku, e) {
  x = c(x);
  y = c(y);
  n = length(x) / 2;
  
  Q <- blockMatrix(kv, ku, R);
  
  eta0 = x0[(length(x0)/2+1):length(x0)];
  cv <- e * exp(eta0);
  b  <- y + e * exp(eta0) * (eta0 - 1);
  
  L <- c(rep(0, n), b);
  
  diagC <- diag.spam(c(rep(0, n), cv));
  Q <- Q + diagC;
  
  return(dmvnorm.canonical(x, L, Q, log=TRUE));
  
  #return(-0.5 * t(eta) %*% (diag(kv + cv)) %*% eta + t(eta) %*% (kv * u + b));
}

etaQ <- function(eta, eta0, y, kv, ku, u, e) {
  eta0 = c(eta0);
  eta = c(eta);
  u = c(u);
  y = c(y)
  
  cv <- e * exp(eta0);
  b <- y + e * exp(eta0) * (eta0 - 1);
  L <- b + kv * u;
  Q <- diag.spam(kv + c(cv));
  
  return(dmvnorm.canonical(eta, L, Q, log=TRUE));
  
  #return(-0.5 * t(eta) %*% (diag(kv + cv)) %*% eta + t(eta) %*% (kv * u + b));
}

## Function for the density function of a multivariate normal distribution
## in canonical representation. 
##
## Arguments
## x : is the vector where we would like to evaluate the density
## b : b vector used in the canonical parameterisation of the normal distribution
## Q : precision matrix (with full rank)
## log : boolean indicating whether the log density should be returned
## memory : argument passed to the spam package (for the exercise this does not need to be changed).

dmvnorm.canonical <- function(x, b, Q, log=TRUE, memory=list(nnzcolindices=6467)){
  
  # some checks
  if (length(x) != NCOL(Q)) {
    stop("x and Q have non-conforming size")
  }
  if (length(b) != NROW(Q)) {
    stop("b and Q have non-conforming size")
  }
  # compute the log determinant
  logdet <- as.numeric(determinant(Q, logarithm=TRUE, memory=memory)$modulus)
  # get the mean
  mu <- solve.spam(Q, b, memory=memory)
  xmu <- (x-mu)
  # get the log-density
  logdens <- (- length(x) * log(2*pi) + logdet - t(xmu)%*%Q%*%xmu)/2
  
  if(log)
    return(logdens)
  exp(logdens)
}
