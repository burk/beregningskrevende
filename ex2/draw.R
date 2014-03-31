#Draws a sample from the full-conditional for kappa_u
drawku <- function(au, bu, u, R){
  u = c(u);
  
  n <- dim(R)[1];
  alpha <- (n-1) / 2 + au;
  beta  <- bu + t(u) %*% R %*% u / 2;
  return(rgamma(1, alpha, rate = beta));
}

drawkv <- function(av, bv, u, eta, R) {
  u = c(u);
  eta = c(eta);
  
  n <- dim(R)[1];
  alpha <- n / 2 + av;
  beta  <- bv + t(eta - u) %*% (eta - u) / 2;
  return(rgamma(1, alpha, rate = beta));
}

blockMatrix <- function(kv, ku, R) {
  n <- dim(R)[1];
  ## Enlarge the dimension, now R is (automatically) only in the
  ## upper left block while the other blocks contain zeros
  Q1 <- R;
  dim(Q1) <- c(2*n, 2*n);
  
  ## Check via
  #par(mfrow=c(1,3));
  #display(R, cex=2);
  #title("Adjacency matrix");
  #display(Q1,  cex=2);
  #title("Q1");
  
  ## Combine different diag matrices
  Q2 <- rbind(cbind(diag.spam(n), -diag.spam(n)), cbind(-diag.spam(n), diag.spam(n)));
  
  ## generate a matrix which will later contain the entries of diag(c)
  ## by now the entries are 1.
  #diagC <- as.spam(diag.spam(c(rep(0,n), rep(1,n))));
  
  ## The entries can be changed (updated) by assigning 
  ## diagC@entries a new vector of length 544. This will be 
  ## necessary at each step of your sampler as the values of
  ## c will change.
  ## For example,
  #diagC@entries <- c(rep(0, n), c(cv));
  
  #diagC <- diag.spam(c(rep(0, n), cv));
  
  Q <- ku * Q1 + kv * Q2;
  #display(Q, cex=2);
  #title("Block matrix");
  
  return(Q);
}

drawx <- function(x0, e, y, kv, ku) {
  n = length(e);
  
  x0 = c(x0);
  eta0 = x0[(length(x0)/2+1):length(x0)];
  
  e = c(e);
  y = c(y);
  
  cv <- e * exp(eta0);
  b <- y + e * exp(eta0) * (eta0 - 1);
  
  Q <- blockMatrix(kv, ku, R);
  diagC <- diag.spam(c(rep(0, n), cv));
  Q <- Q + diagC;
  
  b <- c(rep(0, n), b);
  
  return(c(t(rmvnorm.canonical(n = 1, b, Q, memory=list(nnzcolindices=6467)))));
}

draweta <- function(eta0, e, y, u, kv) {
  eta0 = c(eta0);
  e = c(e);
  y = c(y);
  u = c(u);
  
  cv <- e * exp(eta0);
  b <- y + e * exp(eta0) * (eta0 - 1);
  L <- b + kv * u;
  Q <- diag.spam(kv + cv);
  return(c(t(rmvnorm.canonical(n = 1, L, Q, memory=list(nnzcolindices=6467)))));
}

#Draw a sample from the full conditional for u
drawu <- function(ku, kv, eta, R){
  eta = c(eta);
  
  n <- dim(R)[1];
  #sigma <- solve(diag(kv, n) + ku * R);
  L <- kv * eta;
  Q <- diag(kv, n) + ku * R;
  
  #mu <- rep(0, n);
  return(c(t(rmvnorm.canonical(n = 1, L, Q, memory=list(nnzcolindices=6467)))));
  #return(mvrnorm(1, mu, sigma));
}