# Draw kappa_u from the full conditional
drawku <- function(au, bu, u, R){
  u = c(u);
  
  n <- dim(R)[1];
  
  # Calculate parameters of the gamma distribution
  alpha <- (n-1) / 2 + au;
  # Rate parameter!
  beta  <- bu + t(u) %*% R %*% u / 2;
  return(rgamma(1, alpha, rate = beta));
}

# Draw kappa_v from the full conditional
drawkv <- function(av, bv, u, eta, R) {
  u = c(u);
  eta = c(eta);
  
  n <- dim(R)[1];
  alpha <- n / 2 + av;
  beta  <- bv + t(eta - u) %*% (eta - u) / 2;
  return(rgamma(1, alpha, rate = beta));
}

# Create the block matrix appearing in the full conditional of x.
# Note that before usage in the proposal distribution, the c-vector has to be
# added to the diagonal in the lower right block.
blockMatrix <- function(kv, ku, R) {
  n <- dim(R)[1];
  
  ## Enlarge the dimension, now R is (automatically) only in the
  ## upper left block while the other blocks contain zeros
  Q1 <- R;
  dim(Q1) <- c(2*n, 2*n);
  
  ## Combine different diag matrices
  Q2 <- rbind(cbind(diag.spam(n), -diag.spam(n)), cbind(-diag.spam(n), diag.spam(n)));
   
  Q <- ku * Q1 + kv * Q2;

  return(Q);
}

# Draw from the full conditional of x = u, eta
drawx <- function(x0, e, y, kv, ku) {
  n = length(e);
  
  x0 = c(x0);
  eta0 = x0[(length(x0)/2+1):length(x0)];
  
  e = c(e);
  y = c(y);
  
  cv <- e * exp(eta0);
  b <- y + e * exp(eta0) * (eta0 - 1);
  
  # Create block matrix and add the c-vector
  Q <- blockMatrix(kv, ku, R);
  diagC <- diag.spam(c(rep(0, n), cv));
  Q <- Q + diagC;
  
  b <- c(rep(0, n), b);
  
  return(c(t(rmvnorm.canonical(n = 1, b, Q, memory=list(nnzcolindices=6467)))));
}

# Draw from the proposal distribution of eta
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

# Draw a sample from the full conditional for u
drawu <- function(ku, kv, eta, R){
  eta = c(eta);
  
  n <- dim(R)[1];
  L <- kv * eta;
  Q <- diag(kv, n) + ku * R;
  
  return(c(t(rmvnorm.canonical(n = 1, L, Q, memory=list(nnzcolindices=6467)))));
}