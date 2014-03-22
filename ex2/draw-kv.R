drawkv <- function(alpha_u, beta_v, kv, eta, u, R, n = 544) {
  alpha <- (n-1) / 2 + alpha_u;
  beta  <- beta_u * kv + kv * t(u) %*% R %*% u / 2;
  return(rgamma(1, alpha, beta));
}

draweta <- function(eta0, e, y, u, kv) {
  c <- e * exp(eta0);
  b <- y + e * exp(eta0);
  mu <- diag(1 / (kv + c)) %*% (b + kv * u);
  sigma <- diag(1 / (kv + c));
  return(mvrnorm(n = 1, mu, sigma));
}