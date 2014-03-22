#Draws a sample from the full-conditional for kappa_u
drawKappaU <- function(alpha_u, beta_u,u,R){
  alpha <- (n-1)/2 + alpha_u;
  beta <- beta_u - (1/2)*t(u)%*%R%*%u
  drawKappaU <- rgamma(1, alpha,rate=beta);
  return(drawKappaU);
}