#Draw a sample from the full conditional for u
drawU <- function(kappaV, eta, kappaU,R){
  n=dim(R)[1];
  sigma <-solve(  kappaV*diag(n)+kappaU*R);
 
  mu=vector(mode="numeric",length=n)
  drawU <- mvrnorm(1,mu,sigma); 
}