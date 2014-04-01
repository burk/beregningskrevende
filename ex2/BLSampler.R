BLSampler <- function(R, e, y, nsamples, burnin){
  nsamples <- nsamples + burnin;
  
  accepted <- 0;
  rejected <- 0;
  n  <- dim(R)[1];
  au <- 1;
  av <- 1;
  bu <- 0.01;
  bv <- 0.01;
  
  ku <- matrix(nrow = 1, ncol = nsamples);
  kv <- matrix(nrow = 1, ncol = nsamples);
  
  x <- matrix(data=1, nrow = 2*n, ncol = nsamples);
  
  ku[,1]  <- 0.5;
  kv[,1]  <- 0.5;
  x[,1] <- 1;
  
  for(i in 2:nsamples) {
    ku[,i] <- drawku(au, bu, x[1:n,i-1], R);
    kv[,i] <- drawkv(av, bv, x[1:n,i-1], x[(n+1):(2*n),i-1], R);
    xtmp   <- drawx(x[,i-1], e, y, kv[,i], ku[,i]);
    
    #print(xtmp);
    print(i);
    
    alpha <- min(0,
                 + xProb(xtmp, y, kv[,i], ku[,i], e)
                 - xProb(x[,i-1], y, kv[,i], ku[,i], e)
                 + xQ(x[,i-1], xtmp, y, kv[,i], ku[,i], e)
                 - xQ(xtmp, x[,i-1], y, kv[,i], ku[,i], e)
    );
    
    
    #print(alpha);
    
    if (runif(1) < exp(alpha)) {
      print("ACCEPT");
      accepted <- accepted + 1;
      x[,i] <- xtmp;
    }
    else {
      print("REJECT");
      rejected <- rejected + 1;
      x[,i] <- x[,i-1];
    }
  }
  
  cat("TOTAL ACCEPTED:", accepted, " (burnin:", burnin, ")\n");
  cat("TOTAL REJECTED:", rejected, "\n");
  
  return(list(
    "kv" = kv[(burnin+1):nsamples],
    "ku" = ku[(burnin+1):nsamples],
    "x"  = x[,(burnin+1):nsamples],
    "u"  = x[1:544,(burnin+1):nsamples],
    "eta"= x[545:1088,(burnin+1):nsamples],
    "v"  = x[545:1088,(burnin+1):nsamples] - x[1:544,(burnin+1):nsamples]
  ));
}