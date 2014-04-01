library(MASS);
library(spam);

GISampler <- function(R, e, y, nsamples, burnin){
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
  
  u   <- matrix(data=1, nrow = n, ncol = nsamples);
  eta <- matrix(data=1, nrow = n, ncol = nsamples);
  
  ku[,1]  <- 20;
  kv[,1]  <- 400;
  u[,1]   <- 0;
  eta[,1] <- 1;

  for(i in 2:nsamples) {
    ku[,i] <- drawku(au, bu, u[,i-1], R);
    kv[,i] <- drawkv(av, bv, u[,i-1], eta[,i-1], R);
    u[,i]  <- drawu(ku[,i], kv[,i], eta[,i-1], R);
    etatmp <- draweta(eta[,i-1], e, y, u[,i], kv[,i]);
    
    #print(etatmp);
    
    alpha <- min(0,
                + etaProb(etatmp, y, kv[,i], ku[,i], u[,i], e)
                - etaProb(eta[,i-1], y, kv[,i], ku[,i], u[,i], e)
                + etaQ(eta[,i-1], etatmp, y, kv[,i], ku[,i], u[,i], e)
                - etaQ(etatmp, eta[,i-1], y, kv[,i], ku[,i], u[,i], e)
    );
    
    #print(alpha);
    
    print(i);
    if (runif(1) < exp(alpha)) {
      print("ACCEPT");
      accepted <- accepted + 1;
      eta[,i] <- etatmp;
    }
    else {
      print("REJECT");
      rejected <- rejected + 1;
      eta[,i] <- eta[,i-1];
    }
  }
  
  cat("TOTAL ACCEPTED:", accepted, " (burnin:", burnin, ")\n");
  cat("TOTAL REJECTED:", rejected, "\n");
  
  return(list(
    "kv"  = kv[(burnin+1):nsamples],
    "ku"  = ku[(burnin+1):nsamples],
    "u"   = u[,(burnin+1):nsamples],
    "eta" = eta[,(burnin+1):nsamples],
    "v"   = eta[,(burnin+1):nsamples] - u[,(burnin+1):nsamples],
    "accepted" = accepted,
    "rejected" = rejected
  ));
}