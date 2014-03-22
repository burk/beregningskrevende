GISampler <- function(intialValues,n){
  #initialValues is a vector of size 4 contaiong intial values for y, kappav, eta and u.
    for(i in 1:n){
      kappaUTemp = drawKappaU();
      kappaVtemp = drawKappaV();
      uTemp= drawU();
      etaTemp=drawEta();
    }
}