data(Germany);

g <- system.file("demodata/germany.graph", package="INLA");
source(system.file("demodata/Bym-map.R", package="INLA"));
summary(Germany);

Oral <- cbind(Oral,region.struct=1:544, region=1:544);

formula <- Y ~ f(region.struct, model="besag", graph=g, hyper=list(prec=list(param=c(1,0.01))), constr=F) + f(region, model="iid", hyper=list(prec=list(param=c(1,0.01)))) - 1

## just make a duplicated column
Germany = cbind(Germany, region.struct=Germany$region);

result = inla(formula, family="poisson", data=Germany, E=E);

Bym.map(result$summary.random$region.struct$mean)
