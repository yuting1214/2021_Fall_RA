function(y,yfit,family=binomial(),variance=NULL)
{
  # Calculate the residual for given observed y and its fitted value yfit: 
  # the length between y and yfit along the quardratic variance function:
  #           V(mu) = v2*mu^2+v1*mu+v0
  qvresidual<-function(y,yfit,v2,v1)
  {
    vpa <- 2*v2*yfit+v1
    svpa2 <- sqrt(1+vpa*vpa)
    
    vpb <- 2*v2*y+v1
    svpb2 <- sqrt(1+vpb*vpb)
    
    vr <- (log((vpb+svpb2)/(vpa+svpa2))+vpb*svpb2-vpa*svpa2)/(4*v2)
    vr
  }
  
  if( is.character(family) ) {
    cf <- family
    family <- get(family, mode="function", envir=parent.frame())
  } else
    cf <- family$family

  if( pmatch("Negative Binomial",cf,nomatch=F) )
  {
    theta <- as.numeric(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","", cf, perl=T))
    cf <- "negative.binomial"
  } else if( pmatch("Tweedie",cf,nomatch=F) )
  {
    dv <- Deriv(family$variance,"mu")
    theta <- dv(1)
  }
  
  if( is.null(variance) )
  {
    switch(cf,
      binomial={DFUN<-function(x) qvresidual(x[1],x[2],-1,1)}, # Need modify for Y~Bin(n,p)
      gaussian={DFUN<-function(x) x[1]-x[2]},
      Gamma={DFUN<-function(x) qvresidual(x[1],x[2],1,0)},
      negative.binomial={DFUN<-function(x) qvresidual(x[1],x[2],1/theta,1)},
      poisson={DFUN<-function(x) x[1]-x[2]},
      quasibinomial={DFUN<-function(x) qvresidual(x[1],x[2],-1,1)},
      quasipoisson={DFUN<-function(x) x[1]-x[2]},
      inverse.gaussian={
        DFUN<-function(x) integrate(function(mu){sqrt(1+9*mu^4)},x[1],x[2])$value},
      Tweedie={ # var.power: 0, 1, (1,2), 2, >2
        if( (theta==0)|(theta==1) )
          DFUN<-function(x) x[1]-x[2]
        else if( theta==2 ) 
          DFUN<-function(x) qvresidual(x[1],x[2],1,0)
        else
          DFUN<-function(x) integrate(function(mu){sqrt(1+theta^2*mu^(2*theta-2))},x[1],x[2])$value},
      quasi={  # variance for quasi: "constant","mu(1-mu)","mu","mu^2","mu^3", or other
        if( (family$varfun=="constant")|(family$varfun=="mu") )
          DFUN <- function(x) x[1]-x[2]
        else if( family$varfun=="mu(1-mu)" )
          DFUN<-function(x) qvresidual(x[1],x[2],-1,1)
        else if( family$varfun=="mu^2" )
          DFUN<-function(x) qvresidual(x[1],x[2],1,0)
        else
          DFUN<-function(x) integrate(function(mu){sqrt(1+Deriv(family$variance,"mu")^2)},x[1],x[2])$value})
  }
  else
    DFUN<-function(x) integrate(function(mu){sqrt(1+Deriv(variance,"mu")^2)},x[1],x[2])$value
  
  vresidual <- apply(cbind(y,yfit),1,DFUN)
}
