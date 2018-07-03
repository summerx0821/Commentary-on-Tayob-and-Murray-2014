
####################################################################################################################################
#  PROGRAM NAME   : Function.proportion.events.captured.R
#
#  AUTHOR         : 3JUL2018 by Meng Xia
#
#  EXECUTED UNDER : R 3.4.1 (2017-06-30)
#
#  FUNCTION NAME  : prop.events.captured
#
#  DESCRIPTION    : To calculate the proportion of events captured for a given
#                   s      : overall follow up time
#                   lambda : control group recurrent event rate
#                   a      : space follow-up windows every a unit apart
#
#  EXAMPLE        : prop.events.captured(s      = 48,
#                                        lambda = 1/3, 
#                                        a      = 1.5)
#
#  EXAMPLE OUTPUT : Time difference of 4.705428 mins
#                   [1] "The proportion of events captured is 0.8"
#
####################################################################################################################################

rm(list=ls())

prop.events.captured <- function(s, lambda, a){
  
  time1 <- Sys.time()
  
  #Compute the number of windows
  b=ceiling(s/a)
  
  # Function to compute the average number of events missed
  sumsum <- function(J){ 
    ww <- 1:b
    jj <- 2:J
    
    prob_w <- function (wj){
      w <- wj[1]
      j <- wj[2]
      
      #Function to calculate pdf_Gamma(j,lambda) (r) * cdf_Gamma(l, lambda) (s-r)
      prob_RR <- function(j, l, w){
        f <- function(R_j){dgamma(R_j, shape=j, rate=lambda) * pgamma(s-R_j, shape = l, rate = lambda)}
        integral <- integrate(f, lower = 0, upper = min(a*w, s))$value
        return(integral)
      }
      
      #Function to calculate pdf_Gamma(j-1, lambda) (r) * pdf_Exp(lambda) (g) * cdf_Gamma(l, lambda) (s-r-g)
      prob_RGR <- function(j, l, w){
        InnerFunc <- function(R_j_1, G_j) 
        {dgamma(R_j_1, shape=j-1, rate=lambda) * dexp(G_j, rate = lambda) * pgamma(s-R_j_1-G_j, shape = l, rate = lambda)}
        InnerIntegral <- function(R_j_1) 
        { sapply(R_j_1, function(x) { integrate(function(G_j) InnerFunc(x, G_j), 0, min(a*w, s)-x)$value }) }
        OuterIntegral <- integrate(InnerIntegral, 0, (w-1)*a)$value
        return(OuterIntegral)
      }
      
      #Calculate numerator and denominator
      if (w == 1) {
        top1 <- prob_RR(j=j, l=J-j, w=w)
        top2 <- prob_RR(j=j, l=J-j+1, w=w)
      }
      if (w > 1) {
        top1 <- prob_RR(j=j, l=J-j, w=w)-prob_RGR(j=j, l=J-j, w=w)
        top2 <- prob_RR(j=j, l=J-j+1, w=w)-prob_RGR(j=j, l=J-j+1, w=w)
      }
      bottom <- lambda^J*exp(-lambda*s)*s^J/gamma(J+1)
      
      return((top1-top2)/bottom)
    }
    
    element <- apply(expand.grid(ww,jj),1,prob_w)
    return(sum(element))
  }
  
  prop=array()
  for (J in 1:100){
    if (J==1) {prop[J] <- (s*lambda)^J*exp(-lambda*s)/factorial(J)/(1-(s*lambda)^0*exp(-lambda*s)/factorial(0))}
    if (J>1) {prop[J] <- (J-sumsum(J))/J*(s*lambda)^J*exp(-lambda*s)/factorial(J)/(1-(s*lambda)^0*exp(-lambda*s)/factorial(0))}
  }
  
  #Output
  print(Sys.time()-time1)
  return(paste("The proportion of events captured is",round(sum(prop), digits = 2)))
}

prop.events.captured(s      = 48,
                     lambda = 1/3, 
                     a      = 1.5)


