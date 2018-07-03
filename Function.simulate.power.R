

####################################################################################################################################
#  PROGRAM NAME   : Function.simulate.power.R
#
#  AUTHOR         : 3JUL2018 by Meng Xia
#
#  EXECUTED UNDER : R 3.4.1 (2017-06-30)
#
#  FUNCTION NAME  : simulate.power
#
#  DESCRIPTION    : To calculate the power by simulation
#
#  EXAMPLE        : simulate.power(
#                     n.sim=1000,            #Number of simulations
#                     prop.base=0.7,         #Proportion of patients enrolled at baseline
#                     study.time=48,         #Overall study time 
#                     accrual=24,            #Accrual time
#                     n1=100,                #Sample size in control group
#                     n2=100,                #Sample size in treatment group
#                     recurrent.rate=1/3,    #Recurrent event rate in control group
#                     terminal.haz=1/1000,   #Terminal event hazard in control group; taking very small number is like no terminal events
#                     Tau=12,                #Length of follow-up windows
#                     ratio=0.75,            #ratio = recurrent.rate(terminal.haz) in treatment / recurrent.rate(terminal.haz) in control
#                     rho1=0,                #Correlation between recurrent events
#                     rho2=0,                #Correlation between recurrent events and terminal events
#                     sp=1.5,                #Space between follow-up windows
#                     seeds=21
#                   )
#
#  EXAMPLE OUTPUT : [1] "Power is 0.975"
#
####################################################################################################################################

rm(list=ls())

simulate.power <- function(
  n.sim=1000,            #Number of simulations
  prop.base=0.7,         #Proportion of patients enrolled at baseline
  study.time,            #Overall study time 
  accrual,               #Accrual time
  n1,                    #Sample size in control group
  n2,                    #Sample size in treatment group
  recurrent.rate,        #Recurrent event rate in control group
  terminal.haz,          #Terminal event hazard in control group; taking very small number is like no terminal events
  Tau,                   #Length of follow-up windows
  ratio,                 #ratio = recurrent.rate(terminal.haz) in treatment / recurrent.rate(terminal.haz) in control
  rho1=0,                #Correlation between recurrent events
  rho2=0,                #Correlation between recurrent events and terminal events
  sp,                    #Space between follow-up windows
  seeds=123
){
  #Load packages
  require(survival)
  require(MASS)
  require(gdata)
  
  #A cumulative sum function:
  sum_function=function(X)
  {
    temp=array(X,c(length(X),length(X)))
    lowerTriangle(temp,diag=F)=0
    apply(temp,2,sum)
  }
  
  #Function to generate data of recurrent and terminal events
  generate_data=function(n,lambda_1,lambda_2,rho_1=0,rho_2=0)
  {
    use_dataset=1
    C1=rbinom(n,1,p)
    C2=runif(n,min=A-A_star,max=A)
    C3=C2*C1 + A*(1-C1) #censoring time of recurrent and terminating event- end of study
    
    J=200 #number of recurrent events simulated per subject
    sigma=diag(rep(1,J+1))+ (1-diag(rep(1,J+1)))*rho_1 #set up covariance matrix
    sigma[J+1,1:J]=rho_2
    sigma[1:J,J+1]=rho_2
    
    Y=mvrnorm(n,mu=rep(0,J+1),sigma) #simulate multivariate normal distribution
    U=apply(Y,c(1,2),pnorm) #transform to correlated uniform distributions
    U_Z=U[,1:J]
    U_T=U[,J+1]
    Z=-log(U_Z)/lambda_1 #transform to correlated exponential distributions for recurrent events
    T=-log(U_T)/lambda_2 #transform to correlated exponential distribution for terminating event
    
    X=apply(cbind(C3,T),1,min)
    delta=as.numeric(T<=C3)
    
    Z_star=t(apply(Z,1,sum_function)) #create variable with time from baseline to each recurrent event
    if(sum(Z_star[,J]<X)>0){use_dataset=0}
    Z_star[Z_star>X]=NA
    list(X=X, delta=delta,Z_star=Z_star,use_dataset=use_dataset)
  }
  
  #Function to format data
  format_data_ourmethod=function(data)
  {
    #load data
    X=data$X
    delta=data$delta
    Z_star=data$Z_star
    n=length(X)
    
    X_tk=array(NA,c(n,b))
    delta_tk=array(NA,c(n,b))
    for(k in 1:b)
    {
      Z_star_k=Z_star-t[k]
      Z_star_k[Z_star_k<0]=NA #recurrent events observed before time t[k]
      
      X_k=X-t[k]
      delta_k=delta
      delta_k[X_k<0]=0 #terminating events observed before time t[k]->patient censored at 0
      X_k[X_k<0]=0
      
      X_Z_star_k=apply(cbind(Z_star_k,X_k),1,min,na.rm=T) #time to next recurrent event/terminating event
      X_tk[,k]=X_Z_star_k
      
      delta_k2=as.numeric(X_Z_star_k<X_k) #=1 if event observed was a recurrent event
      delta_k2[X_Z_star_k==X_k]=delta_k[X_Z_star_k==X_k] #=delta for terminating event
      delta_tk[,k]=delta_k2
    }
    t_k=rep(t,times=rep(n,length(t)))
    X_tk=array(X_tk,c(n*b,1))
    delta_tk=array(delta_tk,c(n*b,1))
    ID=rep(seq(1:n),b)
    list(ID=ID,X_tk=X_tk,delta_tk=delta_tk,t_k=t_k)
  }
  
  #Function to calculate restricted means for test
  get_mu_hat_star_tau=function(X_km,delta_km)
  {
    n=length(X_km)/b
    X=array(X_km,c(n,b))
    delta=array(delta_km,c(n,b))
    
    observed_events=sort(X*delta)
    T=unique(c(observed_events[observed_events<=Tau],Tau))
    M=length(T)-1
    T_array=round(t(array(T[1:M],c(M,b))),6)
    
    dN_i=array(NA,c(n,M))
    Y_i=array(NA,c(n,M))
    for(j in 1:n)
    {
      temp1=array(round(X[j,],6),c(b,M))
      temp2=array(delta[j,],c(b,M))
      dN_i[j,]=apply((temp1==T_array & temp2==1),2,sum) 
      Y_i[j,]=apply(temp1>=T_array,2,sum)
    }
    dN=apply(dN_i,2,sum)
    Y=apply(Y_i,2,sum)
    
    time_int=T[2:(M+1)]-T[1:M]
    temp=array(dN/Y,c(M,M))
    lowerTriangle(temp,diag=F)=0
    CH=apply(temp,2,sum)
    S_hat=exp(-CH)
    mean=sum(time_int*S_hat)
    
    q=dN/Y
    z_i_q=t(t(dN_i-t(q*t(Y_i)))*(n/Y))
    temp3=t(apply(z_i_q,1,sum_function))
    time_int_array=t(array(time_int,c(M,n)))
    S_hat_array=t(array(S_hat,c(M,n)))
    z_i_ATau=apply(time_int_array*S_hat_array*temp3,1,sum)
    williams_var=var(z_i_ATau)/n
    
    list(mean=mean, var=williams_var)
  }
  
  #Run simulations
  #Alias
  p=prop.base
  A=study.time
  A_star=accrual
  l1=recurrent.rate
  l2=terminal.haz
  alpha=ratio
  r1=rho1
  r2=rho2
  t=seq(from=0,to=A,by=sp)
  b=length(t)
  
  mean1=array()
  mean2=array()
  test_stat=array()
  test_stat_p=array()
  chose_correct_treatment_tk=array()
  reject=array()
  
  set.seed(seeds)
  
  i=1
  
  while(i<=n.sim)
  {
    #generate datasets for two groups
    data1=generate_data(n=n1,lambda_1=l1,lambda_2=l2,rho_1=r1,rho_2=r2)
    data2=generate_data(n=n2,lambda_1=l1*alpha[1],lambda_2=l2*alpha[1],rho_1=r1,rho_2=r2)
    
    if(data1$use_dataset==1 & data2$use_dataset==1)
    {
      
      time1 <- Sys.time()
      
      #Get restricted mean for each gruop
      data1_tk=format_data_ourmethod(data1)
      results1=get_mu_hat_star_tau(data1_tk$X_tk,data1_tk$delta_tk)
      data2_tk=format_data_ourmethod(data2)
      results2=get_mu_hat_star_tau(data2_tk$X_tk,data2_tk$delta_tk)
      
      #Get Tayob and Murray statistic
      test_stat[i]=(results1$mean-results2$mean)/sqrt(results1$var+results2$var)
      chose_correct_treatment_tk[i]=as.numeric(results1$mean<results2$mean)
      test_stat_p[i]=2*(1-pnorm(abs(test_stat[i])))
      
      time2 <- Sys.time()
      
      reject[i]=as.numeric(test_stat_p[i]<=0.05)
      mean1[i]=results1$mean
      mean2[i]=results2$mean
      
      print(i)
      print(time2-time1)
      
      i=i+1
    }
  }
  
  #Power of correct rejection
  power=mean(chose_correct_treatment_tk*reject, na.rm = T)
  return(paste("Power is",power))
}


simulate.power(
  n.sim=1000,            #Number of simulations
  prop.base=0.7,         #Proportion of patients enrolled at baseline
  study.time=48,         #Overall study time 
  accrual=24,            #Accrual time
  n1=100,                #Sample size in control group
  n2=100,                #Sample size in treatment group
  recurrent.rate=1/3,    #Recurrent event rate in control group
  terminal.haz=1/1000,   #Terminal event hazard in control group; taking very small number is like no terminal events
  Tau=12,                #Length of follow-up windows
  ratio=0.75,            #ratio = recurrent.rate(terminal.haz) in treatment / recurrent.rate(terminal.haz) in control
  rho1=0,                #Correlation between recurrent events
  rho2=0,                #Correlation between recurrent events and terminal events
  sp=1.5,                #Space between follow-up windows
  seeds=21
)






