library(varbvs)
library(dplyr)

DMCG<- function(tau = 50, a = 2, k = 0.8, b = 50, sigma=0,  lambdad = 1/50, ad = 66,mu_c =  0.4, sides = F, bottom = F, runs = 100000){
  resultdf<-data.frame()
  for(i in 1:runs){
    # # Show both automatic and controlled signal 
    
    # DMC N
    # # Model parameters
    A<- 0;                                                                        # Set A=0 for Neutral trials
    dt<- 0.1; tmax<-1000;                                                          # Set change in time and max time
    
    # Compute time-dependent drift rate mu(t)
    t <- seq(from = dt, to = tmax, by = dt)                                       # Set time
    mu<-A*exp(-t/tau)*(exp(1)*t/(a-1)/tau)^(a-1)*((a-1)/t-1/tau)+ mu_c            # calculate slope at each time point
    # Simulate time-dependent Wiener process X(t) for a single trial
    dX<-mu*dt + sigma*sqrt(dt)*randn(1,length(t));                                # calculate change in x positon                                   
    XN<-cumsum(dX);                                                               # cumulative sum corresponds to X(t)
    
    #Calculate Cross of Boundaries - Neutral
    crossN <-min(which(XN > b))                                                   #Check Correct boundary
    IcrossN <- min(which(XN < -b))                                                #Check Inorrect boundary
    
    # DMC C
    # # Model parameters
    A<- 20;                                                                       # Set A=20 for Congruent trials
    
    # Compute time-dependent drift rate mu(t)
    t <- seq(from = dt, to = tmax, by = dt)                                       # Set time
    mu<-A*exp(-t/tau)*(exp(1)*t/(a-1)/tau)^(a-1)*((a-1)/t-1/tau) +  mu_c          # calculate slope at each time point
    # Simulate time-dependent Wiener process X(t) for a single trial
    dX<-mu*dt + sigma*sqrt(dt)*randn(1,length(t));                                # calculate change in x positon                           
    XC<-cumsum(dX);                                                               # cumulative sum corresponds to X(t)
    
    #Calculate Cross of Boundaries - Congruent
    crossC <-min(which(XC > b))                                                   #Check Correct boundary
    IcrossC <- min(which(XC < -b))                                                #Check Inorrect boundary 
    
    # DMC I
    # # Model parameters
    A<- -20;                                                                      # Set A=-20 for Incongruent trials
    
    # Compute time-dependent drift rate mu(t)
    t <- seq(from = dt, to = tmax, by = dt)                                       # Set time
    mu<-A*exp(-t/tau)*(exp(1)*t/(a-1)/tau)^(a-1)*((a-1)/t-1/tau) + mu_c           # calculate slope at each time point
    # Simulate time-dependent Wiener process X(t) for a single trial
    dX<-mu*dt + sigma*sqrt(dt)*randn(1,length(t));                                # calculate change in x positon  
    XI<-cumsum(dX);                                                               # cumulative sum corresponds to X(t)
    
    crossI <-min(which(XI > b))                                                   #Check Correct boundary
    IcrossI <- min(which(XI < -b))                                                #Check Inorrect boundary 
  
  
    N_C <- crossN-crossC                                                          #Calculate Neutral - Congruent Differnece
    IN_N <- crossI - crossN                                                       #Calculate Incongruent - Neutral Differnece
    result <-N_C/IN_N                                                             #Calculate R
    
    
    resultdf <- rbind(resultdf,(c(crossC, crossN, crossI)))
    
  
    if (is.na(result)) {
      result<- 123456789
    }
    else if(result == Inf){                                                            #Check for any time courses that didn't cross the boundary
      result<- 123456789
    }
    else if(result == -Inf){                                                            #Check for any time courses that didn't cross the boundary
      result<- 123456789
    }
    else if((IcrossI <= crossI) | (IcrossN <= crossN) | (IcrossC <= crossC)){          #Check that a time course didn't hit the incorrect boundary first
      result<- 123456789
    }
  }
  
  colnames(resultdf) <- c('Congruent','Neutral','Incongruent')
  bill <- resultdf %>% filter( (Congruent != Inf) &(Neutral != Inf)&(Incongruent != Inf))%>%
    transmute(Congruent = mean(Congruent),
              Neutral = mean(Neutral),
              Incongruent = mean(Incongruent))
  
  
  resultdf <-(bill$Neutral[1] - bill$Congruent[1])/(bill$Incongruent[1] - bill$Neutral[1])
  return(resultdf)

}

mc.DMCa <- replicate(2,DMCG(sigma = 4))



library(future)
plan(multisession)

for (z in 1:20){
  runs <- 100000                                                                               #Set the number of runs
  sig = 4
  
  #Create containers to hold R values, then replicate
  a <- future({
    mc.DMCa <- replicate(runs,DMCG(sigma = sig))
  }, seed=TRUE)
  b <- future({
    mc.DMCb <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  c <- future({
    mc.DMCc <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  d <- future({
    mc.DMCd <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  e <- future({
    mc.DMCe<- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  f <- future({
    mc.DMCf <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  g <- future({
    mc.DMCg <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  h <- future({
    mc.DMCh <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  i <- future({
    mc.DMCi <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  j <- future({
    mc.DMCj <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  k <- future({
    mc.DMCk <- replicate(runs,DMCG(sigma = sig))
  },seed=TRUE)
  system.time(                                                                               #Bind all R value containers to a single matrix
    c <- rbind(value(a), value(b), value(c), value(d), value(e), value(f), value(g), value(h), value(i), value(j), value(k))
  )
  
  storage <- data.frame(t(c))                                                                #Transpose Matrix
  stor <- storage                                                                            #Make another storage just in case
  #Check for errors and remove them
  bill <- stor %>% filter( (X1 != 123456789) & (X2 != 123456789) & (X3 != 123456789) & (X4 != 123456789) & (X5 != 123456789) &( X6 != 123456789) & (X7 != 123456789)&(X8 != 123456789)&(X9 != 123456789)&(X10 != 123456789)&(X11 != 123456789))
  
  
  
  means1<- apply(bill, 2, mean)                                                              #Apply mean to all columns
  write.csv(means1, paste("~/DMCsimulation100Data",z,".csv", sep = ""))
}
means1
mean(means1)
