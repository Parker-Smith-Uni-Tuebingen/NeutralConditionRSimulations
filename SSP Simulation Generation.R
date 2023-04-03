#Models of Conflict Tasks
library(varbvs)
library(dplyr)
library(ggplot2)

#SSP Model
SSP_Process <- function(sda = 1.861, gamma = 0.018, c = 100, runs = 100000){
  resultdf <- data.frame()
  
  
  for(i in 1:runs){
  dt <- 0.01                                                                                  #Difference in time
  t <- seq(from = 0, to = 1000, by = 0.01)                                                    #Set up up time sequence
  sd_of_t <- sda-gamma*t                                                                      #Calculate shrinking of selective attention
  sd_of_t <- replace(sd_of_t, sd_of_t<0, 0)                                                   #Set spread of attention to 0 to avoid negative sd
  pouter<-1-pnorm(1.5, mean = 0, sd = sd_of_t)                                                #Calculate spread for outer flankers
  pinner<-pnorm(1.5, mean = 0, sd = sd_of_t) - pnorm(.5, mean = 0, sd = sd_of_t)              #Calculate spread for inner flankers
  ptarget <- pnorm(.5, mean = 0, sd = sd_of_t)- pnorm(-.5, mean = 0, sd = sd_of_t)            #Calculate spread for target
  
  bias <- 1                                                                                   #Bias for Congruent Condition
  v_t <- bias*2*pouter+bias*2*pinner+ptarget + c*sqrt(dt)*rnorm(length(t), mean = 0, sd = 1)  #Slope at each time point
  totalC <-cumsum(v_t)                                                                        #overall time-course
  CCrossC <- min(which(totalC > 4000))                                                        #Check Correct boundary
  ICrossC <- min(which(totalC < -4000))                                                       #check Incorrect boundary
  
  
  bias <- -1                                                                                   #Bias for Congruent Condition
  v_tIN <- bias*2*pouter+bias*2*pinner+ptarget +c*sqrt(dt)*rnorm(length(t), mean = 0, sd = 1)  #Slope at each time point
  totalIN <-cumsum(v_tIN)                                                                      #overall time-course
  CCrossIN <- min(which(totalIN > 4000))                                                       #Check Correct boundary
  ICrossIN <- min(which(totalIN < -4000))                                                      #check Incorrect boundary
  
  
  bias <- 0                                                                                   #Bias for Congruent Condition
  v_tN <- bias*2*pouter+bias*2*pinner+ptarget + c*sqrt(dt)*rnorm(length(t), mean = 0, sd = 1) #Slope at each time point
  totalN <-cumsum(v_tN)                                                                       #overall time-course
  CCrossN <-min(which(totalN > 4000))                                                         #Check Correct boundary
  ICrossN <- min(which(totalN < -4000))                                                       #check Incorrect boundary
  
  
  N_C <- CCrossN-CCrossC                                                                      #Calculate Neutral - Congruent Differnece
  IN_N <- CCrossIN - CCrossN                                                                  #Calculate Incongruent - Neutral Differnece
  result <-N_C/IN_N                                                                           #Calculate R
  
  
  resultdf <- rbind(resultdf,(c(CCrossC, CCrossN, CCrossIN)))
}
  colnames(resultdf) <- c('Congruent','Neutral','Incongruent')
  bill <- resultdf %>% filter( (Congruent != Inf) &(Neutral != Inf)&(Incongruent != Inf))%>%
    transmute(Congruent = mean(Congruent),
              Neutral = mean(Neutral),
              Incongruent = mean(Incongruent))
  
  
  resultdf <-(bill$Neutral[1] - bill$Congruent[1])/(bill$Incongruent[1] - bill$Neutral[1])
  return(resultdf)
}





library(future)
plan(multisession)

for (z in 1:2){
runs <- 100000                                                                               #Set the number of runs

#Create containers to hold R values, then replicate
  a <- future({
    mc.SSPa <- replicate(runs,SSP_Process(c = 100))
  }, seed=TRUE)
  b <- future({
    mc.SSPb <- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  c <- future({
    mc.SSPc <- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  d <- future({
    mc.SSPd <- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  e <- future({
    mc.SSPe<- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  f <- future({
    mc.SSPf <- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  g <- future({
    mc.SSPg <- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  h <- future({
    mc.SSPh <- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  i <- future({
    mc.SSPi <- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  j <- future({
    mc.SSPj <- replicate(runs,SSP_Process(c = 100))
  },seed=TRUE)
  system.time(                                                                               #Bind all R value containers to a single matrix
    c <- rbind(value(a), value(b), value(c), value(d), value(e), value(f), value(g), value(h), value(i), value(j), value(k))
  )
                                                                                            #Check for errors and remove them

  
  
  means1<- apply(c, 2, mean)                                                              #Apply mean to all columns
  write.csv(means1, paste("~/SSPsimulation100Data",z,".csv", sep = ""))
}


