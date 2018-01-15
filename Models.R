#########################################
## Inputs:
## train and data numeric vector of returns.
## trainRV and dataRV numeric vector of intra-day range.
## trainZ and dataZ numeric vector of observed threshold variable, which can be exogenous or self-exciting (i.e. $Z_{t}=Y_{t}$)
## G arbitary constant
## gamma threshold value
## tau VaR level
#########################################


#############################################################################################################################################
####################################### Adaptive CAViaR #####################################################################################
#CAViaR_{t}(\boldsymbol\beta)=CAViaR_{t-1}(\boldsymbol\beta)+\beta_1\{ [ 1+ \text{exp}(G[Y_{t-1}-CAViaR_{t-1}(\beta_1)])]^{-1} - \beta_{2}\}#
#Old: CAViaR[i-1] + beta1*(((1+exp(G*(train[i-1]-CAViaR[i-1])))^(-1))-beta2)
#############################################################################################################################################
Adaptive<-function(betas){
  beta1 <- betas[1]
  for(i in 2:length(train)){ 
    CAViaR[i] <- CAViaR[i-1] + beta1*(((1+exp(G*(train[i-1]-CAViaR[i-1])))^(-1))-tau) 
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  return(res)
}
AdaptiveForecast<-function(betas,data){
  beta1<-betas[1]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- CAViaR[i-1] + beta1*(((1+exp(G*(data[i-1]-CAViaR[i-1])))^(-1))-tau) 
  }
  return(CAViaR)
}


#############################################################################################################################################
####################################### Symmetric absolute value CAViaR  ####################################################################
#################CAViaR_{t}(\boldsymbol\beta)=\beta_1 + \beta_2 CAViaR_{t-1}(\boldsymbol{\beta}) + \beta_3|Y_{t-1}|##########################
#############################################################################################################################################

SymmetricAbs<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  for(i in 2:length(train)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*abs(train[i-1])
  }
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
SymmetricAbsForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*abs(data[i-1])
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
############################################## Symmetric Absolute Value with mu  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta,\mu) = \mu(1-\beta_2) + \beta_1 + \beta_2CAViaR_{t-1}(\boldsymbol\beta,\mu) + \beta_3 |Y_{t-1}-\mu|
#####################################################################################################################################################################################
symmetricAbsmu<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  mu<-betas[4]
  for(i in 2:length(train)){
    CAViaR[i] <- beta1+mu*(1-beta2)+beta2*CAViaR[i-1]+beta3*abs(train[i-1]-mu)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
symmetricAbsmuForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  mu<-betas[4]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- beta1+mu*(1-beta2)+beta2*CAViaR[i-1]+beta3*abs(data[i-1]-mu)
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
##############################################  Assymetric Slope  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta) = \beta_1 + \beta_2 CAViaR_{t-1}(\boldsymbol{\beta})+\beta_3 Y_{t-1}^+ + \beta4 Y_{t-1}^-
#####################################################################################################################################################################################
assymetricSlope<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  for(i in 2:length(train)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*max(train[i-1],0)+beta4*max(-train[i-1],0)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
assymetricSlopeForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*max(data[i-1],0)+beta4*max(-data[i-1],0)
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
##############################################  Assymetric Slope with mu  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta,\mu) = \beta_1+\mu(1-\beta_2) + \beta_2CAViaR_{t-1}(\boldsymbol\beta,\mu) + \beta_3 |Y_{t-1}-\mu|
#####################################################################################################################################################################################
assymetricSlopemu<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  mu<-betas[5]
  for(i in 2:length(train)){
    CAViaR[i] <- beta1+mu*(1-beta2)+beta2*CAViaR[i-1]+beta3*max((train[i-1]-mu),0)+beta4*max(-(train[i-1]-mu),0)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
assymetricSlopemuForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  mu<-betas[5]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- beta1+mu*(1-beta2)+beta2*CAViaR[i-1]+beta3*max((data[i-1]-mu),0)+beta4*max(-(data[i-1]-mu),0)
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
############################################## Indirect GARCH(1, 1)  ################################################################################################################
############################CAViaR_{t}(\boldsymbol{\beta})=(\beta_{1} + \beta_{2}CAViaR_{t-1}^2(\boldsymbol{\beta}) + \beta_{3} Y_{t-1}^2)^{1/2}#####################################
#####################################################################################################################################################################################
IndirectGARCH<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  for(i in 2:length(train)){
    CAViaR[i] <- -(beta1 + beta2 * (CAViaR[i-1]^2) + beta3 * (train[i-1]^2))^(1/2)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
IndirectGARCHForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- -(beta1+beta2*CAViaR[i-1]^2+beta3*data[i-1]^2)^(1/2)
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
############################################## Indirect GARCH with mu  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta,\mu) =  \mu -(\beta_1 + \beta_2 (CAViaR_{t-1}(\boldsymbol{\beta},\mu)-\mu)^2 + \beta_3(Y_{t-1}-\mu)^2)^{1/2}
#####################################################################################################################################################################################
IndirectGARCHmu<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  mu<-betas[4]
  for(i in 2:length(train)){
    CAViaR[i] <- mu-(beta1+beta2* (CAViaR[i-1]-mu)^2+beta3*(train[i-1]-mu)^2)^(1/2)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
IndirectGARCHmuForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  mu<-betas[4]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- mu-(beta1+beta2* (CAViaR[i-1]-mu)^2+beta3*(data[i-1]-mu)^2)^(1/2)
  }
  return(CAViaR)
}

#####################################################################################################################################################################################
############################################## Linear GARCH  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta)  = \beta_1 + \beta_2 CAViaR_{t-1}(\boldsymbol{\beta})+\beta_3|Y_{t-1}|
#####################################################################################################################################################################################
linearGARCH<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  for(i in 2:length(train)){ # 
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*abs(train[i-1])
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}

linearGARCHForecast<-function(betas, data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*abs(data[i-1])
  }
  return(CAViaR)
}

#####################################################################################################################################################################################
############################################## Linear TGARCH  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta)  = \beta_1 + \beta_2 CAViaR_{t-1}(\boldsymbol{\beta})+\beta_3|Y_{t-1}-\mu|
#####################################################################################################################################################################################
linearTGARCH<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  for(i in 2:length(train)){ # 
    CAViaR[i] <-beta1+beta2*CAViaR[i-1]+beta3*abs(train[i-1])
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}

linearTGARCHForecast<-function(betas, data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]

  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*abs(data[i-1])
  }
  return(CAViaR)
}

#####################################################################################################################################################################################
############################################## Linear TGARCH w Mu  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta)  = \mu(1-\beta_2) + \beta_1 + \beta_2 CAViaR_{t-1}(\boldsymbol{\beta})+\beta_3(Y_{t-1}-\mu)^+ + \beta_4(Y_{t-1}-\mu)^-
#####################################################################################################################################################################################
linearTGARCHmu<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  mu<-betas[5]
  for(i in 2:length(train)){ # 
    CAViaR[i] <-mu*(1-beta2)+beta1+beta2*CAViaR[i-1]+beta3*max(train[i-1]-mu)+beta4*min(train[i-1]-mu)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}

linearTGARCHmuForecast<-function(betas, data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  mu<-betas[5]
  
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- mu*(1-beta2)+beta1+beta2*CAViaR[i-1]+beta3*max(data[i-1]-mu)+beta4*min(data[i-1]-mu)
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
############################################## Linear GARCH with mu ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta)  =\mu (1-\beta_2) + \beta_1 + \beta_2 CAViaR_{t-1}(\boldsymbol{\beta})+\beta_3|Y_{t-1} - \mu|
#####################################################################################################################################################################################

linearGARCHmu<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  mu<-betas[4]
  for(i in 2:length(train)){ # 
    CAViaR[i] <-mu*(1-beta2)+beta1+beta2*CAViaR[i-1]+beta3*abs(train[i-1]-mu)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}

linearGARCHmuForecast<-function(betas, data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  mu<-betas[4]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <-mu*(1-beta2)+beta1+beta2*CAViaR[i-1]+beta3*abs(data[i-1]-mu)
  }
  return(CAViaR)
}

#####################################################################################################################################################################################
############################################## GJR-GARCH  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta) =  -(\beta_1 + \beta_2 CAViaR_{t-1}^2(\boldsymbol{\beta}) + \beta_3(Y_{t-1}^+)^2+\beta_4(Y_{t-1}^-)^2)^{1/2}
#####################################################################################################################################################################################
GJRGARCH<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  for(i in 2:length(train)){
    CAViaR[i] <- -(beta1+beta2*CAViaR[i-1]^2+beta3*max(train[i-1],0)^2+beta4*min(train[i-1],0)^2)^(1/2)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
GJRGARCHForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- -(beta1+beta2*CAViaR[i-1]^2+beta3*max(data[i-1],0)^2+beta4*min(data[i-1],0)^2)^(1/2)
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
############################################## GJR-GARCH with mu  ################################################################################################################
############################CAViaR_{t}(\boldsymbol\beta,\mu) =\mu -(\beta_1 + \beta_2 (CAViaR_{t-1}(\boldsymbol\beta,\mu)-\mu)^2 + \beta_3 \{(Y_{t-1}-\mu)^+ \} ^2 + \beta_4\{(Y_{t-1}-\mu)^{-}\}^2)^{1/2}
#####################################################################################################################################################################################
GJRGARCHmu<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  mu<-betas[5]
  for(i in 2:length(train)){
    CAViaR[i] <- mu-(beta1+beta2*(CAViaR[i-1]-mu)^2+beta3*max(train[i-1]-mu,0)^2+beta4*max(-(train[i-1]-mu),0)^2)^(1/2)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
GJRGARCHmuForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  mu<-betas[5]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- mu-(beta1+beta2*(CAViaR[i-1]-mu)^2+beta3*max(data[i-1]-mu,0)^2+beta4*max(-(data[i-1]-mu),0)^2)^(1/2)
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
##############################################  Threshold CAViaR (TCAV)  ################################################################################################################
#####################################################################################################################################################################################
TCAV<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  beta5<-betas[5]
  beta6<-betas[6]
  for(i in 2:length(train)){
    if(trainZ[i-1]<=gamma){
      CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*abs(train[i-1])
    }
    else{
      CAViaR[i] <- beta4+beta5*CAViaR[i-1]+beta6*abs(train[i-1])
    }
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
TCAVForecast<-function(betas,data,dataZ){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  beta5<-betas[5]
  beta6<-betas[6]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    if(dataZ[i-1]<=gamma){
      CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*abs(data[i-1])
    }
    else{
      CAViaR[i] <- beta4+beta5*CAViaR[i-1]+beta6*abs(data[i-1])
    }
  }
  return(CAViaR)
}

#####################################################################################################################################################################################
##############################################   Range Value (RV)  ################################################################################################################
############################CAViaR_{t}(\boldsymbol{\beta}) =\beta_1 + \beta_2 CAViaR_{t-1}(\boldsymbol{\beta}) + \beta_3 R_{t-1}
#####################################################################################################################################################################################
rangeValue<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  for(i in 2:length(trainRV)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*trainRV[i-1]
  }
  #Objective Function
  res<-sum((tau-(trainRV<CAViaR)) * (trainRV-CAViaR)) / length(trainRV) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
rangeValueForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*data[i-1]
  }
  return(CAViaR)
}


#####################################################################################################################################################################################
############################################## Threshold Range Indirect GARCH(1, 1) (TRIG)  ################################################################################################################
#####################################################################################################################################################################################
trigGARCH<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  beta5<-betas[5]
  beta6<-betas[6]
  for(i in 2:length(trainRV)){
    if(trainRV[i-1]<=gamma){
      CAViaR[i]= -(beta1+beta2*(CAViaR[i-1])^2+beta3*trainRV[i-1]^2)^(1/2)
    }
    else{
      CAViaR[i]= -(beta4+beta5*(CAViaR[i-1])^2+beta6*trainRV[i-1]^2)^(1/2)
    }
  }
  #Objective Function
  res<-sum((tau-(trainRV<CAViaR)) * (trainRV-CAViaR)) / length(trainRV) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
trigGARCHForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  beta5<-betas[5]
  beta6<-betas[6]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    if(data[i-1]<=gamma){ 
      CAViaR[i]= -(beta1+beta2*CAViaR[i-1]^2+beta3*data[i-1]^2)^(1/2)
    }
    else{
      CAViaR[i]= -(beta4+beta5*CAViaR[i-1]^2+beta6*data[i-1]^2)^(1/2)
    }
  }
  return(CAViaR)
}

#####################################################################################################################################################################################
############################################## Threshold Range Value (TRV)  ################################################################################################################
#####################################################################################################################################################################################
trvGARCH<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  beta5<-betas[5]
  beta6<-betas[6]
  for(i in 2:length(trainRV)){
    if(trainRV[i-1]<=gamma){
      CAViaR[i]=beta1+beta2*CAViaR[i-1]+beta3*trainRV[i-1]
    }
    else{
      CAViaR[i]=beta4+beta5*CAViaR[i-1]+beta6*trainRV[i-1]
    }
  }
  #Objective Function
  res<-sum((tau-(trainRV<CAViaR)) * (trainRV-CAViaR)) / length(trainRV) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
trvGARCHForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  beta5<-betas[5]
  beta6<-betas[6]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    if(data[i-1]<=gamma){
      CAViaR[i]=beta1+beta2*CAViaR[i-1]+beta3*data[i-1]
    }
    else{
      CAViaR[i]=beta4+beta5*CAViaR[i-1]+beta6*data[i-1]
    }
  }
  return(CAViaR)
}