
post_draws_1 <- rstan::extract(model.2)

test_2 <- matrix(data=NA, nrow=25000, ncol=61)

x <- seq(0,60,1)
t <- c(0,30)

breaks <- seq(30,60,30)

for(i in 1:2500){
  rate4 <- c( post_draws_1$alpha[i,] * post_draws_1$age_0[i]) #age 0 bednet 0 
  test_2[i,] <- ppexp(x, rate4, t)
}  
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_1[i]) #age 1 bednet 0 
  test_2[i+2500,] <- ppexp(x, rate4, t)  
}
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_2[i]) #age 2 bednet 0 
  test_2[i+5000,] <- ppexp(x, rate4, t)
}
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_3[i]) #age 3 bednet 0 
  test_2[i+7500,] <- ppexp(x, rate4, t)
}
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_4[i]) #age 4 bednet 0 
  test_2[i+10000,] <- ppexp(x, rate4, t)
}
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_0[i] * post_draws_1$llin[i]) #age 0 bednet 1 
  test_2[i+12500,] <- ppexp(x, rate4, t)
}
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_1[i] * post_draws_1$llin[i]) #age 1 bednet 1 
  test_2[i+15000,] <- ppexp(x, rate4, t)
}
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_2[i] * post_draws_1$llin[i]) #age 2 bednet 1 
  test_2[i+17500,] <- ppexp(x, rate4, t)
}
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_3[i] * post_draws_1$llin[i]) #age 3 bednet 1 
  test_2[i+20000,] <- ppexp(x, rate4, t)
}
for(i in 1:2500){
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_4[i] * post_draws_1$llin[i]) #age 4 bednet 1 
  test_2[i+22500,] <- ppexp(x, rate4, t)
}

control_5 = matrix(NA, nrow=3, ncol=length(x))
control_95 = matrix(NA, nrow=3, ncol=length(x))

for(j in 1:length(x)){
  control_5[,j] = quantile( test_2[,j], prob=c(0.25, 0.5, 0.75),na.rm = T )
}

for(j in 1:length(x)){
  control_95[,j] = quantile( test_2[,j], prob=c(0.025, 0.5, 0.975),na.rm = T )
}

control_predicted_2 <- data.frame(time = x, 
                                  median = control_5[2,], 
                                  lower_5 = control_5[1,], 
                                  upper_5 = control_5[3,], 
                                  lower_95 = control_95[1,], 
                                  upper_95 = control_95[3,], 
                                  arm = "Control")

#----------------------------SMC ARM-----------------------------------------------------------------------
p4 <- matrix(data=NA, nrow=25000, ncol=61)

#1
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_0[i]) #age 0 bednet 0 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i,] <- 1-exp(-p2)
  
}
#2
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_1[i]) #age 1 bednet 0  
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+2500,] <- 1-exp(-p2)
  
}
#3
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_2[i]) #age 2 bednet 0 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+5000,] <- 1-exp(-p2)
  
}
#4
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_3[i]) #age 3 bednet 0  
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+7500,] <- 1-exp(-p2)
  
}
#5
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_4[i]) #age 4 bednet 0 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+10000,] <- 1-exp(-p2)
  
}
#6
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_0[i] * post_draws_1$llin[1]) # age 0 bednet 1 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+12500,] <- 1-exp(-p2)
  
}
#7
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_1[i] * post_draws_1$llin[1]) #age 1 bednet 1 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+15000,] <- 1-exp(-p2)
  
}
#8
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_2[i] * post_draws_1$llin[1]) #age 2 bednet 1  
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+17500,] <- 1-exp(-p2)
  
}
#9
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_3[i] * post_draws_1$llin[1]) #age 3 bednet 1  
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+20000,] <- 1-exp(-p2)
  
}
#10
for(i in 1:2500){
  
  lambda <- c(post_draws_1$lambda[i])
  k <- c(post_draws_1$k[i])
  rate4 <- c(post_draws_1$alpha[i,] * post_draws_1$age_4[i] * post_draws_1$llin[1]) #age 4 bednet 1  
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:which(x == breaks[1])] = rate4[1]
  control_haz[which(x == breaks[1]):which(x == breaks[2])] = rate4[2]

  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p4[i+22500,] <- 1-exp(-p2)
  
}



smc_5 = matrix(NA, nrow=3, ncol=length(x))

for(j in 1:length(x)){
  smc_5[,j] = quantile( p4[,j], prob=c(0.25, 0.5, 0.75),na.rm = T )
}

smc_95 = matrix(NA, nrow=3, ncol=length(x))

for(j in 1:length(x)){
  smc_95[,j] = quantile( p4[,j], prob=c(0.025, 0.5, 0.975),na.rm = T )
}

smc_predicted_2 <- data.frame(time = x, 
                              median = smc_5[2,], 
                              lower_5 = smc_5[1,], 
                              upper_5 = smc_5[3,], 
                              lower_95 = smc_95[1,], 
                              upper_95 = smc_95[3,], 
                              arm = "SMC")

model.2_predicted <- bind_rows(control_predicted_2, smc_predicted_2)
