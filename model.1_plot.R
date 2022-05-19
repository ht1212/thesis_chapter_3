
post_draws <- rstan::extract(model.1)

# store outputs 
test <- matrix(data=NA, nrow=25000, ncol=61) 

# time sequence 
x <- seq(0,60,1)

#11
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_0[i] ) #age 0 bednet 0 site 2 
  test[i,] <- pexp(x, rate4)  
}
#12
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_1[i]) #age 1 bednet 0 site 2
  test[i+2500,] <- pexp(x, rate4)
}
#13
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_2[i]) #age 2 bednet 0 site 2 
  test[i+5000,] <- pexp(x, rate4)
}
#14
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_3[i]) #age 3 bednet 0 site 2 
  test[i+7500,] <- pexp(x, rate4)
}
#15
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_4[i]) #age 4 bednet 0 site 2 
  test[i+10000,] <- pexp(x, rate4)
}
#16
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_0[i] * post_draws$llin[i]) #age 0 bednet 1 site 2 
  test[i+12500,] <- pexp(x, rate4)
}
#17
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_1[i] * post_draws$llin[i] ) #age 1 bednet 1 site 2 
  test[i+15000,] <- pexp(x, rate4)
}
#18
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_2[i] * post_draws$llin[i] ) #age 2 bednet 1 site 2 
  test[i+17500,] <- pexp(x, rate4)
}
#19
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_3[i] * post_draws$llin[i] ) #age 3 bednet 1 site 2 
  test[i+20000,] <- pexp(x, rate4)
}
#20
for(i in 1:2500){
  rate4 <- c(post_draws$alpha[i] * post_draws$age_4[i] * post_draws$llin[i] ) #age 4 bednet 1 site 2 
  test[i+22500,] <- pexp(x, rate4)
}

control_5 = matrix(NA, nrow=3, ncol=length(x))
control_95 = matrix(NA, nrow=3, ncol=length(x))

for(j in 1:length(x)){
  control_5[,j] = quantile( test[,j], prob=c(0.25, 0.5, 0.75),na.rm = T )
}

for(j in 1:length(x)){
  control_95[,j] = quantile( test[,j], prob=c(0.025, 0.5, 0.975),na.rm = T )
}

control_predicted <- data.frame(time = x, 
                                median = control_5[2,], 
                                lower_5 = control_5[1,], 
                                upper_5 = control_5[3,], 
                                lower_95 = control_95[1,], 
                                upper_95 = control_95[3,], 
                                arm = "Control")

#----------------------------SMC ARM-----------------------------------------------------------------------
p3 <- matrix(data=NA, nrow=25000, ncol=61)

#1 
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i] * post_draws$age_0[i]) #age 0 bednet 0 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i,] <- 1-exp(-p2)
  
}
#2
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i] * post_draws$age_1[i]) #age 1 bednet 0 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+2500,] <- 1-exp(-p2)
  
}
#3 
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i] * post_draws$age_2[i]) #age 2 bednet 0 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+5000,] <- 1-exp(-p2)
  
}
#4 
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i] * post_draws$age_3[i]) #age 3 bednet 0 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+7500,] <- 1-exp(-p2)
  
}
#5
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i] * post_draws$age_4[i]) #age 4 bednet 0 site 2
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+10000,] <- 1-exp(-p2)
  
}
#6
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i] * post_draws$age_0[i] * post_draws$llin[i]) #age 0 bednet 1 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+12500,] <- 1-exp(-p2)
  
}
#7
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i]  * post_draws$age_1[i] * post_draws$llin[i]) #age 1 bednet 1 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  control_haz[1:length(x)] = rate4[1]
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+15000,] <- 1-exp(-p2)
  
}
#8
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i]  * post_draws$age_2[i] * post_draws$llin[i]) #age 2 bednet 1 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+17500,] <- 1-exp(-p2)
  
}
#9
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i]  * post_draws$age_3[i] * post_draws$llin[i]) #age 3 bednet 1 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+20000,] <- 1-exp(-p2)
  
}
#10
for(i in 1:2500){
  
  lambda <- c(post_draws$lambda[i])
  k <- c(post_draws$k[i])
  rate4 <- c(post_draws$alpha[i]  * post_draws$age_4[i] * post_draws$llin[i]) #age 4 bednet 1 site 2 
  weib <- 1 - (exp(-(x/lambda)^k)) 
  control_haz = rep(NA, length(x))
  
  control_haz[1:length(x)] = rate4[1]
  
  p1 <- (control_haz * weib)
  p2 <- cumsum(p1)
  p3[i+22500,] <- 1-exp(-p2)
  
}

smc_5 = matrix(NA, nrow=3, ncol=length(x))

for(j in 1:length(x)){
  smc_5[,j] = quantile( p3[,j], prob=c(0.25, 0.5, 0.75),na.rm = T )
}

smc_95 = matrix(NA, nrow=3, ncol=length(x))

for(j in 1:length(x)){
  smc_95[,j] = quantile( p3[,j], prob=c(0.025, 0.5, 0.975),na.rm = T )
}

smc_predicted <- data.frame(time = x, 
                            median = smc_5[2,], 
                            lower_5 = smc_5[1,], 
                            upper_5 = smc_5[3,], 
                            lower_95 = smc_95[1,], 
                            upper_95 = smc_95[3,], 
                            arm = "SMC")

model.1_predited <- bind_rows(control_predicted, smc_predicted)
