
# functions to make the plots for each model 

#-------------------raw data and predicted proportion with clinical malaria----------------------------------------------- 
survival_plot <- function(predicted_data_set){
  
  ggplot() + 
    geom_ribbon(data=predicted_data_set, mapping=aes(x=time, ymin=lower_95, ymax=upper_95, fill=arm),  alpha=0.2) + 
    geom_ribbon(data=predicted_data_set, mapping=aes(x=time, ymin=lower_5, ymax=upper_5, fill=arm), alpha=0.5) + 
    geom_step(t.d.p, mapping=aes(x=time, y=surv, group=treatment), col="black",size=0.7)+ 
    geom_line(data=predicted_data_set, mapping=aes(x=time, y=median, col=arm),size=1) +
    theme_pubr(12) + 
    theme(legend.position = "right") + 
    scale_fill_manual(values = c("SMC" = smc_col, "Control" = control_col)) + 
    scale_color_manual(values = c("SMC" = smc_col,"Control" = control_col))  + 
    labs(x="Time since SP+AQ dose (days)", y = "Proportion diagnosed", col = "", fill="")
  
}


#------------------------duration of protection curves--------------------------------------------------------------------

duration_protection <- function(stan_output_data){
  
    post_draws_model <- rstan::extract(stan_output_data)

    output_mcmc = cbind(post_draws_model[["lambda"]], post_draws_model[["k"]])

    # median parameter values
    par_med = c(median(post_draws_model[["lambda"]]), median(post_draws_model[["k"]]))

    # time sequence
    time_seq = seq(0,60,1)

    # protection over time function
    model_PE <- function(time,  pars){
      lambda         <- pars[1]
      k              <- pars[2]
      ## Drug resposne - protection against infection
      drug_eff <-  exp(-(time / lambda)^k)
      }

   # sampled across posterior draws to calculate CrIs
    model_PE = cmpfun(model_PE, options=list(optimize=3))
    median_pred = sapply( time_seq, model_PE, par=par_med)

    # 20,000 samples 
    N_sam = 20000 
    
    # indexes
    sam_seq = round(seq(from=1, to=nrow(output_mcmc), length=N_sam)) 
    
    # empty matrix to store values in
    sam_cs = matrix(NA, nrow=N_sam, ncol=length(time_seq)) 
    
    for(k in 1:N_sam){
      sam_cs[k,] = sapply(time_seq, model_PE,  par=output_mcmc[sam_seq[k],1:2])
    }
    
    quant_cs = matrix(NA, nrow=3, ncol=length(time_seq))
    
    for(j in 1:length(time_seq)){
      quant_cs[,j] = quantile( sam_cs[,j], prob=c(0.025, 0.5, 0.975) )
    }
    
    quant_cs5 = matrix(NA, nrow=3, ncol=length(time_seq))
    
    for(j in 1:length(time_seq)){
      quant_cs5[,j] = quantile( sam_cs[,j], prob=c(0.25, 0.5, 0.75) )
    }
    
    quant_cs9 = matrix(NA, nrow=3, ncol=length(time_seq))
    
    for(j in 1:length(time_seq)){
      quant_cs9[,j] = quantile( sam_cs[,j], prob=c(0.001, 0.5, 1) )
    }
    
    model.2_efficacy <- data.frame(time = time_seq, 
                                   median = quant_cs[2,], 
                                   lower_5 = quant_cs5[1,], 
                                   upper_5 = quant_cs5[3,], 
                                   lower_95 = quant_cs[1,], 
                                   upper_95 = quant_cs[3,], 
                                   lower_99 = quant_cs9[1,], 
                                   upper_99 = quant_cs9[3,])

    ggplot() + 
      #geom_line(data=trial_eff, mapping=aes(x=time, y=eff), size=1) +
      # geom_line(data=pw_solver, mapping=aes(x=time, y=eff), size=0.8)+
      #geom_ribbon(data=model.2_efficacy, mapping=aes(x=time, ymin=lower_99, ymax=upper_99), fill=eff_profile, alpha=0.1) + 
      geom_ribbon(data=model.2_efficacy, mapping=aes(x=time, ymin=lower_95, ymax=upper_95), fill=eff_profile,  alpha=0.2) + 
      geom_ribbon(data=model.2_efficacy, mapping=aes(x=time, ymin=lower_5, ymax=upper_5), fill=eff_profile, alpha=0.5) + 
      geom_line(data=model.2_efficacy, mapping=aes(x=time, y=median, col=arm), col=eff_profile, size=1) +
      theme_pubr(12) + 
      theme(legend.position = "right") + 
      labs(x="Time since SP+AQ dose (days)", y = "Protective efficacy", col = "", fill="")

}

#-----mean duration of protection--------------------------------------------------------

dur_weib_protect <- function(model_output){

  post_draws_model <- rstan::extract(model_output) 
  
  duration_protection.1 = rep(NA,10000)
  for(i in 1:10000){
    duration_protection.1[i] = post_draws_model$lambda[i] * gamma(1+(1/post_draws_model$k[i]))
  }  
  
  print(quantile(duration_protection.1, probs=c(0.025, 0.5, 0.975)))
  
}

over_90 <- function(model_output){
  post_draws_model <- rstan::extract(model_output) 
  
  temp = rep(NA,10000)
  
  for(i in 1:10000){
    temp[i] = stats::qweibull(p = c(0.10), shape = post_draws_model$k[i], scale = post_draws_model$lambda[i])
  }  
  
  print(quantile(temp, probs=c(0.025, 0.5, 0.975)))
}
