functions {
  // Define log probability function model fitting
  real customProb_smc_lpdf(vector RT,
                           vector event_indicator, 
                           vector p_spaq_cum,
                           vector treatment_id,  
                           vector bed_net, 
                           vector age, 
                           real   alpha,  
                           real   lambda, 
                           real   k,
                           real   llin, 
                           real   age_0,
                           real   age_1, 
                           real   age_2, 
                           real   age_3, 
                           real   age_4){
                        
    vector[rows(RT)] lprob;    // per person likelihood  
    real sum_log_prob;         // sum log likelihood 
    sum_log_prob = 0;          // set to 0 
    
    for(a in 1:rows(RT)){
   
   //*******************************************************************************************************************
   //************** control ********************************************************************************************
   //*******************************************************************************************************************
       
    //---------------------------AGE GROUP 0------------------------------------------------------- 
   
      // infected control no bed net
      if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0)        
        lprob[a] = log((alpha * age_0) * exp(-RT[a] * (alpha * age_0)));   
     
      // censored control no bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0)
        lprob[a] = log(exp(-RT[a] * (alpha * age_0))); 
        
      // infected control using bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0)        
        lprob[a] = log( (alpha * llin * age_0) * exp(-RT[a] * (alpha * llin * age_0)));   
     
      // censored control using bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0)
        lprob[a] = log( exp(-RT[a] * (alpha * llin * age_0))); 
        
      
   //---------------------------AGE GROUP 1------------------------------------------------------- 
   
      // infected control no bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 1)        
        lprob[a] = log((alpha * age_1) * exp(-RT[a] * (alpha * age_1)));   
     
      // censored control no bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 1)
        lprob[a] = log(exp(-RT[a] * (alpha * age_1))); 
        
      // infected control using bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 1)        
        lprob[a] = log( (alpha * llin * age_1) * exp(-RT[a] * (alpha * llin * age_1)));   
     
      // censored control using bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 1)
        lprob[a] = log( exp(-RT[a] * (alpha * llin * age_1))); 
        
    //---------------------------AGE GROUP 2-------------------------------------------------------- 
   
      // infected control no bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 2)        
        lprob[a] = log((alpha * age_2) * exp(-RT[a] * (alpha * age_2)));   
     
      // censored control no bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 2)
        lprob[a] = log(exp(-RT[a] * (alpha * age_2))); 
        
      // infected control using bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 2)        
        lprob[a] = log( (alpha * llin * age_2) * exp(-RT[a] * (alpha * llin * age_2)));   
     
      // censored control using bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 2)
        lprob[a] = log( exp(-RT[a] * (alpha * llin * age_2))); 
        
  //---------------------------AGE GROUP 3-------------------------------------------------------- 
   
      // infected control no bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 3)        
        lprob[a] = log((alpha * age_3) * exp(-RT[a] * (alpha * age_3)));   
     
      // censored control no bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 3)
        lprob[a] = log(exp(-RT[a] * (alpha * age_3))); 
        
      // infected control using bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 3)        
        lprob[a] = log( (alpha * llin * age_3) * exp(-RT[a] * (alpha * llin * age_3)));   
     
      // censored control using bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 3)
        lprob[a] = log( exp(-RT[a] * (alpha * llin * age_3))); 
        
  //---------------------------AGE GROUP 4-------------------------------------------------------- 
   
      // infected control no bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 4)        
        lprob[a] = log((alpha * age_4) * exp(-RT[a] * (alpha * age_4)));   
     
      // censored control no bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 4)
        lprob[a] = log(exp(-RT[a] * (alpha * age_4))); 
        
      // infected control using bed net
      else if (treatment_id[a] == 0 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 4)        
        lprob[a] = log( (alpha * llin * age_4) * exp(-RT[a] * (alpha * llin * age_4)));   
     
      // censored control using bed net   
      else if (treatment_id[a] == 0 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 4)
        lprob[a] = log( exp(-RT[a] * (alpha * llin * age_4))); 
        
   //*****************************************************************************************************************    
   //************* SMC ***********************************************************************************************
   //*****************************************************************************************************************
   
       //-----------------------------------AGE GROUP 0--------------------------------------------------------
      // infected smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0)      
        lprob[a] = log( ((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * llin * age_0)) * exp(- p_spaq_cum[a]) ); 
   
      // censored smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0)
        lprob[a] = log( exp(-p_spaq_cum[a]) );  
        
       // infected smc no bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0)      
        lprob[a] = log(((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * age_0)) * exp(- p_spaq_cum[a])); 
   
      // censored smc no bed net site 2 
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0)
        lprob[a] = log(exp(-p_spaq_cum[a]));
        
  //-----------------------------------AGE GROUP 1--------------------------------------------------------
      // infected smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 1)      
        lprob[a] = log( ((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * llin * age_1)) * exp(- p_spaq_cum[a]) ); 
   
      // censored smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 1)
        lprob[a] = log( exp(-p_spaq_cum[a]) );  
        
       // infected smc no bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 1)      
        lprob[a] = log(((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * age_1)) * exp(- p_spaq_cum[a])); 
   
      // censored smc no bed net site 2 
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 1)
        lprob[a] = log(exp(-p_spaq_cum[a]));  
      
  //-----------------------------------AGE GROUP 2--------------------------------------------------------
      // infected smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 2)      
        lprob[a] = log( ((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * llin * age_2)) * exp(- p_spaq_cum[a]) ); 
   
      // censored smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 2)
        lprob[a] = log( exp(-p_spaq_cum[a]) );  
        
       // infected smc no bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 2)      
        lprob[a] = log(((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * age_2)) * exp(- p_spaq_cum[a])); 
   
      // censored smc no bed net site 2 
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 2)
        lprob[a] = log(exp(-p_spaq_cum[a]));  
        
//-----------------------------------AGE GROUP 3--------------------------------------------------------
      // infected smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 3)      
        lprob[a] = log( ((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * llin * age_3)) * exp(- p_spaq_cum[a]) ); 
   
      // censored smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 3)
        lprob[a] = log( exp(-p_spaq_cum[a]) );  
        
       // infected smc no bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 3)      
        lprob[a] = log(((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * age_3)) * exp(- p_spaq_cum[a])); 
   
      // censored smc no bed net site 2 
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 3)
        lprob[a] = log(exp(-p_spaq_cum[a]));  
        
  //-----------------------------------AGE GROUP 4--------------------------------------------------------
      // infected smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 1 && age[a] == 4)      
        lprob[a] = log( ((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * llin * age_4)) * exp(- p_spaq_cum[a]) ); 
   
      // censored smc using bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 1 && age[a] == 4)
        lprob[a] = log( exp(-p_spaq_cum[a]) );  
        
       // infected smc no bed net site 2
      else if (treatment_id[a] == 1 && event_indicator[a] == 1 && bed_net[a] == 0 && age[a] == 4)      
        lprob[a] = log(((1 - ((exp(-(RT[a]/lambda)^k)))) * (alpha * age_4)) * exp(- p_spaq_cum[a])); 
   
      // censored smc no bed net site 2 
      else if (treatment_id[a] == 1 && event_indicator[a] == 0 && bed_net[a] == 0 && age[a] == 4)
        lprob[a] = log(exp(-p_spaq_cum[a]));  
       
   
  
    }
   
    sum_log_prob = sum(lprob);    // sum individuals contributions
    return   sum_log_prob;        // return the log likelihood value
    
  }
 
}

//------------------------------------The input data---------------------------------------------------------------------------------------------- 
data{
    int N;                            // number total individuals   
    vector[N] RT;                     // event times 
    vector[N] event_indicator;        // event indicator 
    vector[N] treatment_id;           // indexes treatment allocation 
    int Ti;                           // number of max event times 
    vector[Ti] time_seq;              // daily time sequence 
    int time_id[N];                   // turn event time into integer 
    vector[N] bed_net;                // bed net variable per person
    vector[N] age;                    // age vector
  
   
}
//-------------------------------------Model parameters--------------------------------------------------------------------------------------------
parameters {
  real<lower = 0> alpha;                // exponential hazard 
  real<lower = 1> lambda;               // weibull scale
  real<lower = 1> k;                    // weibull shape 
  real<lower = 0, upper = 10>  llin;    // bed net multiplier
  real<lower = 0> age_0;                // age modifiers 
  real<lower = 0> age_1; 
  real<lower = 0> age_2; 
  real<lower = 0> age_3; 
  real<lower = 0> age_4; 
}
//-------------------------------------Transformed parameters--------------------------------------------------------------------------------------
transformed parameters{
vector[Ti] p_spaq_t;        // vector of 1-Pspaq eveluated at every time step during follow up (daily)
vector[N]  p_spaq_cum;      // vector for each indivudals cumulative hazard up to observed event time 


// protection of spaq evaluated at each daily interval 
for(i in 1:Ti){
      p_spaq_t[i] =  1 - (exp(-(time_seq[i]/lambda)^k)); 
}

// cumulative sum for each person up to observed event time 
for(j in 1:N){
 
  // site 2 
  if( bed_net[j] == 0 && age[j] == 0){  
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_0);
  }
  else if(bed_net[j] == 1  && age[j] == 0){
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_0 * llin);
  }
  
  else  if( bed_net[j] == 0 && age[j] == 1){  
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_1);
  }
  else if(bed_net[j] == 1  && age[j] == 1){
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_1 * llin);
  }
  else  if( bed_net[j] == 0 && age[j] == 2){
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_2);
  }
  else if(bed_net[j] == 1 && age[j] == 2){
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_2 * llin);
  }

  else  if( bed_net[j] == 0 && age[j] == 3){
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_3);
  }
  else if(bed_net[j] == 1 && age[j] == 3){
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_3 * llin);
  }

  else  if( bed_net[j] == 0 && age[j] == 4){
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_4);
  }
  else if(bed_net[j] == 1 && age[j] == 4){
      p_spaq_cum[j] = sum(p_spaq_t[1:time_id[j]]) * (alpha * age_4 * llin);
  }
  
  
}


 }
 
 
//-----------------------------------Model likelihood-----------------------------------------------------------------------------------------------
model {
  // priors 
    alpha  ~ gamma(0.001,0.001);
    lambda ~ gamma(10,0.3);
    k      ~ gamma(4,0.5);
    llin   ~ uniform(0,10); 
    age_0  ~ lognormal(0,0.8);
    age_1  ~ lognormal(0,0.8);
    age_2  ~ lognormal(0,0.8);
    age_3  ~ lognormal(0,0.8);
    age_4  ~ lognormal(0,0.8);

  // model fit likelihood
  RT ~ customProb_smc(event_indicator,
                      p_spaq_cum,
                      treatment_id,  
                      bed_net, 
                      age,
                      alpha,
                      lambda, 
                      k, 
                      llin, 
                      age_0, 
                      age_1, 
                      age_2, 
                      age_3, 
                      age_4);
}

//-----------------------------Likelihood for loo-cross-validation-----------------------------------------------------------------------------------
generated quantities{
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = customProb_smc_lpdf(rep_vector(RT[n],1)| rep_vector(event_indicator[n],1), rep_vector(p_spaq_cum[n],1), rep_vector(treatment_id[n],1), rep_vector(bed_net[n],1), rep_vector(age[n],1), alpha, lambda, k, llin, age_0, age_1, age_2, age_3, age_4); //  rep_vector(bed_net[n],1), rep_vector(age_1[n],1), rep_vector(age_2[n],1), rep_vector(age_3[n],1), rep_vector(age_4[n],1),  rep_vector(site_2[n],1), rep_vector(site_3[n],1), beta_llin, beta_age, beta_site
  }
}


