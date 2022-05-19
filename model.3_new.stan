functions {
  // Define log probability function for the SMC arm malaria eveny
  real customProb_smc_lpdf(vector RT,
                           vector indicator,
                           vector interval, 
                           vector treatment_id, 
                           vector bed_net, 
                           vector age, 
                           vector p_spaq_cum, 
                           vector alpha,
                           real   llin,
                           real   age_0, 
                           real   age_1, 
                           real   age_2, 
                           real   age_3,
                           real   age_4, 
                           real   lambda, 
                           real   k ){
    
    vector[num_elements(RT)] lprob;    // per person likelihood  
    real sum_log_prob;                // sum log likelihood 
    sum_log_prob = 0;                 // set to 0 
    
   
   for(a in 1:num_elements(RT)){
     
//************************************************************************************************************
//**********************************CONTROL NO BED BETS*******************************************************
//************************************************************************************************************

     // interval 1 infected control 
     if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( ((alpha[1] * age_0) * exp(-(alpha[1] * age_0 * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( (alpha[2] * age_0) * exp(-(alpha[2] * age_0 *  RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( (alpha[3] * age_0) * exp(-(alpha[3] * age_0 * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( (alpha[4] * age_0) * exp(-(alpha[4] * age_0 * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( (alpha[5] * age_0) * exp(-(alpha[5] * age_0 *  RT[a]))); 
     }
     // interval 6 infected control
      else if(interval[a] == 6 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( (alpha[6] * age_0) * exp(-(alpha[6] * age_0 * RT[a]))); 
     } 
     // interval 7 infected control
      else if(interval[a] == 7 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( (alpha[7] * age_0) * exp(-(alpha[7] * age_0 * RT[a]))); 
     }
     // interval 8 infected control
      else if(interval[a] == 8 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( (alpha[8] * age_0) * exp(-(alpha[8] * age_0 * RT[a]))); 
     }
     // interval 9 infected control
      else if(interval[a] == 9 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( (alpha[9] * age_0) * exp(-(alpha[9] * age_0 * RT[a]))); 
     }
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[1] * age_0 * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[2] * age_0 * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[3] * age_0 * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[4] * age_0 * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[5] * age_0 * RT[a]))); 
     }
     // interval 6 censored control
      else if(interval[a] == 6 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[6] * age_0 * RT[a]))); 
     }
     // interval 7 censored control
      else if(interval[a] == 7 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[7] * age_0 * RT[a]))); 
     }
     // interval 8 censored control
     else if(interval[a] == 8 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[8] * age_0 * RT[a]))); 
     }
     // interval 9 censored control
     else if(interval[a] == 9 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[9] * age_0 * RT[a]))); 
     }
 
 //****************************************************************************************************************************   
 //****************************CONTROL YES BED BETS****************************************************************************
 //****************************************************************************************************************************
 
     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( ((alpha[1] * age_0 * llin) * exp(-(alpha[1] * age_0 * llin * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( (alpha[2] * age_0 * llin) * exp(-(alpha[2] * age_0 * llin * RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( (alpha[3] * age_0 * llin) * exp(-(alpha[3] * age_0 * llin * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( (alpha[4] * age_0 * llin) * exp(-(alpha[4] * age_0 * llin * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( (alpha[5] * age_0 * llin) * exp(-(alpha[5] * age_0 * llin * RT[a]))); 
     }
     // interval 6 infected control
      else if(interval[a] == 6 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( (alpha[6] * age_0 * llin) * exp(-(alpha[6] * age_0 * llin * RT[a]))); 
     }
     // interval 7 infected control
      else if(interval[a] == 7 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( (alpha[7] * age_0 * llin) * exp(-(alpha[7] * age_0 * llin * RT[a]))); 
     }
     // interval 8 infected control
      else if(interval[a] == 8 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( (alpha[8] * age_0 * llin) * exp(-(alpha[8] * age_0 * llin * RT[a]))); 
     }
     // interval 9 infected control
      else if(interval[a] == 9 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( (alpha[9] * age_0 * llin) * exp(-(alpha[9] * age_0 * llin * RT[a]))); 
     }
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[1] * age_0 * llin * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[2] * age_0 * llin * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[3] * age_0 * llin * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[4] * age_0 * llin * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[5] * age_0 * llin * RT[a]))); 
     }
     // interval 6 censored control
      else if(interval[a] == 6 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[6] * age_0 * llin * RT[a]))); 
     }
     // interval 7 censored control
      else if(interval[a] == 7 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[7] * age_0 * llin * RT[a]))); 
     }
     // interval 8 censored control
     else if(interval[a] == 8 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[8] * age_0 * llin * RT[a]))); 
     }
     // interval 9 censored control
     else if(interval[a] == 9 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 0){
              lprob[a] = log( exp(-(alpha[9] * age_0 * llin * RT[a]))); 
     }    
     
  //********************************************************************************************************************************   
  //***************************************SMC NO BED NETS SITE 2 ******************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 0){
             lprob[a] = log( (alpha[1] * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 0){
              lprob[a] = log( (alpha[2] * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( (alpha[3]  * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( (alpha[4] * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( (alpha[5] * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 6 infected smc
     else if(interval[a] == 6 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( (alpha[6] * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 7 infected smc 
     else if(interval[a] == 7 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( (alpha[7] * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 8 infected smc 
     else if(interval[a] == 8 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( (alpha[8] * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 9 infected smc
     else if(interval[a] == 9 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( (alpha[9] * age_0 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );
     }
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 6 censored smc
     else if(interval[a] == 6 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 7 censored smc
     else if(interval[a] == 7 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 8 censored smc
     else if(interval[a] == 8 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 9 censored smc
     else if(interval[a] == 9 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     
  //********************************************************************************************************************************   
  //***************************************SMC YES BED NETS site 2 *****************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log( (alpha[1] * age_0 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 0){
              lprob[a] = log( (alpha[2] * age_0 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log( (alpha[3] * age_0 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log( (alpha[4] * age_0 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log( (alpha[5] * age_0 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 0){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }

     
//***********************************************************************************************************************************************************
//***********************************************************************************************************************************************************
//***************************************************AGE GROUP 1*********************************************************************************************
//***********************************************************************************************************************************************************
//*********************************************************************************************************************************************************** 
    
    
    //************************************************************************************************************
   //**********************************CONTROL NO BED BETS*******************************************************
   //************************************************************************************************************

     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( ((alpha[1] * age_1) * exp(-(alpha[1] * age_1 * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( (alpha[2] * age_1) * exp(-(alpha[2] * age_1 *  RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( (alpha[3] * age_1) * exp(-(alpha[3] * age_1 * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( (alpha[4] * age_1) * exp(-(alpha[4] * age_1 * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( (alpha[5] * age_1) * exp(-(alpha[5] * age_1 *  RT[a]))); 
     }
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[1] * age_1 * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[2] * age_1 * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[3] * age_1 * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[4] * age_1 * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[5] * age_1 * RT[a]))); 
     }
 
 //****************************************************************************************************************************   
 //****************************CONTROL YES BED BETS****************************************************************************
 //****************************************************************************************************************************
 
     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( ((alpha[1] * age_1 * llin) * exp(-(alpha[1] * age_1 * llin * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( (alpha[2] * age_1 * llin) * exp(-(alpha[2] * age_1 * llin * RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( (alpha[3] * age_1 * llin) * exp(-(alpha[3] * age_1 * llin * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( (alpha[4] * age_1 * llin) * exp(-(alpha[4] * age_1 * llin * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( (alpha[5] * age_1 * llin) * exp(-(alpha[5] * age_1 * llin * RT[a]))); 
     }
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[1] * age_1 * llin * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[2] * age_1 * llin * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[3] * age_1 * llin * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[4] * age_1 * llin * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 1){
              lprob[a] = log( exp(-(alpha[5] * age_1 * llin * RT[a]))); 
     }
     
  //********************************************************************************************************************************   
  //***************************************SMC NO BED NETS SITE 2 ******************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log( (alpha[1] * age_1 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 1){
              lprob[a] = log( (alpha[2] * age_1 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log( (alpha[3]  * age_1 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log( (alpha[4] * age_1 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log( (alpha[5] * age_1 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 1){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     
  //********************************************************************************************************************************   
  //***************************************SMC YES BED NETS site 2 *****************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log( (alpha[1] * age_1 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 1){
              lprob[a] = log( (alpha[2] * age_1 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log( (alpha[3] * age_1 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log( (alpha[4] * age_1 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log( (alpha[5] * age_1 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 1){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
    
     
//***********************************************************************************************************************************************************
//***********************************************************************************************************************************************************
//***************************************************AGE GROUP 2********************************************************************************************
//***********************************************************************************************************************************************************
//*********************************************************************************************************************************************************** 
    
    
    //************************************************************************************************************
   //**********************************CONTROL NO BED BETS*******************************************************
   //************************************************************************************************************

     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( ((alpha[1] * age_2) * exp(-(alpha[1] * age_2 * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( (alpha[2] * age_2) * exp(-(alpha[2] * age_2 *  RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( (alpha[3] * age_2) * exp(-(alpha[3] * age_2 * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( (alpha[4] * age_2) * exp(-(alpha[4] * age_2 * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( (alpha[5] * age_2) * exp(-(alpha[5] * age_2 *  RT[a]))); 
     }
     
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[1] * age_2 * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[2] * age_2 * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[3] * age_2 * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[4] * age_2 * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[5] * age_2 * RT[a]))); 
     }
    
 
 //****************************************************************************************************************************   
 //****************************CONTROL YES BED BETS****************************************************************************
 //****************************************************************************************************************************
 
     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( ((alpha[1] * age_2 * llin) * exp(-(alpha[1] * age_2 * llin * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( (alpha[2] * age_2 * llin) * exp(-(alpha[2] * age_2 * llin * RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( (alpha[3] * age_2 * llin) * exp(-(alpha[3] * age_2 * llin * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( (alpha[4] * age_2 * llin) * exp(-(alpha[4] * age_2 * llin * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( (alpha[5] * age_2 * llin) * exp(-(alpha[5] * age_2 * llin * RT[a]))); 
     }
    
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[1] * age_2 * llin * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[2] * age_2 * llin * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[3] * age_2 * llin * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[4] * age_2 * llin * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 2){
              lprob[a] = log( exp(-(alpha[5] * age_2 * llin * RT[a]))); 
     }
  
     
  //********************************************************************************************************************************   
  //***************************************SMC NO BED NETS SITE 2 ******************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log( (alpha[1] * age_2 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 2){
              lprob[a] = log( (alpha[2] * age_2 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log( (alpha[3]  * age_2 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log( (alpha[4] * age_2 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log( (alpha[5] * age_2 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 2){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     
  //********************************************************************************************************************************   
  //***************************************SMC YES BED NETS site 2 *****************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log( (alpha[1] * age_2 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 2){
              lprob[a] = log( (alpha[2] * age_2 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log( (alpha[3] * age_2 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log( (alpha[4] * age_2 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log( (alpha[5] * age_2 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 2){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
 
     
//***********************************************************************************************************************************************************
//***********************************************************************************************************************************************************
//***************************************************AGE GROUP 3*********************************************************************************************
//***********************************************************************************************************************************************************
//*********************************************************************************************************************************************************** 
    
    
    //************************************************************************************************************
   //**********************************CONTROL NO BED BETS*******************************************************
   //************************************************************************************************************

     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( ((alpha[1] * age_3) * exp(-(alpha[1] * age_3 * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( (alpha[2] * age_3) * exp(-(alpha[2] * age_3 *  RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( (alpha[3] * age_3) * exp(-(alpha[3] * age_3 * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( (alpha[4] * age_3) * exp(-(alpha[4] * age_3 * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( (alpha[5] * age_3) * exp(-(alpha[5] * age_3 *  RT[a]))); 
     }
     
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[1] * age_3 * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[2] * age_3 * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[3] * age_3 * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[4] * age_3 * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[5] * age_3 * RT[a]))); 
     }
    
 
 //****************************************************************************************************************************   
 //****************************CONTROL YES BED BETS****************************************************************************
 //****************************************************************************************************************************
 
     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( ((alpha[1] * age_3 * llin) * exp(-(alpha[1] * age_3 * llin * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( (alpha[2] * age_3 * llin) * exp(-(alpha[2] * age_3 * llin * RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( (alpha[3] * age_3 * llin) * exp(-(alpha[3] * age_3 * llin * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( (alpha[4] * age_3 * llin) * exp(-(alpha[4] * age_3 * llin * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( (alpha[5] * age_3 * llin) * exp(-(alpha[5] * age_3 * llin * RT[a]))); 
     }
     
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[1] * age_3 * llin * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[2] * age_3 * llin * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[3] * age_3 * llin * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[4] * age_3 * llin * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 3){
              lprob[a] = log( exp(-(alpha[5] * age_3 * llin * RT[a]))); 
     }
    
  //********************************************************************************************************************************   
  //***************************************SMC NO BED NETS SITE 2 ******************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log( (alpha[1] * age_3 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 3){
              lprob[a] = log( (alpha[2] * age_3 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log( (alpha[3]  * age_3 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log( (alpha[4] * age_3 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log( (alpha[5] * age_3 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 3){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     
  //********************************************************************************************************************************   
  //***************************************SMC YES BED NETS site 2 *****************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log( (alpha[1] * age_3 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 3){
              lprob[a] = log( (alpha[2] * age_3 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log( (alpha[3] * age_3 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log( (alpha[4] * age_3 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log( (alpha[5] * age_3 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 3){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     
     
//***********************************************************************************************************************************************************
//***********************************************************************************************************************************************************
//***************************************************AGE GROUP 4*********************************************************************************************
//***********************************************************************************************************************************************************
//*********************************************************************************************************************************************************** 
    
    
    //************************************************************************************************************
   //**********************************CONTROL NO BED BETS*******************************************************
   //************************************************************************************************************

     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( ((alpha[1] * age_4) * exp(-(alpha[1] * age_4 * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( (alpha[2] * age_4) * exp(-(alpha[2] * age_4 *  RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( (alpha[3] * age_4) * exp(-(alpha[3] * age_4 * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( (alpha[4] * age_4) * exp(-(alpha[4] * age_4 * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( (alpha[5] * age_4) * exp(-(alpha[5] * age_4 *  RT[a]))); 
     }
     
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[1] * age_4 * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[2] * age_4 * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[3] * age_4 * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[4] * age_4 * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 0 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[5] * age_4 * RT[a]))); 
     }
  
 
 //****************************************************************************************************************************   
 //****************************CONTROL YES BED BETS****************************************************************************
 //****************************************************************************************************************************
 
     // interval 1 infected control 
     else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( ((alpha[1] * age_4 * llin) * exp(-(alpha[1] * age_4 * llin * RT[a]))) ); 
     }
     // interval 2 infected control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( (alpha[2] * age_4 * llin) * exp(-(alpha[2] * age_4 * llin * RT[a])));  
     }
     // interval 3 infected control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( (alpha[3] * age_4 * llin) * exp(-(alpha[3] * age_4 * llin * RT[a]))); 
     }
     // interval 4 infected control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( (alpha[4] * age_4 * llin) * exp(-(alpha[4] * age_4 * llin * RT[a]))); 
     }
     // interval 5 infected control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 1 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( (alpha[5] * age_4 * llin) * exp(-(alpha[5] * age_4 * llin * RT[a]))); 
     }
     
     // interval 1 censored control 
      else if(interval[a] == 1 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[1] * age_4 * llin * RT[a]))); 
     }
     // interval 2 censored control
      else if(interval[a] == 2 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[2] * age_4 * llin * RT[a])));  
     }
     // interval 3 censored control
      else if(interval[a] == 3 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[3] * age_4 * llin * RT[a]))); 
     }
     // interval 4 censored control
      else if(interval[a] == 4 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[4] * age_4 * llin * RT[a]))); 
     }
     // interval 5 censored control
      else if(interval[a] == 5 && treatment_id[a] == 0 && indicator[a] == 0 && bed_net[a] == 1 && age[a] == 4){
              lprob[a] = log( exp(-(alpha[5] * age_4 * llin * RT[a]))); 
     }
     
     
  //********************************************************************************************************************************   
  //***************************************SMC NO BED NETS SITE 2 ******************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log( (alpha[1] * age_4 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 4){
              lprob[a] = log( (alpha[2] * age_4 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log( (alpha[3]  * age_4 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log( (alpha[4] * age_4 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log( (alpha[5] * age_4 * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 0  && age[a] == 4){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
    
     
  //********************************************************************************************************************************   
  //***************************************SMC YES BED NETS site 2 *****************************************************************
  //********************************************************************************************************************************
  
     // interval 1 infected smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log( (alpha[1] * age_4 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(-( ( p_spaq_cum[a]))) ); 
     }
     // interval 2 infected smc 
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 4){
              lprob[a] = log( (alpha[2] * age_4 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) * exp(- ( ( p_spaq_cum[a])) )); 
     }
     // interval 3 infected smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log( (alpha[3] * age_4 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) ); 
     }
     // interval 4 infected smc
     else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log( (alpha[4] * age_4 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     // interval 5 infected smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 1 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log( (alpha[5] * age_4 * llin * ((1 - ((exp(-(RT[a]/lambda)^k)))))) *  exp(- ( ( p_spaq_cum[a]))) );  
     }
     
     // interval 1 censored smc 
      else if(interval[a] == 1 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log(exp(-( ( p_spaq_cum[a])))); 
     }
     // interval 2 censored smc
      else if(interval[a] == 2 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a])))  );   
     }
     // interval 3 censored smc 
      else if(interval[a] == 3 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) ); 
     }
     // interval 4 censored smc
      else if(interval[a] == 4 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }
     // interval 5 censored smc
     else if(interval[a] == 5 && treatment_id[a] == 1 && indicator[a] == 0 && bed_net[a] == 1  && age[a] == 4){
             lprob[a] = log( exp(- ( ( p_spaq_cum[a]) )) );   
     }

    
     
}
    sum_log_prob = sum(lprob);    // sum individuals contributions
    return   sum_log_prob;        // return the log likelihood value
    
  }
 
}

// The input data is a vector 'y' of length 'N'.
data{
    int N;                            // number total   
    int A; 
    vector[N] RT;                     // data times 
    vector[N] indicator;              // censoring status 
    vector[N] treatment_id;           // index arm of individual 
    vector[N] interval;               // intervals 
    vector[N] bed_net;                // bed net variable per person
    vector[N] age;                    // age vector
    int Ti;                           // number of max event times 
    vector[Ti] time_seq;              // time sequence 
    int time_id[N];                   // turn event time into integer censored
    int break_ids[A-1];
     
}
//
parameters {
  vector<lower = 0>[A] alpha; 
  real<lower = 0>      lambda; 
  real<lower = 0>      k; 
  real<lower = 0>      llin;
  real<lower = 0>      age_0; 
  real<lower = 0>      age_1; 
  real<lower = 0>      age_2; 
  real<lower = 0>      age_3; 
  real<lower = 0>      age_4; 
  

}
//
transformed parameters{
vector[Ti] p_spaq_t;        // vector of 1-Pspaq eveluated at every time step in the season
vector[N]  p_spaq_cum;      // vector for each indivudals cumulative hazard till time to event 

// protection of spaq evaluated at each daily interval
for(i in 1:Ti){
      p_spaq_t[i] = 1 - (exp(-(time_seq[i]/lambda)^k)); //- (1-exp(-(time_seq[i]-1/lambda)^k));
}


// cumulative sum for each person up to censoring time point
for(j in 1:N){
  
  //-------AGE GROUP 0--------------------------------------------------------------
  
  // No bed-net site 2
  if(interval[j] == 1 && bed_net[j] == 0  && age[j] == 0){
      p_spaq_cum[j] = alpha[1] * age_0 * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 0  && age[j] == 0){
      p_spaq_cum[j] = alpha[2] * age_0 * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 0  && age[j] == 0){
     p_spaq_cum[j] = alpha[3] * age_0 * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 0  && age[j] == 0){
     p_spaq_cum[j] = alpha[4] * age_0 * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 0  && age[j] == 0){
     p_spaq_cum[j] = alpha[5] * age_0 * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }
 
 // yes bet net site 2 
  else if(interval[j] == 1 && bed_net[j] == 1  && age[j] == 0){
      p_spaq_cum[j] = alpha[1] * age_0 * llin * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 1  && age[j] == 0){
      p_spaq_cum[j] = alpha[2] * age_0 * llin * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 1  && age[j] == 0){
     p_spaq_cum[j] = alpha[3] * age_0 * llin * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 1  && age[j] == 0){
     p_spaq_cum[j] = alpha[4] * age_0 * llin * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 1  && age[j] == 0){
     p_spaq_cum[j] = alpha[5] * age_0 * llin * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }
 

//------------AGE GROUP 1-------------------------------------------------------------

  // No bed-net site 2
  else if(interval[j] == 1 && bed_net[j] == 0  && age[j] == 1){
      p_spaq_cum[j] = alpha[1] * age_1 * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 0  && age[j] == 1){
      p_spaq_cum[j] = alpha[2] * age_1 * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 0  && age[j] == 1){
     p_spaq_cum[j] = alpha[3] * age_1 * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 0  && age[j] == 1){
     p_spaq_cum[j] = alpha[4] * age_1 * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 0  && age[j] == 1){
     p_spaq_cum[j] = alpha[5] * age_1 * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }
 
 // yes bet net site 2 
  else if(interval[j] == 1 && bed_net[j] == 1  && age[j] == 1){
      p_spaq_cum[j] = alpha[1] * age_1 * llin * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 1  && age[j] == 1){
      p_spaq_cum[j] = alpha[2] * age_1 * llin * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 1  && age[j] == 1){
     p_spaq_cum[j] = alpha[3] * age_1 * llin * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 1  && age[j] == 1){
     p_spaq_cum[j] = alpha[4] * age_1 * llin * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 1  && age[j] == 1){
     p_spaq_cum[j] = alpha[5] * age_1 * llin * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }

//--------------AGE GROUP 2------------------------------------------------------------

  // No bed-net site 2
  else if(interval[j] == 1 && bed_net[j] == 0  && age[j] == 2){
      p_spaq_cum[j] = alpha[1] * age_2 * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 0  && age[j] == 2){
      p_spaq_cum[j] = alpha[2] * age_2 * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 0  && age[j] == 2){
     p_spaq_cum[j] = alpha[3] * age_2 * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 0  && age[j] == 2){
     p_spaq_cum[j] = alpha[4] * age_2 * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 0  && age[j] == 2){
     p_spaq_cum[j] = alpha[5] * age_2 * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }

 // yes bet net site 2 
  else if(interval[j] == 1 && bed_net[j] == 1  && age[j] == 2){
      p_spaq_cum[j] = alpha[1] * age_2 * llin * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 1  && age[j] == 2){
      p_spaq_cum[j] = alpha[2] * age_2 * llin * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 1  && age[j] == 2){
     p_spaq_cum[j] = alpha[3] * age_2 * llin * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 1  && age[j] == 2){
     p_spaq_cum[j] = alpha[4] * age_2 * llin * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 1  && age[j] == 2){
     p_spaq_cum[j] = alpha[5] * age_2 * llin * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }


//-------AGE GROUP 3-------------------------------------------------------------------

  // No bed-net site 2
  else if(interval[j] == 1 && bed_net[j] == 0  && age[j] == 3){
      p_spaq_cum[j] = alpha[1] * age_3 * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 0  && age[j] == 3){
      p_spaq_cum[j] = alpha[2] * age_3 * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 0  && age[j] == 3){
     p_spaq_cum[j] = alpha[3] * age_3 * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 0  && age[j] == 3){
     p_spaq_cum[j] = alpha[4] * age_3 * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 0  && age[j] == 3){
     p_spaq_cum[j] = alpha[5] * age_3 * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }
 
 // yes bet net site 2 
  else if(interval[j] == 1 && bed_net[j] == 1  && age[j] == 3){
      p_spaq_cum[j] = alpha[1] * age_3 * llin * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 1  && age[j] == 3){
      p_spaq_cum[j] = alpha[2] * age_3 * llin * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 1  && age[j] == 3){
     p_spaq_cum[j] = alpha[3] * age_3 * llin * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 1  && age[j] == 3){
     p_spaq_cum[j] = alpha[4] * age_3 * llin * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 1  && age[j] == 3){
     p_spaq_cum[j] = alpha[5] * age_3 * llin * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }
 
 //--------AGE GROUP 4-----------------------------------------------------------------

  // No bed-net site 2
  else if(interval[j] == 1 && bed_net[j] == 0  && age[j] == 4){
      p_spaq_cum[j] = alpha[1] * age_4 * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 0  && age[j] == 4){
      p_spaq_cum[j] = alpha[2] * age_4 * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 0  && age[j] == 4){
     p_spaq_cum[j] = alpha[3] * age_4 * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 0  && age[j] == 4){
     p_spaq_cum[j] = alpha[4] * age_4 * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 0  && age[j] == 4){
     p_spaq_cum[j] = alpha[5] * age_4 * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }
 
 // yes bet net site 2 
  else if(interval[j] == 1 && bed_net[j] == 1  && age[j] == 4){
      p_spaq_cum[j] = alpha[1] * age_4 * llin * sum(p_spaq_t[1:time_id[j]]);
  }
  else if(interval[j] == 2 && bed_net[j] == 1  && age[j] == 4){
      p_spaq_cum[j] = alpha[2] * age_4 * llin * sum(p_spaq_t[break_ids[1]:time_id[j]]);
  }
 else if(interval[j] == 3 && bed_net[j] == 1  && age[j] == 4){
     p_spaq_cum[j] = alpha[3] * age_4 * llin * sum(p_spaq_t[break_ids[2]:time_id[j]]);
 }
 else if(interval[j] == 4 && bed_net[j] == 1  && age[j] == 4){
     p_spaq_cum[j] = alpha[4] * age_4 * llin * sum(p_spaq_t[break_ids[3]:time_id[j]]);
 }
 else if(interval[j] == 5 && bed_net[j] == 1  && age[j] == 4){
     p_spaq_cum[j] = alpha[5] * age_4 * llin * sum(p_spaq_t[break_ids[4]:time_id[j]]);
 }

 

}

 }

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
  
  RT ~ customProb_smc( indicator, 
                       interval,
                       treatment_id,
                       bed_net,
                       age, 
                       p_spaq_cum,
                       alpha, 
                       llin, 
                       age_0, 
                       age_1, 
                       age_2, 
                       age_3, 
                       age_4, 
                       lambda,
                       k); 

}

//log likelihood for WAIC
generated quantities{
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = customProb_smc_lpdf(rep_vector(RT[n],1)| rep_vector(indicator[n],1), rep_vector(interval[n],1),  rep_vector(treatment_id[n],1), rep_vector(bed_net[n],1),rep_vector(age[n],1), rep_vector(p_spaq_cum[n],1), alpha, llin, age_0, age_1, age_2, age_3, age_4, lambda, k); //rep_vector(age_1[n],1), rep_vector(age_2[n],1), rep_vector(age_3[n],1), rep_vector(age_4[n],1), rep_vector(site_2[n],1), rep_vector(site_3[n],1),beta_llin, beta_age, beta_site,
  }
  }

