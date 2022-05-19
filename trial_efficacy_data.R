trial_eff <- data.frame(
  time = c(0,  #3.5,
           7,  #10.5,
           14, #17.5,
           21, #24.5,
           28, #31.5,
           35, #38.5,
           42, #45.5, #46.2,
           49, #50.75,#52.5, #54.25,
           56),
  
  eff = c(1,     #1,
           1,     #0.9925,
           0.965, #0.9175,
           0.8575,#0.79,
           0.72,  #0.6575,
           0.5675,#0.435,
           0.27,  #0.055, #0,
           0,     #0, #0, #0,
           0)
  )

w_scale<- 38.07117
w_shape<- 4.3151847
t<-seq(0,60,0.5)
weib<-exp(-(t/w_scale)^w_shape)   # cumulative weibull
#plot(t,weib, col="black", type="l", lwd=2)

pw_solver <- data.frame(time = t, eff = weib)

w_scale<- 41.48 
w_shape<- 3.67 
t<-seq(0,70,0.5)
weib<-exp(-(t/w_scale)^w_shape)   # cumulative weibull
#plot(t,weib, col="black", type="l", lwd=2)



ht_orig <- data.frame(time = t, eff = weib)
