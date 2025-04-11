
# Rdata contains this dat object
# dat$Ct = as.vector(df_effort$catch_otb)  # catches
# dat$Et = as.vector(df_effort$effort_otb) # effort

# Simplified Schaefer SPM, for learning purposes 
# Bt = Latent Biomass
# Ct = Observed Catch
# Et = Observed Effort
# r, q, K, B0 = fixed parameters

# initialize parameter list with initial estimates
par = list()

# Empty 
par$Bt = as.vector(rep(0,length(dat$Ct)))# init empty vector


par$q = 0.5
par$r = 0.5
par$K = 4*max(dat$Ct)
par$B0 = 0.5*par$K # % of stock depletion at time zero

# What if we wanted a "realistic" starting point for B[t]?
for(i in 1:length(dat$Ct)){
  if(i == 1) {par$Bt[i] = par$B0}
  else{par$Bt[i] = par$Bt[i-1] + par$r * 
    par$Bt[i-1]*(1-(par$Bt[i-1]/par$K)) - dat$Ct[i-1]}}

par$logsdBt = 1
par$logsdCt = 1
par$logsdB0 = 1000

# initialize joint negative loglikelihood function
jnll = function(par){
  getAll(par, dat) # invoke all the parameters and data
  
  # gather and cleanup pars
  Ct = OBS(Ct)
  sdBt = exp(logsdBt)
  sdCt = exp(logsdCt)
  sdB0 = exp(logsdB0)
  
  # assemble model    
  jnll = 0 # init jnll
  
  ## Prior distributions
  ## The 1 in the mean fixes the prior distribution of r
  jnll = jnll - dexp(r,1, log = TRUE)
  jnll = jnll - dexp(q,0.1, log = TRUE)
  jnll = jnll - dnorm(K, 1000000, 1000, log = TRUE)

  ## estimation
  biomass_predictions = c(B0)
  catch_predictions = c()
  
  for(i in 1:length(Ct)){
    if(i == 1) {predBt = B0 + r * B0*(1-(B0/K))}
    else{predBt = predBt + r * predBt*(1-(predBt/K)) - Ct[i-1]}
    jnll = jnll - dnorm(Bt[i], predBt, sdBt, log=TRUE)
    biomass_predictions[i] = predBt}
  
  
  for(i in 1:length(Ct)){
    Ct_hat = q*Bt[i]*Et[i]
    jnll = jnll -dnorm(Ct[i], Ct_hat, sdCt, log = TRUE)
    catch_predictions[i] = Ct_hat
  }
  
  REPORT(biomass_predictions)
  REPORT(catch_predictions)
  REPORT(B0)
  jnll
}

# quick test: do we get a number? This number should be a likelihood.
jnll(par)

obj4 = MakeADFun(jnll, par)
fit4 = nlminb(obj4$par, obj4$fn, obj4$gr)

sdr4 = sdreport(obj4)
pl4 = as.list(sdr4,"Est")
plsd4 = as.list(sdr4,"Std")

# sdr4
obj4$report()$biomass_predictions
plot(obj4$report()$biomass_predictions)

sdr4

# At this point, while biomass predictions have been iterated upon
# In the way I expected, it appears I am failing to get RTMB to take those estimates
# into account when updating B[t]