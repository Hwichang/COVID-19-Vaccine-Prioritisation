########################################################################################################
################################ - Japan Bayesian Inference - #########################################
########################################################################################################
library(optiSolve)
library(quadprog)
library(polynom)
library(logitnorm)
library(rGammaGamma)
library(stats)
library(STAR)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)
library(MCMCpack)
rm(list=ls())
gc()


#data
japan_age_group = as.vector(read.csv('japan_population.csv',header=T)$x)
pre = as.vector(read.csv('japan_pre.csv')$x) 
japan_corona_new = read.csv('japan_corona_daily.csv')
contact_matrix = as.matrix(read.csv('contact_japan_res.csv',header=F))

#######################InCUBATION PERIOD################
Incu_param1 = 4.544
Incu_param2 = 1/0.709


######################TRANSMISSION ONSET#################
tran_dist_mu = -4
tran_param1 = 5.2662158
tran_param2 = 1/0.8709042 


##################Infection to Recover#####################
I_R_param1 = 4
I_R_param2 = 4/5


##################Infection to Quarantine######################
C_param = 1.7


#####################Symptom to Quarantine#####################
symp_q_dist = read.csv('symp_q_dist.csv')$x

q_mean = ceiling(mean(symp_q_dist))
C_mean = ceiling(1/1.7)
incu_mean = ceiling(4.544*0.709)
I_R_mean = 5
tran_mean = ceiling(-4+(tran_param1/tran_param2))


#####################Set seed#################################
set.seed(0814)


####################initial value of W#######################
seiq_matrix = read.csv('initial_seiq_japan.csv')
#seiq_matrix = read.csv('initial_seiq_japan_asym0.04.csv')
#seiq_matrix = read.csv('initial_seiq_japan_asym0.4.csv')

seiq_matrix = seiq_matrix[,2:6]


##################initial susceptible########################
suscept = matrix(0,nrow = 400 , ncol=16)
for(t in 1:nrow(suscept)){
  suscept[t,] = sapply(1:16,function(x) {japan_age_group[x]- pre[x] - length(which((seiq_matrix$E_date<t)&(seiq_matrix$age/5+1==x)))})
}
suscept_new = suscept

###############initial Exposed###############################
exposed = matrix(0,nrow = 400,ncol=16)
for(t in 1:nrow(exposed)){
  exposed[t,] = sapply(1:16,function(x) {length(which((seiq_matrix$E_date==t)&(seiq_matrix$age/5+1==x)))})
}


##############initial Infectious#############################
I_total_sym = matrix(0,nrow = 400,ncol=16)
for(t in 1:400){
  I_total_sym[t,] = sapply(1:16,function(x){length(which((seiq_matrix$I_date<=t)&(seiq_matrix$Q_date>t)&(seiq_matrix$age/5+1==x)&(!is.na(seiq_matrix$Y_date))))})
}

I_total_asym = matrix(0,nrow = 400 ,ncol=16)
for(t in 1:400){
  I_total_asym[t,] = sapply(1:16,function(x){length(which((seiq_matrix$I_date<=t)&(seiq_matrix$Q_date>t)&(seiq_matrix$age/5+1==x)&(is.na(seiq_matrix$Y_date))))})
}


I_total = I_total_sym + 0.5*I_total_asym
I_total_new = I_total

#######################pi_matrix###############################
naive_p = function(t){
  res = - (1/(japan_age_group*100)) * (contact_matrix%*%I_total[t-1,])
  return(res)
}


init_q = matrix(0, nrow =3, ncol=16)
temp_q = (exposed[2:nrow(exposed),]/suscept[1:(nrow(exposed)-1),])/t(sapply(2:nrow(exposed), function(x) {-naive_p(x)}))
temp_q = replace(temp_q , is.infinite(temp_q ),NaN)

q_standard = c(1,289,361,390)
for( i in 1:3){
  temp = temp_q[q_standard[i]:(q_standard[i+1]-1),]
  temp = na.omit(temp)
  init_q[i,] = colMeans(temp)
}


init_q[2,] = sapply(1:16, function(x) { max(rgamma(1,0.001,rate = 0.001),1e-300)})  


q_matrix = matrix(0, nrow = 400,ncol = 16)
q_matrix[1,] = init_q[1,]
for( i in 1:3){
  for(j in q_standard[i]:(q_standard[i+1]-1)){
    q_matrix[j,] = init_q[i,]
  }
}

############################ p(t)#############################
## p_unit : force of infection at time t
## p_unit_new : force of infection at time t divided by q
## tilde_p_unit : force of infection at time t about newly imputed value
## tilde_p_unit_new : force of infection at time t divided by q newly imputed value
p_unit = function(t){
  res = (q_matrix[(t-1),]/(japan_age_group*100)) * (contact_matrix%*%I_total[t-1,])
  return(res)
}

p_unit_new = function(t){
  res = (1/(japan_age_group*100)) * (contact_matrix%*%I_total[t-1,])
  return(res)
}

tilde_p_unit = function(t){
  res = (q_matrix[(t-1),]/(japan_age_group*100)) * (contact_matrix%*%I_total_new[t-1,])
  return(res)
}

tilde_p_unit_new = function(t){
  res = (1/(japan_age_group*100)) * (contact_matrix%*%I_total_new[t-1,])
  return(res)
}

################ For caculate P(E=t) ##########################
p = function(s,t){
  res = -colSums(t(sapply(s:t,function(x) {p_unit(x)})))
  return(res)
}

p_new = function(s,t){
  res = -colSums(t(sapply(s:t,function(x) {p_unit_new(x)})))
  return(res)
}

tilde_p = function(s,t){
  res = -colSums(t(sapply(s:t,function(x) {tilde_p_unit(x)})))
  return(res)
}

tilde_p_new = function(s,t){
  res = -colSums(t(sapply(s:t,function(x) {tilde_p_unit_new(x)})))
  return(res)
}

# 289 #start
# 359 #end
start_date = 289 # 2020-10-15
end_date = 360 # 2020-12-25

start = which(seiq_matrix$Q_date>289)[1]
end = which(seiq_matrix$Q_date<=370)[length(which(seiq_matrix$Q_date<=370))]
seiq_temp = seiq_matrix[start:end,]
q_list = matrix(0, nrow=1100,ncol=16)

k=1

for(iter in 1:1100){
  start_time <- Sys.time()
  ######################################### - W update - ##########################################################
  for( i in 1:(end-start)){
    Q_i = seiq_temp$Q_date[i]
    i_age = seiq_temp$age[i]/5+1
    E_old = seiq_temp$E_date[i]
    I_old = seiq_temp$I_date[i]
    alpha = 0
    
    if(!is.na(seiq_temp$Y_date[i])){
      c = 1
      Y_new = Q_i -  sample(symp_q_dist,1)
      E_new = Y_new -  ceiling(rgamma(1,Incu_param1,Incu_param2))
      I_new = max(Y_new + ceiling(tran_dist_mu + rgamma(1,tran_param1,tran_param2)),E_new)
      
    }else{
      c = 0.5
      Y_new = NA
      I_new = Q_i - ceiling(rexp(1,C_param))
      E_new = I_new - max(ceiling(rgamma(1,Incu_param1,Incu_param2) + tran_dist_mu + rgamma(1,tran_param1,tran_param2)),0)
      
    }
    
    exposed[E_old,i_age] = exposed[E_old,i_age] - 1 #substract i
    
    m_E = min(E_old,E_new)
    M_E = max(E_old,E_new)
    
    m_I = min(I_old,I_new)
    M_I = max(I_old,I_new)
    
    
    if(E_old > E_new){
      suscept_new[E_new:(E_old-1),i_age] = sapply(E_new:(E_old-1), function(x){suscept_new[x,i_age]-1})
    }else if(E_old < E_new){
      suscept_new[E_old:(E_new-1),i_age] = sapply(E_old:(E_new-1), function(x){suscept_new[x,i_age]+1})
    }
    
    if((I_old >= Q_i) & (I_new < Q_i)){
      I_total_new[I_new:(Q_i-1),i_age] = sapply(I_new:(Q_i-1), function(x){I_total_new[x,i_age]+c})
      
    }else if((I_old > I_new) & (I_old < Q_i)){
      I_total_new[I_new:(I_old-1),i_age] = sapply(I_new:(I_old-1), function(x){I_total_new[x,i_age]+c})
      
    }else if((I_old < I_new) & (I_new < Q_i)){
      I_total_new[I_old:(I_new-1),i_age] = sapply(I_old:(I_new-1), function(x){I_total_new[x,i_age]-c})
      
    }else if((I_new >= Q_i) & (I_old < Q_i)){
      I_total_new[I_old:(Q_i-1),i_age] = sapply(I_old:(Q_i-1), function(x){I_total_new[x,i_age]-c})
    }
    
    
    ##################### W_(-i) MCMC ratio ############################
    age_num_B = suscept_new[M_I,]
    temp_B = (tilde_p(m_I+1,M_I) - p(m_I+1,M_I))*age_num_B
    alpha = alpha + sum(temp_B)
    
    age_num_A = exposed[(m_I+1):M_I,]
    temp_A = t(sapply((m_I+1):M_I , function(y) { log(tilde_p_unit(y)) + tilde_p(m_I,y-1) - log(p_unit(y)) - p(m_I,y-1)}))
    alpha = alpha + sum(temp_A*age_num_A)
    
    ######################### W_i MCMC ratio #########################################################################
    alpha = alpha + (log(tilde_p_unit(E_new)) + tilde_p(m_E-1,E_new-1) - log(p_unit(E_old)) - p(m_E-1,E_old-1))[i_age]
    
    
    ########################### accept or reject ###################################################################
    accept = runif(1,0,1)
    if( log(accept)  <= alpha){
      seiq_temp[i,3:5] = c(I_new,Y_new,E_new)
      I_total = I_total_new
      suscept = suscept_new
      exposed[E_new,i_age] = exposed[E_new,i_age]+1
    }else{
      I_total_new = I_total
      suscept_new = suscept
      exposed[E_old,i_age] = exposed[E_old,i_age]+1
    }
    
  }
  ########################################### - pi update - ##########################################################
  data_num = exposed[start_date:end_date,]
  
  data_vec = colSums(t(sapply(start_date:end_date, function(x)  {-p_new(start_date,x)}))*data_num)
  data_vec = data_vec - p_new(start_date,end_date)* suscept[end_date,]
  
  q_new = sapply(1:16, function(x) rgamma( 1 , 0.001 + colSums(data_num)[x] , rate = 0.001 + data_vec[x]))
  q_list[k,] = q_new
  
  q_matrix[start_date:end_date,] = t(sapply(start_date:end_date , function(x) {q_new}))
  print(q_matrix[start_date,])
  
  k = k+1 
  
  end_time <- Sys.time()
  print(end_time - start_time)
}

seiq_matrix[start:end,] = seiq_temp