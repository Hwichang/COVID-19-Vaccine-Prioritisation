########################################################################################################
################################ - Japan Meta Modeling - #########################################
########################################################################################################
library(optiSolve)
library(quadprog)
library(polynom)
library(logitnorm)
library(rGammaGamma)
library(stats)
library(STAR)
library(dplyr)
library(matrixStats)
library(MCMCpack)

rm(list=ls())
gc()

#data
japan_age_group = as.vector(read.csv('japan_population.csv',header=T)$x)
pre = as.vector(read.csv('japan_pre.csv')$x)
japan_corona_new = read.csv('japan_corona_daily.csv')
contact_matrix = as.matrix(read.csv('contact_japan_res.csv',header=F))

####################### Incubation Period ################
Incu_param1 = 4.544
Incu_param2 = 1/0.709


###################### Transmission Onset #################
tran_dist_mu = -4
tran_param1 = 5.2662158
tran_param2 = 1/0.8709042 


################## Infection to Recover #####################
I_R_param1 = 4
I_R_param2 = 4/5


################## Infection to Quarantine ######################
C_param = 1.7


##################### Symptom to Quarantine #####################
symp_q_dist = read.csv('symp_q_dist.csv')$x

q_mean = ceiling(mean(symp_q_dist))
C_mean = ceiling(1/1.7)
incu_mean = ceiling(4.544*0.709)
I_R_mean = 5
tran_mean = ceiling(-4+(tran_param1/tran_param2))



q_list = as.matrix(read.csv('q_list_japan.csv',header=T)[,2:17])
#q_list = as.matrix(read.csv('q_list_japan_asym_0.04.csv',header=T)[,2:17])
#q_list = as.matrix(read.csv('q_list_japan_asym_0.4.csv',header=T)[,2:17])

q_res = colMedians(q_list[100:1100,])/100
#q_res = as.vector(sapply(1:16,function(x) quantile(q_list[100:1100,x],0.025))/100)
#q_res = as.vector(sapply(1:16,function(x) quantile(q_list[100:1100,x],0.975))/100)

seiq_matrix = read.csv('seiq_matrix_japan.csv')[,2:6]
suscept = as.matrix(read.csv('suscept_japan.csv')[,2:17])

seiq_res = seiq_matrix[seiq_matrix$E_date<=360,] # remove people who were exposed after 12/25

I_total_sym_res =  matrix(0,nrow = 15 ,ncol=16) # Infectious people who were exposed before 12/25
I_total_asym_res =  matrix(0,nrow = 15 ,ncol=16) # Infectious people who were exposed before 12/25
for(t in 1:15){
  I_total_sym_res[t,] = sapply(1:16,function(x){length(which((seiq_res$I_date<=t+359)&(seiq_res$Q_date>t+359)&(seiq_res$age/5+1==x)&(!is.na(seiq_res$Y_date))))})
}

for(t in 1:15){
  I_total_asym_res[t,] = sapply(1:16,function(x){length(which((seiq_res$I_date<=t+359)&(seiq_res$Q_date>t+359)&(seiq_res$age/5+1==x)&(is.na(seiq_res$Y_date))))})
}

I_total_res = I_total_sym_res + 0.5*I_total_asym_res


suscept_res = matrix(0, nrow=15, ncol=16)
suscept_res[1,] = suscept[360,] #susceptible at 12/25


exposed_res = matrix(0,nrow=14,ncol=16)
I_total_res_meta = I_total_res
for( i in 1:14){
  exposed_res[i,] = round(sapply(1:16, function(x) {suscept_res[i,x]*(q_res[x]/(japan_age_group[x]))*(contact_matrix[x,]%*%I_total_res_meta[i,])}))
  suscept_res[(i+1),] = suscept_res[i,] - exposed_res[i,]
  
  
  if((i+tran_mean+incu_mean)<=15){
    if((i+q_mean + incu_mean-1)<=15){
      I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + round(0.84*exposed_res[i,])
    }else{
      I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.84*exposed_res[i,])
    }
    if((i+C_mean + incu_mean-1)<=15){
      I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + round(0.5*0.16*exposed_res[i,])
    }else{
      I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.5*0.16*exposed_res[i,])
    }
  }
}

#vaccine efficacy (Baden LR, El Sahly HM, Essink B, et al. Efficacy and Safety of the mRNA-1273 SARS-CoV-2 Vaccine. The New England journal of medicine 2021; 384(5): 403-16)
efficacy_sym = c(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.86,0.86,0.86) #vaccine efficacy
#efficacy_sym = c(0.475,0.475,0.475,0.475,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.86,0.86,0.86) #vaccine efficacy

### We assumed the COVID-19 vaccine efficacy against asymptomatic infection as 70% (Tande AJ, Pollock BD, Shah ND, et al. Impact of the COVID-19 Vaccine on Asymptomatic Infection Among Patients Undergoing Pre-Procedural COVID-19 Molecular Screening. Clinical Infectious Diseases 2021). 
efficacy_asym = c(0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.9*0.7,0.9*0.7,0.9*0.7) #vaccine efficacy
mortality = c(0,0,0,0,0,0,0,0,0.1,0.1,0.3,0.3,1.3,1.3,4.4,9.51)/100 #mortality


sum(colSums(exposed_res)) #Expected number of exposed people in coming 2weeks after the 3rd wave
sum(round(colSums(exposed_res)*mortality)) #Expected number of deaths in coming 2weeks after the 3rd wave


num_infection_1 = matrix(0,nrow=20,ncol=17)
num_death_1 = matrix(0,nrow=20,ncol=17)

num_infection_2 = matrix(0,nrow=20,ncol=17)
num_death_2 = matrix(0,nrow=20,ncol=17)

num_infection_3 = matrix(0,nrow=20,ncol=17)
num_death_3 = matrix(0,nrow=20,ncol=17)

num_infection_4 = matrix(0,nrow=20,ncol=17)
num_death_4 = matrix(0,nrow=20,ncol=17)


K = 16 # The number of age group.
N_vac_list = c(0.2,0.4,0.6,0.8)*sum(japan_age_group)

### Change asymptomatic rate when sensitive analysis.
for( N_vac in N_vac_list){
  if(N_vac == 0.8*sum(japan_age_group)){
    d = 5000
    num_sam = 200000
    end = 20
  }else if(N_vac == 0.6*sum(japan_age_group)){
    d = 10000
    num_sam = 200000
    end = 20
  }else{
    d = 100000
    num_sam = 50000
    end = 20
  }
  
  for( iterate in 1:end){
    sample_p = matrix(0, nrow=num_sam, ncol=K)
    
    for ( i in 1:nrow(sample_p)){
      sample_p[i,] = as.vector(rdirichlet(1,japan_age_group/d))
    }
    
    exposed_num = rep(Inf, nrow(sample_p))
    exposed_mor = rep(Inf, nrow(sample_p))
    
    for( j in 1:nrow(sample_p)){
      vac_num = N_vac*sample_p[j,]
      if(max(vac_num-japan_age_group)<=0){
        Immunity_sym = efficacy_sym*(vac_num)
        Immunity_asym = efficacy_asym*(vac_num)
        
        suscept_res_sym_meta = matrix(0, nrow=15, ncol=16)
        suscept_res_asym_meta = matrix(0, nrow=15, ncol=16)
        suscept_res_sym_meta[1,] = suscept[360,] - Immunity_sym
        suscept_res_asym_meta[1,] = suscept[360,] - Immunity_asym
        
        I_total_res_meta = I_total_res
        
        exposed_res_sym_meta = matrix(0,nrow=14,ncol=16)
        exposed_res_asym_meta = matrix(0,nrow=14,ncol=16)
        
        for( i in 1:14){
          exposed_res_sym_meta[i,] = round(sapply(1:16, function(x) {0.84 * suscept_res_sym_meta[i,x]*(q_res[x]/(japan_age_group[x]))*(contact_matrix[x,]%*%I_total_res_meta[i,])}))
          exposed_res_asym_meta[i,] = round(sapply(1:16, function(x) {0.16 * suscept_res_asym_meta[i,x]*(q_res[x]/(japan_age_group[x]))*(contact_matrix[x,]%*%I_total_res_meta[i,])}))
          suscept_res_sym_meta[(i+1),] = suscept_res_sym_meta[i,] - exposed_res_sym_meta[i,] - exposed_res_asym_meta[i,]
          suscept_res_asym_meta[(i+1),] = suscept_res_asym_meta[i,] - exposed_res_sym_meta[i,] - exposed_res_asym_meta[i,]
          
          if((i+tran_mean+incu_mean)<=15){
            if((i+q_mean + incu_mean-1)<=15){
              I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + exposed_res_sym_meta[i,]
            }else{
              I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + exposed_res_sym_meta[i,]
            }
            if((i+C_mean + incu_mean-1)<=15){
              I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + 0.5*exposed_res_asym_meta[i,]
            }else{
              I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + 0.5*exposed_res_asym_meta[i,]
            }
          }
        }
        exposed_num[j] = sum(colSums(exposed_res_sym_meta + exposed_res_asym_meta))
        exposed_mor[j] = sum(round(colSums(exposed_res_sym_meta + exposed_res_asym_meta)*mortality))
      }
    }
    
    
    res_num = rep(0,16)
    res_mor = rep(0,16)
    
    m=16 #Because sum of p is 1, we remove one element.  
    
    
    #The number of vaccines for each group should not exceed the population size of that group. 
    realsample = sample_p[which(exposed_num!=Inf),-m] 
    realexposed_num = exposed_num[which(exposed_num!=Inf)]/100000
    realexposed_mor = exposed_mor[which(exposed_mor!=Inf)]/10000
    real_num = cbind(realsample,realexposed_num)
    real_mor = cbind(realsample,realexposed_mor)
    real_num = data.frame(real_num)
    real_mor = data.frame(real_mor)
    
    lmfit_num = lm(realexposed_num ~ polym(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,degree=2,raw=TRUE),data=real_num)
    summary(lmfit_num)
    coeff_num = lmfit_num$coefficients
    
    # Setting to solve quadratic programming.
    dvec = coeff_num[sapply(1:15,  function(x)((x^2+x+2)/2)  )]
    Dmat = diag(coeff_num[sapply(1:15,  function(x)((x^2+3*x+2)/2)  )])
    off_diag = matrix(0, nrow=15, ncol=15)
    for( i in 1:14){
      off_diag[(1:i),(i+1)] = coeff_num[seq( from= (i^2+3*i+6)/2  , to= ((i^2+3*i+6)/2)+(i-1) , by=1)]
    }
    
    off_diag = off_diag + t(off_diag)
    
    Dmat = off_diag + 2*Dmat
    bvec = sapply(c(1:16)[-m], function(x) min((japan_age_group/N_vac)[x],1))
    bvec = c(1-(japan_age_group/N_vac)[m],-1,-bvec,rep(0,15))
    Amat = rbind( diag(-1,nrow=15),diag(1,nrow=15))
    Amat = rbind(rep(1,15),rep(-1,15),Amat)
    
    
    # Quadratic programming using cccp
    f= quadfun((1/2)*Dmat, a=dvec, name="quad.fun")
    
    l=lincon(Amat, dir=rep(">=",nrow(Amat)), val=bvec,
             use=rep(TRUE,nrow(Amat)), name=seq(1:32))
    
    op = cop(f, lc=l)
    
    b = solvecop(op,quiet=TRUE)
    
    
    # Only convergence case
    if(b$status == "optimal"){
      
      for ( n in 1:16){
        if(n<m){
          res_num[n] = b$x[n]
        }else if(n>m){
          
          res_num[n] = b$x[n-1]
        }
        
      }
      
      res_num[m] = 1-sum(res_num)
    }  
    
    vac_num_res = N_vac*res_num
    sum(vac_num_res)
    if(max(vac_num_res-japan_age_group)<=0){
      Immunity_sym = efficacy_sym*(vac_num_res)
      Immunity_asym = efficacy_asym*(vac_num_res)
      
      suscept_res_sym_meta = matrix(0, nrow=15, ncol=16)
      suscept_res_asym_meta = matrix(0, nrow=15, ncol=16)
      suscept_res_sym_meta[1,] = suscept[360,] - Immunity_sym
      suscept_res_asym_meta[1,] = suscept[360,] - Immunity_asym
      
      I_total_res_meta = I_total_res
      
      exposed_res_sym_meta = matrix(0,nrow=14,ncol=16)
      exposed_res_asym_meta = matrix(0,nrow=14,ncol=16)
      for( i in 1:14){
        exposed_res_sym_meta[i,] = round(sapply(1:16, function(x) {0.84 * suscept_res_sym_meta[i,x]*(q_res[x]/(japan_age_group[x]))*(contact_matrix[x,]%*%I_total_res_meta[i,])}))
        exposed_res_asym_meta[i,] = round(sapply(1:16, function(x) {0.16 * suscept_res_asym_meta[i,x]*(q_res[x]/(japan_age_group[x]))*(contact_matrix[x,]%*%I_total_res_meta[i,])}))
        suscept_res_sym_meta[(i+1),] = suscept_res_sym_meta[i,] - exposed_res_sym_meta[i,] - exposed_res_asym_meta[i,]
        suscept_res_asym_meta[(i+1),] = suscept_res_asym_meta[i,] - exposed_res_sym_meta[i,] - exposed_res_asym_meta[i,]
        if((i+tran_mean+incu_mean)<=15){
          if((i+q_mean + incu_mean-1)<=15){
            I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + exposed_res_sym_meta[i,]
          }else{
            I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + exposed_res_sym_meta[i,]
          }
          if((i+C_mean + incu_mean-1)<=15){
            I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + 0.5*exposed_res_asym_meta[i,]
          }else{
            I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + 0.5*exposed_res_asym_meta[i,]
          }
        }
      }
      final_num = sum(colSums(exposed_res_sym_meta+exposed_res_asym_meta))
    }
    print(final_num)
    
    if(N_vac == 0.2*sum(japan_age_group)){
      num_infection_1[iterate,1:16] = vac_num_res
      num_infection_1[iterate,17] = final_num
      print(num_infection_1[iterate,])
    }else if(N_vac == 0.4*sum(japan_age_group)){
      num_infection_2[iterate,1:16] = vac_num_res
      num_infection_2[iterate,17] = final_num
      print(num_infection_2[iterate,])
    }else if(N_vac == 0.6*sum(japan_age_group)){
      num_infection_3[iterate,1:16] = vac_num_res
      num_infection_3[iterate,17] = final_num
      print(num_infection_3[iterate,])
    }else{
      num_infection_4[iterate,1:16] = vac_num_res
      num_infection_4[iterate,17] = final_num
      print(num_infection_4[iterate,])
    }
    
    
    
    lmfit_mor = lm(realexposed_mor~ polym(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,degree=2,raw=TRUE),data=real_mor)
    summary(lmfit_mor)
    coeff_mor = lmfit_mor$coefficients
    
    # Setting to solve quadratic programming.
    dvec = coeff_mor[sapply(1:15,  function(x)((x^2+x+2)/2)  )]
    Dmat = diag(coeff_mor[sapply(1:15,  function(x)((x^2+3*x+2)/2)  )])
    off_diag = matrix(0, nrow=15, ncol=15)
    for( i in 1:14){
      off_diag[(1:i),(i+1)] = coeff_mor[seq( from= (i^2+3*i+6)/2  , to= ((i^2+3*i+6)/2)+(i-1) , by=1)]
    }
    
    off_diag = off_diag + t(off_diag)
    
    Dmat = off_diag + 2*Dmat
    bvec = sapply(c(1:16)[-m], function(x) min((japan_age_group/N_vac)[x],1))
    bvec = c(1-(japan_age_group/N_vac)[m],-1,-bvec,rep(0,15))
    Amat = rbind( diag(-1,nrow=15),diag(1,nrow=15))
    Amat = rbind(rep(1,15),rep(-1,15),Amat)
    
    
    #Quadratic programming using cccp
    f= quadfun((1/2)*Dmat, a=dvec, name="quad.fun")
    
    l=lincon(Amat, dir=rep(">=",nrow(Amat)), val=bvec,
             use=rep(TRUE,nrow(Amat)), name=seq(1:32))
    
    op = cop(f, lc=l)
    
    b = solvecop(op,quiet=TRUE)
    
    
    #Only convergence case
    if(b$status == "optimal"){
      
      for ( n in 1:16){
        if(n<m){
          res_mor[n] = b$x[n]
        }else if(n>m){
          
          res_mor[n] = b$x[n-1]
        }
        
      }
      
      res_mor[m] = 1-sum(res_mor)
    }  
    
    res_vac_mor = N_vac*res_mor
    
    if(max(res_vac_mor-japan_age_group)<=0){
      Immunity_sym = efficacy_sym*(res_vac_mor)
      Immunity_asym = efficacy_asym*(res_vac_mor)
      
      suscept_res_sym_meta = matrix(0, nrow=15, ncol=16)
      suscept_res_asym_meta = matrix(0, nrow=15, ncol=16)
      suscept_res_sym_meta[1,] = suscept[360,] - Immunity_sym
      suscept_res_asym_meta[1,] = suscept[360,] - Immunity_asym
      
      I_total_res_meta = I_total_res
      
      exposed_res_sym_meta = matrix(0,nrow=14,ncol=16)
      exposed_res_asym_meta = matrix(0,nrow=14,ncol=16)
      for( i in 1:14){
        exposed_res_sym_meta[i,] = round(sapply(1:16, function(x) {0.84 * suscept_res_sym_meta[i,x]*(q_res[x]/(japan_age_group[x]))*(contact_matrix[x,]%*%I_total_res_meta[i,])}))
        exposed_res_asym_meta[i,] = round(sapply(1:16, function(x) {0.16 * suscept_res_asym_meta[i,x]*(q_res[x]/(japan_age_group[x]))*(contact_matrix[x,]%*%I_total_res_meta[i,])}))
        suscept_res_sym_meta[(i+1),] = suscept_res_sym_meta[i,] - exposed_res_sym_meta[i,] - exposed_res_asym_meta[i,]
        suscept_res_asym_meta[(i+1),] = suscept_res_asym_meta[i,] - exposed_res_sym_meta[i,] - exposed_res_asym_meta[i,]
        
        if((i+tran_mean+incu_mean)<=15){
          if((i+q_mean + incu_mean-1)<=15){
            I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + exposed_res_sym_meta[i,]
          }else{
            I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + exposed_res_sym_meta[i,]
          }
          if((i+C_mean + incu_mean-1)<=15){
            I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + 0.5*exposed_res_asym_meta[i,]
          }else{
            I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + 0.5*exposed_res_asym_meta[i,]
          }
        }
      }
      final_mor = sum(round(colSums(exposed_res_sym_meta+exposed_res_asym_meta)*mortality))
    }
    print(final_mor)
    
    if(N_vac == 0.2*sum(japan_age_group)){
      num_death_1[iterate,1:16] = res_vac_mor
      num_death_1[iterate,17] = final_mor
      print(num_death_1[iterate,])
    }else if(N_vac == 0.4*sum(japan_age_group)){
      num_death_2[iterate,1:16] = res_vac_mor
      num_death_2[iterate,17] = final_mor
      print(num_death_2[iterate,])
    }else if(N_vac == 0.6*sum(japan_age_group)){
      num_death_3[iterate,1:16] = res_vac_mor
      num_death_3[iterate,17] = final_mor
      print(num_death_3[iterate,])
    }else{
      num_death_4[iterate,1:16] = res_vac_mor
      num_death_4[iterate,17] = final_mor
      print(num_death_4[iterate,])
    }
    
  }
}
