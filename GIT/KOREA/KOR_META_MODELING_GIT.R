########################################################################################################
################################ - Korea Meta Modeling - #########################################
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

rm(list=ls())
gc()

load("contact_school.RData") 
load("contact_work.RData") 
load("contact_others.RData") 
load("contact_all (1).RData") 
skage.groups_new = as.vector(read.csv('korea_population.csv',header=T)$x)
corona_new = read.csv('corona_daily_new.csv')

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


#####################Set seed#################################
set.seed(0814)

###################Contact matrix#############################
contact_matrix = function(S){
  if(S>=286 & S<=328){
    res = as.matrix(contact_all$KOR) - (1/3)*as.matrix(contact_school$KOR) - (0.1)*as.matrix(contact_others$KOR)  #Social distancing 1
  }else{
    res = as.matrix(contact_all$KOR) - (2/3)*as.matrix(contact_school$KOR) - (0.23)*as.matrix(contact_others$KOR) #Social distancing 2
  }
  return(res)
}


tran_mean = ceiling(-4+(tran_param1/tran_param2))
q_mean = ceiling(mean(symp_q_dist))
C_mean = ceiling(1/1.7)
incu_mean = ceiling(4.544*0.709)

pi_list = as.matrix(read.csv('pi_list_210227.csv'))
#pi_list = as.matrix(read.csv('pi_list_asym_0.04.csv')[,2:17]) #asymtomatic 0.04
#pi_list = as.matrix(read.csv('pi_list_asym_0.4.csv')[,2:17]) #asymtomatic 0.4

pi_res = colMedians(pi_list[100:1100,])/100
#pi_res = as.vector(sapply(1:16,function(x) quantile(pi_list[100:1100,x],0.975))/100) #upper CI
#pi_res = as.vector(sapply(1:16,function(x) quantile(pi_list[100:1100,x],0.025))/100) #lower CI

seiq_matrix = read.csv('seiq_matrix_210227.csv')[,2:6]
suscept = as.matrix(read.csv('suscept_210227.csv')[,2:17])

seiq_res = seiq_matrix[seiq_matrix$E_date<=357,] # remove people who exposed after 12/22

I_total_sym_res = matrix(0,nrow = 15,ncol=16) # Symptomatic infectious people who exposed before 12/22
for(t in 1:nrow(I_total_sym_res)){
  I_total_sym_res[t,] = sapply(1:16,function(x){length(which((seiq_res$I_date<=t+356)&(seiq_res$Q_date>t+356)&(seiq_res$age/5+1==x)&(!is.na(seiq_res$Y_date))))})
}

I_total_asym_res = matrix(0,nrow = 15,ncol=16) # Asymptomatic infectious people who exposed before 12/22
for(t in 1:nrow(I_total_asym_res)){
  I_total_asym_res[t,] = sapply(1:16,function(x){length(which((seiq_res$I_date<=t+356)&(seiq_res$Q_date>t+356)&(seiq_res$age/5+1==x)&(is.na(seiq_res$Y_date))))})
}

I_total_res = I_total_sym_res + 0.5*I_total_asym_res



suscept_res = matrix(0, nrow=15, ncol=16)
suscept_res[1,] = suscept[357,] #susceptible at 12/22

exposed_res = matrix(0,nrow=14,ncol=16)
I_total_res_meta = I_total_res
for( i in 1:14){
  exposed_res[i,] = sapply(1:16, function(x) {suscept_res[i,x]*(pi_res[x]/(skage.groups_new[x]))*(contact_matrix(357)[x,]%*%I_total_res_meta[i,])})
  suscept_res[(i+1),] = suscept_res[i,] - exposed_res[i,]
  
  if((i+tran_mean+incu_mean)<=15){
    if((i+q_mean + incu_mean-1)<=15){
      I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + round(0.84*exposed_res[i,])
    }else{
      I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.84*exposed_res[i,])
    }
    if((i+C_mean + incu_mean-1)<=15){
      I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + round(0.5 * 0.16*exposed_res[i,])
    }else{
      I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.5 * 0.16*exposed_res[i,])
    }
  }
}


efficacy = c(1,1,1,1,1,1,1,1,0.8,0.8,0.8,0.7,0.7,0.5,0.5,0.5) #vaccine efficacy
#efficacy = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) #vaccine efficacy
mortality = c(0,0,0,0,0,0,0.06,0.06,0.1,0.1,0.31,0.31,1.34,1.34,6.51,15.43)/100 #mortality

sum(round(colSums(exposed_res))) #Expected number of exposed people during 2weeks. 
sum(round(colSums(exposed_res)*mortality)) #Expected number of death during 2weeks.



num_infection_1 = matrix(0,nrow=20,ncol=17)
num_death_1 = matrix(0,nrow=20,ncol=17)

num_infection_2 = matrix(0,nrow=20,ncol=17)
num_death_2 = matrix(0,nrow=20,ncol=17)

num_infection_3 = matrix(0,nrow=20,ncol=17)
num_death_3 = matrix(0,nrow=20,ncol=17)

num_infection_4 = matrix(0,nrow=20,ncol=17)
num_death_4 = matrix(0,nrow=20,ncol=17)


suscept_res_meta = matrix(0, nrow=15, ncol=16)
exposed_res_meta = matrix(0,nrow=14,ncol=16)
contact_matrix_res = contact_matrix(357)
C = pi_res/skage.groups_new
K = 16 # The number of age group.

N_vac_list = round(c(0.2,0.4,0.6,0.8)*sum(skage.groups_new))

for( N_vac in N_vac_list){
  
  if(N_vac == N_vac_list[4]){
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
      sample_p[i,] = as.vector(rdirichlet(1,skage.groups_new/d))
    }
    
    exposed_num = rep(Inf, nrow(sample_p))
    exposed_mor = rep(Inf, nrow(sample_p))
    
    
    for( j in 1:nrow(sample_p)){
      vac_num = N_vac*sample_p[j,]
      if(max(vac_num-skage.groups_new)<=0){
        Immunity = efficacy*(vac_num)
        suscept_res_meta[1,] = suscept[357,] - Immunity
        I_total_res_meta = I_total_res
        vac_num = N_vac*sample_p[j,]
        for( i in 1:14){
          exposed_res_meta[i,] = sapply(1:16, function(x) {suscept_res_meta[i,x]*C[x]*(contact_matrix_res[x,]%*%I_total_res_meta[i,])})
          
          suscept_res_meta[(i+1),] = suscept_res_meta[i,] - exposed_res_meta[i,]
          
          if((i+tran_mean+incu_mean)<=15){
            if((i+q_mean + incu_mean-1)<=15){
              I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + round(0.84*exposed_res_meta[i,])
            }else{
              I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.84*exposed_res_meta[i,])
            }
            if((i+C_mean + incu_mean-1)<=15){
              I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + round(0.5 * 0.16*exposed_res_meta[i,])
            }else{
              I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.5 * 0.16*exposed_res_meta[i,])
            }
          }
          
        }
        
        exposed_num[j] = sum(round(colSums(exposed_res_meta)))
        exposed_mor[j] = sum(round(colSums(exposed_res_meta)*mortality))
      }
    }
    
    res_num = rep(0,16)
    res_mor = rep(0,16)
    
    m=16 #Because sum of p is 1, we remove one element.  
    
    
    #The number of vaccines for each group should not exceed the number of populations in that group.
    realsample = sample_p[which(exposed_num!=Inf),-m] 
    realexposed_num = exposed_num[which(exposed_num!=Inf)]/10000
    realexposed_mor = exposed_mor[which(exposed_mor!=Inf)]/1000
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
    bvec = sapply(c(1:16)[-m], function(x) min((skage.groups_new/N_vac)[x],1))
    bvec = c(1-(skage.groups_new/N_vac)[m],-1,-bvec,rep(0,15))
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
    if(max(vac_num_res-skage.groups_new)<=0){
      Immunity = efficacy*(vac_num_res)
      suscept_res_meta[1,] = suscept[357,] - Immunity
      
      I_total_res_meta = I_total_res
      
      for( i in 1:14){
        exposed_res_meta[i,] = sapply(1:16, function(x) {suscept_res_meta[i,x]*C[x]*(contact_matrix_res[x,]%*%I_total_res_meta[i,])})
        suscept_res_meta[(i+1),] = suscept_res_meta[i,] - exposed_res_meta[i,]
        if((i+tran_mean+incu_mean)<=15){
          if((i+q_mean + incu_mean-1)<=15){
            I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + round(0.84*exposed_res_meta[i,])
          }else{
            I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.84*exposed_res_meta[i,])
          }
          if((i+C_mean + incu_mean-1)<=15){
            I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + round(0.5 * 0.16*exposed_res_meta[i,])
          }else{
            I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.5 * 0.16*exposed_res_meta[i,])
          }
        }
      }
      final_num = sum(round(colSums(exposed_res_meta)))
    }
    print(final_num)
    
    if(N_vac == N_vac_list[1]){
      num_infection_1[iterate,1:16] = vac_num_res
      num_infection_1[iterate,17] = final_num
      print(num_infection_1[iterate,])
    }else if(N_vac == N_vac_list[2]){
      num_infection_2[iterate,1:16] = vac_num_res
      num_infection_2[iterate,17] = final_num
      print(num_infection_2[iterate,])
    }else if(N_vac == N_vac_list[3]){
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
    bvec = sapply(c(1:16)[-m], function(x) min((skage.groups_new/N_vac)[x],1))
    bvec = c(1-(skage.groups_new/N_vac)[m],-1,-bvec,rep(0,15))
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
    
    if(max(res_vac_mor-skage.groups_new)<1){
      Immunity = efficacy*(res_vac_mor)
      suscept_res_meta = matrix(0, nrow=15, ncol=16)
      suscept_res_meta[1,] = suscept[357,] - Immunity
      
      I_total_res_meta = I_total_res
      
      exposed_res_meta = matrix(0,nrow=14,ncol=16)
      for( i in 1:14){
        exposed_res_meta[i,] = round(sapply(1:16, function(x) {suscept_res_meta[i,x]*C[x]*(contact_matrix_res[x,]%*%I_total_res_meta[i,])}))
        suscept_res_meta[(i+1),] = suscept_res_meta[i,] - exposed_res_meta[i,]
        if((i+tran_mean+incu_mean)<=15){
          if((i+q_mean + incu_mean-1)<=15){
            I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + round(0.84*exposed_res_meta[i,])
          }else{
            I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.84*exposed_res_meta[i,])
          }
          if((i+C_mean + incu_mean-1)<=15){
            I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + round(0.5*0.16*exposed_res_meta[i,])
          }else{
            I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.5*0.16*exposed_res_meta[i,])
          }
        }
      }
      final_mor = sum(round(colSums(exposed_res_meta)*mortality))
    }
    print(round(colSums(exposed_res_meta)))
    print(final_mor)
    
    if(N_vac == N_vac_list[1]){
      num_death_1[iterate,1:16] = res_vac_mor
      num_death_1[iterate,17] = final_mor
      print(num_death_1[iterate,])
    }else if(N_vac == N_vac_list[2]){
      num_death_2[iterate,1:16] = res_vac_mor
      num_death_2[iterate,17] = final_mor
      print(num_death_2[iterate,])
    }else if(N_vac == N_vac_list[3]){
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


suscept_res_meta = matrix(0, nrow=15, ncol=16)
suscept_res_meta[1,] = suscept[357,] - efficacy*suscept[357,]

I_total_res_meta = I_total_res

exposed_res_meta = matrix(0,nrow=14,ncol=16)
for( i in 1:14){
  exposed_res_meta[i,] = round(sapply(1:16, function(x) {suscept_res_meta[i,x]*(pi_res[x]/(skage.groups_new[x]))*(contact_matrix_res[x,]%*%I_total_res_meta[i,])}))
  suscept_res_meta[(i+1),] = suscept_res_meta[i,] - exposed_res_meta[i,]
  if((i+tran_mean+incu_mean)<=15){
    if((i+q_mean + incu_mean-1)<=15){
      I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+q_mean + incu_mean-1),] + round(0.84*exposed_res_meta[i,])
    }else{
      I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.84*exposed_res_meta[i,])
    }
    if((i+C_mean + incu_mean-1)<=15){
      I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] = I_total_res_meta[(i+tran_mean+incu_mean):(i+C_mean + tran_mean + incu_mean-1),] + round(0.5*0.16*exposed_res_meta[i,])
    }else{
      I_total_res_meta[(i+tran_mean+incu_mean):15,] = I_total_res_meta[(i+tran_mean+incu_mean):15,] + round(0.5*0.16*exposed_res_meta[i,])
    }
  }
}
mor_5 = sum(round(colSums(exposed_res_meta)*mortality))
inf_5 = sum(round(colSums(exposed_res_meta)))


num_infection_5 = c(skage.groups_new,inf_5)
num_death_5 = c(skage.groups_new,mor_5)

num_infection = rbind(num_infection_1[which.min(num_infection_1[,17]),],num_infection_2[which.min(num_infection_2[,17]),],num_infection_3[which.min(num_infection_3[,17]),],num_infection_4[which.min(num_infection_4[,17]),],num_infection_5)
num_death = rbind(num_death_1[which.min(num_death_1[,17]),],num_death_2[which.min(num_death_2[,17]),],num_death_3[which.min(num_death_3[,17]),],num_death_4[which.min(num_death_4[,17]),],num_death_5)



for(j in 1:5){
  num_infection[j,1:16] = num_infection[j,1:16]/skage.groups_new
  num_death[j,1:16] = num_death[j,1:16]/skage.groups_new
}
# 
# write.csv(t(num_infection),'reduce_infection_KOR_asym0.04_210317.csv')
# write.csv(t(num_death),'reduce_mortality_KOR_asym0.04_210317.csv')



library("ggplot2")
library("tidyr")
library(data.table)
library("lattice")
library(RColorBrewer)

###############################################################
#################### HEAT MAP #################################
###############################################################
data <- t(num_infection[,1:16])
colnames(data) = as.character(c('20%','40%','60%','80%','100%'))
rownames(data) = c('[0,5)','[5,10)','[10,15)','[15,20)','[20,25)','[25,30)',
                   '[30,35)','[35,40)','[40,45)','[45,50)','[50,55)','[55,60)','[60,65)','[65,70)','[70,75)','[75,+)')

data = as.data.frame(data)

setDT(data, keep.rownames = T)[]

colnames(data) = c('age','20%','40%','60%','80%','100%')

data = pivot_longer(data , cols = c(2:6) , names_to = 'Vaccine supply',values_to = 'num')
data$`Vaccine supply` <- factor(data$`Vaccine supply`, names(sort(table(data$`Vaccine supply`)))[c(2:5,1)])

data$age <- factor(data$age, names(sort(table(data$age)))[c(1,10,2:9,11:16)])

mine.heatmap <- ggplot(data , mapping = aes(x = age,
                                            y = `Vaccine supply`,
                                            fill = num)) +
  geom_tile() +
  ggtitle("To reduce cases") +
  xlab(label = "age") +
  scale_fill_gradient2(name = "proportion",
                       low = "#f7fbff",
                       mid = "#6baed6",
                       high = "#08306b",
                       midpoint=0.5)+
  theme(plot.title = element_text(hjust = 0.5))

mine.heatmap
############################################################
############################################################

data <- t(num_death[,1:16])
colnames(data) = as.character(c('20%','40%','60%','80%','100%'))
rownames(data) = c('[0,5)','[5,10)','[10,15)','[15,20)','[20,25)','[25,30)',
                   '[30,35)','[35,40)','[40,45)','[45,50)','[50,55)','[55,60)','[60,65)','[65,70)','[70,75)','[75,+)')

data = as.data.frame(data)

data = setDT(data, keep.rownames = T)[]

colnames(data) = c('age','20%','40%','60%','80%','100%')

data = pivot_longer(data , cols = c(2:6) , names_to = 'Vaccine supply',values_to = 'num')
data$`Vaccine supply` <- factor(data$`Vaccine supply`, names(sort(table(data$`Vaccine supply`)))[c(2:5,1)])

data$age <- factor(data$age, names(sort(table(data$age)))[c(1,10,2:9,11:16)])


mine.heatmap <- ggplot(data , mapping = aes(x = age,
                                            y =`Vaccine supply`,
                                            fill = num)) +
  geom_tile() +
  ggtitle("To reduce deaths") +
  theme( plot.title = element_text(hjust = 0.5 ))+
  xlab(label = "age") +
  scale_fill_gradient2(name = "proportion",
                       low = "#fff7ec",
                       mid = "#fdd49e",
                       high = "#b30000",
                       midpoint = 0.5)
mine.heatmap

