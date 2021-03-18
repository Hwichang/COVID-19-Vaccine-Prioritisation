########################################################################################################
################################ - Korea Initial SEIQ - ################################################
########################################################################################################
library(optiSolve)
library(quadprog)
library(polynom)
library(logitnorm)
library(rGammaGamma)
library(stats)
library(STAR)
rm(list=ls())

corona_daily = read.csv('korea_corona_daily.csv')
skage = read.csv('skage.csv',header=T)
skage = skage[1:101,5]
skage = as.numeric(sapply(skage, function(x) gsub(',','',x)))
age_limits = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)
skage.groups = rep(0,17)
skage.groups[1] = sum(skage[1:5])
for( i in 2:16){
  skage.groups[i] = sum(skage[(age_limits[(i-1)]+1):age_limits[i]])
}
skage.groups[17] = sum(skage[81:101])
age_limits_new = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
skage.groups_new = rep(0,16)
skage.groups_new[1] = sum(skage[1:5])
for( i in 2:15){
  skage.groups_new[i] = sum(skage[(age_limits[(i-1)]+1):age_limits[i]])
}
skage.groups_new[16] = sum(skage[76:101])


###############Divide by age group 5 years old############################
temp = matrix(0,nrow = nrow(corona_daily), ncol = 16)
temp[,1] = corona_daily[,2]

for( i in 1:8){
  weight = skage.groups[-2*i+17]/(skage.groups[-2*i+17]+skage.groups[-2*i+18])
  temp[,((2*i-1):(2*i))] = temp[,((2*i-1):(2*i))] + cbind(round((1-weight)*corona_daily[,i+2]),round(weight*corona_daily[,i+2]))
}

corona_new = data.frame('날짜' =corona_daily[,1])
corona_new[,2:17] = temp
colnames(corona_new) = c('날짜', '[75,+)','[70,75)','[65,70)','[60,65)','[55,60)','[50,55)',
                         '[45,50)','[40,45)','[35,40)','[30,35)','[25,30)','[20,25)','[15,20)','[10,15)','[5,10)','[0,5)')


###########Rearrange in order from age 0 and starting from 01/01/2020#######################
corona_new = corona_new[,c(1,17:2)]
corona_new = corona_new[nrow(corona_new):1,]
row.names(corona_new) <- NULL


#######################InCUBATION PERIOD################
Incu_param1 = 4.544
Incu_param2 = 1/0.709
mean(rgamma(1000000, 4.544 , Incu_param2 ))


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
symp_q = read.csv('symp_q.csv')
symp_q = symp_q[which(as.Date(symp_q$Date)>as.Date('2020-03-01')),]
symp_q_dist = as.integer(as.Date(symp_q$Date) - as.Date(symp_q$Symptom))
symp_q_dist = symp_q_dist[which(symp_q_dist>=0 & symp_q_dist<=10)]


#####################Mean of each distributions#####################
q_mean = ceiling(mean(symp_q_dist))
C_mean = ceiling(1/1.7)
incu_mean = ceiling(4.544*0.709)
I_R_mean = 5
tran_mean = ceiling(-4+(tran_param1/tran_param2))


#####################Set seed#################################
set.seed(0814)


####################asymptomatic or symptomatic#############################
total_I = sum(colSums(corona_new[as.Date(corona_new$날짜)>as.Date('2020-10-14'),2:17]))
total_I_all = sum(colSums(corona_new[,2:17]))
index = total_I_all - total_I + 1
asym = rbinom(total_I_all,1,0.16)
#asym = rbinom(total_I_all,1,0.04) #asymomatic rate 0.04
#asym = rbinom(total_I_all,1,0.4) #asymomatic rate 0.4

n_sym = length(which(asym[index:total_I_all]==0))
n_asym = length(which(asym[index:total_I_all]==1))
seiq_matrix = data.frame('age'=NA,'Q_date'=NA,'I_date'=NA,'Y_date'=NA,'E_date'=NA)


#########################initial sampling####################################
asym_C_list = rexp(n_asym,C_param)
asym_L_list = rgamma(n_sym,Incu_param1,Incu_param2) + tran_dist_mu + rgamma(n_sym,tran_param1,tran_param2)
sym_Y_list = rgamma(n_sym,Incu_param1,Incu_param2)
sym_I_list = tran_dist_mu + rgamma(n_sym,tran_param1,tran_param2)
sym_D_list = sample(symp_q_dist,n_sym,replace=TRUE)

#########################initial value or W####################################
k=1
y=1
l=1

sum(corona_new[,2:17])

for( i in 20:nrow(corona_new)){
  if(i >=289){
    for( j in 2:17){
      if(corona_new[i,j]!=0){
        for( h in 1:corona_new[i,j]){
          if( asym[k]==0){ #symptomatic
            seiq_matrix[k,1] = 5*(j-2) #age
            seiq_matrix[k,2] = i #Qurantined date
            seiq_matrix[k,4] = i - sym_D_list[y] #Symptom onset date
            seiq_matrix[k,5] = i - sym_D_list[y] - ceiling(sym_Y_list[y]) #Exposed date
            seiq_matrix[k,3] = max(i - sym_D_list[y] + ceiling(sym_I_list[y]) , seiq_matrix[k,5])  #transmission onset date
            k = k+1
            y = y+1
          }
          else{
            seiq_matrix[k,1] = 5*(j-2) #age
            seiq_matrix[k,2] = i #Qurantined date
            seiq_matrix[k,5] = i - ceiling(asym_C_list[l]) - max(ceiling(asym_L_list[l]),0) #Exposed date
            seiq_matrix[k,3] = i - ceiling(asym_C_list[l]) #transmission onset date
            seiq_matrix[k,4] = NA
            k = k+1
            l = l+1
          }
        }
      }
    }
  }else{
    for( j in 2:17){
      if(corona_new[i,j]!=0){
        for( h in 1:corona_new[i,j]){
          if( asym[k]==0){ #symptomatic
            seiq_matrix[k,1] = 5*(j-2) #age
            seiq_matrix[k,2] = i #Qurantined date
            seiq_matrix[k,4] = i - q_mean #Symptom onset date
            seiq_matrix[k,5] = i - q_mean - incu_mean #Exposed date
            seiq_matrix[k,3] = i - q_mean + tran_mean #transmission onset date
            k = k+1
          }
          else{
            seiq_matrix[k,1] = 5*(j-2) #age
            seiq_matrix[k,2] = i #Qurantined date
            seiq_matrix[k,5] = i - C_mean - tran_mean - incu_mean #Exposed date
            seiq_matrix[k,3] = i - C_mean #transmission onset date
            seiq_matrix[k,4] = NA
            k = k+1
          }
        }
      }
    }
  }
}

