########################################################################################################
################################ - Japan Initial SEIQ - ################################################
########################################################################################################
library(optiSolve)
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

rm(list=ls())
gc()

#age_data
japan_age = read.csv('japan_population.csv',header=F)
japan_age[,2] = as.numeric(sapply(japan_age[,2], function(x) gsub(',','',x)))
age_limits = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)

japan_age = as.numeric(japan_age$V2)*1000
japan_age_group = japan_age[1:16]
japan_age_group[16] = japan_age_group[16]+ sum(japan_age[17:21])

japan_corona = read.csv('japan.csv',header=T)

weight_matrix = matrix(0, nrow= nrow(japan_corona),ncol=9)
for(i in 1:nrow(weight_matrix)){
  weight_matrix[i,] = as.matrix(japan_corona[i,2:10])/sum(japan_corona[i,2:10])
}

for(i in 1:nrow(japan_corona)){
  japan_corona[i,2:10] =  japan_corona[i,2:10] + round(weight_matrix[i,]*japan_corona[i,11])
}

japan_corona = japan_corona[,1:10]

japan_corona_weekly = japan_corona

N = nrow(japan_corona_weekly)

for(i in 1:(N-1)){
  japan_corona_weekly[(N-i),2:ncol(japan_corona_weekly)] = japan_corona[(N-i),2:ncol(japan_corona_weekly)] - japan_corona[(N-i+1),2:ncol(japan_corona_weekly)]
}


############### Dividing age groups by 5 years ############################
japan_daily = read.csv('pcr_positive_daily.csv',header=T)
japan_daily = japan_daily %>%
  separate(colnames(japan_daily), c('Date','Cases'), '\t')

temp = matrix(0,nrow = nrow(japan_corona_weekly), ncol = 16)
temp[,16] = japan_corona_weekly[,10]

for( i in 1:8){
  weight = japan_age_group[2*i-1]/(japan_age_group[2*i-1]+japan_age_group[2*i])
  temp[,((2*i-1):(2*i))] = temp[,((2*i-1):(2*i))] + cbind(round((1-weight)*japan_corona_weekly[,i+1]),round(weight*japan_corona_weekly[,i+1]))
}

japan_corona_new = data.frame('Date' = japan_corona_weekly[,1])
japan_corona_new[,2:17] = temp
colnames(japan_corona_new) = c('Date', '[0,5)','[5,10)','[10,15)','[15,20)','[20,25)','[25,30)',
                               '[30,35)','[35,40)','[40,45)','[45,50)','[50,55)','[55,60)','[60,65)','[65,70)','[70,75)','[75,+)')


expe = data.frame('Date' = rep(as.Date('2020-01-01'),400))
expe[,2:17] = 0
colnames(expe) = c('Date', '[0,5)','[5,10)','[10,15)','[15,20)','[20,25)','[25,30)',
                   '[30,35)','[35,40)','[40,45)','[45,50)','[50,55)','[55,60)','[60,65)','[65,70)','[70,75)','[75,+)')

k=1
for( i in 1:nrow(japan_corona_new)){
  temp1 = japan_corona_new[i,]
  temp_date1 = as.Date(japan_corona_new[i,1])
  temp_date2 = as.Date(japan_corona_new[i+1,1])
  temp3 = japan_daily %>% filter(as.Date(Date) > temp_date2 & as.Date(Date) <= temp_date1)
  temp4 = as.numeric(temp3[,2])/sum(as.numeric(temp3[,2]))
  for(j in 1:(temp_date1-temp_date2)){
    expe[k,1] = as.Date(japan_corona_new[i,1])+1-j
    expe[k,2:17] = round(japan_corona_new[i,2:17]*(temp4[which(temp3$Date==expe[k,1])]))
    k= k+1
  }
}

expe = expe[196:1,]
pre = colSums(japan_corona_new[25:nrow(japan_corona_new),2:17]) # Accumulated before July 16
japan_corona_new = expe

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
symp_q = read.csv('symp_q.csv')
symp_q = symp_q[which(as.Date(symp_q$Date)>as.Date('2020-03-01')),]
symp_q_dist = as.integer(as.Date(symp_q$Date) - as.Date(symp_q$Symptom))
symp_q_dist = symp_q_dist[which(symp_q_dist>=0 & symp_q_dist<=10)]


q_mean = ceiling(mean(symp_q_dist))
C_mean = ceiling(1/1.7)
incu_mean = ceiling(4.544*0.709)
I_R_mean = 5
tran_mean = ceiling(-4+(tran_param1/tran_param2))

#####################Set seed#################################
set.seed(0814)

####################asymptomatic or symptomatic#############################
total_I = sum(colSums(japan_corona_new[as.Date(japan_corona_new$Date)>as.Date('2020-10-14'),2:17]))
total_I_all = sum(colSums(japan_corona_new[,2:17]))
index = total_I_all - total_I + 1
asym = rbinom(total_I_all,1,0.16)
#asym = rbinom(total_I_all,1,0.04)
#asym = rbinom(total_I_all,1,0.4)
n_sym = length(which(asym[index:total_I_all]==0))
n_asym = length(which(asym[index:total_I_all]==1))
seiq_matrix = data.frame('age'=NA,'Q_date'=NA,'I_date'=NA,'Y_date'=NA,'E_date'=NA)


######################### initial sampling ####################################
asym_C_list = rexp(n_asym,C_param)
asym_L_list = rgamma(n_sym,Incu_param1,Incu_param2) + tran_dist_mu + rgamma(n_sym,tran_param1,tran_param2)
sym_Y_list = rgamma(n_sym,Incu_param1,Incu_param2)
sym_I_list = tran_dist_mu + rgamma(n_sym,tran_param1,tran_param2)
sym_D_list = sample(symp_q_dist,n_sym,replace=TRUE)

######################### initial value of W ####################################
k=1
y=1
l=1
sum(japan_corona_new[,2:17])
move = as.Date('2020-07-16')-as.Date('2020-01-01')
for( i in 1:196){
  q_date = japan_corona_new$Date[i]
  if(as.Date(q_date) - as.Date('2020-01-01')>=288){
    for( j in 2:17){
      if(japan_corona_new[i,j]!=0){
        for( h in 1:japan_corona_new[i,j]){
          if( asym[k]==0){ #symptomatic
            seiq_matrix[k,1] = 5*(j-2) #age
            seiq_matrix[k,2] = i+move #Qurantined date
            seiq_matrix[k,4] = i+move - sym_D_list[y] #Symptom onset date
            seiq_matrix[k,5] = i+move - sym_D_list[y] - ceiling(sym_Y_list[y]) #Exposed date
            seiq_matrix[k,3] = max(as.integer(i+move - sym_D_list[y] + ceiling(sym_I_list[y])) , seiq_matrix[k,5])  #transmission onset date
            k = k+1
            y = y+1
          }
          else{
            seiq_matrix[k,1] = 5*(j-2) #age
            seiq_matrix[k,2] = i+move #Qurantined date
            seiq_matrix[k,5] = i+move - ceiling(asym_C_list[l]) - max(ceiling(asym_L_list[l]),0) #Exposed date
            seiq_matrix[k,3] = i+move - ceiling(asym_C_list[l]) #transmission onset date
            seiq_matrix[k,4] = NA
            k = k+1
            l = l+1
          }
        }
      }
    }
  }else{
    for( j in 2:17){
      if(japan_corona_new[i,j]!=0){
        for( h in 1:japan_corona_new[i,j]){
          if( asym[k]==0){ #asymptomatic
            seiq_matrix[k,1] = 5*(j-2) #age
            seiq_matrix[k,2] = i+move #Qurantined date
            seiq_matrix[k,4] = i+move - q_mean #Symptom onset date
            seiq_matrix[k,5] = i+move - q_mean - incu_mean #Exposed date
            seiq_matrix[k,3] = i+move - q_mean + tran_mean #transmission onset date
            k = k+1
          }
          else{
            seiq_matrix[k,1] = 5*(j-2) #age
            seiq_matrix[k,2] = i+move #Qurantined date
            seiq_matrix[k,5] = i+move - C_mean - tran_mean - incu_mean #Exposed date
            seiq_matrix[k,3] = i+move - C_mean #transmission onset date
            seiq_matrix[k,4] = NA
            k = k+1
          }
        }
      }
    }
  }
}

dim(seiq_matrix)

#write.csv(seiq_matrix,'initial_seiq_japan.csv')

