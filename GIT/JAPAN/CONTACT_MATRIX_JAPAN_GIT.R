##########################################################################################################
########################### - Contact Matrix of Japan ###################################################
##########################################################################################################
library(devtools)
library(optiSolve)
library(quadprog)
library(polynom)
library(fluEvidenceSynthesis)
library(logitnorm)

rm(list=ls())
gc()




#contact_data
ibuka = read.csv('contact.csv')
ibuka = ibuka[,2:ncol(ibuka)]


#age_data
japan_age_group = as.vector(read.csv('japan_population.csv',header=T)$x)
age_limits = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)

mat_week = matrix(0,nrow=16,ncol=16)
mat_holi = matrix(0,nrow=16,ncol=16)

ibuka_new = matrix(0,nrow=nrow(ibuka),ncol=18)

for(i in 1:18){
  if(i<=6){
    ibuka_new[,i] = ibuka[,i]
  }else{
    temp = floor( (i+7)/2 )
    ibuka_new[,i] = 0.5*ibuka[,temp]
  }
}

ibuka_new 
ibuka

##bootstrap
bootibuka_week = matrix(0,nrow=1000000,ncol=ncol(ibuka_new))
bootibuka_holi = matrix(0,nrow=1000000,ncol=ncol(ibuka_new))
ibuka_new_week = ibuka_new[which(ibuka_new[,2]==0),]
ibuka_new_holi = ibuka_new[which(ibuka_new[,2]==1),]


temp_week = sample(1:nrow(ibuka_new_week),1000000,replace=T)
temp_holi = sample(1:nrow(ibuka_new_holi),1000000,replace=T)

for( i in 1:1000000){
  bootibuka_week[i,] = as.matrix(ibuka_new_week[temp_week[i],])
}
for( i in 1:1000000){
  bootibuka_holi[i,] = as.matrix(ibuka_new_holi[temp_holi[i],])
}
bootibuka_week = bootibuka_week[,c(1,3:18)]
bootibuka_holi = bootibuka_holi[,c(1,3:18)]

bootibuka_week = as.data.frame(bootibuka_week)
bootibuka_holi = as.data.frame(bootibuka_holi)
colnames(bootibuka_week) = c('age','0_5','5_10','10_15','15_20','20_25','25_30','30_35','35_40','40_45','45_50','50_55','55_60','60_65','65_70','70_75','75+')
colnames(bootibuka_holi) = c('age','0_5','5_10','10_15','15_20','20_25','25_30','30_35','35_40','40_45','45_50','50_55','55_60','60_65','65_70','70_75','75+')


for(i in 1:16){
  if (i==1){
    temp1 = bootibuka_week %>% filter(0<=age & age<age_limits[1])
    temp2 = bootibuka_holi %>% filter(0<=age & age<age_limits[1])
  }
  else if(i==16){
    temp1 = bootibuka_week %>% filter(age>=age_limits[15])
    temp2 = bootibuka_holi %>% filter(age>=age_limits[15])
  }
  else{
    temp1 = bootibuka_week %>% filter(age>=age_limits[i-1] & age<age_limits[i])
    temp2 = bootibuka_holi %>% filter(age>=age_limits[i-1] & age<age_limits[i])
  }
  for(j in 1:16){
    mat_week[i,j] = sum(temp1[,j+1])/nrow(temp1)
    mat_holi[i,j] = sum(temp2[,j+1])/nrow(temp2)
  }
}
mat_week[is.na(mat_week)]=0

temp_20 = (mat_week[6,] + mat_week[5,])
mat_week[5,] = temp_20
mat_week[6,] = temp_20

temp_20 = (mat_holi[6,] + mat_holi[5,])
mat_holi[5,] = temp_20
mat_holi[6,] = temp_20

#####contactmatrix
mat1_week = mat_week
mat1_holi = mat_holi
for(i in 1:16){
  for(j in 1:16){
    mat1_week[i,j] = (1/2) *(sum(japan_age_group)/japan_age_group[i])*(mat_week[i,j]*(japan_age_group[i]/sum(japan_age_group)) + mat_week[j,i]*(japan_age_group[j]/sum(japan_age_group)))
    mat1_holi[i,j] = (1/2) *(sum(japan_age_group)/japan_age_group[i])*(mat_holi[i,j]*(japan_age_group[i]/sum(japan_age_group)) + mat_holi[j,i]*(japan_age_group[j]/sum(japan_age_group)))
  }
}

mat1_week[is.na(mat1_week)]=0

japan_contact = (5/7)*(mat1_week) + (2/7)*(mat1_holi)


japan_contact = japan_contact - japan_contact*(contact_others$KOR/contact_all$KOR)*(0.12) - japan_contact*(contact_work$KOR/contact_all$KOR)*(0.1)


write.csv(japan_contact,'contact_japan_res.csv')

