#### PREPARATION DATA #######
setwd("/media/liduverger/LILLY_USB/raman_omics")
data_expr=read.table("expr_gene.csv",header=T,row.names=1,sep=",",dec = ",")
data2_data_expr=t(data_expr)

dataraman=read.table("raman_data.csv",header=T,sep=",",dec = ",")
library(stringr)
dataraman$Experiment=str_replace(dataraman$Experiment," ", "_") 

data4mean=dataraman
col_cat=c('Experiment' ,'Cell.line','Antibiotic.Class')
data4mean$Category=apply( dataraman[ , col_cat ] , 1 , paste , collapse = "#" )

#data4mean=data4mean[,-c(1:5)]
# Specify data column

datarama_mean=aggregate(x= data4mean[,6:1116],     
                        
                        # Specify group indicator
                        by = list(data4mean$Cell.line),      
                        
                        # Specify function (i.e. mean)
                        FUN = mean)

datarama_mean$Cat= str_split_fixed(datarama_mean$Group.1, "#", 3)[,3]
datarama_mean$Cell= str_split_fixed(datarama_mean$Group.1, "#", 3)[,2]

write.csv(datarama_mean,file="prepa_raman_data.csv",row.names = TRUE)
