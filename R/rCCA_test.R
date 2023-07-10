### Test rCCA ###
setwd("~/Université/Master/M1/Stage/raman_omics/R")

### Packages needed
library(BiocManager)
BiocManager::install('mixOmics')
library(mixOmics)
library(ade4)
library(factoextra)
library(Cairo)
library(stringr)


data_expr=read.table("/Users/lilly/Documents/Université/Master/M1/Stage/raman_omics/R/Data/expr_gene.csv",header=T,row.names=1,sep=",",dec = ",")
data2_data_expr=t(data_expr)
data2_data_expr
dim(data_expr)# 11 2613
dim(data2_data_expr)# 2613 11

dataraman=read.table("/Users/lilly/Documents/Université/Master/M1/Stage/raman_omics/R/Data/raman_data.csv",header=T,sep=",",dec = ",")
dataraman$Experiment=str_replace(dataraman$Experiment," ", "_")
data4mean=dataraman
col_cat=c('Experiment' ,'Cell.line','Antibiotic.Class')
data4mean$Category=apply( dataraman[ , col_cat ] , 1 , paste , collapse = "#" )

raman=read.table("/Users/lilly/Documents/Université/Master/M1/Stage/raman_omics/R/Data/raman_data_origin.csv",header=T,sep=",",dec = ",")



# Specify data column
dataraman_mean=aggregate(x= data4mean[,6:1116],
                         # Specify group indicator
                         by = list(data4mean$Cell.line),
                         # Specify function (i.e. mean)
                         FUN = mean)

row.names(dataraman_mean)=dataraman_mean[,1]
corres_antibio=unique(dataraman[,c(3,5)])
row.names(corres_antibio)=corres_antibio[,2]
test=dataraman_mean
test[,1]=corres_antibio[test[,1],]$Antibiotic.Class
test
data_raman_2 <- t(test)

boxplot(scale(data_expr))
X <- scale(data_expr) #gene_data
Y <- scale(test[,2:1112]) #raman_data
dim(X); dim(Y)

X11()
imgCor(Y, X, sideColors = c("purple", "green"))

### rcc Shrink method###
shrink.rcc <- rcc(Y,X, method = "shrinkage")

# barplot of shrinkage method rCCA canonical correlations
X11()
plot(shrink.rcc, type = "barplot", main = "Shrinkage") 

X11()
plotIndiv(shrink.rcc , group = corres_antibio[,1] ,
          rep.space = "XY-variate", legend = TRUE,
          ind.names = row.names(dataraman_mean),
          legend.title = 'Antibiotic Type', 
          title = 'Celltype : CCA')

plotArrow(shrink.rcc, group = corres_antibio[,1], 
          col.per.group = color.mixo(1:6),
          legend = TRUE,
          title = 'Gene/raman, shrinkage method')

plotVar(shrink.rcc, var.names = c(TRUE,TRUE),
        cex = c(4, 4), cutoff = 0.9,
        title = ' rCCA shrinkage comp 1 - 2')

# deoC, deoA, zipA, deoD semblent correles
correlated_variables <- c("deoC", "deoA", "zipA", "deoD")
extracted_data <- X[, correlated_variables]
extracted_data

# Afficher les noms des variables corrélées
print(correlated_vars)

#########

correlated_variables <- c("deoC", "deoA", "zipA", "deoD")
extracted_data <- X[, correlated_variables]
extracted_data

#########

X11()
network(shrink.rcc, comp = 1,
        cutoff = 0.88)

X11()
cim(shrink.rcc, comp = 1:2, xlab = "genes", ylab = "raman_data")


#####RIDGE#####
# set grid search values for each regularisation parameter
grid1 <- seq(0.001, 0.2, length = 10) 
grid2 <- seq(0.001, 0.2, length = 10)

# optimise the regularisation parameter values
cv.tune.rcc <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, 
                                   validation = "loo") 

