#### PREPARATION DATA ####
setwd("~/espaces/travail/STAGE/1")
data_expr=read.table("expr_gene.csv",header=T,row.names=1,sep=",",dec=",")
data2_data_expr=t(data_expr) #create a new matrix

data_raman=read.table("raman_data.csv",header=T,sep=",",dec=",")

#### OMICS DATA ####
head(data2_data_expr[,1:5]) #show the first 5 columns
summary(data2_data_expr)

boxplot(data2_data_expr)
data3=scale(data2_data_expr)
boxplot(data3)

plot(data3[,1:11])
plot(data3[,2:10])

### PCA OMICS ###
#PCA.res=dudi.pca(data2_data_expr,scale=TRUE, scannf=FALSE, nf=ncol(data2_data_expr))
library(factoextra)
PCA.res=dudi.pca(data_expr,scale=TRUE, scannf=FALSE, nf=10)

PCA.res$eig
inertia.dudi(PCA.res)
fviz_eig(PCA.res)

res.var <- get_pca_var(PCA.res) #results of  variables
head(res.var$coord) #coordinates
head(res.var$contrib) #axes contributions
head(res.var$cos2) # quality representation

#score(PCA.res)

fviz_pca_ind(PCA.res,
             col.ind = "cos2",
             gradient.cols = c("#ff1493","#ffd700","#00fa9a"),
             repel=TRUE
             )
          
# Graphique des variables. Coloration en fonction de la contribution des variables.
# Les variables corrélées positivement sont du même côté du graphique.
# Les variables corrélées négativement sont sur des côtés opposés du graphique.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#fviz_pca_var(PCA.res,
#             col.var = "contrib",
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE
#

# Graphique Biplot des individus et des variables
#fviz_pca_biplot(PCA.res, repel = TRUE,
#                col.var = "#2E9FDF",
#                col.ind = "#696969"
#                )

data_expr_4_gene=t(data_expr)#inverse colonnes lignes

# Calculate row variances
library(matrixStats)
variances <- rowVars(data_expr_4_gene)

# Set a variance threshold (e.g., keep rows (genes) with variance > 0.01)
threshold <- 0.01


# Filter rows (genes) based on variance threshold
data_expr_variant <- data_expr_4_gene[variances > threshold,]

# Check the dimensions of the filtered data
dim(data_expr_4_gene)

dim(data_expr_variant)

#PCA avec 1703/2613
PCA.res=dudi.pca(data_expr_variant,scale=TRUE, scannf=FALSE, nf=10)
PCA.res$eig
inertia.dudi(PCA.res)
fviz_eig(PCA.res)

res.var <- get_pca_var(PCA.res)
res.ind <- get_pca_ind(PCA.res)

fviz_pca_ind(PCA.res,
             col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE
)

fviz_pca_var(PCA.res,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE
)

PCA.res=dudi.pca(t(data_expr_variant),scale=TRUE, scannf=FALSE, nf=10)
PCA.res$eig
inertia.dudi(PCA.res)
fviz_eig(PCA.res)

fviz_pca_ind(PCA.res,
             col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE
)







