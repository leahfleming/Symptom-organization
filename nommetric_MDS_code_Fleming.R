#Nonmetric MDS code:
#Fleming 10.2020

library(tidyr)
#This is correlation matric (excluding items that do not meet criteria from Exclusion_Criteria for MDS code) which we will use
Browne_corr_matrix

#Now for the MDS
library(MASS)
d <- dist(Browne_corr_matrix) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim (we choose 2 for visualization)
fit # view results


# basic plot of solution 
x <- fit$points[,1] #first dimension
y <- fit$points[,2] #second dimension
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Nonmetric MDS", type="n")
text(x, y, labels = row.names(Browne_corr_matrix), cex=1.1) #this makes it so the points are the 3-letter 
#codes of each item for visualizing




New_fit=fit$points
New_fit_dists=dist(New_fit, diag=TRUE, upper=TRUE)
New_fit_dists
Corr_Matrix_comb
d
r=cor(c(d), c(New_fit_dists))
rsquared <- r * r  #compute
rsquared   #display R^2 


##########################
#Permutation test
library(boot)
library(MASS)
#install.packages('steadyICA')
library(steadyICA)
#install.packages('smacof')
library(smacof)
## permutation test (on the raw data matrix)
#run this code through the MDS you want to test (1st run=Yale unmatched exclusion)
#See paper Mair et al for info on permutation test

set.seed(100)


Browne_fBIRN_UCLA_T10 #raw data


holder=Browne_fBIRN_UCLA_T10
holder_corr=matrix(data=0, nrow=40, ncol=40)

Stress=rep(0, 500)


for (i in 1:500){
  for (symp in 1:ncol(Browne_fBIRN_UCLA_T10))
  {
    holder[,symp]=sample(Browne_fBIRN_UCLA_T10[,symp], size=153, replace=FALSE);
    holder_corr[symp,symp:40]=cor(x=holder[,symp], y=holder[,symp:40], method="pearson", use="complete.obs");
    holder_corr1=t(holder_corr);
    holder_combined=holder_corr+holder_corr1;
    diag(holder_combined)=1;
    distance=dist(holder_combined);
    fit_perm=isoMDS(distance, k=2);
    Stress[i]=fit_perm$stress;
  }
}
stressDf=data.frame(Stress)
fit$stress #stress of original model fof comparison
Stress_plot=ggplot(data=stressDf)+geom_histogram(aes(x=Stress), binwidth=.2)+geom_vline(xintercept=10.62, color="red", linetype="dashed")+geom_hline(yintercept=0, color="gray")+
  theme_bw()+theme(panel.grid.major = element_blank())





#######
#gap statistics to determine number of clusters in data
set.seed(34)

gap_stat <- clusGap(fit_coordinates_Open, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)
print(gap_stat, method = "Tibs2001SEmax")
fviz_gap_stat(gap_stat) #5 clusters is iideal

###################################
#Now look at cluster analysis with Open Access data with own inclusion criteria:
fit_coordinates_Open=data.frame(x_Open, y_Open)

Cluster_Open=kmeans(fit_coordinates_Open, 5, iter.max=10, nstart=2)



###################################

#look at congruence between variables in McGill and Open data
#Procrustes Transformation followed by congruency coefficient calculation with new dimension and comparison dimensions
#from Goodness-of-Fit Assessment in Multidimensional Scaling and Unfolding paper
#Mair, Borg and Rusch, 2016


setwd("/Users/leahfleming/Desktop/Open Access Data/fBIRNPhasell")
Fit_McGill=read.table("fit_dimension_McGill.txt", header=TRUE)


colnames(New_fit)=c("First", "Second")

New_fit=as.matrix(New_fit)
Fit_McGill=as.matrix(Fit_McGill)



Proc=procrustes(X=Fit_McGill, Xstar=New_fit, translation=TRUE, dilation=TRUE)
factor.congruence(Proc$X.new, New_fit, digits=4)
#first factor congruence is .9409, second is .8639
#.9 is considered high degree of factor congruency so this is good! 




