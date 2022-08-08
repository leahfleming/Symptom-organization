library(psych)
##########################################################################################

#remove points where >=2 for <10% of people 
greater_two = c()
for(rows in colnames(raw_dataset)) {
  
  greater_two = c(greater_two,sum(raw_dataset[,rows]>=2, na.rm=T)/nrow(raw_dataset)*100)
}

raw_dataset=raw_dataset%>%mutate_if(is.character, as.numeric) #make numeric

raw_dataset_T10=raw_dataset[,-which(greater_two<10)] #these are the rows with more than 10% reporting 2 or more to be used for MDS

#19=APR ;25=INC ; 29=DST  ;30=CLN removed





Corr_Matrix_comb=cor(na.omit(raw_dataset_T10), method="pearson") #get correlation matrix for dataset



rs_squared_comb=smc(Corr_Matrix_comb, covar=FALSE) #this gives squared multiple correlation --feeds into Browne index formula 
#now use Browne formula found from Drasgow and Dorans 1982 paper in Applied Psychological Measurement
#requires Rp^2, Rp^4 in order to get Rcv^2 (Browne's cross-validity coefficient) --in Minas MDS paper they 
#report Rs^2, Rp^2 and Rcv^2
rp_squared_calculator<-function(rs_squared_comb, N, k){    #N is number subjects (174), k is number items (24), rs_squared_comb is MSC matrix --you input your own values there
  return(1 - (N-1)/(N-k-1)*(1-rs_squared_comb))
}

rp_squared_comb=sapply(rs_squared_comb, function(observation){rp_squared_calculator(observation, 174, 24)})

rp_quad_calculator<-function(rp_squared_comb, N, k){
  return(rp_squared_comb^2-((2*k*(1-rp_squared_comb)^2)/((N-1)*(N-k+1))))
}

rp_quad_comb=sapply(rp_squared_comb, function(obs){rp_quad_calculator(obs, 174, 24)}) 

rp_quad_square=data.frame(rp_quad_comb, rp_squared_comb)


r_cv_calculator<-function(rp_squared_comb, rp_quad, N, k){
  return(((N-k-3)*rp_quad-rp_squared_comb)/((N-2*k-2)*rp_squared_comb+k))
}

r_cv_squared_comb=apply(rp_quad_square, 1, function(x){r_cv_calculator(x[2], x[1], 174,24)}) 
#The output of this shows all the SAPS items and the "Browne index" score for each, the criteria is to exclude items below .15
min(r_cv_squared_comb) #this just shows if the min is below .15 so you know items need to be excluded


Browne_corr_matrix=Corr_Matrix_comb[-which(r_cv_squared_comb<.15), -which(r_cv_squared_comb<.15)] #this is the correlation matrix for only items with
#removed AGR (.08), REP (.05), PCSP (.07), LAT (.06), MST (.04)
Browne_raw_dataset_T10=raw_dataset_T10[,-which(r_cv_squared_comb<.15)] #also remove those items from raw data

#so Browne_raw_dataset_T10 is going to be the dataset you input into the MDS code
