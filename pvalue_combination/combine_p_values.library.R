#Autors: Grace Tiao, Ziao Lin
#R library for p-value integration

#Merge tables from a list of tables
merge_tables=function(list_of_tables){
	master=get(list_of_tables[1])
	for (j in 1:(length(list_of_tables)-1)) {
		master=merge(master, get(list_of_tables[j+1]), all=TRUE, by="ID")
	}
	return(master)
}

#Remove methods with all NAs in covariance matrices
NA_rm=function(cov_matrix){
	na_matrix=is.na(cov_matrix)
	idx=which(colSums(na_matrix)==dim(cov_matrix)[1])	
	remove=rownames(cov_matrix)[idx]
	cov_matrix=subset(cov_matrix, select = !colnames(cov_matrix) %in% remove)
	cov_matrix=cov_matrix[!rownames(cov_matrix) %in% remove,]
	return(cov_matrix)
}


#Remove interior NAs in a square covariance matrix
NA_rm_interior=function(matrix){
	clean_matrix=na.omit(matrix)
	square_clean_matrix=clean_matrix[,colnames(clean_matrix) %in% rownames(clean_matrix)]
	return(square_clean_matrix)
}

#Remove missing values from dataframe
remove_nas=function(df){
	idx=is.na(rowSums(df[,2:dim(df)[2]]))
	return(df[!idx,])
}


#Remove duplicate rows from dataframe
remove_dups=function(df){
	idx=duplicated(df)
	return(df[!idx,])
}


#QQ plot
qq_plot=function(v, filename, dir){
v = v[order(v)]
n = length(!is.na(v))
qs = 1:n/(n+1)

pdf(file.path(dir, paste(filename, "QQ.pdf", sep= ".")))
plot(-log10(qs), -log10(v), cex=0.5, xlab="Expected p-value (-log10)", ylab="Observed p-value (-log10)", main=filename, col="grey", pch=19)
abline(0,1)
dev.off()
}


empiricalBrownsMethod <- function(covar_matrix, p_values, extra_info = FALSE) {
    return(combinePValues(covar_matrix, p_values, extra_info))
}


calculateCovariances_NA=function(data_matrix){
    transformed_data_matrix = apply(data_matrix, MARGIN=1, FUN=transformData_NA)
    covar_matrix = cov(transformed_data_matrix, use="pairwise.complete.obs")
    covar_matrix
  }


empiricalBrownsMethod_NA=function(covar_matrix, p_values, extra_info = FALSE) {
	p_values_clean=p_values[!is.na(p_values)]
	if (length(p_values_clean)<=1) {
		return(rep("NA", times=4))
	} else {
	p_idx=which(p_values %in% p_values_clean)
	covar_matrix_clean=covar_matrix[p_idx,p_idx]
    return(combinePValues(covar_matrix_clean, p_values_clean, extra_info))
    }
}


remove_low_counts=function(df, targets) {
	df=df[df$ID %in% targets,]
	return(df)
}

transformData_NA=function(data_vector) {
	data_vector_clean=data_vector[!is.na(data_vector)]
    dvm = mean(data_vector_clean)
    if (is.na(pop.sd(data_vector_clean))) {
    	return(rep(NA, times=length(data_vector)))
    } else {
    dvsd = pop.sd(data_vector_clean)
    }
    if (dvsd==0) {
    	s = data_vector
    } else {
    s = (data_vector-dvm)/dvsd
    }
    distr = ecdf(s)
    sapply(s, function(a) -2*log(distr(a)))
}


#Function to remove NAs from covariance matrices
#Update function to remove NAs from covariance matrices
outer_edge_NA_rm = function(cov_matrix) {
idx=which(is.na(cov_matrix)) %% dim(cov_matrix)[1]
idx[idx==0]=dim(cov_matrix)[1]
while ((length(idx) > 0) & (1 %in% idx | dim(cov_matrix)[1] %in% idx)) {
      edge_idx=idx[(idx==1 | idx==dim(cov_matrix)[1])]
      remove_idx=as.numeric(names(sort(table(edge_idx), decreasing=T)))  # Sort the indices for NAs by frequency, and remove the most frequent one first
      remove=rownames(cov_matrix)[remove_idx[1]]
      cov_matrix=subset(cov_matrix, select = !colnames(cov_matrix) %in% remove)
      cov_matrix=cov_matrix[!rownames(cov_matrix) %in% remove,]
      idx=which(is.na(cov_matrix)) %% dim(cov_matrix)[1]
      idx[idx==0]=dim(cov_matrix)[1]
      }
return(cov_matrix)
}


BH_adjust=function(pvals){return(p.adjust(pvals, method="BH"))}

#Summarize significant hits per method for automatic method pruning
get_sig=function(pvalues) {
#Remove duplicate levels
pvalues=pvalues[!duplicated(pvalues),]
FDR_vals=apply(pvalues,2,BH_adjust)
FDR_idx=FDR_vals<.1
all_sig_genes=colSums(FDR_idx, na.rm=T)
#df=data.frame(all_sig_genes)
#names(df)=tissue
return(all_sig_genes)
}

