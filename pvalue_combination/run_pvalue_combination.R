# Authors: Grace Tiao, Ziao Lin
# Integration of p-values

# Load libraries and arguments
args=commandArgs(TRUE)

# Source the ebm.R package from the following github repo: https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM

source("pvalue_combination/ebm.R")

source("pvalue_combination/combine_p_values.library.R")
library("reshape")


# Load data
tissue=args[1]
target=args[2]
dir=args[3]
sif_filepath=args[4]

# SIF requirements: method name, filepath, name of gene column, name of p-value column, tissue, target
SIF=read.delim(sif_filepath, header=TRUE, as.is=TRUE)

# Subset to target and tissue
data=SIF[SIF$target==target & SIF$tissue==tissue,]
list_of_tables=lapply(data$filepath, read.delim, na.strings=c("", "NA", " "), check.names=FALSE, header=TRUE, as.is=TRUE)

#--------------
# Clean data
for (i in 1:length(list_of_tables)) {
	list_of_tables[[i]]=list_of_tables[[i]][,c(data$gene_col[i], data$p_col[i])]
	names(list_of_tables[[i]])=c("ID", data$method[i])
}
mapply(assign, paste(data$method, "observed", sep="."), list_of_tables, pos=1)

# Merge all methods into master dataframe
observed=merge_tables(paste(data$method, "observed", sep="."))

# Remove duplicate rows
observed=remove_dups(observed)

# For each data type, set p-values greater than 1 to NA
observed$ID=as.factor(observed$ID)
observed[observed>1]<-NA

#----------------
# Plot QQs of reported p-values for each method
sort_p=function(pvals){return(as.numeric(pvals[order(pvals)]))}
max_colors=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','black','#b15928', 'lightblue', 'grey', '#8c510a', '#35978f')

df=observed
methods=names(df)[names(df)!="ID"]
df=apply(df[,methods], 2, sort_p)
n=length(df[,1][!is.na(df[,1])])
qs=1:n/(n+1)
obs=df[,1][!is.na(df[,1])]

if (length(methods)==1){
	pdf(file.path(dir, paste(tissue, target, "reported_p_values.QQ.pdf", sep=".")))
	plot(-log10(qs), -log10(obs), cex=0.5, xlab="Expected p-value (-log10)", main=paste(tissue, target, "Reported P-values", sep=" "), ylab="Observed p-value (-log10)", col="grey", pch=19)
	abline(0,1)
	legend("topleft", pch=19, col="grey", bty="n", legend=colnames(df), cex=.7)
	dev.off()
} else {
	pdf(file.path(dir, paste(tissue, target, "reported_p_values.QQ.pdf", sep=".")))
	plot(-log10(qs), -log10(obs), cex=0.5, xlab="Expected p-value (-log10)", main=paste(tissue, target, "Reported P-values", sep=" "), ylab="Observed p-value (-log10)", col="grey", pch=19)
	abline(0,1)
	for (i in 2:dim(df)[2]){
		n = length(df[,i][!is.na(df[,i])])
		if (n == 0){
		print (paste('error method', colnames(df)[i], sep = ''))
		}
		qs = 1:n/(n+1)
		obs = df[,i][!is.na(df[,i])]
		points(-log10(qs), -log10(obs), cex=0.5, pch=19, col=max_colors[i])
	}
	legend("topleft", pch=19, bty="n", col=c("grey", max_colors[2:dim(df)[2]]), legend=colnames(df), cex=.7)
	dev.off()
}

#--------------
# Assess covariance

# Cap p-values to 1e-16
observed[observed<1e-16]<-1e-16

# Calculate covariance matrix on observed data (raw p-values are transformed to -2logp)
m=t(data.matrix(observed[,2:dim(observed)[2]]))
observed_cov=calculateCovariances_NA(m)

# Remove methods where all the entries are NA
cov_matrix=observed_cov
cov_matrix=NA_rm(cov_matrix)
#Remove all outer edges with NAs until there are no more NAs on edges; then remove NAs within the matrix with na.omit
n=dim(cov_matrix)[1]
na_sum=sum(c(is.na(cov_matrix[,n]), is.na(cov_matrix[n,])))
while (na_sum > 0) {
	cov_matrix=outer_edge_NA_rm(cov_matrix)
	n=dim(cov_matrix)[1]
	na_sum=sum(c(is.na(cov_matrix[,n]), is.na(cov_matrix[n,])))
}
cov_matrix=NA_rm_interior(cov_matrix)
methods=colnames(cov_matrix)

# Combine p-values
observed_tmp=data.matrix(observed[,methods])
observed_brown=data.frame(matrix(unlist(apply(observed_tmp, 1, empiricalBrownsMethod_NA, covar_matrix=cov_matrix, extra_info=T)), ncol=4, byrow=T))
names(observed_brown)=c("Brown_observed", "Fisher_observed", "Brown_Scale_C_observed", "Brown_DF_observed")

evaluate=cbind(observed[,c("ID", methods)], observed_brown)
write.table(evaluate, file.path(dir, paste(target, tissue, "combined_p_values.txt", sep=".")), col.names=T, row.names=F, sep="\t", quote=F)


#--------------------
# Summarize significant genes per method (non-integrated)

sig_gene_counts=get_sig(observed[,methods])
median_sig=median(sig_gene_counts)
c=4
upper_thresh=max(median_sig,1)*c

# Remove only inflated methods for now
upper_outliers=names(sig_gene_counts)[sig_gene_counts>upper_thresh]
keep=names(sig_gene_counts)[!names(sig_gene_counts) %in% upper_outliers]


# Prepare report
removal_report=cbind(sig_gene_counts, median_sig, c, upper_thresh, sig_gene_counts>upper_thresh)
colnames(removal_report)[5]=c("removed_for_inflation")
write.table(removal_report, file.path(dir, paste(target,tissue, "automatic_method_removal_report.txt", sep=".")), col.names=T, row.names=T, sep="\t", quote=F)

#----------------------
# Automatically remove inflated methods and re-run p-value integration

# Remove methods from observed methods and then subset covariance matrices to match
observed_subset=observed[,keep]
methods=colnames(cov_matrix)[colnames(cov_matrix) %in% names(observed_subset)]

#if (length(methods)<2) {next} #Skip data type if there aren't enough methods to combine

observed_tmp=observed_subset[,methods] # Match (subset) observed data to the methods in the data type
cov_matrix=cov_matrix[methods,methods] # Trim covariance matrix accordingly

# Remove methods where all the entries are NA
cov_matrix=NA_rm(cov_matrix)
n=dim(cov_matrix)[1]
na_sum=sum(c(is.na(cov_matrix[,n]), is.na(cov_matrix[n,])))
while (na_sum > 0) {
	cov_matrix=outer_edge_NA_rm(cov_matrix)
	n=dim(cov_matrix)[1]
	na_sum=sum(c(is.na(cov_matrix[,n]), is.na(cov_matrix[n,])))
}
cov_matrix=NA_rm_interior(cov_matrix)
methods=colnames(cov_matrix)[colnames(cov_matrix) %in% names(observed_tmp)]


# Combine p-values
observed_tmp=observed_subset[,methods] #Trim (subset) observed data to match
combined=data.frame(matrix(unlist(apply(observed_tmp, 1, empiricalBrownsMethod_NA, covar_matrix=cov_matrix, extra_info=T)), ncol=4, byrow=T))
names(combined)=c("Brown_observed", "Fisher_observed", "Brown_Scale_C_observed", "Brown_DF_observed")

evaluate=cbind(observed[,c("ID", methods)], combined)
write.table(evaluate, file.path(dir, paste(target, tissue, "combined_p_values.automatic_method_removal.txt", sep=".")), col.names=T, row.names=F, sep="\t", quote=F)
