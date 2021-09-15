# Maintained and Edited by Henry Claussen
# henry.claussen@emory.edu

# Implementation of JdCov screening on RNAseq, Methylation
# and CNV datasets from the same cohort obtained from 
# TCGA database.  Implementation of the JdCov screening is 
# done in parallel to decrease run time.

# This code assumes that the user has uploaded the RNAseq, 
# Methylation, and CNV datasets from their own directories.

# This code assumes that the user has downloaded the desired 
# data from a TCGA database.  Once the data has been downloaded, 
# it must be cleaned of any columns that contain NA values 
# and matched subject-wise across all three datasets.

# The datasets must be oriented so that the subjects are rows, 
# columns are the corresponding covariates for each subject.

# The cleaned and matched RNAseq, Methylation, and CNV data 
# have been loaded as rna, meth, and cnv respectively.

# The JdCov scripts are downloaded from their respective 
# Github repositories.
source('jdcov_basic_functions.r')
source('jdcov_test.r')

# User-created function to implement JdCov screening
# The RNA response matrix is labeled as Y, and the single column
# of the covariate matrix to be screened is labeled as x.

# This function returns a single p-value
screen_func <- function(Y,x){
	p.value <- jdcov.test(list(Y,as.matrix(x)),cc=0)$p.value
	return(p.value)
}

# Libraries necessary for funning the parallel implementation
# of JdCov screening
library(parallel)
library(foreach)
library(doParallel)

# Number of cores to run in parallel
registerDoParallel(6)

# X in the following code is the desired covariate matrix to be screened. 
# This can be the methylation or cnv matrices seperately, or combined
# column-wise into a single covariate matrix.
X <- cbind(meth,cnv)

# Number of bootstrap samples taken from the total subjects.  All subjects can
# be included but this will significantly increase run-time for large sample
# sizes.
n <- 50

start <- 1
finish <- ncol(X)

time <- Sys.time()
out <- foreach(j=start:finish,combine=cbind) %dopar% {
	screen_func(rna[sample(nrow(rna),n),],X[sample(nrow(X),n),j])
}

# Saves the output of p-values from JdCov screening
save(out,file='JdCov_screened_X.Rdata')




