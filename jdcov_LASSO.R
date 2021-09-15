# Maintained and Edited by Henry Claussen
# henry.claussen@emory.edu

# This code assumes that the user has downloaded the desired 
# data from a TCGA database.  Once the data has been downloaded, 
# it must be cleaned of any columns that contain NA values 
# and matched subject-wise across all three datasets.

# The datasets must be oriented so that the subjects are rows, 
# columns are the corresponding covariates for each subject.

# The cleaned and matched RNAseq, Methylation, and CNV data 
# have been loaded as rna, meth, and cnv respectively.

# This code assumes that JdCov screening has already 
# been applied to the datasets.  The user has also applied 
# their desired p-value cutoff resulting in dimension 
# reduction of the combined covariate matrix. 

# Load packages necessary for LASSO application 
library(glmnet)

# Desired number of subjects, randomly sampled, 
# from the overall number of subjects
n <- 100

# Number of bootstrap iterations.
start <- 1
stop <- 50

for (i in start:stop){
	# Random sample of subjects from the response matrix
	subjects <- sample(nrow(rna),n,replace=FALSE)

	rna <- rna[subjects,]
	cnv <- cnv[subjects,]
	meth <- meth[subjects,]

	# Application of LASSO Model.  Parameters can be changed.  Lambdas 
	# restricted to 20 to prevent integer overflow.  This parameter should 
	# be tuned based on the data used.
	model1 <- glmnet(as.matrix(cbind(cnv,meth)),as.matrix(rna),
		family='mgaussian', standardize=TRUE,
		nlambda=20,trace.it=1)

	# Save lambdas and associated deviances from glm model
	lambdas <- model1$lambda
	dev <- model1$dev

	# Save theta-hat matrix of estimated coefficients from glm model
	model1.coef <- coef(model1,s=lambdas[length(lambdas)])
	thetahat1 <- as.matrix(as.data.frame(lapply(model1.coef,as.matrix)))

	# Calculation of RSS 
	w <- cbind(rep(1,nrow(cnv)),cbind(cnv,meth))
	rss <- sum((rna - (w %*% thetahat1))^2)/(nrow(rna)*ncol(rna))
	
	# Retain indexes of nonzero coefficients 
	kept <- which(thetahat1[,1] != 0)

	# Save all outputs as a dataframe for each model iteration
	results <- c()
	results$lambdas <- lambdas
	results$deviance <- dev
	results$rss <- rss
	results$sig <- kept
	results$thetahat <- thetahat1

	# Save output into current working directory
	save(results,file=paste0('LASSO_model_',i,'_Rdata'))

}
