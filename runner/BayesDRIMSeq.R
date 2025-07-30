library(LaplacesDemon)
library(foreach)
library(doMC)
N <- 666

# Notation: 
#	K denotes the number of transcripts for the given gene
#	n1, n2: number of replicates for the 1st and 2nd condition

# arguments: 
#	x: K x n1 matrix of counts for the first condition
#	y: K x n2 matrix of counts for the second condition

# value:
#	posterior probability of DTU
laplaceApproxHyper <- function(x,y,method,samples,lambdaRate){
	n1 <- dim(x)[2]
	n2 <- dim(y)[2]
	x <- t(x)
	y <- t(y)
	J <- ncol(x) + ncol(y)
	K <- J/2
	lambda = 1
	J <- J + 2
	parm.names <- as.parm.names(list(log.gamma=rep(0,J-2),log.beta1 = 0,log.beta2 = 0))
	mon.names <- c("LP",as.parm.names(list(gamma=rep(1,J-2),beta1 = 1,beta2 = 1)))
	mon.names1 <- mon.names[-1][1:K]
	mon.names2 <- mon.names[-1][(K+1):(2*K)]

	PGF <- function(Data){
		u <- rgamma(Data$J - 2,shape = 1,scale = 1)
		u<-c(log(u),log(rexp(2,rate = lambdaRate)))
		return(u)
	}
	nObs1 <- rowSums(x)
	nObs2 <- rowSums(y)
	multinomialCoefficient1 <- sum(lgamma(nObs1+1) - apply(x,1,function(h)sum(lgamma(h+1))))
	multinomialCoefficient2 <- sum(lgamma(nObs2+1) - apply(y,1,function(h)sum(lgamma(h+1))))

	MyData <- list(J=J, PGF=PGF, X=x, Y = y, mon.names=mon.names,parm.names=parm.names)

	Model <- function(parm, Data){
		### Parameters
		gamma1 <- exp(parm[1:K])/sum(exp(parm[1:K]))
		gamma2 <- exp(parm[(K+1):(2*K)])/sum(exp(parm[(K+1):(2*K)]))
		gamma <- c(gamma1,gamma2)
		beta1 <- exp(parm[J-1])
		beta2 <- exp(parm[J])
		### Log(Prior Densities)
		priorCons <- 2*lgamma(K) - 2*sum(lgamma(rep(1,K))) 
		gamma.prior <- priorCons
		beta.prior1 <- dexp(beta1,rate = lambdaRate, log = TRUE)
		beta.prior2 <- dexp(beta2,rate = lambdaRate, log = TRUE)
		### Log-Likelihood
		gammaSum1 <- beta1*sum(gamma1)
		gammaSum2 <- beta2*sum(gamma2)
		gMat1 <- matrix(beta1*gamma1,nrow = n1,ncol=K,byrow=T)
		gMat2 <- matrix(beta2*gamma2,nrow = n2,ncol=K,byrow=T)
		LL1 <- multinomialCoefficient1 + n1*lgamma(gammaSum1) - sum(lgamma(gammaSum1 + nObs1)) + sum(lgamma(gMat1 + Data$X)) - sum(lgamma(gMat1))
		LL2 <- multinomialCoefficient2 + n2*lgamma(gammaSum2) - sum(lgamma(gammaSum2 + nObs2)) + sum(lgamma(gMat2 + Data$Y)) - sum(lgamma(gMat2))
		LL <- LL1 + LL2
		### Log-Posterior
		LP <- LL1 + LL2 + gamma.prior + beta.prior1 + beta.prior2
		Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,gamma,beta1,beta2),yhat=1, parm=parm)
		return(Modelout)
	}

	Initial.Values <- rep(0,MyData$J)
	for (j in 1:n1){
		Initial.Values[1:K] <- Initial.Values[1:K] + (1+x[j,])/sum(1+x[j,])
	}
	Initial.Values[1:K] <- log(Initial.Values[1:K]/n1)
	for (j in 1:n2){
		Initial.Values[(K+1):(2*K)] <- Initial.Values[(K+1):(2*K)] + (1+y[j,])/sum(1+y[j,])
	}
	Initial.Values[(K+1):(2*K)] <- log(Initial.Values[(K+1):(2*K)]/n2)

	Fit <- LaplaceApproximation(Model,Initial.Values,Data=MyData,Iterations=10000,Method=method,CPUs=1,sir=TRUE, Samples = samples)
	#Fit <- LaplaceApproximation(Model,Initial.Values,Data=MyData,Iterations=10000,CPUs=1,sir=FALSE)

#	Fit <- VariationalBayes(Model, Initial.Values, Data=MyData, Covar=NULL,
#	     Iterations=1000, Method="Salimans2", Stop.Tolerance=1e-3, CPUs=1)

	# NULL model
	x <- rbind(x,y)
	J <- ncol(x)
	parm.names <- as.parm.names(list(log.gamma=rep(0,J),log.beta = 0))
	mon.names <- c("LP",as.parm.names(list(gamma=rep(0,J),beta = 1)))
	mon.names0 <- mon.names[-1][1:K]
	J <- J + 1
	PGF <- function(Data){
		u <- rgamma(Data$J - 1,shape = 1,scale = 1)
		u<-c(log(u),log(rexp(1,rate = lambdaRate)))
		return(u)
	}

	nObs1 <- rowSums(x)
	multinomialCoefficient1 <- sum(lgamma(nObs1+1) - apply(x,1,function(h)sum(lgamma(h+1))))

	MyData <- list(J=J, PGF=PGF, X=x, mon.names=mon.names,parm.names=parm.names)

	Model <- function(parm, Data){
		### Parameters
		gamma <- exp(parm[1:K])/sum(exp(parm[1:K]))
		beta <- exp(parm[J])
		### Log(Prior Densities)
		priorCons <- lgamma(K) - sum(lgamma(rep(1,K))) 
		gamma.prior <- priorCons
		beta.prior <- dexp(beta,rate = lambdaRate, log = TRUE)
		### Log-Likelihood
		gammaSum <- beta*sum(gamma)
		gMat <- matrix(beta*gamma,nrow = n1+n2,ncol=K,byrow=T)
		LL <- multinomialCoefficient1 + (n1+n2)*lgamma(gammaSum) - sum(lgamma(gammaSum + nObs1)) + sum(lgamma(gMat + Data$X)) - sum(lgamma(gMat))
		### Log-Posterior
		LP <- LL + gamma.prior + beta.prior
		Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,gamma,beta),yhat=1, parm=parm)
		return(Modelout)
	}

	Initial.Values <- rep(0,MyData$J)
	for (j in 1:(n1+n2)){
		Initial.Values[1:(J-1)] <- Initial.Values[1:(J-1)] + (1+x[j,])/sum(1+x[j,])
	}
	Initial.Values[1:(J-1)] <- log(Initial.Values[1:(J-1)]/(n1+n2))

	FitNULL <- LaplaceApproximation(Model,Initial.Values,Data=MyData,Iterations=10000,Method=method,CPUs=1,sir=TRUE, Samples = samples)
	#FitNULL <- LaplaceApproximation(Model,Initial.Values,Data=MyData,Iterations=10000,CPUs=1,sir=FALSE)

#	print(exp(Fit$Summary1[1:K,'Mode'])/(sum(exp(Fit$Summary1[1:K,'Mode']))))
#	print(exp(Fit$Summary1[(K+1):(2*K),'Mode'])/(sum(exp(Fit$Summary1[(K+1):(2*K),'Mode']))))
#	print(paste("DTU log Marginal likelihood: ",Fit$LML))


#	print(exp(FitNULL$Summary1[,'Mode'])/(sum(exp(FitNULL$Summary1[,'Mode']))))
#	print(paste("NULL model log Marginal likelihood: ",FitNULL$LML))
	parDTUmodel1 <- parDTUmodel2 <- parEEmodel <- rep(0,K)
	if (Fit$Converged == TRUE){
		parDTUmodel1 <- exp(Fit$Summary1[1:K])/sum(exp(Fit$Summary1[1:K]))
		parDTUmodel2 <- exp(Fit$Summary1[(K+1):(2*K)])/sum(exp(Fit$Summary1[(K+1):(2*K)]))
	}
	if (FitNULL$Converged == TRUE){
		parEEmodel <- exp(FitNULL$Summary1[1:K])/sum(exp(FitNULL$Summary1[1:K]))
	}
	#############################  Bayes Factor  ##############################
	Model.list <- list(M1=Fit, M2=FitNULL)
	b1 <- b2 <- b0 <- 0
	if (is.na(Fit$LML) == TRUE){probDTU <- 0}else{
		if (is.na(FitNULL$LML) == TRUE){probDTU <- 0}else{
			#bf <- BayesFactor(Model.list)
			probDTU <- BayesFactor(Model.list)$Posterior.Probability[1]
			b1 <- exp(Fit$Summary1['log.beta1',1])
			b2 <- exp(Fit$Summary2['log.beta2',1])
			b0 <- exp(FitNULL$Summary1['log.beta',1])
			if(is.na(probDTU) == TRUE){probDTU = 0}
		}
	}
	fdr <- probDTU
	warningCriterion <- max( c(max(abs(parEEmodel - parDTUmodel1)), max(abs(parEEmodel - parDTUmodel2))))
	if (is.na(warningCriterion)==TRUE){warningCriterion <- 0}
	#print(paste("warningCriterion = ", warningCriterion))
#	if ((warningCriterion < 0.05)&(probDTU > 0.95)){probDTU <- fdr <- 0; cat("WARNING: not passed diagnosing test of convergence.","\n")}
	if ((warningCriterion < 0.085)&(probDTU > 0.85)){probDTU <- fdr <- 0; cat("WARNING: not passed diagnosing test of convergence.","\n")}

	warningCriterion <- max(abs(parDTUmodel1 - parDTUmodel2))
	if (is.na(warningCriterion)==TRUE){warningCriterion <- 0}
	if ((warningCriterion > 0.20)&(probDTU < 0.8)){probDTU <- fdr <- 1; cat("WARNING: not passed diagnosing test of convergence2.","\n")}


	results <- list(parDTUmodel1 = parDTUmodel1, parDTUmodel2 = parDTUmodel2, parEEmodel = parEEmodel, probDTU = probDTU, fdr = fdr,b1 = b1, b2 = b2, b0 = b0)

	return(results)



}


####################################################################################################################################
####################################################################################################################################
laplaceDM <- function(count_data, gene_data, grouping, min_reads_filter, nCores, lambdaRate, output_prefix = "tmp"){
	myTMPfolder <- paste0("BayesDRIMSeq_", output_prefix)
	if(dir.exists(myTMPfolder)){stop(paste0("Directory: `", myTMPfolder,"` exists, please give different output_prefix."))}
	system(paste0("mkdir ", myTMPfolder))
	setwd(myTMPfolder)
	
	print('step 1')
	if(missing(min_reads_filter)){min_reads_filter <- 20}
	if(missing(count_data)){stop("`count_data` not provided.")}
	if(missing(lambdaRate)){stop("`lambdaRate` not provided.")}
	if(missing(gene_data)){stop("`gene_data` not provided.")}
	if(missing(grouping)){stop("`grouping` not provided.")}
	if(missing(nCores)){nCores <- 1}
	nRows <- dim(count_data)[1]
	nCols <- dim(count_data)[2]
	groupLevels <- levels(grouping)
	nGroups <- length(groupLevels)
	
	if(length(grouping) != nCols){stop("number of columns in `count_data` is not equal to length of `grouping`.")}
	if(nGroups != 2){stop("number of `grouping` levels is not equal to 2.")}
	if(length(gene_data) != nRows){stop("length of `gene_data` is not equal to number of rows in `count_data`.")}
	index1 <- which(grouping == groupLevels[1])	
	index2 <- which(grouping == groupLevels[2])
	n1 <- length(index1)
	n2 <- length(index2)
	#Filter Genes:
	print('Filter genes')
	nTotal <- rowSums(count_data)
	lowExpressedTranscripts <- which(nTotal < min_reads_filter*nCols)
	if(length(lowExpressedTranscripts) > 0){
		count_data <- count_data[-lowExpressedTranscripts,]
	}
	cat(paste("Raw number of genes: ", length(unique(gene_data))),"\n")
	cat(paste("Raw number of transcripts: ", nRows),"\n")
	cat(paste("low expressed transcripts: ", length(lowExpressedTranscripts)),"\n")
	if(length(lowExpressedTranscripts) > 0){
		genes <- gene_data[-lowExpressedTranscripts]
	}
	cat(paste("Filtered gene list: ", length(unique(genes))),"\n")

	nGenes <- length(unique(genes))
	genePartition <- floor(seq(1,nGenes,length = nCores+1))
	genePartition[nCores+1] <- nGenes
	genePartition[1] <- 0
	laplacePosteriorProbs <- vector("list",nCores)
	for(k in 1:nCores){
		print(nCores)
		subsetOfGenes <- (genePartition[k] + 1):genePartition[k+1]
		laplacePosteriorProbs[[k]] <- data.frame(geneNames = unique(genes)[subsetOfGenes],probDE = numeric(length(subsetOfGenes)),FDR1 = numeric(length(subsetOfGenes)), maxDiff = numeric(length(subsetOfGenes)),FDR2 = numeric(length(subsetOfGenes)),b1 = numeric(length(subsetOfGenes)),b2 = numeric(length(subsetOfGenes)),b0 = numeric(length(subsetOfGenes)),switchCondition = numeric(length(subsetOfGenes)))
	}
	cat(paste('Applying Laplace Approximation With Sampling Importance Resampling'),'\n')
	sink('laplace.log')
	registerDoMC(nCores)	
	foreach(k=1:nCores) %dopar% {
		gIter <- 0
		switchViolations <- 0
		subsetOfGenes <- (genePartition[k] + 1):genePartition[k+1]
		genesNameSubset <- unique(genes)[subsetOfGenes]
		for (g in genesNameSubset){
			gIter <- gIter + 1
			trIndex <-  which(genes == g)
			if(length(trIndex)>1){
			cat(paste("gene:",gIter,", nTr =",length(trIndex)),"\n")		
				x  <- count_data[trIndex,index1]
				y  <- count_data[trIndex,index2]
				lapApprox <- laplaceApproxHyper(x,y,method = 'NM', samples = 1000,lambdaRate = lambdaRate) 
#				probDTU <- lapApprox$probDTU
				MaxProbChange <- max(abs(lapApprox$parDTUmodel1 - lapApprox$parDTUmodel2))
#				MaxProbChange <- max(abs( lapApprox$parEEmodel - (lapApprox$parDTUmodel1 + lapApprox$parDTUmodel2)/2))
				laplacePosteriorProbs[[k]][gIter,2:4] <- c(lapApprox$probDTU,lapApprox$fdr,MaxProbChange)
#				if( k == 6){write(paste('iter =', gIter,": ",c(lapApprox$b1,lapApprox$b2,lapApprox$b0)),stderr())}
				laplacePosteriorProbs[[k]][gIter,6:8] <- c(lapApprox$b1,lapApprox$b2,lapApprox$b0)
				dominantTranscripts <- c(order(lapApprox$parDTUmodel1,decreasing=T)[1],order(lapApprox$parDTUmodel2,decreasing=T)[1])
				switchCondition <- diff(lapApprox$parDTUmodel1[dominantTranscripts])*diff(lapApprox$parDTUmodel2[dominantTranscripts])
				if (switchCondition < 0){laplacePosteriorProbs[[k]][gIter,'switchCondition'] <- 1}
	#				print("---------------------------------------------------------> SWITCH CONDITION  <---------");
#					switchViolations <- switchViolations + 1}
			}
			print(gIter)
			if(gIter%%100 == 0){write(paste("Gene: ",gIter," from parallel thread ",k," done.",sep=""),stderr())}

		}
		write.table(laplacePosteriorProbs[[k]], file = paste("laplaceDM.",k,sep=''),quote=FALSE,row.names=FALSE)
	}


	sink()
	cat('OK','\n')
	for(k in 1:nCores){
		laplacePosteriorProbs[[k]] <- read.table(paste("laplaceDM.",k,sep=''),header=TRUE)
		system(paste("rm laplaceDM.",k,sep=''))
	}
	ll<- laplacePosteriorProbs[[1]]
	for (k in 2:nCores){
		ll <- rbind(ll,laplacePosteriorProbs[[k]])
	}
	setwd("../")
	u <- ll$maxDiff
	perm <- order(u,decreasing=TRUE)
	lapPerm <- ll[perm,]
	naIndex <- which(is.na(lapPerm$probDE)==TRUE)
	if(length(naIndex) > 0){
	lapPerm <- lapPerm[-naIndex,]
	}
	K <- dim(lapPerm)[1]
	sigTranscripts <- numeric(K)
	p <- lapPerm
#	orderedP <- p$probDE
	orderedP <- p$FDR1
	K <- dim(p)[1]
	myList <-  1 - orderedP[1]
	k <- 1 
	sigTranscripts[k] <- myList
	for (k in 2:K){
		myList <- myList + 1.0*(1 - orderedP[k])
		sigTranscripts[k] <- myList/k
	}
	lapPerm[,5] <- 1 - sigTranscripts
	ll <- lapPerm

	ll <- ll[,-3]
	ind <- which(ll$switchCondition == 1)
	pTrust <- fdrTrust <- numeric(dim(ll)[1]) 
	pTrust[ind] <- ll$probDE[ind]
	fdrTrust[ind] <- ll$FDR2[ind]
	ll <- ll[,-8]
	ll <- cbind(ll,pTrust,fdrTrust)
	colnames(ll)[4] <- 'FDRraw'
	return(ll)
	#write.table(ll, file = "laplaceDM.HyperPrior.txt",quote=FALSE,row.names=FALSE)
}



