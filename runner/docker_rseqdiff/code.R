##################functions########################################
log_likelihood_model0=function(A_matrix, N_vector, W_vector, theta)
  #A is matrix; N and W are vectors  
  #theta can be a matrix, then the log_likelihood will be a vector
{
  
  #theta can be a matrix, then the log_likelihood will be a vector
  
  I=length(W_vector)
  
  theta=matrix(as.matrix(theta), ncol=I)
  #this line makes the column of theta to be I, make vector conformative with matrix, so that can do the matrix operation
  
  nrow_theta=dim(theta)[2]
  
  if( any(is.na(A_matrix)) && any(is.na(N_vector)) )
    #in this case, no read is mapped to any isoform of the gene. So the expression is 0.
  {
    loglikelihood=rep(0, times=nrow_theta)
    
    return(loglikelihood)
  }            
  
  else if( any(is.na(A_matrix))==F && any(is.na(N_vector))==F )
    #this is the normal case--neither of A and N is NULL
  {
    epsilon = 1e-10
    
    theta_new = theta + epsilon
    
    loglikelihood = -crossprod(W_vector, t(theta_new)) + crossprod( N_vector, log( crossprod(A_matrix, t(theta_new) ) ) ) 
    
    return(loglikelihood)
  }
  
  else 
    #in this case, only one of the A or N is NULL, which means there's something wrong with either the A or N.
  {
    print("There is something wrong with either the A_matrix or N_vector. One of them is NULL, but the other is not NULL.")
    
    return(NA)
  }
  
}

log_likelihood_model1 = function(A_matrix, N_vector, W_vector, theta, tau)
  #this function calculate the log_likelihood of model 1. N_vector, A_matrix and W_vector are all list
  #theta and tau can be either vector or matrix
{
  
  K_conditions=length(W_vector)
  
  dimention_theta=length(W_vector[[1]])
  
  theta=matrix(as.matrix(theta), ncol=dimention_theta)
  
  tau=matrix(as.matrix(tau), ncol=K_conditions)
  
  nrow_theta=dim(theta)[2]
  
  #catenate A and N to do some data checking
  N_vector_temp = unlist(lapply(N_vector, na.exclude))
  
  A_matrix_temp = do.call(cbind, lapply(A_matrix, na.exclude))
  
  if(is.logical(A_matrix_temp) && is.logical(N_vector_temp))
    #in this case, no reads is mapped to the gene in every condition. Thus theta and ratio are all 0.
  {
    loglikelihood_m1=rep(0, times=nrow_theta)
  }
  
  else if( is.logical(A_matrix_temp)==F && is.logical(N_vector_temp)==F )
    #this is the normal case--neither of A and N is NULL
  {
    epsilon = 1e-10
    
    theta_new = theta + epsilon
    
    tau_new = tau+epsilon
    
    loglikelihood_m1_K_conditions = list() #loglikelihood_m1_K_conditions is a list which will contain K elements. Each element is a vector of length=N_samp
    
    for(k in 1:K_conditions)
    {
      if( any(is.na(A_matrix[[k]]))==F && any(is.na(N_vector[[k]]))==F )
        #this is the normal case
      {
        loglikelihood_m1_K_conditions[[k]]=log(tau[,k])*sum(N_vector[[k]]) + crossprod( N_vector[[k]], log(crossprod(A_matrix[[k]], t(theta_new))) ) - tau[,k]*crossprod(W_vector[[k]], t(theta_new))
      }
      
      else if( (any(is.na(A_matrix[[k]]))==T && any(is.na(N_vector[[k]]))==F) 
               || (any(is.na(A_matrix[[k]]))==F && any(is.na(N_vector[[k]])==T)) )
        #in this case, under the k'th condition, only one of the A or N is NULL, which means there's something wrong with either the A or N.
      {
        
        print(paste("Under condition", k))
        
        print("There is something wrong with either the A_matrix or N_vector. One of them is NULL, but the other is not NULL.")
        
        loglikelihood_m1=NA
        
        return(loglikelihood_m1)
      }
      
      #the third case is under the k'th condition, there is no read mapped to the gene(both A and N are NULL)
      #but we don't need to do anything in this case.
    }
    
    loglikelihood_m1_matrix_form = do.call(rbind, loglikelihood_m1_K_conditions)
    
    loglikelihood_m1=colSums(loglikelihood_m1_matrix_form)
    
    return(loglikelihood_m1)
  }
  
  else 
    #in this case, only one of the A or N is NA(here A and N is the augmented ones), which means there's something wrong with either the A or N.
  {
    print("There is something wrong with either the A_matrix or N_vector. One of them is NULL, but the other is not NULL.")
    
    loglikelihood_m1=NA
    
    return(loglikelihood_m1)
  }
  
  
}

fisher_information_model0 = function(A_matrix, W_vector, theta)
  #A_matrix is a matrix; theta is a vector.
{  
  A_matrix=as.matrix(A_matrix)#in case A will degenerate to a vector
  
  I = length(theta)
  
  J = dim(A_matrix)[2]
  
  epsilon = 1e-10;
  
  theta_new = theta + epsilon
  
  fisher_I = A_matrix %*% diag(as.vector(1 / crossprod(A_matrix, theta_new)), nrow=J, ncol=J) %*% t(A_matrix)
  
  if(det(fisher_I)<1e-6)
  {
    diag(fisher_I) = diag(fisher_I) + (W_vector^2)/(I*sum(W_vector*theta_new))
  }
  
  return(fisher_I)
  
}


fisher_information_model1 = function(A_matrix, W_vector, theta, tau)
  #A_matrix is a list of length K; theta is a vector of length I; tau is a vector of length K
  #The NAs are deleted
{
  K_conditions=length(tau)
  
  dimention_theta = length(theta)
  
  epsilon = 1e-10
  
  theta_new = theta + epsilon
  
  tau_new = tau + epsilon
  
  fisher_I_theta = list()
  
  diag_tau = list()
  
  fisher_I_theta_tau = matrix( rep( 0, times = dimention_theta*(K_conditions-1) ), nrow=dimention_theta )
  
  for (k in 1:K_conditions)
  {
    A_k = as.matrix(A_matrix[[k]])#in case A[[k]] degenerates to a vector or a number
    
    J=dim(A_k)[2]
    
    fisher_I_theta[[k]] = tau_new[k] * (A_k %*% diag(as.vector(1 / crossprod(A_k, theta_new)), nrow=J, ncol=J) %*% t(A_k))
    
    W_k=W_vector[[k]]
    
    diag(fisher_I_theta[[k]]) = diag(fisher_I_theta[[k]]) + (W_k^2)/(dimention_theta*sum(W_k*theta_new))
    #add an value to the diagonal of fisher_I_theta 
    
    if(k<K_conditions)
    {
      #extraterm_diag_tau[[k]] = sum(crossprod(A_matrix[[k]], theta_new)) / tau_new[k]
      diag_tau[[k]] = sum(W_vector[[k]]*theta_new) / tau_new[k] + sum(W_vector[[K_conditions]]*theta_new) / tau_new[K_conditions]
      
      fisher_I_theta_tau[,k] = W_vector[[k]] - W_vector[[K_conditions]]
      #fisher_I_theta_tau[,k] = rowSums(A_matrix[[k]]) - rowSums(A_matrix[[K_conditions]])
    }
  }
  
  fisher_I_theta = Reduce('+', fisher_I_theta) #this is I_theta
  
  if(K_conditions==1)
    #in this case fisher_I_tau is degenarated, only the theta part exist
  {
    fisher_I=fisher_I_theta
  }
  
  else
    #if K>1, then the tau part comes  
  {
    fisher_I_tau = matrix(rep(sum(crossprod(A_matrix[[k]], theta_new)), (K_conditions-1)*(K_conditions-1)), nrow=(K_conditions-1) )
    
    diag(fisher_I_tau) = unlist(diag_tau) #the diaganol of fisher_I_tau will be needed to add one more term
    
    fisher_I_up = cbind(fisher_I_theta, fisher_I_theta_tau)
    
    fisher_I_down = cbind(t(fisher_I_theta_tau), fisher_I_tau)
    
    fisher_I=rbind(fisher_I_up, fisher_I_down)
  }
  
  return(fisher_I)
}

#imp_estimate function, gives the importance sampling mean estimate and 95% CI for each theta 
imp_estimate = function(theta, weight)
  #theta is a vector; weight is a vector
{
  imp_sample=data.frame(theta, weight)
  
  imp_sample=imp_sample[order(imp_sample$theta),]
  
  cum_weight=cumsum(imp_sample$weight)
  
  imp_sample=cbind(imp_sample,cum_weight)
  
  CI_lower=min(imp_sample$theta[ imp_sample$cum_weight >= 0.025 ])
  CI_upper=min(imp_sample$theta[ imp_sample$cum_weight >= 0.975 ])
  imp_mean=crossprod(imp_sample$theta, imp_sample$weight)
  imp_estimation=c(imp_mean, CI_lower, CI_upper)
  return(imp_estimation)
}

calculating_imp_weight=function(imp_samp, imp_density, A_matrix, N_vector, W_vector, loglikelihood_MLE, model_select)
  #this function organize the imp sample by >=0 or <0 and calculate the normalized weight for the importance sample
  #A and N can be eithei a list or a matrix
{
  N_samp=50000
  
  row_min_of_imp_samp = apply(imp_samp, MARGIN=1, FUN=min)
  
  proposal_samp=data.frame(imp_samp, imp_density, row_min_of_imp_samp)
  
  greater_zero_samp=subset(proposal_samp, row_min_of_imp_samp>=0, select=-row_min_of_imp_samp)
  
  if(model_select==0)
  {
    log_imp_samp=log_likelihood_model0(A_matrix, N_vector, W_vector, subset(greater_zero_samp, select=-imp_density))
  }
  
  if(model_select==1)
  {    
    num_tau = length(lapply(A_matrix, na.exclude))
    #this is the number of free tau's here it is >=1
    
    dimention_theta=length(W_vector[[1]])
    
    theta_imp_samp_greaterzero=greater_zero_samp[,1:dimention_theta]
    
    tau_imp_samp_greaterzero=greater_zero_samp[,(dimention_theta+1):(dimention_theta+num_tau)]
    
    log_imp_samp=log_likelihood_model1(A_matrix, N_vector, W_vector, theta_imp_samp_greaterzero, tau_imp_samp_greaterzero)
  }
  
  likelihood_greater_zero = exp(as.vector(log_imp_samp) - loglikelihood_MLE)
  
  smaller_zero_samp=subset(proposal_samp, row_min_of_imp_samp<0, select=-row_min_of_imp_samp)
  
  likelihood_smaller_zero = rep( 0, times=(N_samp-length(likelihood_greater_zero)) )
  
  proposal_samp_new=rbind(greater_zero_samp, smaller_zero_samp)
  
  likelihood_sample = c(likelihood_greater_zero, likelihood_smaller_zero)
  
  weight = likelihood_sample / proposal_samp_new$imp_density
  
  weight_norm = weight / sum(weight)
  
  return ( list(proposal_samp_new=proposal_samp_new, weight_norm=weight_norm) )
}


#this function uses importance sampling to calculate 95% CI for theta
Calculating_CI_for_theta = function(A_matrix, N_vector, W_vector, theta_MLE, loglikelihood_MLE)
  #theta_MLE is a vector of length I; loglikelihood_MLE is a number; A_matrix is a matrix column-binded of K conditions
  #N_vector and W_vector are all vectors
{
  
  dimention_theta=length(theta_MLE)
  
  if( any(is.na(A_matrix)) && any(is.na(N_vector)) )
    #in this case, no read is mapped to any isoform of the gene under this condition. So the expression is 0.
    #then the CI is 0
  {
    theta_result = matrix( rep(0, times=3*dimention_theta), ncol=3, nrow=dimention_theta )
    
    gene_result = c(0,0,0)
    
    colnames(theta_result)=c('imp_mean', 'CI_lower','CI_upper')
    
    names(gene_result)=c('imp_mean', 'CI_lower','CI_upper')
    
    CI_result=list(theta_CI=theta_result, gene_CI=gene_result)
    
    return(CI_result)
  }
  
  else if( any(is.na(A_matrix))==F && any(is.na(N_vector))==F )
    #this is the normal case--neither of A and N is NULL
  {
    
    fisher_I_matrix=fisher_information_model0(A_matrix, W_vector, theta_MLE)
    
    eigen_fisher_I = eigen(fisher_I_matrix, symmetric=T, only.values=T)$values
    
    r_condition = rcond(fisher_I_matrix)
    #the reciprocal of the condition number
    
    if(any(eigen_fisher_I<=0) || r_condition < 1e-9) 
      #fisher_I is not positive definite or the reciprocal of the condition number is too small
    {  
      CI_result=list(theta_CI=NA, gene_CI=NA)
      
      return(CI_result)
    }
    
    C_matrix = solve(fisher_I_matrix)
    
    C_matrix=4*C_matrix
    
    C_matrix=(C_matrix + t(C_matrix))/2
    
    N_samp=50000
    
    #Using the Multivariate normal distribution
    
    #theta_proposal_samp=rmvnorm(N_samp, mean=theta_MLE, sigma=C_matrix) ##note: need to load mvtnorm library
    
    #theta_proposal_density=dmvnorm(theta_proposal_samp, mean=theta_MLE, sigma=C_matrix)
    
    
    #The multivariate t distribution is actually not quite different as the multivariate normal distribution
    
    theta_proposal_samp=rmvt(N_samp, delta=theta_MLE, sigma=C_matrix, df=5) ##note: need to load mvtnorm library
    
    theta_proposal_density=dmvt(theta_proposal_samp, delta=theta_MLE, sigma=C_matrix, df=5, log=FALSE)
    
    theta_samp_weight=calculating_imp_weight(theta_proposal_samp, theta_proposal_density, A_matrix, N_vector, W_vector, loglikelihood_MLE, model_select=0)
    
    theta_proposal_samp_new=theta_samp_weight$proposal_samp_new
    
    weight_norm=theta_samp_weight$weight_norm
    
    theta_result = matrix( rep(0, times=3*dimention_theta), ncol=3, nrow=dimention_theta )
    
    for(i in 1:dimention_theta)
    {
      #calculate imp_point estimate and imp 95% CI through imp_estimate function
      theta_result[i,]=imp_estimate(theta_proposal_samp_new[,i], weight_norm)
    }
    
    gene_imp_sample = rowSums(subset(theta_proposal_samp_new, select=-imp_density))
    
    gene_imp_estimation = imp_estimate(gene_imp_sample, weight_norm)
    
    colnames(theta_result)=c('imp_mean', 'CI_lower','CI_upper')
    
    names(gene_imp_estimation)=c('imp_mean', 'CI_lower','CI_upper')
    
    CI_result=list(theta_CI=theta_result, gene_CI=gene_imp_estimation)
    
    return (CI_result)
  }
  
  else 
    #in this case, only one of the A or N is NULL, which means there's something wrong with either the A or N.
  {
    print("There is something wrong with either the A_matrix or N_vector. One of them is NULL, but the other is not NULL.")
    
    CI_result=list(theta_CI=NA, gene_CI=NA)
    
    return (CI_result)
  }
  
}

Calculating_CI_for_model_0 = function(A_matrix_m0, N_vector_m0, W_vector_m0, theta_MLE_m0, loglikelihood_MLE_m0)
  #theta_MLE_m0 is a vector of length I; loglikelihood_MLE_m0 is a number; A_matrix_m0 is a matrix;
  #N_vector_m0 and W_vector_m0 are vectors
{
  
  Model_0_CI=Calculating_CI_for_theta(A_matrix_m0, N_vector_m0, W_vector_m0, theta_MLE_m0, loglikelihood_MLE_m0)
  
  return(Model_0_CI)
}

Calculating_CI_for_model_1 = function(A_matrix_m1, N_vector_m1, W_vector_m1, theta_MLE_m1, tau_MLE_m1, log_likelihood_MLE_m1)
{  
  K_conditions = length(W_vector_m1) #K
  
  dimention_theta=length(theta_MLE_m1)#i.e. I
  
  #data checking
  N_vector_temp = unlist(lapply(N_vector_m1, na.exclude))
  
  A_matrix_temp = do.call(cbind, lapply(A_matrix_m1, na.exclude))
  
  Model_1_CI=list()
  
  if(is.logical(A_matrix_temp) && is.logical(N_vector_temp))
    #in this case, no reads is mapped to the gene in every condition. Thus the CI for both are all 0.
  {
    theta_result = matrix( rep(0, times=3*dimention_theta), ncol=3, nrow=dimention_theta )
    
    gene_result = c(0,0,0)
    
    colnames(theta_result)=c('imp_mean', 'CI_lower','CI_upper')
    
    names(gene_result)=c('imp_mean', 'CI_lower','CI_upper')
    
    for(k in 1:K_conditions)
    {
      Model_1_CI[[k]]=list(theta_CI=theta_result, gene_CI=gene_result)
    }
    
    return(Model_1_CI)
  }
  
  else if(any(is.na(theta_MLE_m1)))
    #if theta_MLE_m1 is na, then there is an error
  {
    Model_1_CI=list(theta_CI=NA, gene_CI=NA)
    
    return(Model_1_CI)
  }
  
  else
  {
    position_of_zero_counts=NULL
    #this is vector, which will store the number of conditions that have zero read counts
    
    for(k in 1:K_conditions)
      #check and record the A's that are nulls
    {
      
      if(any(is.na(A_matrix_m1[[k]])))
      {
        position_of_zero_counts = append(position_of_zero_counts, k)
      }      
    }
    
    #the next two lines delete the zero read counts conditions in calculating the fisher_I
    
    if(length(position_of_zero_counts)>0)
      #this means there's at least one condition has 0 reads
    {
      A_matrix_m1_new = A_matrix_m1[!is.na(A_matrix_m1)]
      #delete the NA A_matrix    
      N_vector_m1_new = N_vector_m1[!is.na(N_vector_m1)]
      
      tau_MLE_m1_new = tau_MLE_m1[-position_of_zero_counts]  
      #delete the ratio of 0
      
      W_vector_m1_new = W_vector_m1[-position_of_zero_counts]
    }
    
    else
      #there is no condition with 0 reads count  
    {
      A_matrix_m1_new = A_matrix_m1
      
      N_vector_m1_new = N_vector_m1
      
      tau_MLE_m1_new = tau_MLE_m1
      
      W_vector_m1_new = W_vector_m1
      
    }
    
    fisher_I_matrix=fisher_information_model1(A_matrix_m1_new, W_vector_m1_new, theta_MLE_m1, tau_MLE_m1_new)
    
    #check the fisher_I should be computational invertable and positive definite
    eigen_fisher_I = eigen(fisher_I_matrix, symmetric=T, only.values=T)$values
    
    r_condition = rcond(fisher_I_matrix)
    #the reciprocal of the condition number
    
    if(any(eigen_fisher_I<=0) || r_condition < 1e-9) 
      #fisher_I is not positive definite or the reciprocal of the condition number is too small
    {  
      Model_1_CI=list(theta_CI=NA, gene_CI=NA)
      
      return(Model_1_CI)
    }
    
    num_free_tau= length(tau_MLE_m1_new) - 1 #this is the number of free tau's, >=0
    
    dimention_para=dimention_theta + num_free_tau #dimention of free parameters=dimention of fisher_I_matrix, >=I
    
    C_matrix = solve(fisher_I_matrix)
    
    C_matrix=4*C_matrix
    
    C_matrix=(C_matrix + t(C_matrix))/2
    
    #draw 50000 sample from the proposal MVN dist
    N_samp=50000
    
    mean_imp=c(theta_MLE_m1, tau_MLE_m1_new[0:num_free_tau])
    #mean of the mvn or mvt
    
    sigma_imp=C_matrix
    
    #covariance of mvn or mvt
    
    #using the mvn 
    #theta_tau_imp_samp=rmvnorm(N_samp, mean=mean_imp, sigma=sigma_imp) ##note: need to load mvtnorm library
    
    #theta_tau_imp_density=dmvnorm(theta_tau_imp_samp, mean=mean_imp, sigma=sigma_imp)
    
    #The multivariate t distribution is not quite different as the multivariate normal distribution
    theta_tau_imp_samp=rmvt(N_samp, sigma=sigma_imp, delta=mean_imp, df=5) ##note: need to load mvtnorm library
    
    theta_tau_imp_density=dmvt(theta_tau_imp_samp, delta=mean_imp, sigma=sigma_imp, log=FALSE)
    
    #calculate tau[K], the last tau
    if(num_free_tau>=1)
    { 
      tau_imp_lastcol = 1-rowSums( as.matrix(theta_tau_imp_samp[,(dimention_theta+1):(dimention_para)]) )
    }
    
    else
      #in this case, num_free_tau=0
    {
      tau_imp_lastcol = rep(1, times=50000)
    }
    
    #merge the last tau into the imp_samp
    theta_tau_imp_samp = cbind(theta_tau_imp_samp, tau_imp_lastcol)#this is the sample includes tau[K]
    
    proposal_samp_weight=calculating_imp_weight(theta_tau_imp_samp, theta_tau_imp_density, A_matrix_m1, N_vector_m1,  W_vector_m1, log_likelihood_MLE_m1, model_select=1)
    
    proposal_samp_new=proposal_samp_weight$proposal_samp_new
    
    weight_norm=proposal_samp_weight$weight_norm
    
    #theta_samp_K_conditions=list()#this is a list containing K elements, each element is the theta_imp_sample for each k condition
    
    gene_imp_sample=list()
    
    theta_samp=as.matrix(proposal_samp_new[,1:dimention_theta])#number of columns is I
    
    tau_samp=as.matrix(proposal_samp_new[,(dimention_theta+1):(dimention_theta+num_free_tau+1)])#number of columns is K
    
    num_na=0#this will store the count of 0 read conditions in the following loop
    
    for(k in 1:K_conditions)
    {
      theta_result = matrix( rep(0, times=3*dimention_theta), ncol=3, nrow=dimention_theta )
      
      if(!is.na(match(k, position_of_zero_counts)))#this means the k'th condition has 0 reads counts
      {
        gene_result = c(0,0,0) #then CI for both theta and gene is 0 under this condition 
        
        num_na=num_na+1
      }
      
      else
      {
        theta_samp_k_condition=as.matrix( theta_samp*tau_samp[,(k-num_na)] )
        #the theta for the k'th condition is the theta*tau[k]
        
        for(i in 1:dimention_theta)
        {
          #calculate imp_point estimate and imp 95% CI through imp_estimate function
          theta_result[i,]=imp_estimate(theta_samp_k_condition[,i], weight_norm)
        }
        
        gene_imp_sample = rowSums(theta_samp_k_condition)
        
        gene_result = imp_estimate(gene_imp_sample, weight_norm)
      }
      
      colnames(theta_result)=c('imp_mean', 'CI_lower','CI_upper')
      
      names(gene_result)=c('imp_mean', 'CI_lower','CI_upper')
      
      Model_1_CI[[k]]=list(theta_CI=theta_result, gene_CI=gene_result)
      
    }
    
    return (Model_1_CI)
  }
}

#calculating CI for model 2. theta_MLE_m2 is a K*I matrix; log_likelihood_m2 is a K*1 vector; 
#A_matrix_m2, N_vector_m2, W_vector_m2 are all lists that contain K elements for each condition
Calculating_CI_for_model_2 = function(A_matrix_m2, N_vector_m2, W_vector_m2, theta_MLE_m2, log_likelihood_MLE_m2)
  #return a list of K elements, each elements is for each condition, which is a list containing CI_theta and CI_gene under that condition
{
  K_conditions=length(W_vector_m2)
  
  Model_2_CI=list()
  
  for(k in 1:K_conditions)
  {
    Model_2_CI[[k]]=Calculating_CI_for_theta(A_matrix_m2[[k]], N_vector_m2[[k]], W_vector_m2[[k]], theta_MLE_m2[k,], log_likelihood_MLE_m2[k])
  }
  
  return(Model_2_CI)
}

em_model0 = function(A_matrix, N_vector, W_vector)
  #return the MLE of theta under model 0, which is a vector of lenght I
  #A is a matrix; N and W are vectors
{    
  epsilon=1e-10
  
  I=length(W_vector)
  
  if( any(is.na(A_matrix)) && any(is.na(N_vector)) )
    #in this case, no read is mapped to any isoform of the gene in every condition. So the expression is 0.
  {
    theta=rep(0, times=I)
    
    return(theta)
  }
  
  else if( any(is.na(A_matrix))==F && any(is.na(N_vector))==F )
    #this is the normal case--neither of A and N is NA
  {
    theta_old=rep(1/I, times=I)
    
    theta_old_previous=rep(0, times=I)
    
    num_iter=0
    
    while( sum(abs(theta_old-theta_old_previous))>1e-6 )
    {
      num_iter=num_iter+1
      
      n_i_j_k = matrix( c(t(tcrossprod(theta_old, N_vector))) * c(t(A_matrix)) / c(crossprod(theta_old, A_matrix)), nrow=I, byrow=T)
      
      theta_new = rowSums(n_i_j_k) / W_vector + 1e-10
      
      theta_old_previous=theta_old
      
      theta_old=theta_new
    }
    
    theta=theta_new
    
    return(theta)
  }
  
  else 
    #in this case, only one of the A or N is NULL, which means there's something wrong with either the A or N.
  {
    print("There is something wrong with either the A_matrix or N_vector. One of them is NULL, but the other is not NULL.")
    
    return(NA)
  }
  
}

em_model1 = function(A_matrix, N_vector, W_vector)
  #this function return the MLE (contains two part: theta_MLE and tau_MLE) under model 1.
  #A_matrix, N_vector, W_vector are all list of length K.
  #A and N may contain NA's, which index the 0 reads conditions
{
  K_conditions=length(W_vector)
  
  num_isoform=length(W_vector[[1]])
  
  #data checking
  N_vector_temp = unlist(lapply(N_vector, na.exclude))
  
  A_matrix_temp = do.call(cbind, lapply(A_matrix, na.exclude))
  
  if(is.logical(A_matrix_temp) && is.logical(N_vector_temp))
    #in this case, no reads is mapped to the gene in every condition. Thus theta and ratio are all 0.
  {
    theta=rep(0, times=num_isoform)
    
    ratio=rep(0, times=K_conditions)
    
    estimate = list(theta=theta, ratio=ratio)
    
    return(estimate)
  }
  
  else if( is.logical(A_matrix_temp)==F && is.logical(N_vector_temp)==F )
    #this is the normal case, which means in at least one of the K conditions, there is some reads being mapped to the gene.
  {
    theta_old=rep(1/num_isoform, times=num_isoform)
    
    theta_old_previous=rep(0, times=num_isoform)
    
    num_iter=0
    
    ratio=rep(0, times=K_conditions)
    
    n_i_j=list()
    
    while( sum(abs(theta_old-theta_old_previous))>1e-6 )
    {
      num_iter=num_iter+1
      
      for(k in 1:K_conditions)
      { 
        
        
        if( any(is.na(A_matrix[[k]]))==F && any(is.na(N_vector[[k]]))==F )
          #this is the normal case
        {
          ratio[k]=sum(N_vector[[k]])/sum(W_vector[[k]]*theta_old)
          
          n_i_j[[k]] = matrix( c(t(tcrossprod(theta_old, N_vector[[k]]))) * c(t(A_matrix[[k]])) / c(crossprod(theta_old, A_matrix[[k]])), nrow=num_isoform, byrow=T)          
        }
        
        else if( (any(is.na(A_matrix[[k]]))==T && any(is.na(N_vector[[k]]))==F) 
                 || (any(is.na(A_matrix[[k]]))==F && any(is.na(N_vector[[k]]))==T) )
          #in this case, under the k'th condition, only one the A or N is NULL, which means there's something wrong with either the A or N.
        {
          k_error=k
          #record which condition is wrong
          
          estimate = list(theta=NA, ratio=NA, k_error=k_error)
          
          return(estimate)
          
        }
        
        #the third case is under the k'th condition, there is no read mapped to the gene(both A and N are NULL)
        #We don't need to do anything, since ratio[k] and n_i_j are pre-set to 0 and NULL.      
      }
      
      theta_new = rowSums(do.call(cbind, n_i_j)) / crossprod(do.call(rbind, W_vector), ratio)
      
      theta_old_previous=theta_old
      
      theta_old=theta_new
      
    }
    
    #Normalize ratio to sum up to 1 
    theta=theta_new * sum(ratio)
    
    ratio=ratio/sum(ratio)
    
    estimate = list(theta=theta, ratio=ratio)
    return(estimate)
  }
  
  else 
    #in this case, only one of the A or N is NULL(here A and N is the augmented ones), which means there's something wrong with either the A or N.
  {
    print("There is something wrong with either the A_matrix or N_vector. One of them is NA, but the other is not NA.")
    
    estimate = list(theta=NA, ratio=NA)
  }
}

likelihood_ratio_test = function(loglikelihood_MLE_m0, loglikelihood_MLE_m1, loglikelihood_MLE_m2, K_conditions, num_isoform, alpha)
  #return the best model and p-value
{
  LRT_model1_vs_model0 = (-2)*(loglikelihood_MLE_m0-loglikelihood_MLE_m1)
  
  P_model1_vs_model0 = 1-pchisq( LRT_model1_vs_model0, df=(K_conditions-1) )
  
  if(num_isoform==1)
    #for single-isoform gene, only compare model 0 vs model 1
  {
    if(P_model1_vs_model0 >= (2*alpha))
    {
      opt_model='model_0'
      
      return(list(opt_model=opt_model, p_value=P_model1_vs_model0))
    }
    
    else
    {
      opt_model='model_1'
      
      return(list(opt_model=opt_model, p_value=P_model1_vs_model0))
    }
    
  }
  
  LRT_model2_vs_model0 = (-2)*(loglikelihood_MLE_m0-loglikelihood_MLE_m2)
  
  P_model2_vs_model0 = 1 - pchisq( LRT_model2_vs_model0, df=(K_conditions-1)*num_isoform )
  
  p_value=c(P_model1_vs_model0, P_model2_vs_model0)
  
  if(P_model1_vs_model0 >= alpha  &&  P_model2_vs_model0 >= alpha)#in this case, model_0 is the opt_model
  {
    opt_model='model_0'
    
    return(list(opt_model=opt_model, p_value=p_value))
  }
  
  if(P_model1_vs_model0 >= alpha  &&  P_model2_vs_model0 < alpha)#model_2 is the opt model
  {
    
    opt_model='model_2'
    
    return(list(opt_model=opt_model, p_value=p_value))
  }
  
  if(P_model1_vs_model0 < alpha  &&  P_model2_vs_model0 >= alpha)#model_1 is the opt model
  {
    
    opt_model='model_1'
    
    return(list(opt_model=opt_model, p_value=p_value))
  }
  
  if(P_model1_vs_model0 < alpha  &&  P_model2_vs_model0 < alpha)#in this case, need to compare model1 and model2
  {
    LRT_model2_vs_model1 = (-2)*(loglikelihood_MLE_m1-loglikelihood_MLE_m2)
    
    P_model2_vs_model1 = 1 - pchisq( LRT_model2_vs_model1, df=(K_conditions-1)*(num_isoform-1) )
    
    p_value = c(p_value, P_model2_vs_model1)
    
    if(P_model2_vs_model1 >= (alpha*2))
    {
      
      opt_model='model_1'
      
      return(list(opt_model=opt_model, p_value=p_value))
    }
    
    else
    {
      
      opt_model='model_2'
      
      return(list(opt_model=opt_model, p_value=p_value))
    }
    
  }
}

read_file=function(file)
  #this function read one sampling rate file from output of rseq and convert the file to a list 
  #file is the path to the sampling rate file given by the user
{
  
  text=readLines(file)
  
  stringsplit=strsplit(text,"\t")
  
  J_lines=NA #this is a vector record the J of each gene, this is used how many lines should be skipped when start to read next gene block
  
  num_iter=0
  
  lines_skip=0
  
  final_result=list()
  
  while(!is.null(stringsplit[lines_skip+1][[1]]))
    #pass the list to final result
  {
    
    num_iter=num_iter+1
    
    gene_name=stringsplit[[lines_skip+1]][1]
    
    I=as.numeric(stringsplit[[lines_skip+1]][2])
    
    J_lines[num_iter]=as.numeric(stringsplit[[lines_skip+1]][3])
    
    J=J_lines[num_iter]
    
    isoform_names=stringsplit[[lines_skip+2]]
    
    temp_number_vector=as.numeric(unlist(stringsplit[(lines_skip+3):(lines_skip+J+3)]))
    
    temp_number_matrix=matrix(c(temp_number_vector,0), nrow=I+1, ncol=J+1)
    
    if(J==0)
      #in this case, no read is mapped to the gene under this condition
    {
      A_matrix=NA
      
      N_vector=NA
    }
    
    else
    {
      A_matrix=as.matrix(temp_number_matrix[1:I, 1:J])
      
      N_vector=temp_number_matrix[I+1, 1:J]
    }
    
    W_vector=temp_number_matrix[1:I, J+1]#W is always there even though there's zero reads counts mapped to the gene 
    
    gene=list(gene_name=gene_name, isoform_names=isoform_names, I=I, J=J,  A_matrix=A_matrix, N_vector=N_vector, W_vector=W_vector)
    
    final_result[[gene_name]]=gene
    
    lines_skip=sum(J_lines[0:num_iter])+3*num_iter
    
  }
  
  num_genes=num_iter
  
  final_result[['number_genes']]=num_genes
  
  return(final_result)
  
}

merge_file = function(list_K_files, condition_names)
{
  K_conditions=length(list_K_files)
  
  if(missing(condition_names))
  {
    print('The user did not specify condition/sample names. The program will rank the condition numbers according to the input order by the user.')
    
    condition_names=NA
    
    for(k in 1:K_conditions)
    {
      condition_names[k]=paste('condition_', k, sep='')
    }
  }
  
  num_genes=list_K_files[[1]]$number_genes
  
  final_result=list()
  
  for(num_iter in 1:num_genes)#num_iter index gene
  {
    
    #gene_name, isoform_names and I are the same across K conditons
    #J, A_matrix, N_vector and W_vector are different across K conditions
    
    #insert here: input checking
    
    gene_name=list_K_files[[1]][[num_iter]]$gene_name
    
    isoform_names=list_K_files[[1]][[num_iter]]$isoform_names
    
    I=list_K_files[[1]][[num_iter]]$I
    
    J=NA
    
    A_matrix=list()
    
    N_vector=list()
    
    W_vector=list()
    
    for(k in 1:K_conditions)#k indexes condition
    {
      J[k]=list_K_files[[k]][[num_iter]]$J
      
      A_matrix[[ condition_names[k] ]]=list_K_files[[k]][[num_iter]]$A_matrix      
      N_vector[[ condition_names[k] ]]=list_K_files[[k]][[num_iter]]$N_vector
      W_vector[[ condition_names[k] ]]=list_K_files[[k]][[num_iter]]$W_vector      
      
      
    }
    
    gene=list(gene_name=gene_name, isoform_names=isoform_names, I=I, J=J,  A_matrix=A_matrix, N_vector=N_vector, W_vector=W_vector)
    
    final_result[[gene_name]]=gene
    
  }
  
  return(final_result)
  
}

calculate_num_each_model = function(result)
  #result is the result returned by isoform_estimation function
{
  num_genes=length(result)
  
  num_model0_low_reads=NULL
  num_model0=NULL
  num_model1=NULL
  num_model2=NULL
  
  for(i in 1:num_genes)
  {
    if(result[[i]]$OPT_model_estimation$opt_model=='model_0(low_reads)')
    { num_model0_low_reads=append(num_model0_low_reads, i)    }
    
    if(result[[i]]$OPT_model_estimation$opt_model=='model_0')
    { num_model0=append(num_model0, i) }
    
    if(result[[i]]$OPT_model_estimation$opt_model=='model_1')
    { num_model1=append(num_model1, i) }
    
    if(result[[i]]$OPT_model_estimation$opt_model=='model_2')
    { num_model2=append(num_model2, i) }
  }
  
  print("Number of genes that have been successfully processed:")
  print(num_genes)
  
  print("Number of genes with low reads mapped to:")
  print(length(num_model0_low_reads))
  
  print("Number of genes that belongs model 0:")
  print(length(num_model0))
  
  print("Number of genes that belongs model 1:")
  print(length(num_model1))
  
  print("Number of genes that belongs model 2:")
  print(length(num_model2))
  
  count_each_model=list(num_genes=num_genes, num_model0_low_reads=num_model0_low_reads, num_model0=num_model0, num_model1=num_model1, num_model2=num_model2)
  
  return(count_each_model)
}

isoform_estimation = function(list_K_files, condition_names, CI=FALSE, alpha=0.025, read_lim=5)
{
  #check user input of alpha and read_lim
  if(alpha<0 || alpha>1)
    stop('Invalid alpha value!')
  
  if(read_lim<0 || read_lim>50)
    stop('Invalid read_lim value. Valid range: 0<= read_lim <=50')
  
  
  K_conditions=length(list_K_files)
  
  if(missing(condition_names))
  {
    print('The user did not specify condition/sample names. The program will rank the condition numbers according to the input order by the user.')
    
    condition_names=NA
    
    for(k in 1:K_conditions)
    {
      condition_names[k]=paste('condition_', k, sep='')
    }
  }
  
  if(length(condition_names)!=K_conditions)
  {
    stop("The number of condition/sample names given by the user does not match the number of condtions/samples. 
         Please check the input of condition_names.")
  }
  
  if(CI==T)
    #calculating CI will need the package mvtnorm  
  {
    library(mvtnorm)
  }
  
  num_genes=list_K_files[[1]]$number_genes
  
  merged_file = merge_file(list_K_files, condition_names)
  
  final_result=list()#contain the final result
  
  for(num_iter in 1:num_genes) #each iteration calculates one specific gene
  {    
    
    num_isoform = merged_file[[num_iter]]$I #number of isoforms of the gene, I
    
    gene_name=merged_file[[num_iter]]$gene_name
    
    isoform_names=merged_file[[num_iter]]$isoform_names
    
    #A_K, N_K and W_K are all list of length K 
    A_K=merged_file[[num_iter]]$A_matrix
    #note: A_K may contains NA's, where NA indicate the condition with 0 reads count
    
    N_K=merged_file[[num_iter]]$N_vector
    #note: N_K may contains NA's, where NA indicate the condition with 0 reads count
    
    W_K=merged_file[[num_iter]]$W_vector
    #note: W_K should not contains NA     
    
    #first check if total read counts less than read_lim*K, if so, return model 0
    total_read=sum(unlist(lapply(N_K, na.exclude)))
    
    if(total_read <= read_lim*K_conditions)
    {
      opt_model="model_0(low_reads)"
      
      theta_MLE=rep(0, times=num_isoform)

      names(theta_MLE)=isoform_names
      
      loglikelihood_MLE=0
      
      gene_MLE=rep(0, times=K_conditions)
      
      names(gene_MLE)=condition_names
      
      if(CI==FALSE)
      {
        final_result[[gene_name]]=list(OPT_model_estimation=list(gene_name=gene_name,
                                                                 isoform_names=isoform_names,
                                                                 condition_names=condition_names,
                                                                 opt_model=opt_model, 
                                                                 theta_MLE=theta_MLE, 
                                                                 loglikelihood_MLE=loglikelihood_MLE,
                                                                 gene_MLE=gene_MLE))
      }
      
      else
      {
        theta_result = matrix( rep(0, times=3*num_isoform), ncol=3, nrow=num_isoform )
        
        gene_result = c(0,0,0)
        
        colnames(theta_result) = c('imp_mean', 'CI_lower','CI_upper')
        
        names(gene_result) = c('imp_mean', 'CI_lower','CI_upper')
        
        CI_estimation = list(theta_CI=theta_result, gene_CI=gene_result)
        
        final_result[[gene_name]]=list(OPT_model_estimation=list(gene_name=gene_name,
                                                                 isoform_names=isoform_names,
                                                                 condition_names=condition_names,
                                                                 opt_model=opt_model, 
                                                                 theta_MLE=theta_MLE, 
                                                                 loglikelihood_MLE=loglikelihood_MLE,
                                                                 gene_MLE=gene_MLE,
                                                                 CI_estimation=CI_estimation)) 
      }
    }
    
    
    else
      #total reads greater than the lower limit, then do the test and calculations
    {  
      #data checking, record the conditions with 0 reads count  
      condition_zero_read=NULL#this is a vector recording the position of 0 reads conditions
      
      for(k in 1:K_conditions)
        #check and record the condition where A and K are both NA's
      {              
        if( any(is.na(A_K[[k]])) && any(is.na(N_K[[k]])) )
          #zero read counts happened, record it  
        {
          condition_zero_read = append(condition_zero_read, k)
        }
        
        else if( any(is.na(A_K[[k]]))==F && any(is.na(N_K[[k]]))==F )
          #this is the normal case--neither of A and N is NA
        {
          
          if(any(is.na( match(0, W_K[[k]]) ))==F )
            #in case w can contain zero(this happens in paired read case), add a small value to W
          {
            W_K[[k]] = W_K[[k]] + 1e-8
          }
        }
        
        else 
          #in this case, under the k'th condition, only one of the A or N is NA, which means there's something wrong with either the A or N.
        {
          print(paste("When calculating gene:", gene_name)) 
          
          print(paste("There is an error with either the A_matrix or N_vector under condition", k, "One of them is NA, but the other is not NA."))
          
          print("Please check the rseq output files of that gene. The program will continue for the next gene and the estimation for this gene will be marked as NA.")
          
          final_result[[gene_name]]=logical(0)#use logical(0) to mark the error
          
          break
        }
        
        
      }
      
      if(is.logical(final_result[[gene_name]]))
      {
        final_result[[gene_name]] = NA
        next
      }
      
      if(length(condition_zero_read)>0)
        #this means there's at least one condition has 0 reads
      {
        A_matrix_m0 = do.call(cbind, lapply(A_K, na.exclude)) 
        #this will remove NA's and catenate A_matrix of all the non-0-read condition together
        #A_matrix_m0 is a matrix of dimention=I*(J(K-M)), where M is the conditions with 0 reads counts
        #note: if all conditions have 0 reads, then N_vector_m0 = logical(0)
        
        N_vector_m0 = unlist(lapply(N_K, na.exclude))
        #this will remove NA's and catenate N_vector of all the non-0-read condition together
        #N_vector_m0 is a vector of length=J*(K-M), where M is the conditions with 0 reads counts
        #note: if all conditions have 0 reads, then N_vector_m0 = logical(0)
        
        W_vector_m0 = Reduce('+',W_K)
        
      }
      
      else
        #there is no condition with 0 reads count  
      {
        A_matrix_m0 = do.call(cbind, A_K)
        
        N_vector_m0 = unlist(N_K)
        
        W_vector_m0 = Reduce('+',W_K)
        
      }
      
      
      if( is.logical(A_matrix_m0)==F && is.logical(N_vector_m0)==F )
        #this means at least one condition has read counts more than 0
      {            
        
        theta_MLE_m0 = em_model0(A_matrix_m0, N_vector_m0, W_vector_m0)
        #theta_MLE_m0 is a vector of length I
        
        names(theta_MLE_m0)=isoform_names
        
        loglikelihood_MLE_m0 = log_likelihood_model0(A_matrix_m0, N_vector_m0, W_vector_m0, theta_MLE_m0)
        
        gene_MLE_model0 = sum(theta_MLE_m0)
        
        Model0_MLE=list(theta_MLE=theta_MLE_m0, loglikelihood_MLE=loglikelihood_MLE_m0, gene_MLE=gene_MLE_model0)
        
        
        #1-2.if there is only 1 condition, return model 0.
        if(K_conditions==1)
        {
          opt_model="model_0"
          
          theta_MLE=Model0_MLE$theta_MLE
          
          names(theta_MLE)=isoform_names
          
          loglikelihood_MLE=Model0_MLE$loglikelihood_MLE
          
          gene_MLE=Model0_MLE$gene_MLE
          
          names(gene_MLE)=condition_names
          
          if(CI==FALSE)
          {
            final_result[[gene_name]]=list(OPT_model_estimation=list(gene_name=gene_name,
                                                                     isoform_names=isoform_names,
                                                                     condition_names=condition_names,
                                                                     opt_model=opt_model, 
                                                                     theta_MLE=theta_MLE, 
                                                                     loglikelihood_MLE=loglikelihood_MLE,
                                                                     gene_MLE=gene_MLE))
          }
          
          else
          {
            CI_for_model0 = Calculating_CI_for_model_0(A_matrix_m0, N_vector_m0, W_vector_m0, theta_MLE_m0, loglikelihood_MLE_m0)
            
            rownames(CI_for_model0$theta_CI)=isoform_names
            
            final_result[[gene_name]]=list(OPT_model_estimation=list(gene_name=gene_name,
                                                                     isoform_names=isoform_names,
                                                                     condition_names=condition_names,
                                                                     opt_model=opt_model, 
                                                                     theta_MLE=theta_MLE, 
                                                                     loglikelihood_MLE=loglikelihood_MLE,
                                                                     gene_MLE=gene_MLE,
                                                                     CI_estimation=CI_for_model0))
          }
          
        }
        
        
        #if there are K>1 condition, calculate other model
        if(K_conditions>=2)
        {
          #1.3.model_1
          
          model1_em = em_model1(A_K, N_K, W_K)
          
          theta_MLE_m1 = model1_em$theta
          
          if( any(is.na(theta_MLE_m1)) )
            #if theta_MLE_m1 is wrong, then something is wrong with some specific conditions
          {
            k_error = model1_em$k_error
            
            print(paste("Something is wrong with the A_matrix or N_vector of gene:", gene_name, ", under", condition_names[k_error])) 
            
            print("Please check the rseq output files of that gene. The program will continue for the next gene and the estimation for this gene will be marked as NA.")
            
            final_result[[gene_name]]=NA      
            
            next
            #continue to do the next gene
          }
          
          
          rownames(theta_MLE_m1) = isoform_names
          
          ratio_MLE_m1 = model1_em$ratio
          
          names(ratio_MLE_m1) = condition_names
          
          theta_MLE_m1_under_each_condition = outer(ratio_MLE_m1, theta_MLE_m1)#this is the outer product of theta and ratio, which is a K*I matrix, each row is the theta_MLE for each condition
          
          rownames(theta_MLE_m1_under_each_condition)=condition_names
          
          colnames(theta_MLE_m1_under_each_condition)=isoform_names
          
          loglikelihood_MLE_m1 = log_likelihood_model1(A_K, N_K, W_K, theta_MLE_m1, ratio_MLE_m1)
          
          gene_MLE_model1 = rowSums(theta_MLE_m1_under_each_condition)
          
          names(gene_MLE_model1) = condition_names
          
          Model1_MLE=list(basic_theta_vector_MLE=theta_MLE_m1, 
                          tau_MLE=ratio_MLE_m1, 
                          theta_matrix_MLE = theta_MLE_m1_under_each_condition, 
                          loglikelihood_MLE=loglikelihood_MLE_m1, 
                          gene_MLE=gene_MLE_model1)
          
          
          #1.4.Model_2
          theta_MLE_m2 = matrix(rep(0, times=num_isoform*K_conditions), nrow=K_conditions)
          
          loglikelihood_MLE_m2 = rep(0, times=K_conditions)
          
          for(k in 1:K_conditions)
          {
            
            if( any(is.na(A_K[[k]]))==F && any(is.na(N_K[[k]]))==F )
              #this is the normal case--neither of A and N under condition k is NULL
            {        
              theta_MLE_m2[k,] = em_model0(A_K[[k]],N_K[[k]], W_K[[k]])
              
              loglikelihood_MLE_m2[k] = log_likelihood_model0(A_K[[k]], N_K[[k]], W_K[[k]], theta_MLE_m2[k,])
            }
            
            #if the k'th condition has 0 reads counts, we don't need to do anything
            #also, the wrong condition has been checked when calculating model 1
          }
          
          rownames(theta_MLE_m2)=condition_names
          
          colnames(theta_MLE_m2)=isoform_names
          
          loglikelihood_MLE_m2_sum = sum(loglikelihood_MLE_m2)
          
          gene_MLE_model2 = rowSums(theta_MLE_m2)
          
          names(gene_MLE_model2)=condition_names
          
          Model2_MLE=list(theta_matrix_MLE=theta_MLE_m2, loglikelihood_MLE=loglikelihood_MLE_m2_sum, gene_MLE=gene_MLE_model2)
          
          MLE_under_3_models=list(Model0_MLE=Model0_MLE, Model1_MLE=Model1_MLE, Model2_MLE=Model2_MLE)
          
          
          #1.5. Do the LRT test
          
          LRT_result = likelihood_ratio_test(loglikelihood_MLE_m0, loglikelihood_MLE_m1, loglikelihood_MLE_m2_sum, K_conditions, num_isoform, alpha)
          
          opt_model=LRT_result$opt_model
          
          p_value=LRT_result$p_value
          
          if(length(p_value)==2)
          {
            names(p_value)=c('model1_vs_model0', 'model2_vs_model0')
          }
          
          if(length(p_value)==3)
          {
            names(p_value)=c('model1_vs_model0', 'model2_vs_model0', 'model2_vs_model1')
          }
          
          #1.6. Return the result. If CI=T, calculate it
          
          #1.6.1 First calculate p1 (p value from model 1 vs 0), log2(fold_change); p2 (p value from model 2vs 0), T_value
          p1=p_value[1]
          
          log2.fold_change=log2(MLE_under_3_models$Model1_MLE$tau_MLE[2]/MLE_under_3_models$Model1_MLE$tau_MLE[1])
          
          if(log2.fold_change==Inf)
          {
            log2.fold_change=(.Machine$double.xmax)
          }
          
          if(log2.fold_change==-Inf)
          {
            log2.fold_change=(.Machine$double.xmax)
          }
          
          
          p2=p_value[2]
          
          theta_matrix=MLE_under_3_models$Model2_MLE$theta_matrix_MLE
          
          if(sum(theta_matrix[1,])==0 || sum(theta_matrix[2,])==0)
          {
            T_value=0
          }
          
          else
          {
            T_value=0.5*sum(abs(theta_matrix[1,]/sum(theta_matrix[1,])-theta_matrix[2,]/sum(theta_matrix[2,])))
          }
          
          if(opt_model=='model_0')
          {
            if(CI==TRUE)
            {
              CI_for_model0 = Calculating_CI_for_model_0(A_matrix_m0, N_vector_m0, W_vector_m0, theta_MLE_m0, loglikelihood_MLE_m0)        
              
              rownames(CI_for_model0$theta_CI)=isoform_names
              
              OPT_model_estimation=list(opt_model=opt_model,
                                        p_value=p_value,
                                        MLE=Model0_MLE,
                                        CI_estimation=CI_for_model0,
                                        p1=p1,
                                        log2.fold_change=log2.fold_change,
                                        p2=p2,
                                        T_value=T_value)
              
            }
            
            else
            {
              OPT_model_estimation=list(opt_model=opt_model,
                                        p_value=p_value,
                                        MLE=Model0_MLE,
                                        p1=p1,
                                        log2.fold_change=log2.fold_change,
                                        p2=p2,
                                        T_value=T_value)
            }
            
            final_result[[gene_name]]=list(gene_name=gene_name,
                                           isoform_names=isoform_names,
                                           condition_names=condition_names,
                                           MLE_under_3_models=MLE_under_3_models, 
                                           OPT_model_estimation=OPT_model_estimation)
          }
          
          if(opt_model=='model_1')
          {
            if(CI==T)
            {
              CI_for_model1 = Calculating_CI_for_model_1(A_K, N_K, W_K, theta_MLE_m1, ratio_MLE_m1, loglikelihood_MLE_m1)
              
              if(any(is.na(CI_for_model1[[1]]))==F)
              {
                for(k in 1:K_conditions)
                {
                  rownames(CI_for_model1[[k]]$theta_CI)=isoform_names                    
                }
                
                names(CI_for_model1) = condition_names
              }
              
              else
              {
                print(paste('The fisher information matrix for gene:', gene_name, 'is ill-conditioned. The 95% CI for the isoform abundance of that gene is NA.'))
              }
              
              OPT_model_estimation=list(opt_model=opt_model,
                                        p_value=p_value,
                                        MLE=Model1_MLE,
                                        CI_estimation=CI_for_model1,
                                        p1=p1,
                                        log2.fold_change=log2.fold_change,
                                        p2=p2,
                                        T_value=T_value)
              
            }
            
            else
            {
              OPT_model_estimation=list(opt_model=opt_model,
                                        p_value=p_value,
                                        MLE=Model1_MLE,
                                        p1=p1,
                                        log2.fold_change=log2.fold_change,
                                        p2=p2,
                                        T_value=T_value)
            }
            
            final_result[[gene_name]]=list(gene_name=gene_name,
                                           isoform_names=isoform_names,
                                           condition_names=condition_names,
                                           MLE_under_3_models=MLE_under_3_models, 
                                           OPT_model_estimation=OPT_model_estimation)
          }
          
          if(opt_model=='model_2')
          {
            if(CI==T)
            {
              CI_for_model2 = Calculating_CI_for_model_2(A_K, N_K, W_K, theta_MLE_m2, loglikelihood_MLE_m2)
              
              for(k in 1:K_conditions)
              {
                rownames(CI_for_model2[[k]]$theta_CI)=isoform_names                    
              }
              
              names(CI_for_model2) = condition_names
              
              OPT_model_estimation=list(opt_model=opt_model,
                                        p_value=p_value,
                                        MLE=Model2_MLE,
                                        CI_estimation=CI_for_model2,
                                        p1=p1,
                                        log2.fold_change=log2.fold_change,
                                        p2=p2,
                                        T_value=T_value)
              
            }
            
            else
            {
              OPT_model_estimation=list(opt_model=opt_model,
                                        p_value=p_value,
                                        MLE=Model2_MLE,
                                        p1=p1,
                                        log2.fold_change=log2.fold_change,
                                        p2=p2,
                                        T_value=T_value)
            }
            
            final_result[[gene_name]]=list(gene_name=gene_name,
                                           isoform_names=isoform_names,
                                           condition_names=condition_names,
                                           MLE_under_3_models=MLE_under_3_models, 
                                           OPT_model_estimation=OPT_model_estimation)       
          }      
        }    
      }
      
      
      
      else
        #in this case, either the A is all NA's but N is not, or N is all NA's but A is not, which means something is wrong.
        #give the gene name to the user  
      {
        print(paste("Something is wrong with the A_matrix or N_vector of gene:", gene_name)) 
        
        print("Please check the rseq output files of that gene. The program will continue for the next gene and the estimation for this gene will be marked as NA.")
        
        final_result[[gene_name]]=NA    
      }
      
      
      
    } 
    
    #print(paste('gene', num_iter, 'done!'))
  }
  
  count_each_model=calculate_num_each_model(final_result)
  
  return(final_result)
  
  }

output_model0=function(result)
{
  num_genes=length(result)
  condition_names=result[[1]]$condition_names
  
  K=length(condition_names)
  
  gene_name=NULL
  isoform_name=NULL
  model=NULL
  log_likelihood=NULL
  p_value=NULL  
  theta_MLE=NULL
  #theta_MLE is K column matrix, K is the number of conditions
  
  gene_MLE=NULL
  #gene_MLE is K column matrix, K is the number of conditions
  
  p1=NULL
  log2.fold_change=NULL
  p2=NULL
  T_value=NULL
  
  num_model0=0
  
  for(i in 1:num_genes)
  {   
    
    if(result[[i]]$OPT_model_estimation$opt_model=='model_0')
    { 
      num_model0=num_model0+1
      
      gene_name=append(gene_name, result[[i]]$gene_name)
      
      isoform_name_temp=paste(result[[i]]$isoform_names, collapse=",")
      
      isoform_name=append(isoform_name, isoform_name_temp)
      
      model=append(model, result[[i]]$OPT_model_estimation$opt_model)
      
      log_likelihood=append(log_likelihood, result[[i]]$OPT_model_estimation$MLE$loglikelihood_MLE)
      
      p_value_temp=paste(result[[i]]$OPT_model_estimation$p_value, collapse=",")
      
      p_value=append(p_value, p_value_temp)
      
      theta_temp=paste(result[[i]]$OPT_model_estimation$MLE$theta_MLE, collapse=",")
      
      theta_MLE=append(theta_MLE, theta_temp)
      
      gene_MLE=rbind(gene_MLE, result[[i]]$OPT_model_estimation$MLE$gene_MLE)
      
      p1=append(p1, result[[i]]$OPT_model_estimation$p1)
      
      log2.fold_change=append(log2.fold_change, result[[i]]$OPT_model_estimation$log2.fold_change)
      
      p2=append(p2, result[[i]]$OPT_model_estimation$p2)
      
      T_value=append(T_value, result[[i]]$OPT_model_estimation$T_value)      
    }
    
  }
  
  print(paste(num_model0, "genes belong to model 0."))
  
  options(warn=-1)
  genes_model0=data.frame(gene_name, isoform_name, model, log_likelihood, p_value, 
                          theta_MLE, gene_MLE, p1, log2.fold_change, p2, T_value)
  options(warn=0)
  
  return(genes_model0)
  
}

output_model1=function(result)
{
  num_genes=length(result)
  condition_names=result[[1]]$condition_names
  
  K=length(condition_names)
  
  gene_name=NULL
  isoform_name=NULL
  model=NULL
  log_likelihood=NULL
  p_value=NULL
  basic_theta_vector_MLE=NULL
  tau_MLE=NULL
  
  theta_MLE=NULL
  #theta_MLE is K column matrix, K is the number of conditions
  
  gene_MLE=NULL
  #gene_MLE is K column matrix, K is the number of conditions
  
  p1=NULL
  log2.fold_change=NULL
  p2=NULL
  T_value=NULL
  
  num_model1=0
  
  for(i in 1:num_genes)
  {   
    
    if(result[[i]]$OPT_model_estimation$opt_model=='model_1')
    { 
      num_model1=num_model1+1
      
      gene_name=append(gene_name, result[[i]]$gene_name)
      
      isoform_name_temp=paste(result[[i]]$isoform_names, collapse=",")
      
      isoform_name=append(isoform_name, isoform_name_temp)
      
      model=append(model, result[[i]]$OPT_model_estimation$opt_model)
      
      log_likelihood=append(log_likelihood, result[[i]]$OPT_model_estimation$MLE$loglikelihood_MLE)
      
      p_value_temp=paste(result[[i]]$OPT_model_estimation$p_value, collapse=",")
      
      p_value=append(p_value, p_value_temp)
      
      basic_theta_vector_temp=paste(as.vector(result[[i]]$OPT_model_estimation$MLE$basic_theta_vector_MLE), collapse=",")
      
      basic_theta_vector_MLE=append(basic_theta_vector_MLE, basic_theta_vector_temp)
      
      tau_MLE=rbind(tau_MLE, result[[i]]$OPT_model_estimation$MLE$tau_MLE)
      
      theta_matrix=result[[i]]$OPT_model_estimation$MLE$theta_matrix_MLE
      
      theta_temp=apply(theta_matrix, MARGIN=1, paste, collapse=',')
      
      theta_MLE=rbind(theta_MLE, theta_temp)
      
      gene_MLE=rbind(gene_MLE, result[[i]]$OPT_model_estimation$MLE$gene_MLE)
      
      p1=append(p1, result[[i]]$OPT_model_estimation$p1)
      
      log2.fold_change=append(log2.fold_change, result[[i]]$OPT_model_estimation$log2.fold_change)
      
      p2=append(p2, result[[i]]$OPT_model_estimation$p2)
      
      T_value=append(T_value, result[[i]]$OPT_model_estimation$T_value)      
    }
    
  }
  
  print(paste(num_model1, "genes belong to model 1."))
  
  tau_MLE=as.data.frame(tau_MLE)
  
  theta_MLE=as.data.frame(theta_MLE)
  
  gene_MLE=as.data.frame(gene_MLE)
  
  for(k in 1:K)
  {
    
    names(tau_MLE)[k] = paste('tau_MLE_', condition_names[k], sep='')
    
    names(theta_MLE)[k] = paste('theta_MLE_', condition_names[k], sep='')
    
    names(gene_MLE)[k] = paste('gene_MLE_', condition_names[k], sep='')
    
  }
  
  options(warn=-1)
  genes_model1=data.frame(gene_name, isoform_name, model, log_likelihood, p_value, 
                          basic_theta_vector_MLE, tau_MLE,  theta_MLE, gene_MLE, 
                          p1, log2.fold_change, p2, T_value)
  options(warn=0)
  
  model1_rank_by_fold=genes_model1[order(genes_model1$log2.fold_change, decreasing = TRUE),]
  
  return(model1_rank_by_fold)
  
}

output_model2=function(result)
{
  num_genes=length(result)
  condition_names=result[[1]]$condition_names
  
  K=length(condition_names)
  
  gene_name=NULL
  isoform_name=NULL
  model=NULL
  log_likelihood=NULL
  p_value=NULL
  gene_MLE=NULL
  #gene_MLE is K column matrix, K is the number of conditions
  
  theta_MLE=NULL
  #theta_MLE is K column matrix, K is the number of conditions
  
  p1=NULL
  log2.fold_change=NULL
  p2=NULL
  T_value=NULL
  
  num_model2=0
  
  for(i in 1:num_genes)
  {   
    
    if(result[[i]]$OPT_model_estimation$opt_model=='model_2')
    { 
      num_model2=num_model2+1
      
      gene_name=append(gene_name, result[[i]]$gene_name)
      
      isoform_name_temp=paste(result[[i]]$isoform_names, collapse=",")
      
      isoform_name=append(isoform_name, isoform_name_temp)
      
      model=append(model, result[[i]]$OPT_model_estimation$opt_model)
      
      log_likelihood=append(log_likelihood, result[[i]]$OPT_model_estimation$MLE$loglikelihood_MLE)
      
      p_value_temp=paste(result[[i]]$OPT_model_estimation$p_value, collapse=",")
      
      p_value=append(p_value, p_value_temp)
      
      gene_MLE=rbind(gene_MLE, result[[i]]$OPT_model_estimation$MLE$gene_MLE)
      
      theta_matrix=result[[i]]$OPT_model_estimation$MLE$theta_matrix_MLE
      
      theta_temp=apply(theta_matrix, MARGIN=1, paste, collapse=',')
      
      theta_MLE=rbind(theta_MLE, theta_temp)
      
      p1=append(p1, result[[i]]$OPT_model_estimation$p1)
      
      log2.fold_change=append(log2.fold_change, result[[i]]$OPT_model_estimation$log2.fold_change)
      
      p2=append(p2, result[[i]]$OPT_model_estimation$p2)
      
      T_value=append(T_value, result[[i]]$OPT_model_estimation$T_value)            
    }
    
  }
  
  print(paste(num_model2, "genes belong to model 2."))
  
  gene_MLE=as.data.frame(gene_MLE)
  
  theta_MLE=as.data.frame(theta_MLE)
  
  for(k in 1:K)
  {
    names(gene_MLE)[k] = paste('gene_MLE_', condition_names[k], sep='')
    
    names(theta_MLE)[k] = paste('theta_MLE_', condition_names[k], sep='')
  }
  
  options(warn=-1)
  genes_model2=data.frame(gene_name, isoform_name, model, log_likelihood, p_value, 
                          gene_MLE, theta_MLE, p1, log2.fold_change, p2, T_value)
  options(warn=0)
  
  model2_rank_by_T=genes_model2[order(genes_model2$T_value, decreasing = TRUE),]
  
  return(model2_rank_by_T)
  
}

output_each_model=function(result, model)
  #this function output a dataframe for all three models
{
  
  if(model=="model0")
  {
    gene_list=output_model0(result)
  }
  
  if(model=="model1")
  {
    gene_list=output_model1(result)
  }
  
  if(model=="model2")
  {
    gene_list=output_model2(result)
  }
  
  if(model=="all")
  {
    list0=output_model0(result);
    
    list1=output_model1(result);
    
    list2=output_model2(result);
    
    gene_list=list(model0=list0, model1=list1, model2=list2)
  }
  
  return(gene_list)
}

volcano_DF=function(result, pch=20, ...)
{
  p1=NULL
  
  log2.fold_change=NULL
  
  num_genes=length(result)
  
  for(i in 1:num_genes)
  {
    if(result[[i]]$OPT_model_estimation$opt_model != "model_0(low_reads)")
    {
      if(result[[i]]$OPT_model_estimation$p1==0)
      {
        p1_temp=.Machine$double.neg.eps
      }
      
      else
      {
        p1_temp=result[[i]]$OPT_model_estimation$p1
      }
      
      p1=append(p1, p1_temp)
      
      log2.fold_change=append(log2.fold_change, result[[i]]$OPT_model_estimation$log2.fold_change)
    }
  }
  
  lg10p1=log10(p1)
  
  plot(log2.fold_change, -lg10p1, pch=pch, ...)
  
}

volcano_DS=function(result, range=c(0, 1), pch=20, ...)
{
  p2=NULL
  
  T_value=NULL
  
  num_genes=length(result)
  
  for(i in 1:num_genes)
  {
    if(result[[i]]$OPT_model_estimation$opt_model != "model_0(low_reads)" && (!is.na(result[[i]]$OPT_model_estimation$p2)))
    {
      if(result[[i]]$OPT_model_estimation$p2==0)
      {
        p2_temp=.Machine$double.neg.eps
      }
      
      else
      {
        p2_temp=result[[i]]$OPT_model_estimation$p2
      }
      
      p2=append(p2, p2_temp)
      
      T_value=append(T_value, result[[i]]$OPT_model_estimation$T_value)
    }
  }
  
  lg10p2=log10(p2)
  
  plot(T_value, -lg10p2, xlim=range, pch=pch, ...)
  
}
