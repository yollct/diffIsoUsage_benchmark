source("/code.R")
library(dplyr)

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
    
    isoform_names=merged_file[[num_iter]]$isoform_names[4]
    
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
      
      ### change
      num_isoform = ifelse(num_isoform==0, 1, num_isoform)

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
      print(K_conditions)
      print(A_K)
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

args <- commandArgs(trailingOnly = T)
outdir <- args[1]
con1 <- args[2]
con2 <- args[3]

metadata <- read.csv("/MOUNT/meta.txt", sep="\t")
group1 <- metadata %>% dplyr::filter(group==con1)
group2 <- metadata %>% dplyr::filter(group==con2)

rseqdir <- "/MOUNT/rseqdiff"

exp.list <- lapply(c(group1$sample_id[1], group2$sample_id[2]), function(y){paste0(rseqdir, "/", y, "_out.txt.3.b.exp")}) %>% unlist
print(exp.list)

list2files <- list()
i<-1
for (x in exp.list){
    tmp <- read_file(x)
    list2files[[i]] <- tmp
    i <<- i+1
}

print("check")
print(length(list2files))



expression = isoform_estimation(list2files, unique(metadata$group))
print("done")
print(head(expression))