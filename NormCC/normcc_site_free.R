#Normalize metagenomic Hi-C data and detect spurious contacts using zero-inflated Negative Binominal regression frameworks
#Auther and maintainer: Yuxuan Du <yuxuandu@usc.edu>
#HiCzin R script depends on 'glmmTMB' package
library('glmmTMB')


normcc = function(contig_info_file , row_sum_file)
{
  contig_info = read.csv(contig_info_file , header = F , sep = ',' )
  contig_info = as.data.frame(contig_info)
  colnames(contig_info) = c('contig_name' , 'length' , 'covcc' , 'signal')
  
  data_sample = cbind(log(contig_info$length) , log(contig_info$covcc) , contig_info$signal)
  data_sample = as.data.frame(data_sample)
  colnames(data_sample) = c('sample_len' , 'sample_covcc', 'sampleCon')
  
  tryCatch(
    {
      
      fit1 = glmmTMB(sampleCon~sample_len+sample_covcc, data = data_sample,
                     ziformula=~0 , family=nbinom2)
        
      },
      error = function(e){
        message(e)
        message(paste("\nskip",  sep=" "))
      },
      warning = function(w){
        message(w)
        message(paste("\nskip",  sep=" "))
      }
    )
  coeff = as.numeric(fit1$fit$par)
  result = c(coeff[1:3])
  return(result)
}













