#Normalize metagenomic Hi-C data and detect spurious contacts using zero-inflated Negative Binominal regression frameworks
#Auther and maintainer: Yuxuan Du <yuxuandu@usc.edu>
#HiCzin R script depends on 'glmmTMB' package
library('glmmTMB')


normcc = function(contig_info_file , row_sum_file)
{
  row_sum_contact = read.csv(row_sum_file , sep = ',' , header  = F)
  sampleCon = as.numeric(row_sum_contact[[1]])
  contig_info = read.csv(contig_info_file , header = F , sep = ',' )
  contig_info = as.data.frame(contig_info)
  colnames(contig_info) = c('contig_name' , 'site' , 'length' , 'coverage')
    
  data_sample = cbind(log(contig_info$site+1) , log(contig_info$length) , 
                      log(contig_info$cov) , sampleCon)
  data_sample = as.data.frame(data_sample)
  colnames(data_sample) = c('sample_site' , 'sample_len' , 'sample_cov', 'sampleCon')
    
  tryCatch(
      {
        
        fit1 = glmmTMB(sampleCon~sample_site+sample_len+sample_cov, data = data_sample,
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
  result = c(coeff[1:4])
  return(result)
}













