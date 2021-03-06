
#open function
max.depth <- function(vcfR, maxdepth=NULL){
  
  #extract depth from the vcf
  dp.matrix<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)
  
  #calculate vector of depth per SNP
  snpdepth<-rowSums(dp.matrix, na.rm = TRUE)/rowSums(is.na(dp.matrix) == FALSE)
  
  if (is.null(maxdepth)){
    par(mfrow=c(2,1))
    #plot histogram of depth
    hist(snpdepth, xlab = "mean of the depth of all samples successfully genotyped at a given SNP",
         main="histogram showing the depth of all called SNPs")
    abline(v=mean(snpdepth, na.rm = TRUE), col="red", lty="dashed")
    
    #zoomed in histogram
    hist(snpdepth[snpdepth < 200], xlab = "mean of the depth of all samples successfully genotyped at a given SNP",
         main ="distribution of SNPs below a depth of 200")
    abline(v=mean(snpdepth, na.rm = TRUE), col="red", lty="dashed")
    
    #print
    print(paste0("dashed line indicates a mean depth across all SNPs of ",round(mean(snpdepth, na.rm = TRUE),1)))
    
  }
  
  else {
    #plot the maxdepth cutoff
    hist(snpdepth, xlab = "mean of the depth of all samples successfully genotyped at a given SNP", main ="max depth cutoff")
    abline(v=maxdepth, col="red")
    
    #calculate % of snps that fail the max depth filter
    i<-round(sum(snpdepth > maxdepth, na.rm = TRUE)/length(snpdepth)*100, 2)
    
    print(paste0(i, "% of SNPs were above a mean depth of ", maxdepth, " and were removed from the vcf"))
    
    #filter vcf
    vcfR<-vcfR[snpdepth < maxdepth,]
    
    return(vcfR)
  }
}

