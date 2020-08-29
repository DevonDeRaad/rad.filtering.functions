library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap) 
library(ggplot2)


#This only needs to be done once at the beginning of the script to avoid repeatedly reading in a large vcf
#read in vcf as vcfR
vcfR <- read.vcfR(x)
### convert vcfR to genlight
genlight <- vcfR2genlight(vcfR)
#turn genlight into snp matrix
gen.mat<-as.matrix(genlight)


missing.by.snp <- function(gen.mat){
  
  #calculate the proportion of individuals successfully genotyped at each SNP
  miss<-colSums(is.na(gen.mat) == TRUE)/nrow(gen.mat)
  #list snp names
  snp<-colnames(gen.mat)
  #create dataframe of missing data by snp
  miss.by.snp<-data.frame(snp=snp, proportion.missing=as.numeric(miss), row.names = NULL)

  #make histogram of the data
  print(
  ggplot(miss.by.snp, aes(x=proportion.missing))+
    geom_histogram(color="black", fill="white")+
    theme_classic()+
    geom_vline(aes(xintercept=mean(proportion.missing)),
               color="red", linetype="dashed", size=1)
  )
  
  #return the dataframe
  return(miss.by.snp)
}

miss<-missing.by.snp(gen.mat=gen.mat)


  