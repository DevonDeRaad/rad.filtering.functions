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

#overlayed dot plots of the number of total snps retained in one color
#and the mean missingness percentage of the remaining snps in another color

compare.snps.retained <- function(x){
  
  #calculate the proportion of individuals successfully genotyped at each SNP
  miss<-colSums(is.na(gen.mat) == TRUE)/nrow(gen.mat)
  
  #loop that stores a vector of # non-missing loci retained for each individual
  #looped over all possible completeness filter values
  #initialize df.x
  df.x<- data.frame(filt=numeric(), missingness=numeric(), snps.retained=numeric())
  #loop
  for (i in c(0,.1,.2,.3,.4,.5,.6,.65,.7,.75,.8,.85,.9,.95,1)){
    #subset the matrix to only snps passing a given filter level
    gen.z.mat<-gen.mat[,1-miss >= i]
    #calc the total number of snps retained at the given filter level
    snps.retained<-ncol(gen.z.mat)
    #calculate the total missing data % for the given cutoff
    missingness<-sum(is.na(gen.z.mat))/(ncol(gen.z.mat)*nrow(gen.z.mat))
    #calculate the completeness cutoff for each snp to be retained
    filt<-i
    df.x<-rbind(df.x, as.data.frame(cbind(filt, missingness, snps.retained)))
  }
  
  #make columns numeric for plotting
  df.x$filt<-as.numeric(as.character(df.x$filt))
  df.x$missingness<-as.numeric(as.character(df.x$missingness))
  df.x$snps.retained<-as.numeric(as.character(df.x$snps.retained))
  
  #visualize as a dotplot for each filter level
  plot1<-ggplot(df.x, aes(x=as.factor(filt))) + 
    geom_point(aes(y=snps.retained)) +
    theme_bw() +
    labs(x = "SNP completeness cutoff", y = "total loci retained")
  
  plot2<-ggplot(df.x, aes(x=as.factor(filt))) + 
    geom_point(aes(y=missingness)) +
    theme_bw() +
    labs(x = "SNP completeness cutoff", y = "total % missing data")
  
  print(grid.arrange(plot1,plot2))
}


#run filter function
compare.snps.retained(gen.mat)






