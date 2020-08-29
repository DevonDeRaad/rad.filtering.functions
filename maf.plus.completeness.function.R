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


#pseudocode
# block that will maf filter the original genlight object
#block that will loop over genlight and create filtered genlight at each possible level of completeness
#block will pull the distibution of NAs (missing genotypes) from the genlight object and hist it
#put all hists together in a single graph and label by filt scheme

#function takes a genlight object as input
#turns genlight into matrix
#filters by maf if applicable
#uses for loop to calulate:
#(number of loci retained at each possible filtering level and store in a vector)
#(corresponding filtration percentage)
#sticks those two vectors together with cbind
#function returns you this dataframe to make your own plots from if you wish
#and autoprints to your plot window a ggplot of (x = matrix completeness, y = # of loci retained)
#fin

maf.filter.completeness <- function(x){
  #turn genlight into snp matrix
  gen.mat<-as.matrix(x)
  
  #filter genlight to maf = .01, .03, .05
  #calc maf: divide the number of total alternate alleles by 2 * the number of called genotypes
  colSums(gen.mat, na.rm = T) / (2*(colSums(is.na(gen.mat)=="FALSE")))
  gen.mat.01<-gen.mat[,(colSums(gen.mat, na.rm = T) / (2*(colSums(is.na(gen.mat)=="FALSE"))) > .01)]
  gen.mat.03<-gen.mat[,(colSums(gen.mat, na.rm = T) / (2*(colSums(is.na(gen.mat)=="FALSE"))) > .03)]
  gen.mat.05<-gen.mat[,(colSums(gen.mat, na.rm = T) / (2*(colSums(is.na(gen.mat)=="FALSE"))) > .05)]
  
  #loop that stores total loci retained for each possible completeness filter value no maf
  loci.retained <- vector("numeric", length = nrow(gen.mat))
  filt<- vector("numeric", length = nrow(gen.mat))
  maf<- rep(0, times =nrow(gen.mat))
  for (i in nrow(gen.mat):1){
    #calculate the number of loci retained
    loci.retained[i]<-ncol(gen.mat[,(colSums(is.na(gen.mat)) < i)])
    #calculate the completeness cutoff for each snp to be retained
    filt[i]<-1-((i-1)/nrow(gen.mat))
  }
  #build df
  df.x<-as.data.frame(cbind(maf, filt, loci.retained))
  
  #loop that stores total loci retained for each possible completeness filter value maf=.01
  loci.retained <- vector("numeric", length = nrow(gen.mat.01))
  filt<- vector("numeric", length = nrow(gen.mat.01))
  maf<- rep(.01, times =nrow(gen.mat.01))
  for (i in nrow(gen.mat.01):1){
    #calculate the number of loci retained
    loci.retained[i]<-ncol(gen.mat[,(colSums(is.na(gen.mat.01)) < i)])
    #calculate the completeness cutoff for each snp to be retained
    filt[i]<-1-((i-1)/nrow(gen.mat.01))
  }
  #append new df to old df
  df.x<-rbind(df.x, as.data.frame(cbind(maf, filt, loci.retained)))
  
  #loop that stores total loci retained for each possible completeness filter value maf=.03
  loci.retained <- vector("numeric", length = nrow(gen.mat.03))
  filt<- vector("numeric", length = nrow(gen.mat.03))
  maf<- rep(.03, times =nrow(gen.mat.03))
  for (i in nrow(gen.mat.03):1){
    #calculate the number of loci retained
    loci.retained[i]<-ncol(gen.mat[,(colSums(is.na(gen.mat.03)) < i)])
    #calculate the completeness cutoff for each snp to be retained
    filt[i]<-1-((i-1)/nrow(gen.mat.03))
  }
  #append new df to old df
  df.x<-rbind(df.x, as.data.frame(cbind(maf, filt, loci.retained)))
  
  #loop that stores total loci retained for each possible completeness filter value maf=.05
  loci.retained <- vector("numeric", length = nrow(gen.mat.05))
  filt<- vector("numeric", length = nrow(gen.mat.05))
  maf<- rep(.05, times =nrow(gen.mat.05))
  for (i in nrow(gen.mat.05):1){
    #calculate the number of loci retained
    loci.retained[i]<-ncol(gen.mat[,(colSums(is.na(gen.mat.05)) < i)])
    #calculate the completeness cutoff for each snp to be retained
    filt[i]<-1-((i-1)/nrow(gen.mat.05))
  }
  #append new df to old df
  df.x<-rbind(df.x, as.data.frame(cbind(maf, filt, loci.retained)))
  
  
  #turn maf column into categorical for plotting
  df.x$maf<-as.character(df.x$maf)
  #return plot of df and df itself
  print( 
    ggplot(df.x, aes(x=filt, y=loci.retained, color=maf))+
      geom_point()+
      ggtitle("SNPs retained by filtering scheme") +
      xlab("fraction of non-missing genotypes required to retain each SNP (0-1)") + ylab("total SNPs retained")+
      theme_light()
  )
  return(df.x)
}



#give it a whirl
todi.vcf <- read.vcfR("~/Downloads/todi.allsnps.vcf") #read in all data
head(todi.vcf) #check the vcf object
### convert to genlight
todi.full.genlight <- vcfR2genlight(todi.vcf)

#run filter function
maf.filter.completeness(todi.full.genlight)
#success


