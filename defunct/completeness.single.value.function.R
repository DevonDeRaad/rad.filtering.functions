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




#functionize this shit

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

filter.completeness <- function(x){
  #turn genlight into snp matrix
  gen.mat<-as.matrix(x)
  
  #loop that stores total loci retained for each possible completeness filter value
  loci.retained <- vector("numeric", length = nrow(gen.mat))
  filt<- vector("numeric", length = nrow(gen.mat))
  for (i in nrow(gen.mat):1){
    #calculate the number of loci retained
    loci.retained[i]<-ncol(gen.mat[,(colSums(is.na(gen.mat)) < i)])
    #calculate the completeness cutoff for each snp to be retained
    filt[i]<-1-((i-1)/nrow(gen.mat))
  }
  
  #build df
  df.x<-as.data.frame(cbind(filt, loci.retained))
  #return plot of df and df itself
  print( 
    ggplot(df.x, aes(x=filt, y=loci.retained))+
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
df.x<-filter.completeness(todi.full.genlight)
#success






