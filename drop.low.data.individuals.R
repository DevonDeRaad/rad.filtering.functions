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


#this function should be used in conjunction with missing.data.by.individual()
#output must end in ".vcf.gz"

drop.low.data.individuals <- function(input.vcf, missing.cutoff=NULL, output){
  #turn genlight into snp matrix
  vcf.r<- read.vcfR(input.vcf) #read in all data
  #convert to genlight
  genlight<- vcfR2genlight(vcf.r)
  #convert genlight to matrix
  gen.mat<- as.matrix(genlight)
  #get a vector that tells you miss % of each sample
  miss.dat<-(rowSums(is.na(gen.mat)))/(ncol(gen.mat))
  #filter genlight to drop samples (matrix rows) above the missing data cutoff defined by user
  filtered.genlight <- new("genlight", gen.mat[miss.dat < missing.cutoff,])
  #print message to user
  print(paste("filtered vcf contains all input SNPs and input samples with less than", missing.cutoff, "% missing genotypes"))
  #convert genlight to vcfR
  filtered.vcf.r<-
  #write filtered genlight out as an output vcf
  write.vcf(filtered.vcf.r, file= out )
}
  
  
  drop.low.data.individuals(input.vcf="~/Downloads/todi.vcf", missing.cutoff=.7, out="~/Downloads/todi.filt.vcf.gz" )
  
  
  
  #initialize dataframe
  df.x<- data.frame(indiv=character(), loci.retained=numeric(), filt=numeric())
  #loop
  for (i in nrow(gen.mat):1){
    #get vector of individuals we are looping over
    indiv<-rownames(gen.mat)
    #calculate the completeness cutoff for each snp to be retained in this iteration
    filt<-rep(1-((i-1)/nrow(gen.mat)), times =nrow(gen.mat))
    #calculate the number of loci successfully genotyped in each sample in this iteration
    snps.retained<-rowSums(is.na(gen.mat[,(colSums(is.na(gen.mat)) < i)]) == "FALSE")
    #append the data from this iteration to existing df
    df.x<-rbind(df.x, as.data.frame(cbind(indiv, filt, snps.retained)))
    #close for loop
  }
  
  #make columns numeric for plotting
  df.x$filt<-as.numeric(df.x$filt)
  df.x$snps.retained<-as.numeric(df.x$snps.retained)
  #visualize color coded by individual
  print( 
    ggplot(df.x, aes(x=filt, y=snps.retained, color=indiv))+
      geom_point()+
      ggtitle("SNPs retained by filtering scheme") +
      xlab("fraction of non-missing genotypes required to retain each SNP (0-1)") + ylab("non-missing SNP genotypes")+
      theme_light()
  )
  
  #calculate average missing data per individual so that we can figure out which individuals are pulling distribution down
  print("missing data % in each individual in input vcf")
  miss.dat<-(rowSums(is.na(gen.mat)))/(ncol(gen.mat))
  #
  dotchart(sort(miss.dat),labels=row.names(miss.dat), cex=.7, xlab="proportion of missing genotypes")
  
  return(df.x)
}



#give it a whirl

missing.data.by.individual("~/Downloads/todi.allsnps.vcf")
#success






