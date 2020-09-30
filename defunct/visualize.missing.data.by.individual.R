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

#missing.data.by.individual("~/Downloads/todi.vcf")


missing.data.by.individual <- function(input.vcf){
  #turn genlight into snp matrix
  vcf.r<- read.vcfR(input.vcf) #read in all data
  #convert to genlight
  genlight<- vcfR2genlight(vcf.r)
  #convvert genlight to matrix
  gen.mat<- as.matrix(genlight)

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






