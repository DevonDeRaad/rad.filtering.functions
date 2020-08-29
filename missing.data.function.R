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

#function takes a vcfobject as input
#turns genlight into matrix
#uses for loop to calulate:
#(number of SNPs retained at each possible filtering level and store in a vector)
#(corresponding filtration percentage)
#(%missing at each SNP)
#sticks those three vectors together with cbind
#function returns you this dataframe to make your own plots from if you wish
#and autoprints to your plot window ggplot stacked histograms of (x = missing %, y = # of loci retained)
#fin

#This only needs to be done once at the beginning of the script to avoid repeatedly reading in a large vcf
#read in vcf as vcfR
vcfR <- read.vcfR(x)

### convert vcfR to genlight
genlight <- vcfR2genlight(vcfR)


filter.completeness <- function(x){

  #turn genlight into snp matrix
  gen.mat<-as.matrix(genlight)
  
  #calculate the proportion of individuals successfully genotyped at each SNP
  completeness<-colSums(is.na(gen.mat) == FALSE)/nrow(gen.mat)
  #list snp names
  snp<-colnames(gen.mat)
  #create dataframe of missing data by snp
  df.x<-data.frame(snp=snp, completeness=completeness)
  
  
  
  #loop that stores a vector of # non-missing loci retained for each individual
  #looped over all possible completeness filter values
  #initialize df.x
  df.x<- data.frame(indiv=character(), loci.retained=numeric(), filt=numeric())
  #loop
  for (i in c(0,.1,.2,.3,.4,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1)){
    #calculate the number of individuals successfully genotyped at each SNP
    loci.retained<-colSums(is.na(gen.mat[,(colSums(is.na(gen.mat)) >= i*nrow(gen.mat))]) == "FALSE")
    #calculate the completeness cutoff for each snp to be retained
    filt<-rep(i, times =ncol(gen.mat))
    snp<-colnames(gen.mat)
    df.x<-rbind(df.x, as.data.frame(cbind(snp, filt, loci.retained)))
    
  }
  #make columns numeric for plotting
  df.x$filt<-as.numeric(as.character(df.x$filt))
  df.x$loci.retained<-as.numeric(as.character(df.x$loci.retained))

  #visualize as a histogram for each filter level
  library(ggridges)
  print(
    ggplot(df.x, aes(x = loci.retained, y = as.character(round(filt, digits = 2)), fill = filt, color = filt)) +
      geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .15, cex=.5) +
      theme_classic() +
      labs(x = "# non-missing SNP genotypes retained", y = "% missing data floor") +
      theme(legend.position = "none")
  )
  
  #return plot of df and df itself
  print( 
    ggplot(df.x, aes(x=as.factor(filt), y=loci.retained))+
      geom_violin()+
      ggtitle("SNPs retained by filtering scheme") +
      xlab("fraction of non-missing genotypes required to retain each SNP (0-1)") + ylab("total SNPs retained")+
      theme_light()
  )
  return(df.x)
}


#run filter function
df.x<-filter.completeness("~/Downloads/todi.allsnps.vcf")






