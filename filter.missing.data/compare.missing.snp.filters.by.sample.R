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
#turn genlight into snp matrix
gen.mat<-as.matrix(genlight)


compare.missing.snp.filters.by.sample <- function(gen.mat){

  #calculate the proportion of individuals successfully genotyped at each SNP
  miss<-colSums(is.na(gen.mat) == TRUE)/nrow(gen.mat)
  
  #loop that stores a vector of # non-missing loci retained for each individual
  #looped over all possible completeness filter values
  #initialize df.x
  df.x<- data.frame(sample=character(), filt=numeric(), missingness=numeric())
  #loop
  for (i in c(0,.1,.2,.3,.4,.5,.6,.65,.7,.75,.8,.85,.9,.95,1)){
    #subset the matrix to only snps passing a given filter level
    gen.z.mat<-gen.mat[,1-miss >= i]
    #calculate the missing % in each sample for the given cutoff
    missingness<-rowSums(is.na(gen.z.mat))/ncol(gen.z.mat)
    #calculate the completeness cutoff for each snp to be retained
    filt<-rep(i, times= length(missingness))
    #list samples
    sample<-rownames(gen.z.mat)
    df.x<-rbind(df.x, as.data.frame(cbind(sample, filt, missingness)))
  }
  
  #make columns numeric for plotting
  df.x$filt<-as.numeric(as.character(df.x$filt))
  df.x$missingness<-as.numeric(as.character(df.x$missingness))
  rownames(df.x)<-NULL
  
  #visualize as a histogram for each filter level
  library(ggridges)
  print(
    ggplot(df.x, aes(x = missingness, y = as.character(round(filt, digits = 2)), fill = filt, color = filt)) +
      geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .15, cex=.5) +
      theme_classic() +
      labs(x = "missing proportion of genotypes in each sample", y = "% missing data floor") +
      theme(legend.position = "none")
  )
  
  #plot colored by sample
  print( 
    ggplot(df.x, aes(x=as.factor(filt), y=missingness, color=sample))+
      geom_point()+
      xlab("fraction of non-missing genotypes required to retain each SNP (0-1)")+
      ylab("proportion of SNPs missing in each sample")+
      theme_light()
  )
  
  return(df.x)
}


#run filter function
compare.missing.snp.filters.by.sample(gen.mat = gen.mat)





