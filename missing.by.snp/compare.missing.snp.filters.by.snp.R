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


compare.missing.snp.filters.by.snp <- function(gen.mat){
  
  #calculate the proportion of individuals successfully genotyped at each SNP
  miss<-colSums(is.na(gen.mat) == TRUE)/nrow(gen.mat)

  #loop that stores a vector of # non-missing loci retained for each individual
  #looped over all possible completeness filter values
  #initialize df.x
  df.x<- data.frame(filt=numeric(), missingness=numeric())
  #loop
  for (i in c(0,.1,.2,.3,.4,.5,.6,.65,.7,.75,.8,.85,.9,.95,1)){
    #subset the matrix to only snps passing a given filter level
    missingness<-miss[1-miss >= i]
    #calculate the completeness cutoff for each snp to be retained
    filt<-rep(i, times= length(missingness))
    df.x<-rbind(df.x, as.data.frame(cbind(filt, missingness)))
  }
  
  #make columns numeric for plotting
  df.x$filt<-as.numeric(as.character(df.x$filt))
  df.x$missingness<-as.numeric(as.character(df.x$missingness))
  
  #visualize as a histogram for each filter level
  library(ggridges)
  print(
    ggplot(df.x, aes(x = missingness, y = as.character(round(filt, digits = 2)), fill = filt, color = filt)) +
      geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .15, cex=.5) +
      theme_classic() +
      labs(x = "missing data proportion in retained SNPs", y = "% missing data floor") +
      theme(legend.position = "none")
  )
  
}


#run filter function
compare.missing.snp.filters.by.snp(gen.mat=gen.mat)






