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
  
  #loop that stores a vector of # non-missing loci retained for each individual
  #looped over all possible completeness filter values
  
  #plot a histogram of each distribution and label it by filtering scheme
  #initialize df.x
  df.x<- data.frame(indiv=character(), loci.retained=numeric(), filt=numeric())
  #loop
  for (i in nrow(gen.mat):1){
    #calculate the number of snps successfully genotyped in each sample
    loci.retained<-rowSums(is.na(gen.mat[,(colSums(is.na(gen.mat)) < i)]) == "FALSE")
    #calc num of polymorphic loci retained at given filt level
    #calculate the completeness cutoff for each snp to be retained
    filt<-rep(1-((i-1)/nrow(gen.mat)), times =nrow(gen.mat))
    indiv<-rownames(gen.mat)
    df.x<-rbind(df.x, as.data.frame(cbind(indiv, filt, loci.retained)))
    
  }
  #make columns numeric for plotting
  df.x$filt<-as.numeric(df.x$filt)
  df.x$loci.retained<-as.numeric(df.x$loci.retained)
  #visualize color coded by individual
  print( 
    ggplot(df.x, aes(x=filt, y=loci.retained, color=indiv))+
      geom_point()+
      ggtitle("SNPs retained by filtering scheme") +
      xlab("fraction of non-missing genotypes required to retain each SNP (0-1)") + ylab("non-missing SNP genotypes")+
      theme_light()
  )
  #print average missing data per individual so that we can figure out which individuals are pulling distribution down
  print(
    
  )
  
  
  #visualize as a histogram for each filter level
  library(ggridges)
  print(
    ggplot(df.x, aes(x = loci.retained, y = as.character(round(filt, digits = 2)), fill = filt, color = filt)) +
      geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .15, cex=.5) +
      theme_classic() +
      labs(x = "# non-missing SNP genotypes", y = "lower completeness limit to retain SNP") +
      theme(legend.position = "none")
  )
  
  #calculate variance btwn samples within each filter level
  snp.var<-vector("numeric",length = length(levels(as.factor(df.x$filt))))
  snp.range<-vector("numeric",length = length(levels(as.factor(df.x$filt))))
  filtr<-levels(as.factor(df.x$filt))
  n=1
  for (i in levels(as.factor(df.x$filt))){
    snp.var[n]<-var(df.x$loci.retained[df.x$filt == i])
    snp.range[n]<-diff(range(df.x$loci.retained[df.x$filt == i]))
    n=n+1
  }
  
  #calculate number of loci retained at each filter level
  for (i in nrow(gen.mat):1){
    #calculate the number of loci retained
    loci.retained[i]<-ncol(gen.mat[,(colSums(is.na(gen.mat)) < i)])
    poly.loci.retained[i]<-length(unique(str_extract(colnames(gen.mat[,(colSums(is.na(gen.mat)) < i)]), pattern = "[0-9]+")))
    #calculate the completeness cutoff for each snp to be retained
    filt[i]<-1-((i-1)/nrow(gen.mat))
  }
  
  #plot these side by side
  plot(x=filtr,y=snp.var/1000, col="red")
  points(x=filt, y = loci.retained)
  points(x=filt, y = poly.loci.retained, col="green")
  points(x=filtr, y=snp.range, col="blue")
  
  
  
  
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



#
snps<-colnames(gen.mat)
str_replace(snps, ":[0-9]{1,2}:[+-]", "")
#calc number of polymorphic loci
length(unique(str_replace(colnames(gen.mat[,(colSums(is.na(gen.mat)) < 5)]), ":[0-9]{1,2}:[+-]", "")))


str_replace(colnames(gen.mat[,(colSums(is.na(gen.mat)) < 5)]), ":[0-9]{1,2}:[+-]", "")

str_extract(colnames(gen.mat[,(colSums(is.na(gen.mat)) < 5)]), pattern = "[0-9]+")

colnames(gen.mat[,(colSums(is.na(gen.mat)) < 5)])

