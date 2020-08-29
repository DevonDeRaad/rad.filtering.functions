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
library(ggridges)


#function requires snps to be named as "uniquelocusname*any non [0-9] character*snpname" in order for poly loci calculation to work
#give the function full path to m 1-7 

#eg function(m1="~/Downloads/todi.allsnps.vcf", m2= _, m3= _, m4= _, m5= _, m6= _, m7 =_)

#open function:
optimize_m <- function(m1,m2,m3,m4,m5,m6,m7){

#initialize depth.df
depth.df<- data.frame(m=character(), avg.depth=numeric())
#initialize m.df
m.df<- data.frame(m=character(), snps.retained=numeric(), poly.loci.retained=numeric())

if(missing(m1)){} else{
##read in vcfR
vcf.r<- read.vcfR(m1) #read in all data
###calc avg depth of each individual
dep<- (colSums(extract.gt(vcf.r, element='DP', as.numeric=TRUE), na.rm = T)) / (colSums(is.na(extract.gt(vcf.r, element='DP', as.numeric=TRUE)) == "FALSE"))
###rep m1, times = length of avg depth vector
m<- rep("m1", times = length(dep))
###cbind these vectors into depth df
depth.df<- rbind(depth.df, as.data.frame(cbind(m,dep)))

##convert vcfR to genlight
genlight<- vcfR2genlight(vcf.r)
##convert m1 genlight to matrix
gen.mat<- as.matrix(genlight)
#initialize vectors to hold filt level, snps retained, poly loci retained
filt<- vector("numeric", length = nrow(gen.mat))
snps<- vector("numeric", length = nrow(gen.mat))
poly.loci<- vector("numeric", length = nrow(gen.mat))
##run loop to fill up vectors with a value for each filter level
for (i in nrow(gen.mat):1){
  #calculate the completeness cutoff for each snp to be retained
  filt[i]<-1-((i-1)/nrow(gen.mat))
  #calculate the number of snps retained at this cutoff
  snps[i]<-ncol(gen.mat[,(colSums(is.na(gen.mat)) < i)])
  #calculate number of polymorphic loci retained at this cutoff
  poly.loci[i]<-length(unique(str_extract(colnames(gen.mat[,(colSums(is.na(gen.mat)) < i)]), pattern = "[0-9]+")))
}
##cbind these three vectors with m identifier and append it all to the m df
m.df<- rbind(m.df, as.data.frame(cbind(m,filt,snps,poly.loci)))
}



#initialize depth.df
depth.df<- data.frame(m=character(), avg.depth=numeric())
#initialize m.df
m.df<- data.frame(m=character(), snps.retained=numeric(), poly.loci.retained=numeric())
#set m identifiers
ms<-c(m1,m2,m3,m4,m5,m6,m7)
j=1
#for loop this entire thing
for(x in c(m1,m2,m3,m4,m5,m6,m7)){
if(missing(x)){j=j+1} else{
  ##read in vcfR
  vcf.r<- read.vcfR(x) #read in all data
  ###calc avg depth of each individual
  dep<- (colSums(extract.gt(vcf.r, element='DP', as.numeric=TRUE), na.rm = T)) / (colSums(is.na(extract.gt(vcf.r, element='DP', as.numeric=TRUE)) == "FALSE"))
  ###rep m identifier, times = number of samples in the vcf
  m<- rep(ms[j], times = length(dep))
  #set j for the next m vlaue
  j=j+1
  ###cbind these vectors into depth df
  depth.df<- rbind(depth.df, as.data.frame(cbind(m,dep)))
  
  ##convert vcfR to genlight
  genlight<- vcfR2genlight(vcf.r)
  ##convert m1 genlight to matrix
  gen.mat<- as.matrix(genlight)
  #initialize vectors to hold filt level, snps retained, poly loci retained
  filt<- vector("numeric", length = nrow(gen.mat))
  snps<- vector("numeric", length = nrow(gen.mat))
  poly.loci<- vector("numeric", length = nrow(gen.mat))
  ##run loop to fill up vectors with a value for each filter level
  for (i in nrow(gen.mat):1){
    #calculate the completeness cutoff for each snp to be retained
    filt[i]<-1-((i-1)/nrow(gen.mat))
    #calculate the number of snps retained at this cutoff
    snps[i]<-ncol(gen.mat[,(colSums(is.na(gen.mat)) < i)])
    #calculate number of polymorphic loci retained at this cutoff
    poly.loci[i]<-length(unique(str_extract(colnames(gen.mat[,(colSums(is.na(gen.mat)) < i)]), pattern = "[0-9]+")))
  }
  ##cbind these three vectors with m identifier and append it all to the m df
  m.df<- rbind(m.df, as.data.frame(cbind(m,filt,snps,poly.loci)))
}

}








#take depth df output from all of these possibilities
#plot hist of depth at each m value on same plot

#take m df output from all these possibilities
#plot number of SNPs retained colored by m at each filt level, as open circles
#plot the number of polymorphic loci retained colored by m at each filt level, as closed circles
#plot a vertical line at x=.8
#return a message as output telling the user to pick the m value with the most polyloci retained at r80 (.8)


#close function
}


#do I want to plot depth of sequencing at each m too?
#this takes vcfr, calcs total reads in each individual
colSums(extract.gt(todi.vcf, element='DP', as.numeric=TRUE), na.rm = T)
#now calc number of GTs successfully called in each individual
colSums(is.na(extract.gt(todi.vcf, element='DP', as.numeric=TRUE)) == "FALSE")
#divide total reads by total GTs to get a vector of avg read depth in each individual
(colSums(extract.gt(todi.vcf, element='DP', as.numeric=TRUE), na.rm = T)) / (colSums(is.na(extract.gt(todi.vcf, element='DP', as.numeric=TRUE)) == "FALSE"))





