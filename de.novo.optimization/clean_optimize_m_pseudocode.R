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



#read in m 1-7

#give the function (m1= _, m2= _, m3= _, m4= _, m5= _, m6= _, m7 =_)

#open function:
#if(missing(m1)){}
#else{
##read in vcfR
###calc avg depth of each individual
###rep m1, times = length of avg depth vector
###cbind these vectors into depth df
##convert vcfR to genlight
##convert m1 genlight to matrix
##rep m1, times = # of rows in m1
##calculate num snps retained at each filt level
##calculate num poly loci retained at each filt level
##cbind these three vectors into m df
#}

#if(missing(m2)){}
#else{
##read in vcfR
###calc avg depth of each individual
###rep m2, times = length of avg depth vector
###cbind these vectors into depth df
##convert vcfR to genlight
##convert m2 genlight to matrix
##rep m2, times = # of rows in m2
##calculate num snps retained at each filt level
##calculate num poly loci retained at each filt level
##cbind these three vectors into m df
#}

#...




#take depth df output from all of these possibilities
#plot hist of depth at each m value on same plot

#take m df output from all these possibilities
#plot number of SNPs retained colored by m at each filt level, as open circles
#plot the number of polymorphic loci retained colored by m at each filt level, as closed circles
#plot a vertical line at x=.8
#return a message as output telling the user to pick the m value with the most polyloci retained at r80 (.8)


#close function


#do I want to plot depth of sequencing at each m too?
#this takes vcfr, calcs total reads in each individual
colSums(extract.gt(todi.vcf, element='DP', as.numeric=TRUE), na.rm = T)
#now calc number of GTs successfully called in each individual
colSums(is.na(extract.gt(todi.vcf, element='DP', as.numeric=TRUE)) == "FALSE")
#divide total reads by total GTs to get a vector of avg read depth in each individual
(colSums(extract.gt(todi.vcf, element='DP', as.numeric=TRUE), na.rm = T)) / (colSums(is.na(extract.gt(todi.vcf, element='DP', as.numeric=TRUE)) == "FALSE"))





