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


#write a function that will make a heatmap of depth of sequencing, with sample on y axis, snp on x axis
#going to lean heavily on the tutorial from Knaus:
#https://knausb.github.io/vcfR_documentation/sequence_coverage.html


#This only needs to be done once at the beginning of the script to avoid repeatedly reading in a large vcf
#read in vcf as vcfR
vcfR <- read.vcfR("~/Downloads/todi.allsnps.vcf")

### convert vcfR to genlight
genlight <- vcfR2genlight(vcfR)


vis.depth.and.missing <- function(vcfR, popmap=NULL){
  
      #extract depth from the vcf
      dp<- extract.gt(vcfR, element='DP', as.numeric=TRUE)
      
            if (is.null(popmap)) {
                print("No popmap provided")
                } else {
                #calculate missingness by pop here and make dotplot
                  #popmap must be a two column dataframe with 'id' and 'pop' columns
                  #id's must match the ids in the vcf file
                  #calculate missingness by individual
                  miss<-colSums(is.na(dp))/nrow(dp)
                  #calculate avg depth by individual
                  avg.depth<-colMeans(dp, na.rm = TRUE)
                  #store ordered column names
                  samples<-colnames(dp)
                  #create df
                  df.x<-data.frame(id=samples,missingness=miss,avg.depth=avg.depth, row.names = NULL)
                  #bring in popmap
                  df.f<-merge(df.x, popmap, by="id")
                  #plot missingness and depth by pop
                  plot1<-ggplot(df.f, aes(x= reorder(pop, -missingness), y=missingness)) + 
                    geom_violin()+ theme_classic()+
                    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))+
                    ylab("proportion of SNPs missing data")+
                    geom_dotplot(binaxis='y', stackdir='center', dotsize=.8)
                  #dotplot avg depth of sequencing       
                  plot2<-ggplot(df.f, aes(x= reorder(pop, -missingness), y=avg.depth)) + 
                    geom_violin()+ theme_classic()+
                    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))+
                    ylab("average depth of sequencing for called SNPs")+
                    geom_dotplot(binaxis='y', stackdir='center', dotsize=.8)
                  
                  #
                  grid.arrange(plot1,plot2)
                  
                }
      
    #calculate missingness by individual
    miss<-colSums(is.na(dp))/nrow(dp)
    #calculate avg depth by individual
    avg.depth<-colMeans(dp, na.rm = TRUE)
    #store ordered column names
    samples<-colnames(dp)
    #create df
    df.x<-data.frame(id=samples,missingness=miss,avg.depth=avg.depth, row.names = NULL)
    #dotplot missing data
    plot1<-ggplot(df.x, aes(x= reorder(id, -missingness), y=missingness)) + 
      geom_point()+ theme_light()+
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))+
      ylab("proportion of SNPs missing data")
    #dotplot avg depth of sequencing       
    plot2<-ggplot(df.x, aes(x= reorder(id, -missingness), y=avg.depth)) + 
      geom_point()+ theme_light()+
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))+
      ylab("average depth of sequencing for called SNPs")
    #
    grid.arrange(plot1,plot2)
  #return the dataframe showing the missingness and avg depth per individual
  return(df.x)
}


#run filter function
df.x<-vis.depth.and.missing(vcfR=vcfR, popmap = popmap)






