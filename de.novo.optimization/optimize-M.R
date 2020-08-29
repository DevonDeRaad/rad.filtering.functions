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


#open function:
optimize_M <- function(M1=NULL,M2=NULL,M3=NULL,M4=NULL,M5=NULL,M6=NULL,M7=NULL){
  #initialize empty depth.df
  depth.df<- data.frame(M=character(), avg.depth=numeric())
  #initialize empty M.df
  M.df<- data.frame(M=character(), filt=as.numeric(), snps=numeric(), V4=character())
  #set vector of m identifiers
  Ms<-c("M1","M2","M3","M4","M5","M6","M7")
  #start on first position in vector of m identifiers
  j=1
  
  #open for loop for each m identifier
  for(x in list(M1,M2,M3,M4,M5,M6,M7)){
    #open if else statement, if no m of given value, move j up to next m identifier, else calculate snps/loci retained
    if(is.null(x)){j=j+1} else{
      #calculate depth first
      ##read in vcfR
      vcf.r<- read.vcfR(x) #read in all data
      ###calc avg depth of each individual
      dep<- (colSums(extract.gt(vcf.r, element='DP', as.numeric=TRUE), na.rm = T)) / (colSums(is.na(extract.gt(vcf.r, element='DP', as.numeric=TRUE)) == "FALSE"))
      ###rep m identifier, times = number of samples in the vcf
      M<- rep(Ms[j], times = length(dep))
      #set j for the next m identifier for next time we go through this loop
      j=j+1
      ###cbind depth and m identifier into depth df
      depth.df<- rbind(depth.df, as.data.frame(cbind(M,dep)))
      
      #start calculating polymorphic loci and snps retained at all possible filter levels for this m identifier
      ##convert vcfR to genlight
      genlight<- vcfR2genlight(vcf.r)
      ##convert genlight to matrix
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
        #close for loop
      }
      ##cbind these three vectors with m identifier and append it all to the m df
      snpsubset<-as.data.frame(cbind(M, filt, snps, rep("snp", times = nrow(gen.mat))))
      locisubset<-as.data.frame(cbind(M, filt, poly.loci, rep("loci", times = nrow(gen.mat))))
      #match colnames so you can rbind these together in tidy format
      colnames(locisubset)[3]<-"snps"
      #append to existing df
      M.df<- rbind(M.df, as.data.frame(rbind(snpsubset,locisubset)))
      #close if else statement
    }
    #close for loop
  }
  
  
  
  #take depth df output from all of these possibilities
  #plot hist of depth at each m value on same plot
  depth.df$M<-as.factor(depth.df$M)
  depth.df$dep<-as.numeric(depth.df$dep)
  print(
    ggplot(depth.df, aes(x = dep, y = M, fill = M, color = M)) +
      geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .35, cex=.5) +
      theme_classic() +
      labs(x = "average depth in each sample", y = "M value (tolerated mismatches in a stack)") +
      theme(legend.position = "none")
  )
  #take m df output from all these possibilities
  #plot number of SNPs retained colored by m at each filt level, as open circles
  #plot the number of polymorphic loci retained colored by m at each filt level, as closed circles
  #plot a vertical line at x=.8
  #return a message as output telling the user to pick the m value with the most polyloci retained at r80 (.8)
  #rename columns
  colnames(M.df)<-c("M","filt","retained","snp.locus")
  m.df$M<-as.character(m.df$M)
  m.df$filt<-as.numeric(m.df$filt)
  m.df$retained<-as.numeric(m.df$retained)
  m.df$snp.locus<-as.character(m.df$snp.locus)
  print(
    ggplot(M.df, aes(x=filt, y=retained, col = M, shape=snp.locus))+
      geom_point()+
      ggtitle("total SNPs and polymorphic loci retained by filtering scheme") +
      xlab("fraction of non-missing genotypes required to retain each SNP (0-1)") + ylab("# SNPs/loci")+
      theme_light()+
      geom_vline(xintercept=.8)+
      labs(col = "mismatches allowed", shape="")
  )
  
  #return the depth and snp/loci dataframes in case you want to do your own visualizations
  out <- list()
  out$depth<-depth.df
  out$m.comparisons<-M.df
  return(out)
  
  print("Check out how different values of M affect depth of coverage")
  print("Optimal M value returns the most polymorphic loci in the 80% complete matrix (Paris et al. 2017)")
  
  #close function
}
