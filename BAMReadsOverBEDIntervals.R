# Meer Mustafa
# 11.1.17



# script/function to compute BAM reads over a set of window intervals (commonly provided by BED file)
# usage: Rscript Rscriptfilename.R


# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite('Biostrings')
library(Biostrings)
#biocLite('Rsamtools')
library(Rsamtools)
#biocLite('GenomicRanges')
library(GenomicRanges)


# function that converts any file (BED or otherwise columnar file) into a GRange object (RSamtools needs a GRange interval in order to work)
bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  library("GenomicRanges")
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}



# input bed interval file to compute bam reads over intervals, convert to Grange object
gr = bed_to_granges('FusionsOncogenesAndOncogeneFusionGenes.bed')

# if bed files are centered at the same site (e.g. chr1 100 100), extend this interval by n bp
flankedgr = flank(gr, width = 1000000, both = T)
head(flankedgr)



# turn the below for loop into a function to change the bin sizes
getBAMReadsWithinBinnedEnhancers = function (inputBEDfile, enhancerBinSizebp) { # topNEnhancersToSelect
  
  # load in the BED file - 
  #require(genomation)
  
  #read in the cancer-wide compiled enhancer BED
  # fullBEDGRange = readGeneric(inputBEDfile,
  #                             #meta.cols = list(score = 5),
  #                             header = F, keep.all.metadata = T)
  # cat('successfully read in BED file and converted to GRanges object')
  # cat('The number of duplicated ranges in the GRange is:',length(which(!duplicated(fullBEDGRange))),'\n')
  # 
  # # subset GRange only to those that are not duplicated ranges
  # fullBEDGRange = fullBEDGRange[ which(!duplicated(fullBEDGRange)) , ]
  # 
  # # # subset enhancer to top N enhancers, 
  # # # likely will change this to include equal representation of each cancer's enhancers
  # # n = topNEnhancersToSelect
  # top_fullBEDGRange = fullBEDGRange[c(200:400,600:800,1000:1600,1800:2000,2200:2800)] # pick the A549 enhancers
  # cat('successfully subsetted the GRange to the top',topNEnhancersToSelect,'enhancers\n')
  
  # make sure the enhancers have a uniform width
  # if ( length(which(!duplicated(width(top_fullBEDGRange)))) != 1 ) { 
  #   print('hello')
  #   
  #   # if they don't set it equal to something
  #   width(top_fullBEDGRange) = 1001
  #   hist(width(top_fullBEDGRange), breaks = c(100))
  # }
  
  
  # BED to GRANGE version ------
  top_fullBEDGRange = inputBEDfile
  
  
  # loop through the bam files. make sure each .bam has a .bai index
  for (BAM in list.files(pattern = ".bam$")) { #list.files(path = '../',pattern = '.bam$')
    
    cat('Now working on:',BAM,'\n')
    
    # assign bam file name
    bf = BamFile(file = BAM) 
    
    # initialise matrix to hold all enhancers
    currentCancerMatrix = matrix(
      nrow = length(top_fullBEDGRange), # number of enhancers
      ncol = length(seq(start(top_fullBEDGRange[1,]), end(top_fullBEDGRange[1,]) , enhancerBinSizebp)) - 1 ) # number of bins
    
    cat('The dimensions of the current cancers df is:',dim(currentCancerMatrix),'\n')
    
    # loop through the enhancers in the BED file 
    for (BEDRowIndex in 1:length(top_fullBEDGRange)) {
      
      # set chr
      binPositionChr = as.character(seqnames(top_fullBEDGRange[BEDRowIndex, ]))
      
      # set the edge positions of the bins in the enhancer
      binPositionEdges = seq(start(top_fullBEDGRange[BEDRowIndex, ]), end(top_fullBEDGRange[BEDRowIndex, ]) , enhancerBinSizebp)
      
      # set the singular enhancer region as a variable to be used later in the ScanBamParam
      singularEnhancer = top_fullBEDGRange[BEDRowIndex, ]
      #cat('This is the singular enhancer that is about to be analysed:') ; print(singularEnhancer)
      
      # initialise an empty vector to hold the MAX read pileup
      #### NAME THIS E1, E2, ETC... IN PARALLEL WITH THE ROW NUMBER BEING PROCESSED
      binPileupForCurrentEnhancer = c()
      
      # segment enhancers and compute pileups in bins of the enhancers
      # use length of loop to stop at length(bins - 1)
      for (binIndex in 1: (length( seq(start(singularEnhancer), end(singularEnhancer) , enhancerBinSizebp)) - 1))  {
        
        #cat('This is the bin index that pileup will be computed for:',binIndex,'\n')
        
        # use the i and i + 1 bin. Starting with bin 1 and 2. Stop at length(bins - 1)
        # set this GRange as the ScanBamParam for RSamTools to find the pileups per bp in this genomic range
        param = ScanBamParam(which = GRanges(seqnames = binPositionChr, IRanges(start = binPositionEdges[binIndex], end = binPositionEdges[binIndex + 1])))
        
        # compute pileup in the current GRange
        binPileup = pileup(bf, scanBamParam = param)
        
        
        # Take the max pileup in the binned genomic region
        maxCountInBin = max(binPileup$count)
        #maxCountInBin = round(max (binPileup$count), digits = 1)
        cat('the max read count in the bin is:',maxCountInBin,'\n')
        
        # if the maxCountInBin is a -Inf or a non-numeric value then convert it to 0
        if (maxCountInBin == -Inf) {
          cat('Found a Inf value, replacing it with a 0 read count!\n')
          maxCountInBin = 0
          cat('the NEW max read count in the bin is:',maxCountInBin,'\n')
        }
        
        # check for Not a Number cases, convert to 0
        if (is.nan(maxCountInBin) == TRUE) {
          cat('Found a NaN value, replacing it with a 0 read count!\n')
          maxCountInBin = 0
          cat('the NEW max read count in the bin is:',maxCountInBin,'\n')
        }
        
        
        # append this bin's pileup to the running vector of max read count pileups for this enhancer
        # add a value of the next binned pilup to the vector
        binPileupForCurrentEnhancer = append(binPileupForCurrentEnhancer, maxCountInBin)
        
      } # end of bins
      
      cat('about the append this enhancers bins to the full matrix \n')
      cat('the size of the binned',BEDRowIndex,'enhancer vector is',length(binPileupForCurrentEnhancer),'and the # of columns of the matrix is', ncol(currentCancerMatrix),'\n')
      
      # check to make sure that the length of bin values should == # of cols of the matrix
      if (length(binPileupForCurrentEnhancer == ncol(currentCancerMatrix))) {
        
        # append the most recent enhancer's binned pileups () to the full cancer matrix
        currentCancerMatrix[BEDRowIndex, ] = binPileupForCurrentEnhancer
        cat('The new dimensions of the cancers matrix after adding the',BEDRowIndex,'enhancer is:',dim(currentCancerMatrix),'\n')
      }
      
    } # end of all enhancers
    
    cat('About to compute library size of',BAM)
    
    # normlize the read count to the library size
    CurrentBAMLibrarySizeInMillionsOfReads =  ( as.numeric(system(paste("samtools view", 
                                                                        BAM, 
                                                                        " | wc -l",sep=" "), 
                                                                  intern = T)) 
                                                / 1000000 ) # divide by 1e6 to get Reads Per Million
    cat('The RPM value for',BAM,'is',print(CurrentBAMLibrarySizeInMillionsOfReads),'\n')
    
    # divide the matrix by this RPM value
    currentCancerMatrix = currentCancerMatrix / CurrentBAMLibrarySizeInMillionsOfReads
    cat('The dimensions of the current cancers matrix after normalizing to library size is:',dim(currentCancerMatrix),'\n')
    # after looping through all the enhancer GRanges change the name of the currentCancerMatrix 
    # to a name that has the BAM cancer cell line name
    assign( paste0( sapply( strsplit( BAM, split='.', fixed=TRUE ), # pull name from BAM file name currently loaded in
                            function(x) (x[1])) # split the bam file name to keep only cell line
                    , 'binnedEnhancerPileups') # add to the name that this is the binnedEnhancerPileups
            , value = currentCancerMatrix # should be empty df on 1st iteration, should have
            , envir = globalenv()) # save to environment outside the function
    cat('successfully assigned matrix to an R object')
    
    # print statement to notify completion of the loop for each BAM
    cat('Successfully computed binned read pileup for',
        paste0( sapply( strsplit( BAM, split='.', fixed=TRUE ), # pull name from BAM file name currently loaded in
                        function(x) (x[1]))),'BAM file. \n')
    
    # troubleshoot the NaNs and Infs 
    save.image(file = paste0('individualCancers_COSMICOncogenesBinnedEnhancerMatrices',Sys.time(),'.RDS', collapse = "") )  
    
  } # end of each BAM file
  
  
} # end of function


# call the function
ptm = proc.time() # calc time it takes to run on system
getBAMReadsWithinBinnedEnhancers(inputBEDfile = flankedgr,
                                 # topNEnhancersToSelect = 50,
                                 enhancerBinSizebp = 1000)
(proc.time() - ptm)[1]/60/60


# save all of the individual DFs for each cancer's read pileup
ls(pattern = '*binnedEnhancer*')
# option to save 1 image
#save(list = ls(pattern = '*binnedEnhancer*'), file = paste0('COSMICOncogenesBinnedEnhancerMatrices.RDS',Sys.Date()))


save.image(file = paste0('COSMICOncogenesBinnedEnhancerMatrices',Sys.time(),'.RDS'))





