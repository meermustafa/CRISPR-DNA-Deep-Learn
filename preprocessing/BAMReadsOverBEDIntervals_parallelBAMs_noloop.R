# Meer Mustafa
# 11.1.17


# script/function to compute BAM reads over a set of window intervals (commonly provided by BED file), and outputs a matrix per BAM that has has as rows each window interval that are divided by columns into binned read pileups 
# usage: Rscript Rscriptfilename.R <InputBEDFileName> <NameOfFinalRDSfile e.g. h3k27ac_oncogenes> <name of BAM file from shell ls loop> <flanking size to expand out the BED file intervals>
# **** requires BAM file named with extension .bam inside the working directory. May also use '.' separators in BAM name for naming purposes
# **** requires BED file inside the working directory


### passing in command line arguments ----
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# print out arguments
print(args[1])
print(args[2])
print(args[3])
print(args[4])

#biocLite('Biostrings')
#biocLite('Rsamtools')
#biocLite('GenomicRanges')
require(Biostrings)
require(Rsamtools)
require(GenomicRanges)


# function that converts any file (BED or otherwise columnar file) into a GRange object (RSamtools needs a GRange interval in order to work)
bed_to_granges = function(file){
  df = read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df = df[,-c(7:length(df))]
  }
  
  if(length(df) < 3){
    stop("File has less than 3 columns")
  }
  
  header = c('chr','start','end','id','score','strand')
  names(df) = header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand = gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  # load GR
  require(GenomicRanges)
  
  # if there are more columns i nthe BED, add as metadata
  if(length(df) == 3){
    gr = with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df )== 4){
    gr = with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df) == 5){
    gr = with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df) == 6){
    gr = with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  
  return(gr)
  
}



# input bed interval file to compute bam reads over intervals, convert to Grange object
gr = bed_to_granges(args[1])

# save the input command line arg as the BAM file name
BAM = args[3]

# save the input command line arg as the flanking size
flankingSize = as.numeric(args[4])

# if bed files are centered at the same site (e.g. chr1 100 100), extend this interval by n bp
flankedgr = flank(gr, width = flankingSize, both = T)
head(flankedgr)


# main function to get BAM read pileups for each window interval that is binned to size n
getBAMReadsWithinBinnedEnhancers = function (inputBEDfile, enhancerBinSizebp, topNEnhancersToSelect = nrow(inputBEDfile)) {
  
  # load in the BED file - 

  # assign GR to name ------
  top_fullBEDGRange = inputBEDfile
  
  # subset GR if desired
  top_fullBEDGRange = top_fullBEDGRange[1:topNEnhancersToSelect]
  
  # make sure each .bam has a .bai index
  system(paste("samtools index ", 
               BAM, 
               sep=" "), 
         intern = T)
  
  
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
  CurrentBAMLibrarySizeInMillionsOfReads =  ( as.numeric(system(paste("samtools view -F 0x04 -c", 
                                                                      BAM, 
                                                                      sep=" "), 
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
  #save.image(file = paste0( args[2], args[3], Sys.time(),'.RDS', collapse = ".") )  
  
 # end of each BAM file
  
} # end of function


# call the function
ptm = proc.time() # calc time it takes to run on system
getBAMReadsWithinBinnedEnhancers(inputBEDfile = flankedgr,
                                 enhancerBinSizebp = 1000,
                                 topNEnhancersToSelect = 500
                                 )
(proc.time() - ptm)[1]/60/60

# save a new name for the output file
BAMName = sapply( strsplit( args[3], split='.', fixed=TRUE ), # pull name from BAM file name currently loaded in
                  function(x) (x[1]))

# save all R objects to dir
save.image(file = paste0(BAMName, args[2], Sys.time(),'.RDS', collapse = ".") )  