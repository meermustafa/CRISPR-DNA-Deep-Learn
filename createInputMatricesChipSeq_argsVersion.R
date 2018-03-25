# 3.20.18
# Meer Mustafa

# process signal data into 
# e.g. pull out values along some index and turn this into a matrix



### passing in command line arguments ----
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)




# create arrays of x entries in rows and y vectors in columns (vectors are l length (200 bp))

# what is the output size of one signal per x input
desiredInputLength = 1
desiredInputWidth = 200
flanksize = desiredInputWidth/2



# file that provides index
indexFile = read.table('../../../CRISPRiScreen_A549/CRISPRiMYCScreen_chrAndCenterGuidePosition_flank200bp.bed', sep = '\t')


### args -----
# args 1 should be the input bedgraph file
# args 2 should be the output CSV name

# create vector of all files the pull signal from
#vectorOfFiles = args[1]

# read file in
readFile = read.table(args[1], sep='\t')


mat = matrix(data = NA, 
       nrow = (nrow(indexFile) * desiredInputLength),
       ncol = (desiredInputWidth)
       )

# fill out the matrix
for (guideIndex in 1:nrow(indexFile)) {
  mat [guideIndex, ] = readFile$V3[indexFile$V2[guideIndex]] : readFile$V3[indexFile$V3[guideIndex]] # start to end
  
}

# write to a casv
write.table(mat, paste0(args[1],'.',desiredInputWidth,'flanking_pileups.csv',collapse = ""),row.names = F,col.names = F,quote = F,na = '0')  




