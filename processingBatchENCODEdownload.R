# Meer Mustafa
# 3.6.18

# processing batch downloaded ENCODE files (ENCF_______.extension) using their metadata.tsv file


# read metadata
m = read.table('metadata.tsv',header=T, stringsAsFactors=F,sep='\t')

# criteria ----
# entries must be certain conditions: Assembly == 'hg19', Output.type == 'alignments' (i.e. not unfiltered alignments), 
# ONLY ONE OF EACH Experiment.target (antibody experiment type e.g. H3K27ac) SHOULD BE USED, NO DUPLICATED
# USE Biological.replicate.s. == '1' TO PICK THE FIRST BIOREP; SOME EXPERIMENTS ONLY HAVE 1 BIOREP, OTHERS HAVE > 1
# USE Biosample.treatments == "" TO PICK ONLY SAMPLES WITH NO TREATMENTS E.G. ETHANOL
# File.accession is used to search the listed files names

# which experiments come from unique experiment sets
length(which(!duplicated(m$Experiment.accession)))

m[!duplicated(m$Experiment.target),]$Experiment.target


# this still has duplicate experiment sets (i.e. multiple expts coming from each separate labs)
subset = m[ m$Assembly == 'hg19' &  m$Output.type == 'alignments' & m$Biological.replicate.s. == '1'  & m$Biosample.treatments == "", ]$Experiment.target

# this removes too many samples
m[ m$Assembly == 'hg19' &  m$Output.type == 'alignments' & m$Biological.replicate.s. == '1'& !duplicated(m$Experiment.accession) , ]$Experiment.target

which(m$Assembly == 'hg19' &  m$Output.type == 'alignments' & m$Biological.replicate.s. == '1'& !duplicated(m$Experiment.accession) )



# Experiment.target, File.accession, and extension will name the ultimate file

# now need to take the File.accession (file name given by ENCODE and rename them) and split them to get just the file name and no extension
listedFiles = sapply(strsplit(list.files(pattern = '.bam$'), split='.', fixed=TRUE), function(x) (x[1]))

# which desired files are in listed files -- should be all of them
subset$File.accession [subset$File.accession %in% listedFiles]

# now move all of these files into a subdir called subset
# make new subdir
system(command = 'mkdir subset', intern = F)

# add the file extension name back to the desired subset files, pick the 1st one as all should be the same
subsetExtensionName = sapply(strsplit(list.files(pattern = '.bam$'), split='.', fixed=TRUE), function(x) (x[2]))[1]
subsetExtensionName

# LOOP THROUG THE SUBSET ROOT FILE NAMES
for (subsetname in subset$File.accession) {
  print(subsetname)
  
  # concatenate the subset file root name and the extension name -- ENCF_____.bam
  fullRawFilename = paste0(subsetname, subsetExtensionName, collapse = "")
  
  # or concatenate the Experiment.target (H3K27ac), Biosample.term.name i.e. cell line (A549), subset root name (ENCF____), the extension name (.bam)
  # pull the row in the df where the entry is located in order to get names of other
  idx = which(subset$File.accession == subsetname)
  
  # create the new file name based on this new string
  NewFullRawFilename = paste0(subset$Experiment.target[idx] , subset$Biosample.term.name[idx], subsetname, subsetExtensionName, collapse = "")
  
  # rename a copy of this file as this string using cp
  system(command = sprintf('cp %s %s', fullRawFilename, NewFullRawFilename), intern = F)
  
  # then mv this file name into the newly create subdir
  system(command = sprintf('mv %s subset/', NewFullRawFilename), intern = F)
  
  
}



# now run analysis on them
# can run multiple qsub/sbatch

# analysis #1 - compute per normalized read counts per bp site (e.g. 1,2,3,4,5,5,5,4,3,2,1 for position x:y)




