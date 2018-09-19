from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import pandas as pd
import numpy as np
from Bio import SeqIO
import csv

# generator function to pull a line at a time
def get_data(filename):
    with open(filename,"r") as csvfile:
        datareader = csv.reader(csvfile)
        count = 0
        for row in datareader:
            if count > 0:
                yield count, row
            count += 1
            if count % 1000 == 0:
                print(count)

flank_size = 100

# master file
crispr_file = "/Users/jason/PycharmProjects/Tsirigos/deeplearning/A549_Rep1_MiSeqRun_finalCRISPRiScores2018-03-01.csv"
crispr_df = pd.read_csv(crispr_file)

label_encoder = LabelEncoder()
label_encoder.fit(["G","C","A","T"])

# chr8 reference genome hg19
seqRec = SeqIO.read(open("chr8.fa"),"fasta")

for count, row in get_data(crispr_file):
# count = 0
    unique_id = row[0]
    start = int(row[2])
    end = int(row[3])
    strand = row[4]
    guide_seq = row[5]
    # generating indexes + flanking region size
    mid_idx = int((start+end)/2)
    s_idx = mid_idx-flank_size
    e_idx = mid_idx+flank_size
    seq = seqRec.seq[s_idx:e_idx].upper()
    # if on neg strand, take reverse complement of sequence
    if strand == "-":
        seq = seq.reverse_complement()
    # data verification
    # note: seq[90:110] returns the guide sequence
    # print(s_idx)
    # print(e_idx)
    print(seq)
    print(seq[90:110])
    print(len(seq))

    # integer encode
    integer_encoded = label_encoder.transform(list(seq))
    # binary encode
    onehot_encoder = OneHotEncoder(sparse=False)
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.fit_transform(integer_encoded)

    output = np.array(onehot_encoded.T)
    # remove comment when running full
    # saving each DNA file as .npy
    # np.save("DNA_arrays/"+unique_id+"_DNA.npy",output)
    print(output)
    # print(output.shape)

    # comment when running in full
    if count == 3:
        break
    count += 1

# first 3 guide sequences (from master file)
# AATGTAGTGTACTCACTGAT
# GCTTGTAGGAGGGAAGAAGG
# GTGCTTGTAGGAGGGAAGAA
