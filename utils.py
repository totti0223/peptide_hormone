from Bio import SeqIO
import csv
import pandas as pd

def read_fasta(path):
    fasta = SeqIO.parse(path, 'fasta')
    return fasta

# a function to convert gff3 format file to pandas dataframe

# convert gff3 format file to dataframe
def gff3_to_df(gff3_path):
    data = []
    with open(gff3_path, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue
            l = l.strip().split('\t')
            data.append(l)
    df = pd.DataFrame(data)

    return df
