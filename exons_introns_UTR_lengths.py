import pandas as pd
import numpy as np
from multiprocessing import Pool


def len_from_fasta(fastafile):
    """function to create a dataframe of seq_name, parent_gene, and length of exons, introns, and UTRs, from FTP fasta files"""

    #opening fasta
    with open(fastafile, "r") as fasta:
        fa_lines = fasta.readlines()

    print("opened", fastafile)

    #list of info lines
    fa_info_lines = [i for i in fa_lines if ">" in i]

    #lists to append to
    seq_name = []
    parent_name = []
    seq_length = []

    print("looping...")

    #looping through info lines
    for line in fa_info_lines:
        line_ls = line.split("; ")

        seq = str(line_ls[0]).replace(">", "")
        seq_name.append(str(seq.split(" ")[0]))

        for i in line_ls:
            
            if "parent" in str(i):
                parents = str(str(i).split('=')[1])
                parent_gene = str(parents.split(",")[0])
                parent_name.append(parent_gene)
            elif "length" in str(i):
                length = int(str(i).split("=")[1])
                seq_length.append(length)
            else:
                continue

    #making dictionary
    dictionary = {"Sequence_name":seq_name, "Gene_name":parent_name, "Length":seq_length}

    #making dataframe
    df = pd.DataFrame(dictionary)

    print("df made for", fastafile)

    #csv filename
    if "exon" in fastafile:
        csv_name = "exon_lengths.csv"
    elif "intron" in fastafile:
        csv_name = "intron_lengths.csv"
    elif "five" in fastafile:
        csv_name = "five_prime_UTR_lengths.csv"
    elif "three" in fastafile:
        csv_name = "three_prime_UTR_lengths.csv"
    else:
        print("something is wrong")


    #saving csv
    df.to_csv(csv_name, index=False)

    print("saved", csv_name)








#multiprocessing
fasta_ls = ['dmel-all-exon-r6.48.fasta', 'dmel-all-five_prime_UTR-r6.48.fasta', 'dmel-all-intron-r6.48.fasta', 'dmel-all-three_prime_UTR-r6.48.fasta']


def run_in_parallel():
    pool = Pool(processes=4)
    pool.map(len_from_fasta, fasta_ls)


if __name__ == '__main__':
    run_in_parallel()

