import pandas as pd
import numpy as np



def site_counter(seq, length, codons):
    """Function to loop thoruhg a CDS codon by codon and calculate the number of potential nonsynoymous sites if there is one mutation per codon"""
    
    #lists for summing
    SS = []
    NS = []

    #looping through sequence
    for i in range(int(length/3)):
        codon = seq[(0+(3*i)):(3+(3*i))]

        SS.append(3-(codons[codon]))
        NS.append(codons[codon])

    return sum(SS), sum(NS)




codon_NS_conversion = {
    "TTT":2.5, "TCT":2, "TAT":2.5, "TGT":2.5,
    "TTC":2.5, "TCC":2, "TAC":2.5, "TGC":2.5,
    "TTA":2, "TCA":2, "TAA":2, "TGA":2.5,
    "TTG":2, "TCG":2, "TAG":2.5, "TGG":3,

    "CTT":2, "CCT":2, "CAT":2.5, "CGT":2,
    "CTC":2, "CCC":2, "CAC":2.5, "CGC":2,
    "CTA":1.5, "CCA":2, "CAA":2.5, "CGA":1.5,
    "CTG":1.5, "CCG":2, "CAG":2.5, "CGG":1.5,

    "ATT":2.25, "ACT":2, "AAT":2.5, "AGT":2.5,
    "ATC":2.25, "ACC":2, "AAC":2.5, "AGC":2.5,
    "ATA":2.25, "ACA":2, "AAA":2.5, "AGA":2,
    "ATG":3, "ACG":2, "AAG":2.5, "AGG":2,

    "GTT":2, "GCT":2, "GAT":2.5, "GGT":2,
    "GTC":2, "GCC":2, "GAC":2.5, "GGC":2,
    "GTA":2, "GCA":2, "GAA":2.5, "GGA":2,
    "GTG":2, "GCG":2, "GAG":2.5, "GGG":2
}


#opneing fasta
with open("dmel-all-CDS-r6.48.fasta", 'r') as fastafile:
    fasta = fastafile.read()

print("opened fasta file")

#fasta list
fasta_ls = fasta.split(">")
del fasta_ls[0]


#lists for df
gene = []
length_CDS = []
num_NS = []
num_SS = []


print("looping...")

#looping through fasta
for record in fasta_ls:

    line_ls = record.split("; ")

    for i in line_ls:
        
        if "parent" in str(i):
            parents = str(str(i).split('=')[1])
            parent_gene = str(parents.split(",")[0])
            gene.append(parent_gene)
        elif "length" in str(i):
            CDS_length = float(str(i).split("=")[1])
            length_CDS.append(CDS_length)
        elif "=" in str(i):
            continue
        else:
            sequence = str(i).replace('\n', '')

            sumSS, sumNS = site_counter(sequence, CDS_length, codon_NS_conversion)

            num_NS.append(sumNS)
            num_SS.append(sumSS)


print("all sites counted")


#making df
dictionary = {"FBgn":gene, "Length_of_CDS":length_CDS, "Number_NS_sites":num_NS, "Number_S_sites":num_SS}

df = pd.DataFrame(dictionary)
df = df.sort_values(by='Length_of_CDS', ascending=False)
df = df.drop_duplicates(subset=['FBgn'])

df["NS_sites_per_bp"] = df["Number_NS_sites"] / df["Length_of_CDS"]

print("dataframe made")

#saving df
df.to_csv('potential_NS_site_counting.csv', index=False)


print("dataframe saved as csv")















""" 

s = "ATGTCAGACGCTGTGGACAAAATGTTGGCGGGAATGGCGGCCAATCGGCAGACCATGAACCGCCAGCTGGCTAAGATAGACGAGATAATGGAGCGGTCGAACAACACATTGCTCCACATCGAGTCCAACAGCAAGGCATTCAGTCAGAATGTGGCTTTGTCGGAGACCCAGAAGATGTACAACCTCAGGCCGGAGGCTGAGATGACCCTGTCCAAGATCCTGGAGAACTTCAAGCTGCTCATGTCTAGCAGTGACCAGCGGGAGGAGACATACAGCGCTCTGGAGGGCTGCCTGGCCTACAGGCATCGCGTGGAGCACCTTGGCAGCTCCGTGCGTAAGCTGGTGGCGCTCTACGACACCGTGGGCCAGATGAAGAACTCACAGGAAGAGCAGTATGCCAGCGAGGATTCCCCCTAG"

print(site_counter(s, 417, codon_NS_conversion))


 """