import pandas as pd
import numpy as np
import os


#openning GO genes and lengths csv
df = pd.read_csv("GO_genes_lengths.csv")

#getting list of GO terms
terms = list(set(df["Gene.goAnnotation.ontologyTerm.parents.name"].to_list()))

#calculating coverage
GOterm = []
genelength = []
CDSlength = []

for term in terms:

    GOterm.append(term)

    filtered_df = df[df['Gene.goAnnotation.ontologyTerm.parents.name'] == term]

    genelen = sum(filtered_df["Gene.length"].to_list())
    CDSlen = sum(filtered_df["CDS.length"].to_list())

    genelength.append(genelen)
    CDSlength.append(CDSlen)


#making coverage df
dictionary = {"GO_term":GOterm, "Gene_coverage":genelength, "CDS_coverage":CDSlength}
final_df = pd.DataFrame(dictionary)

#saving file
final_df.to_csv("GO_coverage.csv", index=False)