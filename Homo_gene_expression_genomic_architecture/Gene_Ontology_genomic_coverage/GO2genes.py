import pandas as pd
from intermine.webservice import Service



def humanmineGO2df(GOterm):

    service = Service("https://www.humanmine.org/humanmine/service")

    # Get a new query on the class (table) you will be querying:
    query = service.new_query("Gene")

    # The view specifies the output columns
    query.add_view(
        "secondaryIdentifier", "symbol", "goAnnotation.ontologyTerm.parents.name",
        "goAnnotation.ontologyTerm.parents.identifier", "length"
    )

    # This query's custom sort order is specified below:
    query.add_sort_order("Gene.symbol", "ASC")

    # You can edit the constraint values below
    query.add_constraint("goAnnotation.ontologyTerm.parents.name", "=", GOterm, code="C")

    #getting df
    df = query.dataframe()

    #extracting columns
    df = df[["Gene.secondaryIdentifier", "Gene.symbol", "Gene.goAnnotation.ontologyTerm.parents.name", "Gene.goAnnotation.ontologyTerm.parents.identifier", "Gene.length"]]

    return df




#open GOterm list
GOterms_df = pd.read_csv("GO_terms.tsv", sep="\t")
GOterms = GOterms_df["term"].to_list()


#creating GO dataframes
GO_df_ls = []

for term in GOterms:
    dfGO = humanmineGO2df(term)
    GO_df_ls.append(dfGO)

#concatenate dfs
full_df = pd.concat(GO_df_ls)

#opening CDS lengths file
CDS_df = pd.read_csv("CDS_lengths.csv")
CDS_df = CDS_df[CDS_df['CDS Length'].notna()]
CDS_df = CDS_df.sort_values("CDS Length")

CDS_ls = CDS_df["Gene stable ID"].to_list()
CDS_len_ls = CDS_df["CDS Length"].to_list()



#CDS length dictionary
CDS_len_dict = dict(zip(CDS_ls, CDS_len_ls))

#matching CDS lengths to full df
genes = full_df["Gene.secondaryIdentifier"].to_list()
full_CDS_lens = []
missing_CDSs = []

for gene in genes:
    try:
        full_CDS_lens.append(CDS_len_dict[gene])
    except KeyError:
        full_CDS_lens.append("")
        missing_CDSs.append(gene)

full_df["CDS.length"] = full_CDS_lens


#saving df
full_df.to_csv("GO_genes_lengths.csv", index=False)

#print missing
missing_CDSs = list(set(missing_CDSs))
for i in missing_CDSs:
    print(i)