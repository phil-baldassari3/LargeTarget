"""
This programs takes the file downloaded from ChIP-Atlas Peak Browser ChIP: TFs and others threshold 200
Filters for only TF peaks using a list of TF probes downloaded from ENCODE Experiment Report: DNA binding,
TF ChIP-seq, transcription factor, Homo sapiens, cell line and outputs 2 csv files. One file convert the bed
file to a csv with columns "Binding_site_length" and "Cell_type". The other csv file has columns "Cell_Type",
"Number_of_binding_sites", and "Total_length_of_binding_sites".
"""

#import modules
import pandas as pd


bed_file = "ChIP_Atlas_TF_peaks_human_thresh200.bed"
TF_target_file = "ENCODE_ChIPseq_TF_targets.tsv"


#getting target names
targets = pd.read_csv(TF_target_file, skiprows=1, sep="\t")
targets_ls = targets["Target of assay"].to_list()
targets_ls = list(set(targets_ls))

#loading bed file
bed = pd.read_csv(bed_file, skiprows=1, header=None, usecols=[0,1,2,3], sep="\t")
cols = ["Chrom", "Start", "Stop", "Metadata"]
bed.columns = cols

#filtering
bed = bed[bed["Metadata"].str.contains("Name=" + "|".join(targets_ls))]

#bed to df to csv
tf_df = pd.DataFrame(columns=["Binding_site_length", "Cell_type"])

tf_df["Binding_site_length"] = bed["Stop"] - bed["Start"]
tf_df["Cell_type"] = bed["Metadata"].str.split(r"Cell%20group=").str[1].str.split(";").str[0].str.replace(r"%20", " ")




#grouping by cell types
cell_types = tf_df["Cell_type"].to_list()
cell_types = list(set(cell_types))

celltypes = []
num_sites = []
length = []

for ctype in cell_types:
    filtered = tf_df[tf_df["Cell_type"] == ctype]
    celltypes.append(ctype)
    num_sites.append(len(filtered))
    length.append(filtered["Binding_site_length"].sum())


dictionary = {"Cell_Type":celltypes, "Number_of_binding_sites":num_sites, "Total_length_of_binding_sites":length}

cell_type_df = pd.DataFrame(dictionary)


#Saving files
tf_df.to_csv("ChIP_seq_TF.csv", index=False)
cell_type_df.to_csv("ChIP_seq_TF_by_cell_type.csv", index=False)