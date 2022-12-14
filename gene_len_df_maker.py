import pandas as pd
import numpy as np
from statistics import mean
from statistics import stdev

#opening dataframes
flymine_neurogenesis = pd.read_csv('neurogenesis.csv')
flymine_circ_sys_dev = pd.read_csv('circulatory_system_development.csv')
flymine_circ_sys_proc = pd.read_csv('circulatory_system_process.csv')
flymine_digest_sys_dev = pd.read_csv('digestive_system_development.csv')
flymine_digest_sys_proc = pd.read_csv('digestive_system_process.csv')
flymine_endo_proc = pd.read_csv('endocrine_process.csv')
flymine_endo_sys_dev = pd.read_csv('endocrine_system_development.csv')
flymine_excretion = pd.read_csv('excretion.csv')
flymine_exo_sys_dev = pd.read_csv('exocrine_system_development.csv')
flymine_immune_sys_dev = pd.read_csv('immune_system_development.csv')
flymine_immune_sys_proc = pd.read_csv('immune_system_process.csv')
flymine_muscle = pd.read_csv('muscle_cell_development.csv')
flymine_renal_sys_dev = pd.read_csv('renal_system_development.csv')
flymine_renal_sys_proc = pd.read_csv('renal_system_process.csv')
flymine_reporduction = pd.read_csv('reproduction.csv')
flymine_respiratory = pd.read_csv('respiratory_system_development.csv')
flymine_urogenital = pd.read_csv('urogenital_system_development.csv')

print("dataframes all loaded\n")


#getting rid of duplicates
df_neurogenesis = flymine_neurogenesis.drop_duplicates(subset=['FBgn'])
df_circ_sys_dev = flymine_circ_sys_dev.drop_duplicates(subset=['FBgn'])
df_circ_sys_proc = flymine_circ_sys_proc.drop_duplicates(subset=['FBgn'])
df_digest_sys_dev = flymine_digest_sys_dev.drop_duplicates(subset=['FBgn'])
df_digest_sys_proc = flymine_digest_sys_proc.drop_duplicates(subset=['FBgn'])
df_endo_proc = flymine_endo_proc.drop_duplicates(subset=['FBgn'])
df_endo_sys_dev = flymine_endo_sys_dev.drop_duplicates(subset=['FBgn'])
df_excretion = flymine_excretion.drop_duplicates(subset=['FBgn'])
df_exo_sys_dev = flymine_exo_sys_dev.drop_duplicates(subset=['FBgn'])
df_immune_sys_dev = flymine_immune_sys_dev.drop_duplicates(subset=['FBgn'])
df_immune_sys_proc = flymine_immune_sys_proc.drop_duplicates(subset=['FBgn'])
df_muscle = flymine_muscle.drop_duplicates(subset=['FBgn'])
df_renal_sys_dev = flymine_renal_sys_dev.drop_duplicates(subset=['FBgn'])
df_renal_sys_proc = flymine_renal_sys_proc.drop_duplicates(subset=['FBgn'])
df_reporduction = flymine_reporduction.drop_duplicates(subset=['FBgn'])
df_respiratory = flymine_respiratory.drop_duplicates(subset=['FBgn'])
df_urogenital = flymine_urogenital.drop_duplicates(subset=['FBgn'])


#combination dataframes
df_circ = pd.concat([df_circ_sys_dev, df_circ_sys_proc], ignore_index=True, sort=False)
df_circ = df_circ.drop_duplicates(subset=['FBgn'])

df_digest = pd.concat([df_digest_sys_dev, df_digest_sys_proc], ignore_index=True, sort=False)
df_digest = df_digest.drop_duplicates(subset=['FBgn'])

df_endo = pd.concat([df_endo_proc, df_endo_sys_dev], ignore_index=True, sort=False)
df_endo = df_endo.drop_duplicates(subset=['FBgn'])

df_immune = pd.concat([df_immune_sys_dev, df_immune_sys_proc], ignore_index=True, sort=False)
df_immune = df_immune.drop_duplicates(subset=['FBgn'])

df_renal = pd.concat([df_renal_sys_dev, df_renal_sys_proc], ignore_index=True, sort=False)
df_renal = df_renal.drop_duplicates(subset=['FBgn'])

df_nonneural = pd.concat([df_circ, df_digest, df_endo, df_excretion, df_exo_sys_dev, df_immune, df_muscle, df_renal, df_reporduction, df_respiratory, df_urogenital], ignore_index=True, sort=False)
df_nonneural = df_nonneural.drop_duplicates(subset=['FBgn'])


print("duplicates removed and combination dataframes made\n")


#genes lists
genes_neurogenesis = df_neurogenesis['FBgn'].to_list()
genes_circ_sys_dev = df_circ_sys_dev['FBgn'].to_list()
genes_circ_sys_proc = df_circ_sys_proc['FBgn'].to_list()
genes_digest_sys_dev = df_digest_sys_dev['FBgn'].to_list()
genes_digest_sys_proc = df_digest_sys_proc['FBgn'].to_list()
genes_endo_proc = df_endo_proc['FBgn'].to_list()
genes_endo_sys_dev = df_endo_sys_dev['FBgn'].to_list()
genes_excretion = df_excretion['FBgn'].to_list()
genes_exo_sys_dev = df_exo_sys_dev['FBgn'].to_list()
genes_immune_sys_dev = df_immune_sys_dev['FBgn'].to_list()
genes_immune_sys_proc = df_immune_sys_proc['FBgn'].to_list()
genes_muscle = df_muscle['FBgn'].to_list()
genes_renal_sys_dev = df_renal_sys_dev['FBgn'].to_list()
genes_renal_sys_proc = df_renal_sys_proc['FBgn'].to_list()
genes_reporduction = df_reporduction['FBgn'].to_list()
genes_respiratory = df_respiratory['FBgn'].to_list()
genes_urogenital = df_urogenital['FBgn'].to_list()
genes_circ = df_circ['FBgn'].to_list()
genes_digest = df_digest['FBgn'].to_list()
genes_endo = df_endo['FBgn'].to_list()
genes_immune = df_immune['FBgn'].to_list()
genes_renal = df_renal['FBgn'].to_list()
genes_nonneural = df_nonneural['FBgn'].to_list()


#gene list text files
with open("neurogenesis_genes.txt", 'w') as file:
    file.write(", ".join(genes_neurogenesis))

with open("circulatory_system_development_genes.txt", 'w') as file:
    file.write(", ".join(genes_circ_sys_dev))

with open("circulatory_system_process_genes.txt", 'w') as file:
    file.write(", ".join(genes_circ_sys_proc))

with open("digestive_system_development_genes.txt", 'w') as file:
    file.write(", ".join(genes_digest_sys_dev))

with open("digestive_system_process_genes.txt", 'w') as file:
    file.write(", ".join(genes_digest_sys_proc))

with open("endocrine_process_genes.txt", 'w') as file:
    file.write(", ".join(genes_endo_proc))

with open("endocrine_system_development_genes.txt", 'w') as file:
    file.write(", ".join(genes_endo_sys_dev))

with open("excretion_genes.txt", 'w') as file:
    file.write(", ".join(genes_excretion))

with open("exocrine_system_development_genes.txt", 'w') as file:
    file.write(", ".join(genes_exo_sys_dev))

with open("immune_system_development_genes.txt", 'w') as file:
    file.write(", ".join(genes_immune_sys_dev))

with open("immune_system_process_genes.txt", 'w') as file:
    file.write(", ".join(genes_immune_sys_proc))

with open("muscle_genes.txt", 'w') as file:
    file.write(", ".join(genes_muscle))

with open("renal_system_development_genes.txt", 'w') as file:
    file.write(", ".join(genes_renal_sys_dev))

with open("renal_system_process_genes.txt", 'w') as file:
    file.write(", ".join(genes_renal_sys_proc))

with open("reproduction_genes.txt", 'w') as file:
    file.write(", ".join(genes_reporduction))

with open("respiratory_system_development_genes.txt", 'w') as file:
    file.write(", ".join(genes_respiratory))

with open("urogenital_system_development_genes.txt", 'w') as file:
    file.write(", ".join(genes_urogenital))

with open("circulatory_genes.txt", 'w') as file:
    file.write(", ".join(genes_circ))

with open("digestive_genes.txt", 'w') as file:
    file.write(", ".join(genes_digest))

with open("endocrine_genes.txt", 'w') as file:
    file.write(", ".join(genes_endo))

with open("immune_genes.txt", 'w') as file:
    file.write(", ".join(genes_immune))

with open("renal_genes.txt", 'w') as file:
    file.write(", ".join(genes_renal))

with open("nonneural_genes.txt", 'w') as file:
    file.write(", ".join(genes_nonneural))


print("gene lists written to files\n")



#length distribution lists
ls_neurogenesis = df_neurogenesis['gene_length'].to_list()
ls_circ_sys_dev = df_circ_sys_dev['gene_length'].to_list()
ls_circ_sys_proc = df_circ_sys_proc['gene_length'].to_list()
ls_digest_sys_dev = df_digest_sys_dev['gene_length'].to_list()
ls_digest_sys_proc = df_digest_sys_proc['gene_length'].to_list()
ls_endo_proc = df_endo_proc['gene_length'].to_list()
ls_endo_sys_dev = df_endo_sys_dev['gene_length'].to_list()
ls_excretion = df_excretion['gene_length'].to_list()
ls_exo_sys_dev = df_exo_sys_dev['gene_length'].to_list()
ls_immune_sys_dev = df_immune_sys_dev['gene_length'].to_list()
ls_immune_sys_proc = df_immune_sys_proc['gene_length'].to_list()
ls_muscle = df_muscle['gene_length'].to_list()
ls_renal_sys_dev = df_renal_sys_dev['gene_length'].to_list()
ls_renal_sys_proc = df_renal_sys_proc['gene_length'].to_list()
ls_reporduction = df_reporduction['gene_length'].to_list()
ls_respiratory = df_respiratory['gene_length'].to_list()
ls_urogenital = df_urogenital['gene_length'].to_list()
ls_circ = df_circ['gene_length'].to_list()
ls_digest = df_digest['gene_length'].to_list()
ls_endo = df_endo['gene_length'].to_list()
ls_immune = df_immune['gene_length'].to_list()
ls_renal = df_renal['gene_length'].to_list()
ls_nonneural = df_nonneural['gene_length'].to_list()




#making dictionaries
length_dict = {
    "Neurogenesis":np.array(ls_neurogenesis), "Circulatory_system_development":np.array(ls_circ_sys_dev), "Circulatory_system_process":np.array(ls_circ_sys_proc), 
    "Digestive_system_development":np.array(ls_digest_sys_dev), "Digestive_system_process":np.array(ls_digest_sys_proc), "Endocrine_process":np.array(ls_endo_proc), 
    "Endocrine_system_development":np.array(ls_endo_sys_dev), "Excretion":np.array(ls_excretion), "Exocrine_system_development":np.array(ls_exo_sys_dev), 
    "Immune_system_development":np.array(ls_immune_sys_dev), "Immune_system_process":np.array(ls_immune_sys_proc), "Muscle_cell_development":np.array(ls_muscle), 
    "Renal_system_development":np.array(ls_renal_sys_dev), "Renal_system_process":np.array(ls_renal_sys_proc), "Reproduction":np.array(ls_reporduction), 
    "Respiratory_system_development":np.array(ls_respiratory), "Urogenital_system_development":np.array(ls_urogenital), "Circulatory":np.array(ls_circ),
    "Digestive":np.array(ls_digest), "Endocrine":np.array(ls_endo), "Immune":np.array(ls_immune), "Renal":np.array(ls_renal), "Nonneural":np.array(ls_nonneural)
}



summary_dict = {
    "Class":[
        "Neurogenesis", "Circulatory_system_development", "Circulatory_system_process", 
        "Digestive_system_development", "Digestive_system_process", "Endocrine_process", 
        "Endocrine_system_development", "Excretion", "Exocrine_system_development", 
        "Immune_system_development", "Immune_system_process", "Muscle_cell_development", 
        "Renal_system_development", "Renal_system_process", "Reproduction", 
        "Respiratory_system_development", "Urogenital_system_development", "Circulatory", 
        "Digestive", "Endocrine", "Immune", "Renal", "Nonneural"
    ],

    "Total":[
        sum(ls_neurogenesis), sum(ls_circ_sys_dev), sum(ls_circ_sys_proc), 
        sum(ls_digest_sys_dev), sum(ls_digest_sys_proc), sum(ls_endo_proc), 
        sum(ls_endo_sys_dev), sum(ls_excretion), sum(ls_exo_sys_dev), 
        sum(ls_immune_sys_dev), sum(ls_immune_sys_proc), sum(ls_muscle), 
        sum(ls_renal_sys_dev), sum(ls_renal_sys_proc), sum(ls_reporduction), 
        sum(ls_respiratory), sum(ls_urogenital), sum(ls_circ), 
        sum(ls_digest), sum(ls_endo), sum(ls_immune), sum(ls_renal), sum(ls_nonneural)
    ],

    "Mean":[
        mean(ls_neurogenesis), mean(ls_circ_sys_dev), mean(ls_circ_sys_proc), 
        mean(ls_digest_sys_dev), mean(ls_digest_sys_proc), mean(ls_endo_proc), 
        mean(ls_endo_sys_dev), mean(ls_excretion), mean(ls_exo_sys_dev), 
        mean(ls_immune_sys_dev), mean(ls_immune_sys_proc), mean(ls_muscle), 
        mean(ls_renal_sys_dev), mean(ls_renal_sys_proc), mean(ls_reporduction), 
        mean(ls_respiratory), mean(ls_urogenital), mean(ls_circ), 
        mean(ls_digest), mean(ls_endo), mean(ls_immune), mean(ls_renal), mean(ls_nonneural)  
    ],

    "SD":[
        stdev(ls_neurogenesis), stdev(ls_circ_sys_dev), stdev(ls_circ_sys_proc), 
        stdev(ls_digest_sys_dev), stdev(ls_digest_sys_proc), stdev(ls_endo_proc), 
        stdev(ls_endo_sys_dev), stdev(ls_excretion), stdev(ls_exo_sys_dev), 
        stdev(ls_immune_sys_dev), stdev(ls_immune_sys_proc), stdev(ls_muscle), 
        stdev(ls_renal_sys_dev), stdev(ls_renal_sys_proc), stdev(ls_reporduction), 
        stdev(ls_respiratory), stdev(ls_urogenital), stdev(ls_circ), 
        stdev(ls_digest), stdev(ls_endo), stdev(ls_immune), stdev(ls_renal), stdev(ls_nonneural)
    ]
}



#saving dataframes
lengths_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in length_dict.items() ]))

summary_df = pd.DataFrame(summary_dict)

print("new dataframes made\n")

lengths_df.to_csv("gene_lengths_GO.csv", index=False)

summary_df.to_csv("total_avg_gene_lengths.csv", index=False)

print("DONE!")


