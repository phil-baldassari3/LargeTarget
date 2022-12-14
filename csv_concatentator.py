#importing modules
import numpy as np
import pandas as pd
from multiprocessing import Pool



#opening dataframes
df_neurogenesis = pd.read_csv('neurogenesis.csv')
df_exo_sys_dev = pd.read_csv('exocrine_system_development.csv')
df_excretion = pd.read_csv('excretion.csv')
df_muscle = pd.read_csv('muscle_cell_development.csv')
df_reproduction = pd.read_csv('reproduction.csv')
df_respiratory = pd.read_csv('respiratory_system_development.csv')
df_urogenital = pd.read_csv('urogenital_system_development.csv')
df_circulatory = pd.read_csv('circulatory.csv')
df_digestive = pd.read_csv('digestive.csv')
df_endocrine = pd.read_csv('endocrine.csv')
df_immune = pd.read_csv('immune.csv')
df_renal = pd.read_csv('renal.csv')
df_nonneural = pd.read_csv('nonneural.csv')


df_neurogenesis["GO"] = "Neurogenesis"
df_exo_sys_dev["GO"] = "Exocrine_system_development"
df_excretion["GO"] = "Excretion"
df_muscle["GO"] = "Muscle"
df_reproduction["GO"] = "Reproduction"
df_respiratory["GO"] = "Respiratory"
df_urogenital["GO"] = "Urogenital"
df_circulatory["GO"] = "Circulatory"
df_digestive["GO"] = "Digestive"
df_endocrine["GO"] = "Endocrine"
df_immune["GO"] = "Immune"
df_renal["GO"] = "Renal"
df_nonneural["GO"] = "Non-neural"






all_df = pd.concat([df_neurogenesis, df_exo_sys_dev, df_excretion, df_muscle, df_reproduction, df_respiratory, df_urogenital, df_circulatory, df_digestive, df_endocrine, df_immune, df_renal], axis=0)

neurononneuro_df = pd.concat([df_neurogenesis, df_nonneural], axis=0)


all_df.to_csv('all.csv', index=False)
neurononneuro_df.to_csv('neurononneuro.csv', index=False)
