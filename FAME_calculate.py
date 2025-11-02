import pandas as pd


df = pd.read_csv("/Users/zhao/Downloads/Phenotype——fig/FAME_ALE100.csv", header=None)  



samples = df.iloc[0, 1:].tolist()      
sample_mass = df.iloc[1, 1:].tolist()  


fa_names = df.iloc[2:, 0].tolist()     
data = df.iloc[2:, 1:].astype(float)   
data.columns = samples
data.index = fa_names


df_tidy = data.reset_index().melt(id_vars="index", 
                                  var_name="Sample", 
                                  value_name="Value_mg")
df_tidy.rename(columns={"index": "FattyAcid"}, inplace=True)


df_tidy["Strain"] = df_tidy["Sample"].str.extract(r"(ALE100|WT)")
df_tidy["Salinity"] = df_tidy["Sample"].str.extract(r"(\d+ppt)")
df_tidy["Replicate"] = df_tidy["Sample"].str.extract(r"-(\d+)$")


mass_dict = dict(zip(samples, sample_mass))
df_tidy["Sample_mg"] = df_tidy["Sample"].map(mass_dict)


df_tidy["Value_mg_gDW"] = df_tidy["Value_mg"] / df_tidy["Sample_mg"] * 1000


df_tidy.to_csv("/Users/zhao/Downloads/Phenotype——fig/FAME_tidy.csv", index=False)
print(df_tidy.head())