
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# === 1. tidy table ===
df = pd.read_csv("/Users/zhao/Downloads/Phenotype——fig/FAME_tidy.csv")

# === 2. cal mg/g DW ===
df["Value_mg_gDW"] = df["Value_mg"] / df["Sample_mg"] * 1000

# === 3. matrix ===
heatmap_data = df.pivot(index="FattyAcid", columns="Sample", values="Value_mg_gDW")

# --- 1. class (SFA / MUFA / PUFA) ---
SFA = ["Methyl tetradecanoate","Methyl pentadecanoate IM","Methyl palmitate",
       "Methyl heptadecanoate ","Methyl octadecanoate (methyl stearate)",
       "Methyl arachidate","Methyl heneicosanoate ","Methyl docosanoate",
       "Methyl tricosanoate","Methyl lignocarate"]
MUFA = ["Methyl palmitoleate","Methyl cis-11-eicosenoate","Methyl erucate"]
PUFA = ["Cis-11,14-Eicosadienoic acid methyl ester"]

fa_class_map = {}
for fa in heatmap_data.index:
    if fa in SFA:
        fa_class_map[fa] = "SFA"
    elif fa in MUFA:
        fa_class_map[fa] = "MUFA"
    elif fa in PUFA:
        fa_class_map[fa] = "PUFA"
    else:
        fa_class_map[fa] = "Other"




row_colors = pd.Series(fa_class_map, index=heatmap_data.index).map(class_colors)


g = sns.clustermap(heatmap_data,
                   metric="correlation",
                   method="average",
                   z_score=0,
                   cmap="viridis",
                   row_colors=row_colors,
                   figsize=(12,8))


for label, color in class_colors.items():
    g.ax_heatmap.bar(0, 0, color=color, label=label, linewidth=0)

g.ax_heatmap.legend(loc="upper center",
                    bbox_to_anchor=(0.5, 1.15),
                    ncol=4, fontsize=10, frameon=False)


plt.show()