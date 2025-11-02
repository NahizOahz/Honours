import pandas as pd


df = pd.read_csv("/Applications/Document/denovo_blast_trinity_own/DEG_list_copy/DEG_ALE35_saline_annot.csv")


df_sig = df[(df["padj"] < 0.05) & (~df["gene_id"].isna())]


wt_up = set(df_sig[df_sig["log2FoldChange"] > 0]["gene_id"].unique())
wt_down = set(df_sig[df_sig["log2FoldChange"] < 0]["gene_id"].unique())


df = pd.read_csv("/Applications/Document/denovo_blast_trinity_own/DEG_list_copy/DEG_ALE100_saline_annot.csv")


df_sig = df[(df["padj"] < 0.05) & (~df["gene_id"].isna())]


ale100_up = set(df_sig[df_sig["log2FoldChange"] > 0]["gene_id"].unique())
ale100_down = set(df_sig[df_sig["log2FoldChange"] < 0]["gene_id"].unique())


df = pd.read_csv("/Applications/Document/denovo_blast_trinity_own/DEG_list_copy/DEG_strain_35ppt_annot.csv")


df_sig = df[(df["padj"] < 0.05) & (~df["gene_id"].isna())]


salt35_up = set(df_sig[df_sig["log2FoldChange"] > 0]["gene_id"].unique())
salt35_down = set(df_sig[df_sig["log2FoldChange"] < 0]["gene_id"].unique())


df = pd.read_csv("/Applications/Document/denovo_blast_trinity_own/DEG_list_copy/DEG_strain_100ppt_annot.csv")


df_sig = df[(df["padj"] < 0.05) & (~df["gene_id"].isna())]


salt100_up = set(df_sig[df_sig["log2FoldChange"] > 0]["gene_id"].unique())
salt100_down = set(df_sig[df_sig["log2FoldChange"] < 0]["gene_id"].unique())


from venn import venn
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'


sets1 = {
    "WT Up": wt_up,
    "WT Down": wt_down,
    "ALE100 Up": ale100_up,
    "ALE100 Down": ale100_down
}
venn(sets1, cmap="tab10")
plt.title("Salt response (WT vs ALE100)  Up & Down", fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig("/Users/zhao/Desktop/OWN_fig/sets1_colored.png", dpi=600, bbox_inches="tight")
plt.close()


sets2 = {
    "35 ppt Up": salt35_up,
    "35 ppt Down": salt35_down,
    "100 ppt Up": salt100_up,
    "100 ppt Down": salt100_down
}
venn(sets2, cmap="Dark2_r")
plt.title("Evolution response (35 ppt vs 100 ppt)  Up & Down", fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig("/Users/zhao/Desktop/OWN_fig/sets2_colored.png", dpi=600, bbox_inches="tight")
plt.close()


sets3 = {
    "WT": wt_up.union(wt_down),
    "ALE100": ale100_up.union(ale100_down),
    "35 ppt": salt35_up.union(salt35_down),
    "100 ppt": salt100_up.union(salt100_down)
}
venn(sets3, cmap="RdBu")
plt.title("Four-condition overlap (All DEGs)", fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig("/Users/zhao/Desktop/OWN_fig/sets3_colored.png", dpi=600, bbox_inches="tight")
plt.close()