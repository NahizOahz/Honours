import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os


wt_path     = "/Applications/Document/denovo_blast_trinity_own/degslist/DEG_ALE35_saline_annot.csv"  # ALE35=WT
ale100_path = "/Applications/Document/denovo_blast_trinity_own/degslist/DEG_ALE100_saline_annot.csv"
outdir      = "/Applications/Document/denovo_blast_trinity_own/degslist"
os.makedirs(outdir, exist_ok=True)


padj_sig   = 0.05
lfc_up     = 1.0
lfc_down   = -1.0
top_n      = 20
vmin, vmax = -12, 12
metric, method = "correlation", "average"
use_zscore = False


wt_raw     = pd.read_csv(wt_path)
ale100_raw = pd.read_csv(ale100_path)


wt_dedup = (wt_raw
            .dropna(subset=["uniprot_acc","log2FoldChange","padj"])
            .sort_values(by="log2FoldChange", key=lambda s: s.abs(), ascending=False)
            .drop_duplicates(subset="uniprot_acc", keep="first"))

ale100_dedup = (ale100_raw
                .dropna(subset=["uniprot_acc","log2FoldChange","padj"])
                .sort_values(by="log2FoldChange", key=lambda s: s.abs(), ascending=False)
                .drop_duplicates(subset="uniprot_acc", keep="first"))


wt = wt_dedup.rename(columns={"log2FoldChange":"WT_100vs35","padj":"WT_padj"})[["uniprot_acc","WT_100vs35","WT_padj"]]
ale100 = ale100_dedup.rename(columns={"log2FoldChange":"ALE100_100vs35","padj":"ALE100_padj"})[["uniprot_acc","ALE100_100vs35","ALE100_padj"]]
df = pd.merge(wt, ale100, on="uniprot_acc", how="inner").dropna()


WT_up       = (df["WT_padj"]    < padj_sig) & (df["WT_100vs35"]    >  lfc_up)
WT_down     = (df["WT_padj"]    < padj_sig) & (df["WT_100vs35"]    <  lfc_down)
WT_neutral  = (df["WT_padj"]    >= padj_sig) & (df["WT_100vs35"].abs() <= 1)

ALE100_up   = (df["ALE100_padj"] < padj_sig) & (df["ALE100_100vs35"] >  lfc_up)
ALE100_down = (df["ALE100_padj"] < padj_sig) & (df["ALE100_100vs35"] <  lfc_down)



mask_A = ALE100_up & (WT_down | WT_neutral)

mask_B = ALE100_down & (WT_up   | WT_neutral)

sel_A = df.loc[mask_A].copy()
sel_A["category"] = "ALE↑ & WT(↓/neutral)"
sel_B = df.loc[mask_B].copy()
sel_B["category"] = "ALE↓ & WT(↑/neutral)"

sel = pd.concat([sel_A, sel_B], axis=0, ignore_index=True)

if sel.empty:
    raise SystemExit("none")


for d in (sel_A, sel_B):
    d["max_abs"] = d[["WT_100vs35","ALE100_100vs35"]].abs().max(axis=1)

if top_n is not None:
    sel_A = sel_A.sort_values("max_abs", ascending=False).head(top_n)
    sel_B = sel_B.sort_values("max_abs", ascending=False).head(top_n)
else:
    sel_A = sel_A.sort_values("max_abs", ascending=False)
    sel_B = sel_B.sort_values("max_abs", ascending=False)


sel = pd.concat([sel_A, sel_B], axis=0, ignore_index=True)


csv_out = os.path.join(outdir, "ALE_both_patterns_top.csv")
sel[["uniprot_acc","WT_100vs35","ALE100_100vs35","WT_padj","ALE100_padj","category"]].to_csv(csv_out, index=False)
print(f"save: {csv_out}")


data = sel.set_index("uniprot_acc")[["WT_100vs35","ALE100_100vs35"]].astype(float)
if use_zscore:
    data = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1).replace(0,1), axis=0)
    cbar_label = "log2FC (Z-score by gene)"
else:
    cbar_label = "log2FC"

sns.set(context="notebook", style="white")
g = sns.clustermap(
    data,
    cmap="RdBu_r", vmin=vmin, vmax=vmax,
    metric=metric, method=method,
    row_cluster=True, col_cluster=False,
    figsize=(10, max(6, 0.35*data.shape[0])),
    cbar_kws={"label": cbar_label},
    dendrogram_ratio=(0.15, 0),
    xticklabels=True, yticklabels=True
)

g.cax.set_position([1, 0.4, 0.05, 0.2])


g.ax_heatmap.set_title("ALE↑ & WT(↓/neutral)  +  ALE↓ & WT(↑/neutral) — Row clustering only", pad=10, fontsize=11)
g.ax_heatmap.set_xticklabels(["WT 100 vs 35", "ALE100 100 vs 35"], rotation=0)

png = os.path.join(outdir, "ALE_both_patterns_clustermap.png")
pdf = os.path.join(outdir, "ALE_both_patterns_clustermap.pdf")
plt.savefig(png, dpi=600); plt.savefig(pdf)
plt.show()
print("save:", png, pdf)


sel_A = df.loc[mask_A].copy()
sel_A["category"] = "ALE↑ & WT(↓/neutral)"
sel_B = df.loc[mask_B].copy()
sel_B["category"] = "ALE↓ & WT(↑/neutral)"

print(f"🔹 ALE↑ & WT(↓/neutral): {sel_A.shape[0]} ")
print(f"🔸 ALE↓ & WT(↑/neutral): {sel_B.shape[0]} ")
sel = pd.concat([sel_A, sel_B], axis=0, ignore_index=True)
print(f": {sel.shape[0]} ")
