#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


vsd_csv = "/Applications/Document/denovo_blast_trinity_own/DEG_list/VST_expression_matrix.csv"
outdir  = "/Applications/Document/denovo_blast_trinity_own/DEG_list"
out_png = f"{outdir}/Heatmap_VST_top500_mean_noAnn.png"
out_pdf = f"{outdir}/Heatmap_VST_top500_mean_noAnn.pdf"


df = pd.read_csv(vsd_csv, index_col=0)


def infer_group(name: str) -> str:
    strain = "WT" if re.search(r"^WT", name) else "ALE100"
    salinity = "35ppt" if re.search(r"_35_", name) else "100ppt"
    return f"{strain}_{salinity}"


group_map = {col: infer_group(col) for col in df.columns}
df_mean = df.groupby(group_map, axis=1).mean()


topn = min(500, df_mean.shape[0])
var_series = df_mean.var(axis=1)
df_top = df_mean.loc[var_series.nlargest(topn).index]


df_z = df_top.sub(df_top.mean(axis=1), axis=0)
df_z = df_z.div(df_top.std(axis=1).replace(0, np.nan), axis=0).fillna(0)


sns.set_context("notebook", font_scale=1.0)
sns.set_style("white")


cmap = sns.diverging_palette(240, 10, as_cmap=True)

g = sns.clustermap(
    df_z,
    cmap=cmap,
    row_cluster=True,
    col_cluster=True,
    metric="correlation",
    method="average",
    xticklabels=True,
    yticklabels=False,
    figsize=(8, 8),
    cbar_pos=(0.02, 0.8, 0.02, 0.15),
    dendrogram_ratio=(0.15, 0.1)
)

# g.ax_heatmap.set_title("Top 500 Most Variable Genes (WT vs ALE100, averaged replicates)", pad=12)


plt.savefig(out_png, dpi=600, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
