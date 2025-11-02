import pandas as pd
import numpy as np
from pathlib import Path

path = Path("/Applications/Document/denovo_blast/de_blastp_trinity_sprot.tsv")


cols = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","slen","qlen"]


df = pd.read_csv(path, sep="\t", header=None, names=cols)


for c in ["pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","slen","qlen"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")


df["qcov"] = df["length"] / df["qlen"]


df["pass_filter"] = (df["evalue"] < 1e-5) & (df["pident"] > 30) & (df["qcov"] > 0.5)


all_passed = df[df["pass_filter"]].copy()
all_passed.to_csv("all_hits_filtered.csv", index=False)


best = (
    df.sort_values(by=["evalue","bitscore","pident"], ascending=[True,False,False])
      .groupby("qseqid", as_index=False)
      .head(1)
)
best_passed = best[best["pass_filter"]].copy()
best_passed.to_csv("best_hits_filtered.csv", index=False)


total_queries = df["qseqid"].nunique()
queries_with_hits = best_passed["qseqid"].nunique()
print(f"total query : {total_queries}")
print(f" 1 query : {queries_with_hits} ({queries_with_hits/total_queries:.1%})")
print(f"no query : {total_queries - queries_with_hits}")