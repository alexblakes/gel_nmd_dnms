"""Split annotated DNMs into dominant, recessive, and non-morbid genes."""

# Import modules
import numpy as np
import pandas as pd

# Load data
df = pd.read_csv("../outputs/dnms_annotated_clinical.tsv", sep="\t").fillna(
    value={"omim": ""}
)

# Masks for disease categories
## Dominant
m1 = df["omim"].str.lower().str.contains("autosomal dominant")
m2 = df["omim"].str.lower().str.contains("x-linked")
m3 = m1 | m2

## Non morbid
m4 = df["omim"] == ""

## Recessive
m5 = df["omim"].str.lower().str.contains("autosomal recessive")
m6 = ~m3 & m5  # Exclude genes with both AD and AR phenotypes

## Other (often somatic or "Inheritance not available")
m7 = ~(m3 | m4 | m6)

# Define masks and names
masks = [m3, m4, m6, m7]
disease = ["dominant", "non_morbid", "recessive", "other"]

# Check that each variant belongs to a distinct category
assert sum([x.sum() for x in masks]) == len(df)

# Apply masks and write to output
for m, d in zip(masks, disease):
    filtered_df = df[m]
    outpath = f"../outputs/dnms_annotated_clinical_{d}.tsv"
    filtered_df.to_csv(outpath, sep="\t", index=False)