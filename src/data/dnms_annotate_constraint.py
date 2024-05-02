"""Annotate DNMs with regional constraint and OMIM data."""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger
from src.data.dnms_tidy_logging import read_tidy_dnms


# Logging
logger = setup_logger(Path(__file__).stem)


# Module constants
_USECOLS_REGIONAL_CONSTRAINT = [
    "enst",
    "region",
    "constraint",
    "pli",
    "loeuf",
    "n_obs",
    "n_exp",
    "oe",
    "p",
    "fdr_p",
]


# Functions
def get_nonsense_constraint(path):
    """Get regional constraint annotations."""

    return pd.read_csv(path, sep="\t", usecols=_USECOLS_REGIONAL_CONSTRAINT)[
        _USECOLS_REGIONAL_CONSTRAINT
    ]


def get_nmd(path):
    """Get NMD annotations."""

    logger.info("Getting NMD annotations.")

    return pd.read_csv(
        path, sep="\t", usecols=["chr", "pos", "transcript_id", "nmd_definitive"]
    ).set_axis(["chr", "pos", "enst", "region"], axis=1)


def get_genes(path):
    """Read gene IDs."""

    genes = pd.read_csv(
        path,
        sep="\t",
    ).set_axis(["ensg", "enst", "symbol"], axis=1)

    logger.info(f"Unique transcript IDs: {genes.enst.nunique()}")
    logger.info(f"Unique gene synbols: {genes.symbol.nunique()}")
    
    return genes


def condense_phenotypes_and_inheritance(omim):
    """Condense OMIM annotations, so each gene is represented only once."""

    omim = (
        omim.groupby("ensg")[["phenotype", "inheritance"]]
        .agg(lambda x: " | ".join(x))
        .drop_duplicates()
        .reset_index(drop=False)
    )

    logger.info(
        f"OMIM entries after condensing phenotype and inheritance annotations: {len(omim)}"
    )
    logger.info(f"Unique ensg IDs: {omim.ensg.nunique()}")

    return omim


def simplify_inheritance(df):
    """Simplify inheritance annotation for OMIM morbid genes."""

    # Conditions
    dominant = df["inheritance"].str.lower().str.contains("dominant")
    recessive = df["inheritance"].str.lower().str.contains("recessive")
    x_linked = df["inheritance"].str.lower().str.contains("x-linked")

    df.loc[dominant, "inheritance_simple"] = "AD"
    df.loc[recessive, "inheritance_simple"] = "AR"
    df.loc[dominant & recessive, "inheritance_simple"] = "AD_AR"
    df.loc[x_linked, "inheritance_simple"] = "XL"
    df = df.fillna({"inheritance_simple": "Other"})

    logger.info(
        f"Simplified inheritance value counts:\n{df.inheritance_simple.value_counts(dropna=False)}"
    )

    return df


def get_omim(path, gene_ids_path):

    gene_ids = get_genes(gene_ids_path)

    omim = (
            pd.read_csv(path, sep="\t")
            .pipe(condense_phenotypes_and_inheritance)
            .pipe(simplify_inheritance)
            .merge(gene_ids[["ensg","enst"]])
            .drop("ensg", axis=1)
        )
        
    return omim

def main():
    """Run as script."""

    dnms = read_tidy_dnms(C.DNMS_VEP_TIDY)
    cst = get_nonsense_constraint(C.REGIONAL_NONSENSE_CONSTRAINT)
    nmd = get_nmd(C.NMD_ANNOTATION)
    gene_ids = get_genes(C.GENE_IDS)
    omim = get_omim(C.GENEMAP2_SIMPLE, C.GENE_IDS)


    logger.info(f"OMIM entries after merging gene IDs: {len(omim)}")
    logger.info(f"Unique transcript IDs: {omim.enst.nunique()}")

    # Find gene symbols for DNMs
    df = dnms.merge(gene_ids[["enst","symbol"]], how="left")
    logger.info(f"DNMs before merging with gene symbol data: {len(dnms)}")
    logger.info(f"DNMs after merging with gene symbols: {len(df)}")
    logger.info(f"DNMs with a gene symbol: {(~df.symbol.isna()).sum()}")

    # Merge NMD annotations
    df = df.merge(nmd, how="left")

    logger.info(f"DNMs after merging NMD annotations: {len(df)}")
    logger.info(f"DNMs without an NMD region: {df.region.isna().sum()}")
    logger.info(
        "The DNMs without an NMD region annotation are all indels. "
        "The most 3' position of these indels is outside of the CDS, "
        "and they therefore lack an NMD annotation. "
        "There is no easy way to rectify this, but we keep them for completeness."
    )

    # Merge constraint data
    df = df.merge(cst, how="left")

    logger.info(f"Transcripts in DNM data: {dnms.enst.nunique()}")
    logger.info(f"Transcripts in constraint data: {cst.enst.nunique()}")
    logger.info(
        f"DNM transcripts missing from constraint data: {(~dnms.enst.isin(cst.enst)).sum()}"
    )
    logger.info(f"DNMs after merging constraint data: {len(df)}")
    logger.info(f"DNMs with a constraint annotation: {df.n_obs.count()}")

    # Merge OMIM annotations
    df = df.merge(omim, how="left")

    logger.info(f"DNMs in OMIM morbid genes: {(~df.phenotype.isna()).sum()}")
    logger.info(
        f"DNMs in morbid genes by mode of inheritance:\n{df.inheritance_simple.value_counts(dropna=False)}"
    )

    # Write to output
    logger.info("Writing to output.")
    df.to_csv(C.DNMS_ANNOTATED, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
