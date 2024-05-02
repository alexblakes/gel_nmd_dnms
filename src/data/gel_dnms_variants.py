"""Process the GEL DNMs.

Filter for DNMs passing "stringent" filters, in our cleaned list of trio offspring. 
Where a trio has DNMs in both GRCh37 and GRCh38, keep only those in GRCh38. Variants in 
GRCh37 and GRCh38 are written out to separate files. They will be recombined after 
liftover.
"""


# Imports
from pathlib import Path

import pandas as pd

from src.data.ddd_dnms_to_vcf import write_vcf
from src import setup_logger
from src import constants as C


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_dnms(path):
    """Get DNM variant data."""

    logger.info("Reading GEL de novo variant data.")

    dnms = pd.read_csv(
        path,
        sep="\t",
        usecols=[
            "Trio Id",
            "Family Id",
            "Assembly",
            "Chrom",
            "Position",
            "Reference",
            "Alternate",
            "Base Filter",
            "Stringent Filter",
        ],
        low_memory=False,
    ).rename(
        columns={
            "Chrom": "chr",
            "Position": "pos",
            "Reference": "ref",
            "Alternate": "alt",
            "Base Filter": "base",
            "Stringent Filter": "stringent",
        }
    )

    logger.info(f"DNMs: {len(dnms)}")

    return dnms


def stringent_dnms(dnms):
    """Filter for stringent DNMs."""

    dnms = dnms[dnms["stringent"] == 1]

    logger.info(f"Stringent DNMs: {len(dnms)}")

    return dnms


def dnms_in_cleaned_trios(dnms, offspring):
    """Find DNMs in cleaned trios."""

    logger.info(f"Trio IDs in offspring data: {offspring['Trio Id'].nunique()}")
    logger.info(f"Trio IDs in variant data (before filter): {dnms['Trio Id'].nunique()}")

    dnms = dnms[dnms["Trio Id"].isin(offspring["Trio Id"])]

    logger.info(f"Trio IDs in variant data (after filter): {dnms['Trio Id'].nunique()}")
    logger.info(f"Stringent DNMs in slected trios: {len(dnms)}")

    return dnms


def prioritise_grch38_variants(dnms):
    """Where a trio has DNMs in both GRCh37 and GRCh38, keep only those in GRCh38."""

    # Find Trio Ids with variants in both GRCh37 and GRCh38
    trios = dnms[["Trio Id", "Assembly"]].drop_duplicates()
    trios = trios["Trio Id"]
    dup_trios = trios[trios.duplicated()]
    
    logger.info(f"Trios in both GRCh37 and GRCh38: {trios.duplicated().sum()}")

    # Drop GRCh37 variants, where GRCh38 variants are available for that trio
    m1 = dnms["Trio Id"].isin(dup_trios)
    m2 = dnms["Assembly"] == "GRCh37"

    logger.info(f"DNMs in GRCh37, where GRCh38 is available for the trio: {(m1 & m2).sum()}")

    dnms = dnms[~(m1 & m2)]

    logger.info(f"DNMs after dropping redundant GRCh37 variants: {len(dnms)}")

    return dnms


def subset_by_assembly(df, assembly):
    """Subset to variants in the given assembly."""
    
    df = df[df["Assembly"] == assembly]
    
    logger.info(f"DNMs in {assembly}: {len(df)}")
    
    return df


def reformat_to_vcf(df):
    """Reformat to VCF."""

    vcf = df.assign(QUAL=".", FILTER=".", INFO=".")
    vcf = vcf[["chr", "pos", "Trio Id", "ref", "alt", "QUAL", "FILTER", "INFO"]]

    return vcf


def main():
    """Run as script."""

    # Get offpsring annotations
    offspring = pd.read_csv(C.DE_NOVO_OFFSPRING, sep="\t")

    # Tidy DNM variant data
    dnms = (
        get_dnms(C.GEL_DNMS)
        .pipe(stringent_dnms).pipe(dnms_in_cleaned_trios, offspring)
        .pipe(prioritise_grch38_variants)
    )

    # Process GRCh37 and GRCh38 variants separately
    for assembly, path in zip(["GRCh37", "GRCh38"], [C.GEL_DNMS_VCF_37, C.GEL_DNMS_VCF_38]):
        variants = dnms.pipe(subset_by_assembly, assembly).pipe(reformat_to_vcf)

        logger.info(f"Writing {assembly} variants to output.")
        write_vcf(variants, path)

    return None


if __name__ == "__main__":
    main()
