""" Combine LabKey and DNM data for clinical review.

This script merges the annotated DNMs with participant-level clinical data from 
LabKey. The resulting DataFrame is written to output.
"""

import logging

import pandas as pd

import src

_DE_NOVO_OFFSPRING = "data/interim/gel_dnm_offspring_clean.tsv"
_DNMS_ANNOTATED = "data/interim/dnms_annotated.tsv"
_LABKEY_CLINICAL = "data/interim/labkey_participant_clinical.tsv"
_FILE_OUT = "data/interim/dnms_annotated_clinical.tsv"
_SPREADSHEET_COLS = [
    "symbol",
    "enst",
    "chr",
    "pos",
    "ref",
    "alt",
    "phenotype",
    "inheritance",
    "inheritance_simple",
    "region",
    "constraint",
    "oe_ci_hi",
    "pli",
    "loeuf",
    "n_obs",
    "n_exp",
    "oe",
    "p",
    "fdr_p",
    "n_trunc_region",
    "n_trunc_transcript",
    "cohort",
    "participant_id",
    "id",
    "csq",
    "case_solved",
    "max_tier",
    "gmc",
    "birth_year",
    "hpo",
]

logger = logging.getLogger(__name__)


def read_offspring_ids(path):
    """Read the offpsring IDs"""

    ids = pd.read_csv(
        path,
        sep="\t",
        usecols=["Trio Id", "Participant Id"],
    ).set_axis(["id", "participant_id"], axis=1)

    logger.info(f"GEL participant IDs: {len(ids)}")
    logger.info(f"Unique GEL participant IDs: {ids.participant_id.nunique()}")

    return ids


def tidy_dnm_data_for_manual_review(df):
    """Tidy the combined dataset for manual review"""

    # Find the number of truncating variants per region
    n_trunc = lambda x: x.isin(["frameshift_variant", "stop_gained"]).sum()

    df["n_trunc_region"] = df.groupby(["enst", "region"], dropna=False)[
        "csq"
    ].transform(n_trunc)
    df["n_trunc_transcript"] = df.groupby("enst")["csq"].transform(n_trunc)

    # Reorder columns for ease of use in spreadsheet
    df = df[_SPREADSHEET_COLS].copy()

    # Rename columns
    df = df.rename(
        columns={
            "phentype": "omim_phenotype",
            "inheritance": "omim_inheritance",
            "inheritance_simple": "omim_inheritance_simple",
        }
    )

    df["comment"] = "."  # Forces wrapping of HPO text in the spreadsheet

    return df


def main():
    """Run the script"""

    # Read data sets to memory
    ids = read_offspring_ids(_DE_NOVO_OFFSPRING)
    dnms = pd.read_csv(_DNMS_ANNOTATED, sep="\t")
    clinical = pd.read_csv(_LABKEY_CLINICAL, sep="\t")

    # Merge and tidy the data
    df = (
        dnms.merge(ids, how="left")
        .merge(clinical, how="left")
        .pipe(tidy_dnm_data_for_manual_review)
    )

    # Logging
    logger.info(f"Starting DNMs: {len(dnms)}")
    logger.info(f"DNMs after merging annotations: {len(df)}")
    logger.info(f"DNMs in GEL participants: {(df.cohort == 'GEL').sum()}")
    logger.info(f"Unique GEL participants: {df.participant_id.nunique()}")
    logger.info(
        f"GEL participants lacking an ID: "
        f"{((df.cohort == 'GEL') & (df.participant_id.isna())).sum()}"
    )
    logger.info(
        f"GEL participants with no HPO terms: "
        f"{df[(df.cohort=='GEL') & (df.hpo.isna())].participant_id.nunique()}. "
        "These individuals are (presumably unaffected) relatives of a proband."
    )
    logger.info(
        f"DNMs in GEL participants with no HPO terms: "
        f"{((df.cohort=='GEL') & (df.hpo.isna())).sum()}"
    )

    # Write to output
    logger.info("Writing to output.")
    df.to_csv(_FILE_OUT, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
