"""Find transcripts with >= 3 de novo PTVs in a constrained region."""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/interim/dnms_annotated_clinical.tsv"
_FILE_OUT = "data/final/recurrent_dnms.tsv"

logger = logging.getLogger(__name__)


def read_dnm_ptvs(path):
    ptvs = (
        pd.read_csv(
            path,
            sep="\t",
            usecols=[
                "enst",
                "symbol",
                "omim_inheritance_simple",
                "region",
                "n_trunc_region",
                "constraint",
                "pli",
                "loeuf",
                "csq",
            ],
        )
        .query("csq == 'frameshift_variant' | csq == 'stop_gained'")
        .query("n_trunc_region >= 3")
        .query("constraint == 'constrained'")
        .drop(["csq", "constraint"], axis=1)
        .drop_duplicates()
    )
    m1 = ~(ptvs.omim_inheritance_simple.str.contains("AD|XL").fillna(False))
    ptvs = ptvs[m1].sort_values("n_trunc_region", ascending=False)

    logger.info(
        f"Non-morbid genes with >= 3 PTVs in a constrained region:"
        f"{ptvs.symbol.nunique()}"
    )
    return ptvs


def main():
    """Run as script."""

    ptvs = read_dnm_ptvs(_FILE_IN)

    ptvs.to_csv(_FILE_OUT, sep="\t", index=False)

    return ptvs


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
