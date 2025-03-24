import logging
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from config import Config
from data_services import PhenotypeService

# Keep the global constants??? WHY?
# TODO: WHY?
GWAS_PHENO_DIR = Config.GWAS_PHENO_DIR
UKBB_PHENO_DIR = Config.UKBB_PHENO_DIR

SOURCE_PLINK_GENOME = Config.SOURCE_PLINK_GENOME
PLINK_BINARY = Config.PLINK_BINARY
NOMALY_VARIANTS_PATH = (
    Config.NOMALY_VARIANTS_PATH
)  # TODO: Replace with nomaly_data service
logger = logging.getLogger(__name__)


def run_gwas(
    phecode: str,
    ancestry: str,
    phenotype_service: PhenotypeService,
    no_cache: bool = False,
) -> pd.DataFrame:
    """Run GWAS for a phecode."""

    # TODO: Refactor to use the nomaly_data service instead of direct file access

    output_suffix = f"{GWAS_PHENO_DIR}/phecode_{phecode}_ancestry_{ancestry}"
    assoc_file = GWAS_PHENO_DIR / f"{output_suffix}.assoc"
    nomaly_file = GWAS_PHENO_DIR / f"{output_suffix}.assoc_nomaly.tsv"
    fam_file = GWAS_PHENO_DIR / f"{output_suffix}.fam"

    # Return cached results if they exist and no_cache is False
    if nomaly_file.exists() and not no_cache:
        logger.info(f"Loading cached GWAS results for {phecode}")
        return pd.read_csv(nomaly_file, sep="\t")

    # Create FAM file if needed
    if not fam_file.exists() or no_cache:
        create_fam_file(
            fam_file=fam_file,
            phecode=phecode,
            ancestry=ancestry,
            phenotype_service=phenotype_service,
        )

    # Run PLINK if needed
    cmd = f"""{PLINK_BINARY} \
        --bed {SOURCE_PLINK_GENOME}-{ancestry}.bed \
        --bim {SOURCE_PLINK_GENOME}-{ancestry}.bim \
        --fam {fam_file} \
        --assoc fisher \
        --out {output_suffix}"""

    try:
        process = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )

        # Print output in real-time
        if process.stdout:
            while True:
                output = process.stdout.readline()
                if output == "" and process.poll() is not None:
                    break
                if output:
                    logger.info(output.strip())

        # Check for any errors
        if process.returncode != 0 and process.stderr:
            stderr_output = process.stderr.read()
            logger.error(f"PLINK process failed: {stderr_output}")
            raise subprocess.CalledProcessError(process.returncode, cmd)

    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK command failed with error: {str(e)}")
        raise

    # Process results
    assoc = pd.read_csv(f"{output_suffix}.assoc.fisher", sep=r"\s+", dtype={"CHR": str})
    assoc["CHR"] = assoc["CHR"].replace({"23": "X", "24": "Y", "25": "XY", "26": "MT"})
    assoc["CHR_BP_A1_A2"] = assoc.apply(
        lambda x: f"{x.CHR}:{x.BP}_{x.A1}/{x.A2}", axis=1
    )

    # Add variant annotations
    # TODO: Replace with services.nomaly_data.df
    nomaly_variants = pd.read_csv(NOMALY_VARIANTS_PATH, sep="\t")
    assoc = assoc.merge(
        nomaly_variants[["CHR_BP_A1_A2", "gene_id", "nomaly_variant", "RSID"]],
        on="CHR_BP_A1_A2",
        how="left",
    )

    # Save and return results
    result_cols = [
        "nomaly_variant",
        "gene_id",
        "RSID",
        "CHR_BP_A1_A2",
        "F_A",
        "F_U",
        "OR",
        "P",
    ]
    assoc = assoc[result_cols].sort_values("P")
    assoc.to_csv(nomaly_file, sep="\t", index=False)

    return assoc


def create_fam_file(
    fam_file: Path,
    phecode: str,
    ancestry: str,
    phenotype_service: PhenotypeService,
) -> None:
    """Create FAM file for GWAS."""

    # Get phenotype data
    disease_sex = phenotype_service.get_disease_sex_for_phecode(phecode)
    phenotype_data = phenotype_service.get_cases_for_phecode(phecode)

    # Load case information
    logger.info(f"Creating FAM file for {phecode} ({ancestry})")

    # Double check that we don't have cases for the 'wrong' sex...
    if disease_sex == "Both":
        pass
    elif disease_sex == "Female":
        assert np.all(phenotype_data[phenotype_data["sex"] == "M"]["phenotype"] == 9)
    elif disease_sex == "Male":
        assert np.all(phenotype_data[phenotype_data["sex"] == "F"]["phenotype"] == 9)

    # Reformat sex and phenotype columns for PLINK... (1=Female, 2=Male /
    # (1=control, 2=case, -9=missing)
    phenotype_data["sex"] = phenotype_data["sex"].map({"M": 1, "F": 2})
    phenotype_data["phenotype"] = phenotype_data["phenotype"].map({0: 1, 1: 2, 9: -9})

    # Add garbage columns
    phenotype_data["FID"] = phenotype_data["eid"]
    phenotype_data["Father"] = 0
    phenotype_data["Mother"] = 0

    # Create FAM file
    phenotype_data.to_csv(
        fam_file,
        sep=" ",
        header=False,
        index=False,
        columns=["eid", "FID", "Father", "Mother", "sex", "phenotype"],
    )


def format_gwas_results(
    assoc_df: pd.DataFrame, significance_threshold: float = 0.05
) -> list:
    """Format GWAS results for JSON response."""
    if assoc_df.empty:
        return []

    # Filter significant results first
    sig_results = assoc_df[assoc_df["P"] < significance_threshold].copy()

    # Basic column renaming
    sig_results = sig_results.rename(
        columns={"CHR_BP_A1_A2": "Variant", "gene_id": "Gene"}
    )

    # Format RSID as HTML link
    sig_results["RSID"] = sig_results["RSID"].apply(
        lambda x: f'<a href="https://www.ncbi.nlm.nih.gov/snp/{x}">{x}</a>'
        if pd.notna(x)
        else None
    )

    # Replace NaN with None
    sig_results = sig_results.replace({np.nan: None})

    # Convert to records, explicitly replacing NaN with None
    return sig_results.where(pd.notna(sig_results), None).to_dict(orient="records")
