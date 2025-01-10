import os
import pickle
import subprocess
import pandas as pd
import logging
import numpy as np

from config import Config

# Keep the global constants
GWAS_PHENO_DIR = Config.GWAS_PHENO_DIR
UKBB_PHENO_DIR = Config.UKBB_PHENO_DIR

SOURCE_PLINK_GENOME = Config.SOURCE_PLINK_GENOME
PLINK_BINARY = Config.PLINK_BINARY
NOMALY_VARIANTS_PATH = Config.NOMALY_VARIANTS_PATH
logger = logging.getLogger(__name__)


def run_gwas(phecode: str) -> pd.DataFrame:
    """Run GWAS for a phecode if not already done and return results."""

    output_suffix = f"{GWAS_PHENO_DIR}/phecode_{phecode}"
    assoc_file = GWAS_PHENO_DIR / f"{output_suffix}.assoc"
    nomaly_file = GWAS_PHENO_DIR / f"{output_suffix}.assoc_nomaly.tsv"
    fam_file = GWAS_PHENO_DIR / f"{output_suffix}.fam"

    # Return cached results if they exist
    if nomaly_file.exists():
        logger.info(f"Loading cached GWAS results for {phecode}")
        return pd.read_csv(nomaly_file, sep="\t")

    # Create FAM file if needed
    if not fam_file.exists():
        # Load case information
        logger.info(f"Running new GWAS for {phecode}")
        with open(
            UKBB_PHENO_DIR / "phecode_cases_excludes" / f"phecode_{phecode}.pkl", "rb"
        ) as f:
            cases = pickle.load(f)

        fam = pd.read_csv(f"{SOURCE_PLINK_GENOME}.fam", header=None, sep=r"\s+")
        fam.columns = ["FID", "IID", "Father", "Mother", "sex", "phenotype"]

        # Set phenotypes (1=control, 2=case, -9=missing)
        fam["phenotype"] = 1
        fam.loc[fam["IID"].isin(cases["cases"]), "phenotype"] = 2

        # Handle sex-specific cases
        if cases["Sex"] == "Female":
            fam.loc[fam["sex"] == 1, "phenotype"] = -9
        elif cases["Sex"] == "Male":
            fam.loc[fam["sex"] == 2, "phenotype"] = -9

        # Handle exclusions
        if cases["exclude"]:
            fam.loc[fam["IID"].isin(cases["exclude"]), "phenotype"] = -9

        fam.to_csv(fam_file, sep=" ", header=False, index=False)

    # Run PLINK if needed
    cmd = f"{PLINK_BINARY} --allow-no-sex --bed {SOURCE_PLINK_GENOME}.bed --bim {SOURCE_PLINK_GENOME}.bim --fam {fam_file} --assoc --out {output_suffix}"

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
    assoc = pd.read_csv(f"{output_suffix}.assoc", sep=r"\s+", dtype={"CHR": str})
    assoc["CHR"] = assoc["CHR"].replace({"23": "X", "24": "Y", "25": "XY", "26": "MT"})
    assoc["CHR_BP_A1_A2"] = assoc.apply(
        lambda x: f"{x.CHR}:{x.BP}_{x.A1}/{x.A2}", axis=1
    )

    # Add variant annotations
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
