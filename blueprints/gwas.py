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

    output_prefix = f"phecode_{phecode}"
    output_path = GWAS_PHENO_DIR / output_prefix
    assoc_path = output_path / f"{output_path}.assoc"
    nomaly_path = output_path / f"{assoc_path}_nomaly.tsv"

    # Return cached results if they exist
    if os.path.exists(nomaly_path):
        logger.info(f"Loading cached GWAS results for {phecode}")
        return pd.read_csv(nomaly_path, sep="\t")

    # Load case information
    logger.info(f"Running new GWAS for {phecode}")
    with open(
        UKBB_PHENO_DIR / "phecode_cases_excludes" / f"phecode_{phecode}.pkl", "rb"
    ) as f:
        cases = pickle.load(f)

    # Create FAM file if needed
    if not os.path.exists(f"{output_path}.fam"):
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

        fam.to_csv(f"{output_path}.fam", sep=" ", header=False, index=False)

    # Run PLINK if needed
    if not os.path.exists(assoc_path):
        cmd = f"{PLINK_BINARY} --allow-no-sex --bed {SOURCE_PLINK_GENOME}.bed --bim {SOURCE_PLINK_GENOME}.bim --fam {output_path}.fam --assoc --out {output_path} --silent"
        subprocess.run(cmd, shell=True, check=True)

    # Process results
    assoc = pd.read_csv(assoc_path, sep=r"\s+", dtype={"CHR": str})
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
    assoc.to_csv(nomaly_path, sep="\t", index=False)

    return assoc


def format_gwas_results(
    assoc_df: pd.DataFrame, significance_threshold: float = 0.05
) -> list:
    """Format GWAS results for JSON response.

    This function handles all formatting of GWAS data including:
    - Converting numeric columns to proper types
    - Handling missing values
    - Formatting RSID links
    - Renaming columns for display
    """
    if assoc_df.empty:
        return []

    # Make a copy to avoid modifying the original
    formatted_df = assoc_df.copy()

    # Ensure numeric columns are float type
    numeric_cols = ["P", "OR", "F_A", "F_U"]
    for col in numeric_cols:
        formatted_df[col] = pd.to_numeric(formatted_df[col], errors="coerce")

    # Format for display
    formatted_df = formatted_df.rename(
        columns={
            "CHR_BP_A1_A2": "Variant",
            "gene_id": "Gene",
        }
    )

    # Handle RSID links
    formatted_df["RSID"] = formatted_df["RSID"].apply(
        lambda x: f'<a href="https://www.ncbi.nlm.nih.gov/snp/{x}">{x}</a>'
        if pd.notna(x)
        else None
    )

    # Convert numeric columns to float and replace NaN with None
    for col in numeric_cols:
        formatted_df[col] = formatted_df[col].astype(float).replace({np.nan: None})

    sig_results = formatted_df[formatted_df["P"] < significance_threshold].copy()

    return sig_results.to_dict(orient="records")
