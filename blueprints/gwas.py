import os
import pickle
import subprocess
import pandas as pd
import logging
import numpy as np

# Keep the global constants
GWAS_PHENO_DIR = "/data/clu/ukbb/by_pheno/"
UKBB_PHENO_DIR = "/data/general/UKBB/Phenotypes/"
source_plink_genome = "/data/clu/ukbb/genotypes_nomaly"
plink_binary = "/data/clu/ukbb/plink"

logger = logging.getLogger(__name__)


def run_gwas(phecode: str) -> pd.DataFrame:
    """Run GWAS for a phecode if not already done and return results."""

    output_prefix = f"phecode_{phecode}"
    output_path = f"{GWAS_PHENO_DIR}{output_prefix}"
    assoc_path = f"{output_path}.assoc"
    nomaly_path = f"{assoc_path}_nomaly.tsv"

    # Return cached results if they exist
    if os.path.exists(nomaly_path):
        logger.info(f"Loading cached GWAS results for {phecode}")
        return pd.read_csv(nomaly_path, sep="\t")

    # Load case information
    logger.info(f"Running new GWAS for {phecode}")
    with open(
        f"{UKBB_PHENO_DIR}phecode_cases_excludes/phecode_{phecode}.pkl", "rb"
    ) as f:
        cases = pickle.load(f)

    # Create FAM file if needed
    if not os.path.exists(f"{output_path}.fam"):
        fam = pd.read_csv(f"{source_plink_genome}.fam", header=None, sep=r"\s+")
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
        cmd = f"{plink_binary} --bed {source_plink_genome}.bed --bim {source_plink_genome}.bim --fam {output_path}.fam --assoc --out {output_path} --silent"
        subprocess.run(cmd, shell=True, check=True)

    # Process results
    assoc = pd.read_csv(assoc_path, sep=r"\s+", dtype={"CHR": str})
    assoc["CHR"] = assoc["CHR"].replace({"23": "X", "24": "Y", "25": "XY", "26": "MT"})
    assoc["CHR_BP_A1_A2"] = assoc.apply(
        lambda x: f"{x.CHR}:{x.BP}_{x.A1}/{x.A2}", axis=1
    )

    # Add variant annotations
    nomaly_variants = pd.read_csv("/data/clu/ukbb/nomaly_variants.tsv", sep="\t")
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


def format_gwas_results(assoc_df: pd.DataFrame) -> list:
    """Format GWAS results for JSON response."""
    if assoc_df.empty:
        return []

    # Filter significant results
    sig_results = assoc_df[assoc_df["P"] < 0.05].copy()

    # Format for display
    sig_results = sig_results.rename(
        columns={"CHR_BP_A1_A2": "Variant", "gene_id": "Gene"}
    )

    # Handle RSID links
    sig_results["RSID"] = sig_results["RSID"].apply(
        lambda x: f'<a href="https://www.ncbi.nlm.nih.gov/snp/{x}">{x}</a>'
        if pd.notna(x)
        else None
    )

    # Convert numeric columns to float and replace NaN with None
    numeric_cols = ["F_A", "F_U", "OR", "P"]
    for col in numeric_cols:
        # First convert to float, then replace NaN with None
        sig_results[col] = sig_results[col].astype(float).replace({np.nan: None})

    return sig_results.to_dict(orient="records")
