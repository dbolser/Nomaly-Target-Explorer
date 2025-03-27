import logging
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from config import Config
from data_services import PhenotypeService, NomalyDataService

# Keep the global constants??? WHY?
# TODO: WHY?
GWAS_PHENO_DIR = Config.GWAS_PHENO_DIR
UKBB_PHENO_DIR = Config.UKBB_PHENO_DIR

logger = logging.getLogger(__name__)


def run_gwas(
    phecode: str,
    ancestry: str,
    phenotype_service: PhenotypeService,
    nomaly_data_service: NomalyDataService,
    no_cache: bool = False,
) -> pd.DataFrame:
    """Run GWAS for a phecode."""

    file_stem = f"phecode_{phecode}_ancestry_{ancestry}"
    nom_file = GWAS_PHENO_DIR / f"{file_stem}.assoc.nom.tsv"
    fam_file = GWAS_PHENO_DIR / f"{file_stem}.fam"

    # Return cached results if they exist and no_cache is False
    if nom_file.exists() and not no_cache:
        logger.info(
            f"Loading cached GWAS results for {phecode} ({ancestry}) from {nom_file}"
        )
        return pd.read_csv(nom_file, sep="\t")

    # Create FAM file if needed
    if not fam_file.exists() or no_cache:
        create_fam_file(fam_file, phecode, ancestry, phenotype_service)

    # Run PLINK
    bed = Config.GENOTYPES_BED
    genotypes_bed = Path(f"{bed.parent}/{bed.stem}-{ancestry}.bed")
    genotypes_bim = Path(f"{bed.parent}/{bed.stem}-{ancestry}.bim")

    cmd = f"""{Config.PLINK_BINARY} \
        --bed {genotypes_bed} \
        --bim {genotypes_bim} \
        --fam {fam_file} \
        --assoc fisher \
        --out {fam_file.with_suffix("")}"""

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
    assoc = pd.read_csv(
        f"{fam_file.with_suffix('.assoc.fisher')}", sep=r"\s+", dtype={"CHR": str}
    )
    assoc["CHR"] = assoc["CHR"].replace({"23": "X", "24": "Y", "25": "XY", "26": "MT"})
    assoc["CHR_BP"] = assoc.apply(lambda x: f"{x.CHR}:{x.BP}", axis=1)
    assoc["CHR_BP_A1_A2"] = assoc.apply(lambda x: f"{x.CHR_BP}_{x.A1}/{x.A2}", axis=1)

    # Add variant annotations
    print(f"NOMALY VARIANTS: {nomaly_data_service.colapsed_df.shape}")
    nomaly_variants = nomaly_data_service.colapsed_df

    merged = assoc.merge(
        nomaly_variants[["CHR_BP", "gene_id", "nomaly_variant"]], how="left"
    )

    # NOTE: We get 80k assoiations from GWAS, but we only have 30k
    #       variants in the 'nomaly_variants' dataframe
    assert len(merged) == len(assoc)

    # Save and return results
    result_cols = [
        "CHR",
        "BP",
        "A1",
        "A2",
        "CHR_BP",
        "CHR_BP_A1_A2",
        "SNP",
        "F_A",
        "F_U",
        "P",
        "OR",
        "gene_id",
        "nomaly_variant",
    ]
    merged = merged[result_cols].sort_values("P")
    merged.to_csv(nom_file, sep="\t", index=False)

    # TODO: NO CLUE WHY RETURNING merged HERE FAILS WITH failed to json serialsei.set...
    return pd.read_csv(nom_file, sep="\t")


def create_fam_file(
    fam_file: Path,
    phecode: str,
    ancestry: str,
    phenotype_service: PhenotypeService,
) -> None:
    """Create FAM file for GWAS."""

    # Get phenotype data
    phenotype_data = phenotype_service.get_cases_for_phecode(phecode, ancestry)

    # Load case information
    logger.info(f"Creating FAM file for {phecode} ({ancestry}) at {fam_file}")

    # Reformat sex and phenotype columns for PLINK...

    # 1=Male, 2=Female
    phenotype_data["sex"] = phenotype_data["sex"].map({"M": 1, "F": 2})

    # 0=control, 1=case, 9=missing
    phenotype_data["phenotype"] = phenotype_data["phenotype"].map({0: 1, 1: 2, 9: -9})

    # Add FID, Father, Mother columns
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
    assoc_df: pd.DataFrame, significance_threshold: float = 1
) -> list:
    """Format GWAS results for JSON response."""

    # Filter significant results first
    sig_results = assoc_df[assoc_df["P"] < significance_threshold].copy()

    # Basic column renaming
    sig_results = sig_results.rename(
        columns={"CHR_BP_A1_A2": "Variant", "gene_id": "Gene"}
    )

    # Format RSID as HTML link
    # TODO: This should probably be done in the frontend
    sig_results["RSID"] = sig_results["SNP"].apply(
        lambda x: f'<a href="https://www.ncbi.nlm.nih.gov/snp/{x}">{x}</a>'
        if pd.notna(x)
        else None
    )

    # Replace NaN with None for JSON serialization
    sig_results = sig_results.replace({np.nan: None})
    sig_results = sig_results.where(pd.notna(sig_results), None)

    sig_results["Gene"] = sig_results["Gene"].apply(
        lambda x: list(x) if isinstance(x, list) else x
    )

    # Convert to records, explicitly replacing NaN with None
    return sig_results.to_dict(orient="records")


def main():
    phenotype_service = PhenotypeService(Config.PHENOTYPES_HDF)
    nomaly_data_service = NomalyDataService(Config.NOMALY_VARIANTS_PATH)

    phecode = "290"
    phecode = "290.11"
    phecode = "722.7"

    print(f"RUNNING GWAS FOR {phecode}")
    # result = run_gwas(
    #     phecode, "EUR", phenotype_service, nomaly_data_service, no_cache=True
    # )

    result = run_gwas("290", "EUR", phenotype_service, nomaly_data_service)

    result = format_gwas_results(result)

    import json

    print(json.dumps(result[:6], indent=4))

    x = json.dumps(result, indent=4)

    print("Done")


if __name__ == "__main__":
    main()