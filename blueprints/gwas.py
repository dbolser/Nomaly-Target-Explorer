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

    # TODO: Migrate to using plink2?

    bfile = f"{bed.parent}/{bed.stem}-{ancestry}"
    cmd2 = f"""plink2 \
        --glm sex allow-no-covars --1 \
        --bfile {bfile} \
        --pheno {fam_file} \
        --pheno-col-nums 6 \
        --require-pheno \
        --out {fam_file.with_suffix("")}-plink2
    """
    logger.info(f"We could be running PLINK2: {cmd2}")

    try:
        run_subprocess(cmd)
        # run_subprocess(cmd2)
    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK command failed with error: {str(e)}")
        raise
    except ValueError:
        raise
    except Exception as e:
        logger.error(f"PLINK command failed with error: {str(e)}")
        raise

    # Process results
    assoc = pd.read_csv(
        f"{fam_file.with_suffix('.assoc.fisher')}", sep=r"\s+", dtype={"CHR": str}
    )
    assoc["CHR"] = assoc["CHR"].replace({"23": "X", "24": "Y", "25": "XY", "26": "MT"})
    assoc["CHR_BP"] = assoc["CHR"].astype(str) + ":" + assoc["BP"].astype(str)
    assoc["CHR_BP_A1_A2"] = (
        assoc["CHR_BP"].astype(str)
        + "_"
        + assoc["A1"].astype(str)
        + "/"
        + assoc["A2"].astype(str)
    )

    # Add variant annotations
    nomaly_variant_df = nomaly_data_service.colapsed_df
    logger.info(f"NOMALY VARIANTS: {nomaly_variant_df.shape}")

    merged = assoc.merge(nomaly_variant_df, on="CHR_BP", how="left")

    # NOTE: We get 80k assoiations from GWAS, but we only have 30k
    #       variants in the 'nomaly_variants' dataframe
    assert len(merged) == len(assoc)

    # Squash gene_id to a single string here, before we save to cache
    merged["gene_id"] = merged["gene_id"].apply(
        lambda x: ", ".join(x) if isinstance(x, set) else x
    )

    merged = merged.sort_values("P")
    merged.to_csv(nom_file, sep="\t", index=False)

    return merged


def run_subprocess(cmd: str) -> None:
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
        if process.returncode == 7:
            logger.info(f"Hopefully forgivable... {process.returncode}")
            raise ValueError(f"Hopefully forgivable... {process.returncode}")
        else:
            raise subprocess.CalledProcessError(process.returncode, cmd)


def create_fam_file(
    fam_file: Path,
    phecode: str,
    ancestry: str,
    phenotype_service: PhenotypeService,
) -> None:
    """Create FAM file for GWAS."""

    logger.info(f"Creating FAM file for {phecode} ({ancestry}) at {fam_file}")

    # Get phenotype data
    phenotype_df = phenotype_service.get_cases_for_phecode(phecode, ancestry)
    phenotype_df = phenotype_df[["eid", "phenotype"]]

    # Map the phenotype to the correct values for PLINK
    # 0=control, 1=case, 9=missing
    phenotype_df["phenotype"] = phenotype_df["phenotype"].map({0: 1, 1: 2, 9: -9})

    # Read the original FAM file to get the correct order of the samples
    genotypes_fam = Config.GENOTYPES_FAM

    # NOTE: The ancestry specficific FAM files were created 'manually' using plink
    ancestry_fam = (
        genotypes_fam.parent / f"{genotypes_fam.stem}-{ancestry}{genotypes_fam.suffix}"
    )
    logger.info(f"Reading FAM file for ({ancestry}) from {ancestry_fam}")

    # NOTE: None of these columns 'conflict' prior to merging
    fam_file_columns = ["FID", "IID", "Father", "Mother", "sex", "phenotype"]
    fam_df = pd.read_csv(
        ancestry_fam,
        sep=r"\s+",
        header=None,
        names=fam_file_columns,
        # Drop the original 'random' phenotype column
        usecols=["FID", "IID", "Father", "Mother", "sex"],
    )

    # Merge in the phecode specific phenotype
    merged_df = fam_df.merge(phenotype_df, left_on="IID", right_on="eid", how="left")

    assert len(merged_df) == len(fam_df)
    assert np.all(merged_df["IID"] == fam_df["IID"])

    if np.any(merged_df.phenotype.isna()):
        logger.warning(
            f"Some samples in {ancestry_fam} are not returned from the phenotype service"
            f"{len(merged_df[merged_df.phenotype.isna()])} samples of {len(fam_df)} are missing"
        )
        merged_df = merged_df[merged_df.phenotype.notna()]
        merged_df["phenotype"] = merged_df["phenotype"].astype(int)

    # Write the merged FAM file
    merged_df.to_csv(
        fam_file,
        sep=" ",
        header=False,
        index=False,
        columns=fam_file_columns,
    )


def format_gwas_results(
    assoc_df: pd.DataFrame, significance_threshold: float = 1
) -> pd.DataFrame:
    """Format GWAS results for JSON response."""

    # Filter significant results first
    if significance_threshold < 1:
        assoc_df = assoc_df[assoc_df["P"] < significance_threshold]

    # Basic column renaming
    assoc_df = assoc_df.rename(columns={"CHR_BP_A1_A2": "Variant", "gene_id": "Gene"})

    # Format RSID as HTML link. TODO: This should probably be done in the frontend
    assoc_df["RSID"] = np.where(
        pd.notna(assoc_df["SNP"]),
        '<a href="https://www.ncbi.nlm.nih.gov/snp/'
        + assoc_df["SNP"]
        + '">'
        + assoc_df["SNP"]
        + "</a>",
        "",
    )

    # Replace NaN with 1 for OR
    assoc_df["OR"] = assoc_df["OR"].replace(np.nan, 1)

    # Replace 'isna' with None for nomaly_variant and Gene
    assoc_df["nomaly_variant"] = assoc_df["nomaly_variant"].replace(np.nan, None)
    assoc_df["Gene"] = assoc_df["Gene"].replace(np.nan, None)

    return assoc_df


def main():
    phenotype_service = PhenotypeService(Config.PHENOTYPES_HDF)
    nomaly_data_service = NomalyDataService(Config.NOMALY_VARIANT_MAPPING_PATH)

    # import sys

    # phecode = sys.argv[1]
    # ancestry = sys.argv[2]

    phecode = "290.11"
    ancestry = "EUR"

    print(f"Running GWAS for {phecode} ({ancestry})")

    run_gwas(phecode, ancestry, phenotype_service, nomaly_data_service, no_cache=True)

    exit(0)

    # Run some random GWAS for fun...
    ancestries = np.array(["AFR", "EAS", "EUR", "SAS"])
    phecodes = phenotype_service.phecodes

    done = set()
    for _ in range(10000):
        phecode = np.random.choice(phecodes)
        ancestry = np.random.choice(ancestries)
        if (phecode, ancestry) in done:
            continue
        done.add((phecode, ancestry))
        print(f"Running GWAS for {phecode} ({ancestry})")
        try:
            run_gwas(
                phecode, ancestry, phenotype_service, nomaly_data_service, no_cache=True
            )
            print(f"Successfully ran GWAS for {phecode} ({ancestry})")
        except Exception as e:
            print(f"Error running GWAS for {phecode} ({ancestry}): {e}")
            continue


if __name__ == "__main__":
    main()