import pandas as pd

from config import Config
from data_services.genotype import GenotypeService


def main():
    gs = GenotypeService(hdf5_file=Config.GENOTYPES_HDF)

    print(gs.get_genotypes().shape)  # (83011, 487950)

    # NOTE: REF/ALT is defined relative to the order given in the
    # 'plink_variant_id'. Check the corresponding 'nomaly_variant_id' to see if
    # the alleles are flipped in the genotype file relative to the nominal
    # variant ID. What you do with that information is apparently so obvious
    # only an idiot would consider it.

    nomaly_variant_id = gs.nomaly_variant_ids
    plink_variant_id = gs.plink_variant_ids

    flip_me = nomaly_variant_id != plink_variant_id

    genotype_counts = gs.get_genotype_counts_and_freqs()

    het_count = genotype_counts["het_count"]
    ref_count = genotype_counts["ref_count"]
    alt_count = genotype_counts["alt_count"]

    ref_allele_count = (2 * ref_count) + het_count
    alt_allele_count = (2 * alt_count) + het_count

    freq_01 = genotype_counts["het_freq"]
    freq_00 = genotype_counts["ref_freq"]
    freq_11 = genotype_counts["alt_freq"]

    # Create a dataframe
    df = pd.DataFrame(
        {
            "genotype_variant_id": plink_variant_id,
            "nomaly_variant_id": nomaly_variant_id,
            "plink_variant_id": plink_variant_id,
            "htrz": het_count,
            "ref_hmoz": ref_count,
            "alt_hmoz": alt_count,
            "ref_allele_count": ref_allele_count,
            "alt_allele_count": alt_allele_count,
            "freq_01": freq_01,
            "freq_00": freq_00,
            "freq_11": freq_11,
            "flip_me": flip_me,
        }
    )

    df.to_csv("genotype_counts.tsv", sep="\t", index=False)

    # For reference, here is a MySQL table to store the above:
    #
    """
    CREATE TABLE genotype_counts (
        genotype_variant_id VARCHAR(255) NOT NULL,
        nomaly_variant_id VARCHAR(255) NOT NULL,
        plink_variant_id VARCHAR(255) NOT NULL,
        htrz INT UNSIGNED NOT NULL,
        ref_hmoz INT UNSIGNED NOT NULL,
        alt_hmoz INT UNSIGNED NOT NULL,
        ref_allele_count INT UNSIGNED NOT NULL,
        alt_allele_count INT UNSIGNED NOT NULL,
        freq_01 FLOAT NOT NULL,
        freq_00 FLOAT NOT NULL,
        freq_11 FLOAT NOT NULL,
        flip_me BOOLEAN NULL,
        PRIMARY KEY (genotype_variant_id),
        INDEX (nomaly_variant_id),
        INDEX (plink_variant_id)
    );
    """

    # And here is the SQL command to load it:
    """
    LOAD DATA LOCAL INFILE '/home/danbolser/Work/nomaly-disease-browser/genotype_counts.tsv'
    INTO TABLE genotype_counts
    FIELDS TERMINATED BY '\t'
    LINES TERMINATED BY '\n'
    IGNORE 1 LINES
    """

    # Now you can do something like this if you want to:
    """
    query="
    SELECT
      nomaly_variant_id,
      hmm_score,
      htrz,
      ref_hmoz,
      alt_hmoz,
      freq_01, freq_00, freq_11,
      CAST(hmm_score*hmm_score * freq_01 + hmm_score*hmm_score * 4 * freq_11 AS FLOAT) AS vs00,
      CAST(hmm_score*hmm_score * freq_01 + hmm_score*hmm_score * 4 * freq_00 AS FLOAT) AS vs11,
      CAST(hmm_score*hmm_score * (freq_00 + freq_11)                         AS FLOAT) AS vs01
    FROM
      hmm_score
    INNER JOIN
      variants USING (variant_id)
    INNER JOIN
      genotype_counts USING (plink_variant_id)
    #LIMIT
    #  10
    "
    mysql ukbb -e "$query" > genotype_counts_with_vs.tsv
    
    """


if __name__ == "__main__":
    main()
