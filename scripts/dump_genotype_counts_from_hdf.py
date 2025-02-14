from config import Config

from data_services.genotype import GenotypesHDF5

import pandas as pd
import numpy as np


def main():
    genotypes = GenotypesHDF5(Config.GENOTYPES_H5)
    genotypes.allele_flipped_in_genotype_file_relative_to_nomaly_variant_id

    print(genotypes.genotype_matrix.shape)  # (83011, 487950)

    print(genotypes.genotype_counts.shape)  # (83011, 3)

    # Somehow we just know that
    #  - col 0 = genotype 0 (REF/REF or 00)
    #  - col 1 = genotype 1 (REF/ALT or 01)
    #  - col 2 = genotype 2 (ALT/ALT or 11)

    # NOTE: REF/ALT is defined relative to the order given in the
    # 'genotype_variant_id'. Check the conveniently named
    # 'allele_flipped_in_genotype_file_relative_to_nomaly_variant_id' array to
    # see if the alleles are flipped in the genotype file relative to the
    # nominal variant ID. What you do with that information is apparently so
    # obvious only an idiot would consider it.

    # Lets always present things in the order determined by the NOMALY VARIANT
    # ID. It at least has the advantage of history, and using the reference
    # genome would be too easy.

    genotype_variant_id = genotypes.genotype_variant_id
    nomaly_variant_id = genotypes.nomaly_variant_id
    plink_variant_id = genotypes.plink_variant_id
    genotype_counts = genotypes.genotype_counts

    htrz = genotype_counts[:, 1]
    ref_hmoz = genotype_counts[:, 0]
    alt_hmoz = genotype_counts[:, 2]
    total = genotype_counts.sum(axis=1)

    ref_allele_count = (2 * ref_hmoz) + htrz
    alt_allele_count = (2 * alt_hmoz) + htrz

    # NOTE: We ignore missing alleles when calculating frequencies
    freq_01 = htrz / total
    freq_00 = ref_hmoz / total
    freq_11 = alt_hmoz / total

    assert np.allclose(freq_01 + freq_00 + freq_11, 1)

    # Create a dataframe
    df = pd.DataFrame(
        {
            "genotype_variant_id": genotype_variant_id,
            "nomaly_variant_id": nomaly_variant_id,
            "plink_variant_id": plink_variant_id,
            "htrz": htrz,
            "ref_hmoz": ref_hmoz,
            "alt_hmoz": alt_hmoz,
            "ref_allele_count": ref_allele_count,
            "alt_allele_count": alt_allele_count,
            "freq_01": freq_01,
            "freq_00": freq_00,
            "freq_11": freq_11,
            "flip_me": genotypes.allele_flipped_in_genotype_file_relative_to_nomaly_variant_id,
        }
    )

    # Convert None to \N in flip_me column for MySQL compatibility
    df["flip_me"] = df["flip_me"].replace({None: r"\N", True: 1, False: 0})

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
