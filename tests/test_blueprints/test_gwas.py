import pandas as pd
import pytest
from unittest.mock import MagicMock, patch

from blueprints.gwas import create_fam_file

# Define mock phenotype data
MOCK_PHENOTYPE_DATA = pd.DataFrame(
    {
        "eid": [101, 102, 103, 104, 105],
        "phenotype": [1, 0, 1, 9, 0],  # 1=case, 0=control, 9=missing
    }
)

# Define sample source FAM content

SOURCE_FAM_CONTENT = """\
101 101 0 0 1 0
102 102 0 0 2 0
103 103 0 0 1 0
104 104 0 0 2 0
105 105 0 0 1 0
"""

SOURCE_FAM_CONTENT_WITH_EXTRA_SAMPLES = """\
101 101 0 0 1 2
102 102 0 0 2 1
103 103 0 0 1 2
104 104 0 0 2 -9
105 105 0 0 1 1
106 106 0 0 2 -9
"""

# Define expected output FAM content for both of the above
# phenotype_PLINK maps {0: 1, 1: 2, 9: -9}. Missing eids get -9.
EXPECTED_FAM_CONTENT = """\
101 101 0 0 1 2
102 102 0 0 2 1
103 103 0 0 1 2
104 104 0 0 2 -9
105 105 0 0 1 1
"""

SOURCE_FAM_CONTENT_WITH_MISSING_SAMPLES = """\
102 102 0 0 2 1
103 103 0 0 1 2
105 105 0 0 1 1
"""

EXPECTED_FAM_CONTENT_WITH_MISSING_SAMPLES = """\
102 102 0 0 2 1
103 103 0 0 1 2
105 105 0 0 1 1
"""


@pytest.mark.parametrize(
    "source_and_expected_fam_content",
    [
        (SOURCE_FAM_CONTENT, EXPECTED_FAM_CONTENT),
        (SOURCE_FAM_CONTENT_WITH_EXTRA_SAMPLES, EXPECTED_FAM_CONTENT),
        (
            SOURCE_FAM_CONTENT_WITH_MISSING_SAMPLES,
            EXPECTED_FAM_CONTENT_WITH_MISSING_SAMPLES,
        ),
    ],
)
def test_create_fam_file(tmp_path, source_and_expected_fam_content):
    """Test the creation of the FAM file for PLINK input."""

    source_fam_content, expected_fam_content = source_and_expected_fam_content

    # Both of thse are ignored in our mock setup
    phecode = "123.4"
    ancestry = "XYZ"

    output_fam_file = tmp_path / f"phecode_{phecode}_ancestry_{ancestry}.fam"

    # The source fam file needs to follow the naming convention defined in
    # create_fam_file. It derives the name from Config.GENOTYPES_FAM and adds
    # the ancestry suffix. We mock Config.GENOTYPES_FAM to control this base
    # path within tmp_path!
    base_genotypes_path = tmp_path / "dummy_genotypes.fam"
    source_fam_file = (
        base_genotypes_path.parent
        / f"{base_genotypes_path.stem}-{ancestry}{base_genotypes_path.suffix}"
    )

    # Write the dummy source FAM file
    source_fam_file.parent.mkdir(parents=True, exist_ok=True)
    source_fam_file.write_text(source_fam_content)

    # Mock PhenotypeService
    mock_phenotype_service = MagicMock()
    mock_phenotype_service.get_cases_for_phecode.return_value = MOCK_PHENOTYPE_DATA

    # Mock Config.GENOTYPES_FAM path using patch. We patch
    # 'blueprints.gwas.Config' as it's imported and used there
    with patch("blueprints.gwas.Config") as mock_config:
        # Point Config.GENOTYPES_FAM to the base path within our tmp_path
        mock_config.GENOTYPES_FAM = base_genotypes_path

        # --- Execution ---
        create_fam_file(
            fam_file=output_fam_file,
            phecode=phecode,
            ancestry=ancestry,
            phenotype_service=mock_phenotype_service,
        )

    # Check if PhenotypeService was called correctly
    mock_phenotype_service.get_cases_for_phecode.assert_called_once_with(
        phecode, ancestry
    )

    # Check if the output file was created
    assert output_fam_file.exists()

    # Check the content of the output file matches what we expect. Read both
    # expected and created, normalizing whitespace for comparison...
    created_content = output_fam_file.read_text().strip()
    normalized_expected = "\n".join(
        " ".join(line.split()) for line in expected_fam_content.strip().splitlines()
    )

    # For easier debugging if it fails, compare line by line
    created_lines = created_content.splitlines()
    expected_lines = normalized_expected.splitlines()

    assert len(created_lines) == len(expected_lines), (
        f"Number of lines differ. Got {len(created_lines)}, Expected {len(expected_lines)}"
    )

    for i, (c_line, e_line) in enumerate(zip(created_lines, expected_lines)):
        assert c_line == e_line, (
            f"Line {i + 1} differs. Got '{c_line}', Expected '{e_line}'"
        )

    # Final check on the full content
    assert created_content == normalized_expected
