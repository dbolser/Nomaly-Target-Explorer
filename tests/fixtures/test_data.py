import pandas as pd
import numpy as np


def create_test_phecode_data():
    """Create sample phecode data for testing."""
    return pd.DataFrame(
        {
            "phecode": ["250.2", "401.1", "296.2"],
            "description": ["Type 2 diabetes", "Hypertension", "Depression"],
            "sex": ["Both", "Both", "Both"],
            "phecode_group": ["endocrine", "circulatory", "psychiatric"],
            "affected": [1000, 2000, 1500],
            "excluded": [100, 200, 150],
            "phecode_exclude": ["250.1", "401.2", "296.1"],
        }
    )


def create_test_variant_data():
    """Create sample variant data for testing."""
    return pd.DataFrame(
        {
            "variant_id": ["11:69083946:T:C", "9:35066710:A:G"],
            "gene_id": ["GENE1", "GENE2"],
            "consequence": ["missense", "missense"],
            "P": [0.001, 0.005],
            "OR": [1.5, 0.8],
        }
    )


def create_test_icd10_data():
    """Create sample ICD10 coding data for testing."""
    return pd.DataFrame(
        {
            "icd10": ["E11", "I10", "F32"],
            "meaning": ["Type 2 diabetes", "Hypertension", "Depression"],
            "icd10_count": [500, 1000, 750],
        }
    )


MOCK_GWAS_RESULTS = pd.DataFrame(
    {
        "P": [0.01, 0.02, 0.06],
        "gene_id": ["GENE1", "GENE2", "GENE3"],
        "RSID": ["rs1", "rs2", "rs3"],
        "CHR_BP_A1_A2": ["1:100_A/T", "1:200_C/G", "1:300_G/A"],
        "F_A": [0.1, 0.2, 0.3],
        "F_U": [0.2, 0.3, 0.4],
        "OR": [1.5, 0.8, 1.2],
    }
)
