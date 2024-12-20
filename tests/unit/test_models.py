import pandas as pd
from db import get_all_phecodes

# TODO: Think about what we actually want to test here.


def test_get_all_phecodes():
    phecodes = get_all_phecodes()
    assert len(phecodes) > 0
    assert isinstance(phecodes, pd.DataFrame)
    assert 'phecode' in phecodes.columns
    assert 'description' in phecodes.columns
