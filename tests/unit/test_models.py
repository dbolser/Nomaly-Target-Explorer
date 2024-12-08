import pytest
import pandas as pd
import numpy as np
from db import get_all_phecodes, get_db_connection

def test_get_all_phecodes(mocker):
    # Mock data that matches the expected database return
    mock_data = [
        ('123.1', 'Test Description 1', 'Both', 'Test Group 1'),
        ('456.2', 'Test Description 2', 'Female', 'Test Group 2'),
    ]
    
    # Mock the database connection and cursor
    mock_conn = mocker.patch('mysql.connector.connect')
    mock_cursor = mocker.Mock()
    mock_conn.return_value.cursor.return_value = mock_cursor
    
    # Set up the cursor mock to return our test data
    mock_cursor.fetchall.return_value = mock_data
    mock_cursor.description = [
        ('phecode', None), 
        ('description', None),
        ('sex', None),
        ('phecode_group', None)
    ]

    # Call the function
    result = get_all_phecodes()

    # Assert the result is a pandas DataFrame with expected columns and data
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ['phecode', 'description', 'sex', 'phecode_group']
    assert len(result) == 2
    assert result.iloc[0]['phecode'] == '123.1' 