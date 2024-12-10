import json
from flask import url_for

def test_index_route(client):
    """Test the index route redirects to login when not authenticated."""
    response = client.get('/')
    assert response.status_code == 302  # Redirect to login
    assert '/login' in response.location

def test_login_route(client):
    """Test login functionality."""
    # Test GET request
    response = client.get('/login')
    assert response.status_code == 200
    assert b'login' in response.data.lower()

    # Test POST request with valid credentials
    response = client.post('/login', data={
        'username': 'testuser',
        'password': 'testpass'
    })
    assert response.status_code == 302  # Redirect to index
    assert '/' == response.location

def test_search_route(client):
    """Test the disease search functionality."""
    response = client.get('/diseasesearch?q=diabetes')
    assert response.status_code == 200
    data = json.loads(response.data)
    assert isinstance(data, list)
    # Verify the structure of returned data
    if len(data) > 0:
        assert all(key in data[0] for key in [
            'phecode', 'description', 'sex', 'affected', 
            'excluded', 'phecode_exclude', 'phecode_group'
        ])

def test_phecode_route(client):
    """Test the phecode detail route."""
    test_phecode = '250.2'  # Example phecode
    response = client.get(f'/phecode/{test_phecode}')
    assert response.status_code == 200
    assert bytes(test_phecode, 'utf-8') in response.data 