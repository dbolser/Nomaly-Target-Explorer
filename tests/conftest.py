import os
from unittest.mock import patch
import pytest
from werkzeug.security import generate_password_hash

from app import create_app
from blueprints.nomaly_services import services
from db import get_db_connection


@pytest.fixture(scope="session")
def app():
    """Create and configure a Flask app instance for unit tests."""
    _app = create_app(config_name="testing")  # Uses UnitTestConfig
    yield _app


@pytest.fixture(scope="session")
def integration_app():
    """Create and configure a Flask app instance for integration tests."""
    _app = create_app(config_name="integration")  # Uses IntegrationTestConfig
    yield _app


@pytest.fixture(scope="session")
def hdf5_integration(integration_app):
    """Fixture for tests that need real HDF5 services.
    Use this fixture ONLY in tests that absolutely need real HDF5 files."""
    yield services


@pytest.fixture
def client(app):
    """Provide a test client for unit tests."""
    return app.test_client()


@pytest.fixture
def integration_client(integration_app):
    """Provide a test client for integration tests."""
    return integration_app.test_client()


@pytest.fixture
def test_admin():
    """Create a test admin user with full permissions."""
    conn = get_db_connection()
    cursor = conn.cursor()

    # Hash the password
    hashed_password = generate_password_hash("test_password")

    # Create test admin user
    cursor.execute(
        """
        INSERT INTO users2 (username, password, email, is_active)
        VALUES ('test_admin', %s, 'test_admin@example.com', TRUE)
    """,
        (hashed_password,),
    )
    user_id = cursor.lastrowid

    # Grant admin permissions
    cursor.execute(
        """
        INSERT INTO user_permissions (user_id, allowed_paths)
        VALUES (%s, '*')
    """,
        (user_id,),
    )

    conn.commit()
    cursor.close()
    conn.close()

    yield {"id": user_id, "username": "test_admin", "password": "test_password"}

    # Cleanup
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute("DELETE FROM user_permissions WHERE user_id = %s", (user_id,))
    cursor.execute("DELETE FROM users2 WHERE id = %s", (user_id,))
    conn.commit()
    cursor.close()
    conn.close()


def cleanup_test_admin_after_test_timeout():
    """Cleanup the test admin user."""
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute(
        """
        DELETE
            user_permissions, users2
        FROM
            user_permissions 
        INNER JOIN
            users2
        ON
          user_permissions.user_id = users2.id
        WHERE
            username = "test_admin"
        """
    )
    conn.commit()
    cursor.close()
    conn.close()


@pytest.fixture
def auth_client(client, test_admin):
    """Authenticate the client as the test admin user."""
    response = client.post(
        "/login",
        data={"username": test_admin["username"], "password": "test_password"},
        follow_redirects=True,
    )
    assert response.status_code == 200
    assert b"welcome" in response.data.lower()
    return client


@pytest.fixture
def test_limited_user():
    """Create a test user with limited permissions."""
    conn = get_db_connection()
    cursor = conn.cursor()

    # Hash the password
    hashed_password = generate_password_hash("test_password")

    # Create test user
    cursor.execute(
        """
        INSERT INTO users2 (username, password, email, is_active)
        VALUES ('test_limited', %s, 'test_limited@example.com', TRUE)
    """,
        (hashed_password,),
    )
    user_id = cursor.lastrowid

    # Grant limited permissions
    cursor.execute(
        """
        INSERT INTO user_permissions (user_id, allowed_paths)
        VALUES (%s, '250.2,649.1,561')
    """,
        (user_id,),
    )

    conn.commit()
    cursor.close()
    conn.close()

    yield {
        "id": user_id,
        "username": "test_limited",
        "password": "test_password",
        "allowed_paths": ["250.2", "649.1", "561"],
    }

    # Cleanup
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute("DELETE FROM user_permissions WHERE user_id = %s", (user_id,))
    cursor.execute("DELETE FROM users2 WHERE id = %s", (user_id,))
    conn.commit()
    cursor.close()
    conn.close()


@pytest.fixture
def mock_config(tmp_path):
    """Create a temporary directory for cache testing"""
    with patch("config.Config") as mock_config:
        mock_config.VARIANT_SCORES_DIR = str(tmp_path)
        yield mock_config


@pytest.fixture(scope="module")
def db_connection(integration_app):
    """Provide a database connection for the test module.
    This fixture automatically provides an app context."""
    with integration_app.app_context():
        conn = get_db_connection()
        yield conn
        conn.close()  # Connection is closed after the test module finishes


@pytest.fixture
def auth_integration_client(integration_client, test_admin):
    """Authenticate the integration client as the test admin user."""
    response = integration_client.post(
        "/login",
        data={"username": test_admin["username"], "password": "test_password"},
        follow_redirects=True,
    )
    assert response.status_code == 200
    assert b"welcome" in response.data.lower()
    return integration_client


def main():
    # Disaster recovery...
    cleanup_test_admin_after_test_timeout()
    exit(0)


if __name__ == "__main__":
    main()
