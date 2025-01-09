import pytest
from app import app as flask_app
from db import get_db_connection


@pytest.fixture
def app():
    """Basic Flask app configured for testing."""
    flask_app.config["TESTING"] = True
    flask_app.config["WTF_CSRF_ENABLED"] = False  # Disable CSRF for testing
    return flask_app


@pytest.fixture
def client(app):
    """Basic test client without authentication."""
    return app.test_client()


@pytest.fixture
def test_admin():
    """Create a test admin user.
    This fixture manages the test user lifecycle (create/cleanup).
    """
    # Set up
    conn = get_db_connection()
    cursor = conn.cursor()

    # Create test user
    cursor.execute("""
        INSERT INTO users2 (username, password, email, is_active)
        VALUES ('test_admin', 'test_password', 'test@example.com', TRUE)
    """)
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


@pytest.fixture
def auth_client(client, test_admin):
    """Client with admin authentication by performing a login."""
    # Perform login via POST request
    response = client.post(
        "/login",
        data={"username": test_admin["username"], "password": test_admin["password"]},
        follow_redirects=True,
    )
    assert response.status_code == 200
    assert b"welcome" in response.data.lower()
    return client


@pytest.fixture
def test_limited_user():
    """Create a test user with limited permissions.
    This fixture manages the test user lifecycle (create/cleanup).
    """
    conn = get_db_connection()
    cursor = conn.cursor()

    # Create test user
    cursor.execute("""
        INSERT INTO users2 (username, password, email, is_active)
        VALUES ('test_limited', 'test_password', 'limited@example.com', TRUE)
    """)
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
