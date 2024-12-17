import pytest
from app import app as flask_app
from db import get_db_connection
# from blueprints.gwas import GWASTaskManager, GWASTaskStatus


@pytest.fixture
def app():
    flask_app.config["TESTING"] = True
    flask_app.config["LOGIN_DISABLED"] = True  # If using Flask-Login
    return flask_app


@pytest.fixture
def client(app):
    return app.test_client()


@pytest.fixture
def auth_client(client, test_user):
    """Client with authentication."""
    client.post("/login", data={"username": "test_admin", "password": "test_password"})
    return client


@pytest.fixture
def test_user(app):
    """Create a test admin user with full permissions."""
    conn = get_db_connection()
    cursor = conn.cursor()

    # Create test user
    cursor.execute("""
        INSERT INTO users2 (username, password, email, is_active)
        VALUES ('test_admin', 'test_password', 'test@example.com', TRUE)
    """)
    user_id = cursor.lastrowid

    # Create root page entry
    cursor.execute("""
        INSERT IGNORE INTO pages (path, description)
        VALUES ('/*', 'Root access for testing')
    """)
    page_id = cursor.lastrowid or 1

    # Grant admin permissions
    cursor.execute(
        """
        INSERT INTO user_permissions (user_id, page_id, wildcard_path)
        VALUES (%s, %s, '/*')
    """,
        (user_id, page_id),
    )

    conn.commit()
    cursor.close()
    conn.close()

    yield

    # Cleanup after tests
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute("DELETE FROM user_permissions WHERE user_id = %s", (user_id,))
    #cursor.execute("DELETE FROM pages WHERE id = %s", (page_id,))
    cursor.execute("DELETE FROM users2 WHERE id = %s", (user_id,))
    conn.commit()
    cursor.close()
    conn.close()


@pytest.fixture
def limited_user(app):
    """Create a test user with limited permissions (e.g., only phecode access)."""
    conn = get_db_connection()
    cursor = conn.cursor()

    # Create test user
    cursor.execute("""
        INSERT INTO users2 (username, password, email, is_active)
        VALUES ('test_limited', 'test_password', 'limited@example.com', TRUE)
    """)
    user_id = cursor.lastrowid

    # Create phecode page entry
    cursor.execute("""
        INSERT INTO pages (path, description)
        VALUES ('/phecode/*', 'Phecode access for testing')
    """)
    page_id = cursor.lastrowid

    # Grant limited permissions
    cursor.execute(
        """
        INSERT INTO user_permissions (user_id, page_id, wildcard_path)
        VALUES (%s, %s, '/phecode/*')
    """,
        (user_id, page_id),
    )

    conn.commit()
    cursor.close()
    conn.close()

    yield user_id  # You can use the user_id in tests if needed

    # Cleanup after tests
    # ... similar cleanup code ...
