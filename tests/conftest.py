import pytest
from app import app as flask_app
# from blueprints.gwas import GWASTaskManager, GWASTaskStatus


@pytest.fixture
def app():
    flask_app.config["TESTING"] = True
    flask_app.config["LOGIN_DISABLED"] = True  # If using Flask-Login
    return flask_app


@pytest.fixture
def client(app):
    return app.test_client()


# @pytest.fixture
# def mock_gwas_manager(mocker):
#     manager = GWASTaskManager()
#     mocker.patch("blueprints.gwas.variant_level_assoc", return_value=None)
#     return manager
