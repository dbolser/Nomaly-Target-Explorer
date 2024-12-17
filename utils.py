from contextlib import contextmanager
from flask import current_app


@contextmanager
def app_context():
    """Context manager that provides app context"""
    if current_app:
        # If we're already in an app context, just yield
        yield
    else:
        # Otherwise, create a new app context
        from app import app  # Import here to avoid circular imports

        with app.app_context():
            yield
