from contextlib import contextmanager
from flask import current_app, has_app_context


@contextmanager
def app_context():
    """Context manager that provides app context"""
    if has_app_context():
        # If we're already in an app context, just yield
        yield
    else:
        # If no app context exists, raise an error
        raise RuntimeError(
            "No application context found. "
            "This operation requires an application context. "
            "Please ensure this is called within a Flask app context."
        )
