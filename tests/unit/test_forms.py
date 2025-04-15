from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField
from wtforms.validators import DataRequired


# Test form class (if you have forms in your application)
class LoginForm(FlaskForm):
    username = StringField("Username", validators=[DataRequired()])
    password = PasswordField("Password", validators=[DataRequired()])


def test_login_form_validation(unit_test_app):
    """Test form validation with CSRF enabled."""
    with unit_test_app.test_request_context():
        form = LoginForm()
        assert not form.validate()

        # Test with valid data
        form = LoginForm(
            data={
                "username": "testuser",
                "password": "testpass",
                # Er... not test config or something?
                #"csrf_token": form.csrf_token.current_token,
            }
        )
        assert form.validate(), f"Form validation failed: {form.errors}"

        # Test with missing data
        form = LoginForm(data={"username": "", "password": "testpass"})
        assert not form.validate()
        assert "This field is required." in form.username.errors
