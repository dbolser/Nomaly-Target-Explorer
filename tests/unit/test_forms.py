from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField
from wtforms.validators import DataRequired

# Test form class (if you have forms in your application)
class LoginForm(FlaskForm):
    username = StringField('Username', validators=[DataRequired()])
    password = PasswordField('Password', validators=[DataRequired()])

def test_login_form_validation(app):
    # Add secret key configuration for CSRF
    app.config['SECRET_KEY'] = 'test-secret-key'
    app.config['WTF_CSRF_SECRET_KEY'] = 'test-csrf-secret-key'  # Optional but can be used for specific CSRF secret

    with app.test_request_context():
        form = LoginForm()
        assert not form.validate()
        
        # Test with data
        form = LoginForm(data={
            'username': 'testuser',
            'password': 'testpass',
            'csrf_token': form.csrf_token.current_token
        })
        
        # Print validation errors for debugging
        if not form.validate():
            print("Form errors:", form.errors)
            
        assert form.validate()
        
        # Test with missing data
        form = LoginForm(data={
            'username': '',
            'password': 'testpass'
        })
        assert not form.validate()
        assert 'This field is required.' in form.username.errors 