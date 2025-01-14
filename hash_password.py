from werkzeug.security import generate_password_hash
from db import get_db_connection
from getpass import getpass


def update_password(username, new_password):
    """Update user's password with a hashed version."""
    # Generate hashed password
    hashed_password = generate_password_hash(new_password)

    # Update in database
    conn = get_db_connection()
    cursor = conn.cursor()

    try:
        cursor.execute(
            "UPDATE users2 SET password = %s WHERE username = %s",
            (hashed_password, username),
        )
        conn.commit()
        print(f"Successfully updated password for user: {username}")
    except Exception as e:
        print(f"Error updating password: {e}")
    finally:
        cursor.close()
        conn.close()


if __name__ == "__main__":
    username = input("Enter username: ")
    new_password = getpass("Enter new password: ")
    update_password(username, new_password)
