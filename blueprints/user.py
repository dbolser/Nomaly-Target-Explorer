from flask import Blueprint, render_template, abort
from flask_login import current_user, login_required
from db import get_db_connection

user_bp = Blueprint("user", __name__, template_folder="../templates")


@user_bp.route("/my-phecodes")
@login_required
def my_phecodes():
    """Show phecodes that the current user has permission to view."""
    if not current_user.is_authenticated:
        abort(401)

    conn = get_db_connection()
    categories = {}

    try:
        # First check if user is admin
        cursor = conn.cursor()
        cursor.execute("SELECT is_admin FROM users2 WHERE id = %s", (current_user.id,))
        user_result = cursor.fetchone()
        is_admin = user_result[0] if user_result else False

        # Get user's permissions
        cursor.execute(
            "SELECT allowed_paths FROM user_permissions WHERE user_id = %s",
            (current_user.id,),
        )
        perm_result = cursor.fetchone()
        allowed_paths = perm_result[0] if perm_result else ""
        cursor.close()

        # Get phecodes with dictionary results
        cursor = conn.cursor(dictionary=True)
        allowed_phecodes = []

        if is_admin:
            # Admin users can see all phecodes
            cursor.execute(
                """
                SELECT DISTINCT phecode, description, phecode_group, sex 
                FROM phecode_definition
                ORDER BY phecode
                """
            )
            allowed_phecodes = cursor.fetchall()
        elif allowed_paths:
            # Non-admin users can only see phecodes they have permission for
            phecode_list = allowed_paths.split(",")

            # Get all the base codes (like '250')
            base_codes = [code for code in phecode_list if "." not in code]

            # Get all the specific codes (like '250.1')
            specific_codes = [code for code in phecode_list if "." in code]

            # Build the query conditions
            conditions = []
            params = []

            # Add explicit matches for specific codes
            if specific_codes:
                placeholders = ", ".join(["%s"] * len(specific_codes))
                conditions.append(f"phecode IN ({placeholders})")
                params.extend(specific_codes)

            # Add prefix matches for base codes
            for base_code in base_codes:
                conditions.append("phecode = %s OR phecode LIKE %s")
                params.extend([base_code, f"{base_code}.%"])

            if conditions:
                query = f"""
                    SELECT DISTINCT phecode, description, phecode_group, sex 
                    FROM phecode_definition
                    WHERE {" OR ".join(conditions)}
                    ORDER BY phecode
                """
                cursor.execute(query, params)
                allowed_phecodes = cursor.fetchall()

        # Organize phecodes by category
        for phecode_data in allowed_phecodes:
            category = phecode_data["phecode_group"]
            if category not in categories:
                categories[category] = []
            categories[category].append(phecode_data)

        cursor.close()

    except Exception as e:
        print(f"Error retrieving user's phecodes: {str(e)}")
    finally:
        conn.close()

    return render_template("user/my_phecodes.html", categories=categories)
