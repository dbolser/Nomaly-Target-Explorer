from flask import render_template, Blueprint

disease_select = Blueprint("disease", __name__, template_folder="../templates")


# Route for page1 with tables
@disease_select.route("/page1")
def render_page1():
    search_categories = {
        "Skin Conditions": [
            "Hidradenitis",
            "Rosacea",
            "Ehlers",
        ],
        "Women's Health": [
            "Polycystic Ovarian",
            "Endometriosis",
        ],
    }
    return render_template("page1.html", search_categories=search_categories)
