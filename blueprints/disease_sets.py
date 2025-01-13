from flask import render_template, Blueprint

disease_select = Blueprint("disease", __name__, template_folder="../templates")


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


@disease_select.route("/page-diabetes")
def render_diabetes():
    search_categories = {
        "Diabetes Primary Conditions": [
            "Diabetes",
            "Type 1 diabetes",
            "Type 2 diabetes",
        ],
        "Diabetes Complications": [
            "Retinopathy",
            "Nephropathy",
            "Neuropathy",
        ],
    }
    return render_template("page1.html", search_categories=search_categories)


@disease_select.route("/page-dementia")
def render_dementia():
    search_categories = {
        "Dementia Cognitive Disorders": [
            "Dementia",
            "Alzheimer",
            "Cognitive decline",
        ],
        "Dementia Related Conditions": [
            "Depression",
            "Anxiety",
        ],
    }
    return render_template("page1.html", search_categories=search_categories)
