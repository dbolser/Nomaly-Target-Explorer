from flask import render_template, Blueprint, abort

disease_sets_bp = Blueprint("disease_sets", __name__, template_folder="../templates")

CATEGORY_SETS = {
    "set1": {
        "Skin Conditions": [
            "Hidradenitis",
            "Rosacea",
            "Ehlers",
        ],
    },
    "set2": {
        "Women's Health": [
            "Polycystic Ovarian",
            "Endometriosis",
        ],
    },
    "set3": {
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
    },
    "set4": {
        "Dementia Cognitive Disorders": [
            "Dementia",
            "Alzheimer",
            "Cognitive decline",
        ],
        "Dementia Related Conditions": [
            "Depression",
            "Anxiety",
        ],
    },
}


@disease_sets_bp.route("/disease-sets/<set_id>")
def show_set(set_id):
    if set_id not in CATEGORY_SETS:
        abort(404)
    return render_template("disease_sets.html", search_categories=CATEGORY_SETS[set_id])
