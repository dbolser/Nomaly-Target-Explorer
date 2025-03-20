from flask import Blueprint, current_app, jsonify

# Create a Blueprint for stats-related routes
stats_bp = Blueprint("stats", __name__, url_prefix="/api/stats")


@stats_bp.route("/<run_version>/<ancestry>/<term>/<phecode>")
def get_stats(run_version, ancestry, term, phecode):
    print(f"Getting stats for {run_version}, {ancestry}, {term}, {phecode}")
    try:
        registry = current_app.extensions["nomaly_services"]

        # Get the stats service for the requested run_version and ancestry
        stats_service = registry.stats_registry.get(run_version, ancestry)

        # Use the service to fetch data
        stats_data = stats_service.get_stats_by_term_phecode(term, phecode)

        return jsonify(stats_data)
    except ValueError as e:
        return jsonify({"error": str(e)}), 404
