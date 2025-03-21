from flask import Blueprint, current_app, jsonify


from data_services import StatsService

# Create a Blueprint for stats-related routes
stats_bp = Blueprint("stats", __name__, url_prefix="/api/stats")


@stats_bp.route("/<run_version>/<ancestry>/<term>/<phecode>")
def get_stats(run_version, ancestry, term, phecode):
    print(f"Getting stats for {run_version}, {ancestry}, {term}, {phecode}")

    try:
        registry = current_app.extensions["nomaly_services"]

        # Get the stats service for the requested run_version and ancestry
        stats_service: StatsService = registry.stats_registry.get(run_version, ancestry)

        # Use the service to fetch data

        # You can use either of these two methods, they're functionally
        # equivalent when only passing one term and one phecode
        # stats_data = stats_service.get_term_stats(term, phecode=phecode)
        stats_data = stats_service.get_phecode_stats(phecode, term=term)

        return jsonify(stats_data.to_dict(orient="records"))

    except ValueError as e:
        return jsonify({"error": str(e)}), 404


# Defining a main function for simple debugging
def main():
    from config import Config
    from data_services import StatsRegistry

    # Is this really the best we can do?
    stats_registry = StatsRegistry(Config.STATS_SELECTOR)
    stats_service = stats_registry.get("Run-v1", "EUR")

    old = stats_service.get_stats_by_phecode("332")
    new = stats_service.get_phecode_stats("332")
    assert old.equals(new)
    print(new)

    old = stats_service.get_stats_by_phecode("332", stats_type=["num_rp"])
    new = stats_service.get_phecode_stats("332", stats_types=["num_rp"])
    assert old.equals(new)
    print(new)

    old = stats_service.get_stats_by_phecode("332", stats_type=["num_rp", "num_rn"])
    new = stats_service.get_phecode_stats("332", stats_types=["num_rp", "num_rn"])
    assert old.equals(new)
    print(new)

    # The old method doesn't support this...
    new = stats_service.get_phecode_stats("332", term="MP:0005179")
    print(new)

    new = stats_service.get_phecode_stats("332", terms=["MP:0005179", "GO:0006915"])
    print(new)

    stats_data = stats_service.get_term_stats("MP:0005179")
    print(stats_data)

    stats_data = stats_service.get_term_stats("MP:0005179", phecode="332")
    print(stats_data)

    stats_data = stats_service.get_term_stats("MP:0005179", phecodes=["290", "332"])
    print(stats_data)

    print("horribly done")


if __name__ == "__main__":
    main()