"""Production CLI for PheWAS analysis."""


# TODO: Don't use the app, get the services via config!
import click
from app import create_app


@click.command()
@click.option("--variant", help="Variant ID to analyze")
@click.option("--phecode", help="Optional phecode to analyze")
@click.option("--full-run", is_flag=True, help="Run analysis on all variants")
def run_phewas(variant, phecode, full_run):
    """Run PheWAS analysis from command line."""
    app = create_app("production")
    with app.app_context():
        if full_run:
            from blueprints.phewas import run_full_analysis

            run_full_analysis()
        else:
            from blueprints.phewas import get_phewas_results

            results = get_phewas_results(variant, phecode)
            print(results)


if __name__ == "__main__":
    run_phewas()
