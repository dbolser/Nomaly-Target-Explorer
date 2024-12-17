import logging
import logging.handlers
import os
from datetime import datetime


def setup_logging(app):
    # Create logs directory if it doesn't exist
    if not os.path.exists("logs"):
        os.makedirs("logs")

    # Configure file handler
    log_file = f'logs/nomaly_{datetime.now().strftime("%Y%m%d")}.log'
    file_handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=10485760, backupCount=10
    )

    # Configure formatters
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)

    # Set up the Flask logger
    app.logger.addHandler(file_handler)
    app.logger.setLevel(logging.INFO)

    # Set up the GWAS logger
    gwas_logger = logging.getLogger("blueprints.gwas")
    gwas_logger.addHandler(file_handler)
    gwas_logger.setLevel(logging.INFO)

    # Set up the Phecode logger
    phecode_logger = logging.getLogger("blueprints.phecode")
    phecode_logger.addHandler(file_handler)
    phecode_logger.setLevel(logging.INFO)

    # Set up the Phecode Term logger
    phecode_term_logger = logging.getLogger("blueprints.phecode_term")
    phecode_term_logger.addHandler(file_handler)
    phecode_term_logger.setLevel(logging.INFO)

    # Return the configured logger
    return app.logger
