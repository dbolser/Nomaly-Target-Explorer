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

    # Configure formatters - use a more detailed one for errors
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)

    # Add a filter to handle exception info for errors and critical logs
    class ErrorStackTraceFilter(logging.Filter):
        def filter(self, record):
            # If this is an ERROR or CRITICAL log and no exc_info is explicitly set
            if record.levelno >= logging.ERROR and record.exc_info is None:
                import sys

                record.exc_info = sys.exc_info()
                # Only add if there's an actual exception
                if record.exc_info == (None, None, None):
                    record.exc_info = None
            return True

    file_handler.addFilter(ErrorStackTraceFilter())

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

    # Also configure the root logger to catch everything else
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)
    root_logger.setLevel(logging.INFO)

    # Return the configured logger
    return app.logger
