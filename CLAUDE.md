# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Setup and Environment
```bash
# Initial setup
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env

# Alternative with UV (preferred)
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

### Running the Application
```bash
export FLASK_APP=app.py
export FLASK_ENV=development
flask run
```

### Testing
```bash
# Run all tests
pytest

# Run specific test categories
pytest tests/unit/
pytest tests/integration/
pytest tests/test_blueprints/

# Run with coverage
pytest --cov

# Run a single test file
pytest tests/test_blueprints/test_variant_scores.py
```

## Architecture Overview

### Flask Application Structure
- **Application Factory Pattern**: `create_app()` in `app.py` with environment-based configuration
- **Blueprint Organization**: Modular structure with 11 blueprints handling different functional areas
- **Service Layer Pattern**: All data access through dedicated service classes in `data_services/`
- **ServiceRegistry**: Central registry managing all data services with dependency injection

### Key Blueprints
- `prioritisation_by_nomaly_scores.py`: Variant prioritization (core functionality)
- `network_analysis.py`: Causal analysis using NetworkX/bnlearn/pgmpy
- `phecode.py` & `phecode_term.py`: Disease phenotype browsing
- `gwas.py`: GWAS analysis functionality
- `variant.py`: Individual variant analysis

### Data Services Architecture
The application uses a service registry pattern with these key services:
- **GenotypeService**: HDF5-based genotype data (488K × 83K, 2.2GB)
- **NomalyScoreService**: Nomaly scores HDF5 data (486K × 13K, 3.5GB)
- **PhenotypeService**: Phenotype and phecode data management
- **StatsRegistry**: Population-specific statistical analysis services
- **NomalyDataService**: Variant-to-gene mapping data

### Database Schema
- **Authentication**: `users2` table with path-based permissions system
- **Genomic Data**: Variants, genes, individuals with population and sex metadata
- **Phenotype Data**: ICD10, phecodes, and phenotype mappings

## Data File Structure
Large genomic datasets are stored in HDF5 format under `/data/general/UKBB/`:
- Genotypes: `Genotypes/GRCh38/genotypes-ukbb.h5`
- Nomaly Scores: `Run-v1/DatabaseInputs/float16_scores.h5`
- Statistics: Population-specific files (EUR, AFR, SAS, EAS, ALL)
- Cached results: `/data/danbolser/ukbb/cache/` with organized subdirectories

## Configuration Management
- Environment configs in `config.py` (development, testing, production)
- Database credentials and data paths in `.env` file
- HDF5 file paths configured per environment in config classes

## Authentication & Permissions
- Flask-Login for session management
- Path-based permissions system: users granted access to specific phecodes or `*` for admin
- Admin user creation via SQL commands in README.md

## Error Handling
- Custom exceptions: `DatabaseConnectionError`, `DataNotFoundError`
- Global error handlers for JSON/HTML responses
- Structured error templates in `templates/errors/`

## Testing Strategy
- 173 tests across unit, integration, and blueprint categories
- Mock-based unit testing with fixtures in `tests/fixtures/`
- Integration tests for HDF5 data services
- Performance benchmarks in `tests/benchmarks/`

## Performance Considerations
- HDF5 for efficient large dataset access
- Structured caching system for computed results
- Service registry for connection pooling and resource management
- Pandas/NumPy optimized data processing pipelines