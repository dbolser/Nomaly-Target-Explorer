# Service Usage Review

## High Priority Updates

### Phenotype Data Access
1. `blueprints/gwas.py`:
   - Currently may be loading phenotype data directly
   - Should use `services.phenotype` with appropriate population filtering
   - GWAS results should consider population context

2. `blueprints/phewas.py`:
   - Similar to GWAS, should use phenotype service
   - PheWAS analysis should respect population settings
   - Consider version compatibility with stats service

3. `db.py`:
   - Functions like `get_phecode_info()` should potentially use phenotype service
   - Population metadata might be available through the service

### Stats Service Integration
1. `blueprints/phecode_term.py`:
   - Should use `services.get_stats_service()` with version/population
   - Term analysis should be consistent with phecode view

2. `blueprints/variant.py`:
   - Check if using latest stats service interface
   - Ensure population context is maintained

### Nomaly Score Service
1. `blueprints/nomaly.py`:
   - Verify using `services.get_nomaly_service()` with correct version
   - QQ plot generation should use consistent data version

## Data Consistency Checks

### Population Handling
1. Review all places using phenotype data:
   - `get_case_counts_for_phecode()`
   - `get_phenotype_data()`
   - Ensure population parameter is passed through

2. Stats calculations:
   - Verify population filtering is applied consistently
   - Check version/population compatibility

### Version Management
1. Review all stats service usage:
   - Replace direct access with `get_stats_service()`
   - Ensure version is passed from UI when needed

2. Check version compatibility:
   - Between phenotype and stats data
   - Between nomaly scores and stats

## Technical Debt

1. Cache directories in config.py:
   - Consider version/population specific caching
   - Review cache invalidation strategy

2. Service initialization:
   - Add validation for version/population combinations
   - Consider logging invalid configurations

3. Error handling:
   - Add specific errors for invalid version/population
   - Improve user feedback for configuration issues

## UI/UX Considerations

1. Population selector:
   - Add to other relevant views
   - Maintain consistency across navigation

2. Version switching:
   - Review impact on cached results
   - Consider adding version indicator in UI

## Testing Needs

1. Add tests for:
   - Population filtering
   - Version switching
   - Service caching
   - Invalid configurations

2. Integration tests:
   - Cross-service data consistency
   - Population/version compatibility

## Documentation Updates Needed

1. Add documentation for:
   - Available populations
   - Version differences
   - Service configuration
   - Population/version compatibility matrix