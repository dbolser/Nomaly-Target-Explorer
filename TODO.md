


## Migration plan



1. **Data Services First**
```python
class VariantDataService:
    def __init__(self, db_connection, cache_service=None):
        self.db = db_connection
        self.cache = cache_service

    def get_term_variants(self, term):
        # Currently might read from files/HDF5/etc
        # Later will read from DB
        pass

class PrioritizationService:
    def __init__(self, variant_service, scoring_service):
        self.variants = variant_service
        self.scoring = scoring_service
```

Reasons:
- Creates clean boundaries between data access and business logic
- Makes future DB migration easier (just update service implementations)
- Makes Celery integration simpler (services can be passed to tasks)
- Easier to test and refactor in small steps
- Reduces risk by keeping core functionality stable

2. **Database Migration**
```python
class DatabaseVariantService(VariantDataService):
    def get_term_variants(self, term):
        # New DB-based implementation
        return self.db.execute("SELECT ...")
```

Reasons:
- Services already abstract data access
- Can migrate one data source at a time
- Can easily A/B test old vs new implementations
- No impact on task execution logic

3. **Celery Integration**
```python
@celery.task(bind=True)
def prioritize_variants(self, disease_code, term):
    variant_service = get_service()  # Get from DI container
    return variant_service.prioritize(disease_code, term)
```

Reasons:
- Services are already task-ready
- Cache logic is already abstracted
- Can gradually move tasks to background
- Minimal risk to data integrity

This order minimizes disruption because:
1. Each step is independent
2. Changes are isolated
3. Can roll back easily
4. Testing remains straightforward
5. Users see minimal disruption


### More on Task Services...



No, the PrioritizationService is actually crucial! I should have shown its role more clearly. Here's how it all fits together:

```python
# 1. Core Service Layer
class PrioritizationService:
    def __init__(self, variant_service, scoring_service, cache_service=None):
        self.variants = variant_service
        self.scoring = scoring_service
        self.cache = cache_service

    def prioritize_variants(self, disease_code, term):
        # Business logic stays here
        if self.cache and (cached := self.cache.get(disease_code, term)):
            return cached
            
        variants = self.variants.get_term_variants(term)
        scores = self.scoring.get_variant_scores(variants)
        result = self._calculate_priorities(variants, scores)
        
        if self.cache:
            self.cache.set(disease_code, term, result)
        return result

# 2. Celery Task Layer
@celery.task(bind=True)
def prioritize_variants_task(self, disease_code, term):
    # Task layer is thin, just coordinates services
    service = get_prioritization_service()  # Get from DI container
    return service.prioritize_variants(disease_code, term)

# 3. Flask Route Layer
@app.route("/prioritize/<disease>/<term>")
def start_prioritization(disease, term):
    # API layer starts task and returns task ID
    task = prioritize_variants_task.delay(disease, term)
    return {"task_id": task.id}
```

The benefits of keeping the service:
1. Business logic stays in one place
2. Task layer just handles Celery-specific concerns
3. Services are easier to test than tasks
4. Can run logic directly or via Celery
5. Clear separation of concerns

Without the service layer, you'd end up with business logic in your tasks, making them harder to test and maintain.


## Celery Cache vs. Task Cache



Celery's result backend (caching) operates at the task level, which is a different layer than our business logic caching. Here's how they can work together:

```python
# 1. Service Layer (Business Logic Cache)
class PrioritizationService:
    def __init__(self, variant_service, scoring_service, cache_service=None):
        self.variants = variant_service
        self.scoring = scoring_service
        self.cache = cache_service  # Redis/Memcached for computed results

    def prioritize_variants(self, disease_code, term):
        cache_key = f"variants:{disease_code}:{term}"
        if self.cache and (cached := self.cache.get(cache_key)):
            return cached
        
        # ... computation ...
        return result

# 2. Celery Task (Task-Level Cache)
@celery.task(bind=True)
@app.cache.memoize(timeout=3600)  # Celery result backend cache
def prioritize_variants_task(self, disease_code, term):
    """
    - Celery cache stores task results (success/failure)
    - Useful for task status checks
    - Handles retry logic
    - Stores task metadata
    """
    service = get_prioritization_service()
    return service.prioritize_variants(disease_code, term)

# 3. API Layer
@app.route("/prioritize/<disease>/<term>")
def start_prioritization(disease, term):
    # Check if task result exists in Celery's cache
    task = prioritize_variants_task.delay(disease, term)
    if task.ready():  # Result in Celery cache
        return task.get()
    return {"task_id": task.id}
```

The two caching layers serve different purposes:
- Service cache: Fast access to computed results
- Celery cache: Task state management and results



## About Celery

Celery's result backend is about task management rather than data caching. It helps with:

1. Task Pipeline Management:
```python
@celery.task
def prioritize_variants_task(disease_code, term):
    # Task status tracking
    # Retry logic
    # Error handling
    return result

@celery.task
def notify_user_task(result):
    # Can check if previous task succeeded
    # Access task metadata
    # Handle failures

# Chain tasks
workflow = (
    prioritize_variants_task.s(disease_code, term) |
    notify_user_task.s()
)
```

2. Task State Inspection:
```python
@app.route("/status/<task_id>")
def check_status(task_id):
    task = prioritize_variants_task.AsyncResult(task_id)
    return {
        "state": task.state,          # PENDING, STARTED, SUCCESS, FAILURE
        "progress": task.info,        # Custom status updates
        "result": task.get() if task.ready() else None
    }
```

Your actual data caching (Redis/Memcached) remains in the service layer where it belongs.


## The Flask 'factory pattern' for configuration

The Flask Factory Pattern is about creating the application object dynamically rather than globally. Here's how it works:

```python
# Without Factory (current approach)
from flask import Flask
app = Flask(__name__)
app.config['DATABASE_URI'] = 'sqlite:///example.db'

# With Factory Pattern
def create_app(config_name='default'):
    app = Flask(__name__)
    
    # Load different configs based on environment
    if config_name == 'development':
        app.config.from_object('config.DevelopmentConfig')
    elif config_name == 'testing':
        app.config.from_object('config.TestingConfig')
    elif config_name == 'production':
        app.config.from_object('config.ProductionConfig')

    # Initialize extensions
    db.init_app(app)
    cache.init_app(app)
    
    # Register blueprints
    from .blueprints.prioritisation import prioritisation_bp
    app.register_blueprint(prioritisation_bp)
    
    return app
```

The config classes might look like:
```python
# config.py
class Config:
    """Base config."""
    VARIANT_SCORES_DIR = '/data/scores'
    CACHE_TYPE = 'redis'

class DevelopmentConfig(Config):
    DEBUG = True
    CACHE_DIR = '/tmp/dev-cache'

class TestingConfig(Config):
    TESTING = True
    CACHE_DIR = '/tmp/test-cache'
    # Use mock services
    VARIANT_SERVICE = 'mock'

class ProductionConfig(Config):
    CACHE_DIR = '/data/cache'
    # Use real services
    VARIANT_SERVICE = 'real'
```

This makes testing much cleaner:
```python
# tests/conftest.py
@pytest.fixture
def app():
    app = create_app('testing')
    return app

# Can also test different configs
def test_development_config():
    app = create_app('development')
    assert app.config['DEBUG']
    assert not app.config['TESTING']
```

Benefits:
1. Different configs for different environments
2. Easier testing
3. Better dependency management
4. More modular application structure

