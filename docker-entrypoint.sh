#!/bin/sh

if [ "$FLASK_ENV" = "production" ]; then
    exec gunicorn wsgi:app --bind 0.0.0.0:8000 --workers 3 --threads 10
else
    exec flask run --host=0.0.0.0 --port=5000
fi
