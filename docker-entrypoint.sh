#!/bin/sh

if [ "$FLASK_ENV" = "development" ]; then
    exec gunicorn app:app --bind 0.0.0.0:5000 --threads 10
else
    exec flask run --host=0.0.0.0 --port=5000
fi
