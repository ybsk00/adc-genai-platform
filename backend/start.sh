#!/bin/sh
set -e

# Cloud Run supplies the PORT environment variable.
# Default to 8080 if not set.
PORT=${PORT:-8080}

echo "Starting uvicorn on port $PORT..."
exec uvicorn app.main:app --host 0.0.0.0 --port $PORT
