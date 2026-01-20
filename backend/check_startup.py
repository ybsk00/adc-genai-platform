import sys
import os

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

print("Checking imports...")
try:
    from app import main
    print("SUCCESS: 'app.main' imported successfully.")
except ImportError as e:
    print(f"FAILURE: Import failed. {e}")
    sys.exit(1)
except Exception as e:
    print(f"FAILURE: Startup crashed. {e}")
    sys.exit(1)
