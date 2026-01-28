"""
Sandbox Library Version Verification Script
Constraint 7: Version mismatch causes build failure
"""
import sys
import importlib.metadata
import json

# Required versions (sync with requirements.lock)
REQUIRED_VERSIONS = {
    "numpy": "1.26.4",
    "scipy": "1.12.0",
    "pandas": "2.1.4",
    "scikit-learn": "1.4.0",
}


def verify_versions():
    """Verify all library versions match requirements"""
    errors = []
    installed = {}

    # Check RDKit separately
    try:
        from rdkit import __version__ as rdkit_version
        installed["rdkit"] = rdkit_version
        # RDKit version format differs, so we check prefix
        if not rdkit_version.startswith("2024"):
            errors.append({
                "package": "rdkit",
                "required": "2024.x.x",
                "actual": rdkit_version
            })
    except ImportError:
        errors.append({
            "package": "rdkit",
            "required": "2024.x.x",
            "actual": "NOT INSTALLED"
        })

    # Check other packages
    for package, required_version in REQUIRED_VERSIONS.items():
        try:
            actual_version = importlib.metadata.version(package)
            installed[package] = actual_version

            if actual_version != required_version:
                errors.append({
                    "package": package,
                    "required": required_version,
                    "actual": actual_version
                })

        except importlib.metadata.PackageNotFoundError:
            errors.append({
                "package": package,
                "required": required_version,
                "actual": "NOT INSTALLED"
            })

    if errors:
        print("=" * 60)
        print("VERSION MISMATCH DETECTED - BUILD BLOCKED")
        print("=" * 60)
        for err in errors:
            print(f"  {err['package']}: required={err['required']}, actual={err['actual']}")
        print("=" * 60)
        print("Fix: Update requirements.lock and rebuild")
        sys.exit(1)

    print("All library versions verified successfully")
    print(json.dumps(installed, indent=2))
    return True


if __name__ == "__main__":
    verify_versions()
