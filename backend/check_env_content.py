import os

def check_env():
    env_path = ".env"
    if not os.path.exists(env_path):
        print(".env not found in current directory")
        return

    print(f"Reading {env_path}...")
    with open(env_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key = line.split("=")[0]
                print(f"Key found: {key}")

if __name__ == "__main__":
    check_env()
