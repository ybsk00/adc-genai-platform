import re

def test_smiles_regex():
    print("\n--- [3] Ambeed SMILES Regex Test ---")
    mock_body_text = """
    Product Information
    Catalog No. : A123456
    CAS No. : 12345-67-8
    SMILES Code : CC1=C(C(=O)N(C1=O)C)C2=CC=C(C=C2)C3=CC=C(C=C3)C4=CC=C(C=C4)S(=O)(=O)N
    Molecular Formula : C25H21N3O4S
    """
    
    # 정규식 패턴 (크롤러에 적용한 것과 동일)
    pattern = r"SMILES Code\s*:\s*([A-Za-z0-9@#\(\)\[\]\/\ \=+-]+)"
    
    match = re.search(pattern, mock_body_text)
    if match:
        smiles = match.group(1).strip()
        print(f"   ✅ Success: Found SMILES via Regex: {smiles}")
    else:
        print("   ❌ Error: SMILES not found via Regex")

if __name__ == "__main__":
    test_smiles_regex()
