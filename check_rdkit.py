import sys

def check_rdkit():
    print("--- Checking RDKit Installation & Functionality ---")
    try:
        from rdkit import Chem
        print("✅ RDKit imported successfully.")
        
        # Test Valid SMILES (Aspirin)
        valid_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
        mol = Chem.MolFromSmiles(valid_smiles)
        if mol:
            print(f"✅ Valid SMILES check passed: {valid_smiles}")
        else:
            print(f"❌ Valid SMILES check FAILED: {valid_smiles}")

        # Test Invalid SMILES
        invalid_smiles = "InvalidStringHere"
        mol_invalid = Chem.MolFromSmiles(invalid_smiles)
        if mol_invalid is None:
            print(f"✅ Invalid SMILES check passed (returned None): {invalid_smiles}")
        else:
            print(f"❌ Invalid SMILES check FAILED (returned object): {invalid_smiles}")
            
    except ImportError:
        print("❌ RDKit is NOT installed.")
    except Exception as e:
        print(f"❌ RDKit check error: {e}")

if __name__ == "__main__":
    check_rdkit()
