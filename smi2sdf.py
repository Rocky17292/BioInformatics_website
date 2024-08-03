from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import PyMol
import py3Dmol

# Step 1: Generate 3D coordinates from SMILES
def smiles_to_3D(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # Generate 3D coordinates
    return mol

# Step 2: Convert to PDB format
def mol_to_pdb(mol, filename):
    Chem.MolToPDBFile(mol, filename)

# Example usage
if __name__ == "__main__":
    smiles = "CCO"
    mol = smiles_to_3D(smiles)
    mol_to_pdb(mol, "output.pdb")
