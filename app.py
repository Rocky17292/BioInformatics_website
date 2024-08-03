import os
from flask import Flask, render_template, request, send_from_directory
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

app = Flask(__name__)

# Directory for saving images
STATIC_DIR = os.path.join(os.path.dirname(__file__), 'static')

def formula_to_2d_structure(formula, filename):
    try:
        mol = Chem.MolFromSmiles(formula)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        AllChem.Compute2DCoords(mol)
        Draw.MolToFile(mol, os.path.join(STATIC_DIR, filename))
        return True
    except Exception as e:
        return False

def formula_to_3d_structure(formula, filename):
    try:
        mol = Chem.MolFromSmiles(formula)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        # Add explicit hydrogen atoms
        mol = Chem.AddHs(mol)

        # Embed the molecule in 3D space
        AllChem.EmbedMolecule(mol)

        # Optimize the 3D structure
        AllChem.UFFOptimizeMolecule(mol)

        # Draw 2D image with ball and stick representation
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
        opts = drawer.drawOptions()
        opts.atomRadius = 0.3  # Set atom radius
        opts.bondLineWidth = 1.5  # Set bond line width
        opts.fillHighlights = True  # Fill the atoms with color

        # Draw the molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        # Save the image
        img = drawer.GetDrawingText()
        with open(os.path.join(STATIC_DIR, filename), 'wb') as f:
            f.write(img)

        return True
    except Exception as e:
        print("Error generating 3D ball and stick structure:", e)
        return False

@app.route('/')
def index():
    return render_template('index.html')

def smiles_to_mol_structure(smiles, filename):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        # Convert molecule to Mol object and save it
        mol_block = Chem.MolToMolBlock(mol)
        with open(os.path.join(STATIC_DIR, filename), 'w') as f:
            f.write(mol_block)
        return True
    except Exception as e:
        print("Error generating .mol structure:", e)
        return False

def smiles_to_pdb_structure(smiles, filename):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        # Generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

        # Convert to PDB format
        pdb_block = Chem.MolToPDBBlock(mol)

        # Save to file
        with open(os.path.join(STATIC_DIR, filename), 'w') as f:
            f.write(pdb_block)

        return True
    except Exception as e:
        print("Error generating .pdb structure:", e)
        return False

def smiles_to_sdf_structure(smiles, filename):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        # Generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

        # Convert to SDF format
        sdf_block = Chem.MolToMolBlock(mol)

        # Save to file
        with open(os.path.join(STATIC_DIR, filename), 'w') as f:
            f.write(sdf_block)

        return True
    except Exception as e:
        print("Error generating .sdf structure:", e)
        return False

@app.route('/submit', methods=['POST'])
def submit():
    if 'formula_2d' in request.form:
        formula = request.form['formula_2d']
        success = formula_to_2d_structure(formula, '2d_structure.png')
        if success:
            return render_template('result.html')
    elif 'formula_3d' in request.form:
        formula = request.form['formula_3d']
        success = formula_to_3d_structure(formula, '3d_structure.png')
        if success:
            return render_template('result3D.html', filename='3d_structure.png')
    elif 'smiles_to_mol' in request.form:  # Check if form is for .mol structure
        smiles = request.form['smiles_to_mol']
        success = smiles_to_mol_structure(smiles, 'mol_structure.mol')
        if success:
            return render_template('result_mol.html', filename='mol_structure.mol')
    elif 'smiles_to_pdb' in request.form:  # Check if form is for .pdb structure
        smiles = request.form['smiles_to_pdb']
        success = smiles_to_pdb_structure(smiles, 'pdb_structure.pdb')
        if success:
            return render_template('result_pdb.html', filename='pdb_structure.pdb')
    elif 'smiles_to_sdf' in request.form:  # Check if form is for .sdf structure
        smiles = request.form['smiles_to_sdf']
        success = smiles_to_sdf_structure(smiles, 'sdf_structure.sdf')
        if success:
            return render_template('result_sdf.html', filename='sdf_structure.sdf')
    else:
        return "Error: Invalid form submission"

    return "Error generating structure."

@app.route('/upload_pdb', methods=['POST'])
def upload_pdb():
    if 'pdb_file' not in request.files:
        return "No file uploaded"

    pdb_file = request.files['pdb_file']

    # Check if the file is empty
    if pdb_file.filename == '':
        return "No file selected"

    # Read the content of the file
    pdb_data = pdb_file.read()

    # Get the selected visualization type
    visualization_type = request.form['visualization_type']

    # Render the corresponding HTML template based on the visualization type
    if visualization_type == 'wireframe':
        return render_template('wireframe.html', pdb_data=pdb_data)
    elif visualization_type == 'ball_and_stick':
        return render_template('ball_and_stick.html', pdb_data=pdb_data)
    elif visualization_type == 'stick':
        return render_template('stick.html', pdb_data=pdb_data)
    elif visualization_type == 'space_filling':
        return render_template('space_filling.html', pdb_data=pdb_data)
    else:
        return "Error: Unsupported visualization type"

@app.route('/submitMol', methods=['POST'])
def submits():
    smiles = request.form.get('smiles_converter')
    conversion_type = request.form.get('conversion_type')

    if not smiles or not conversion_type:
        return "Error: Invalid form submission"

    filename = f"{conversion_type}_structure.{conversion_type}"

    if conversion_type == 'sdf':
        success = smiles_to_sdf_structure(smiles, filename)
        if success:
            return render_template('result_sdf.html', filename=filename)
    elif conversion_type == 'pdb':
        success = smiles_to_pdb_structure(smiles, filename)
        if success:
            return render_template('result_pdb.html', filename=filename)
    elif conversion_type == 'mol':
        success = smiles_to_mol_structure(smiles, filename)
        if success:
            return render_template('result_mol.html', filename=filename)
    else:
        return "Error: Unsupported conversion type"

    return "Error: Failed to perform conversion"



@app.route('/ngl/ngl.js')
def serve_ngl_script():
    return send_from_directory('static', 'ngl/ngl.js')

if __name__ == '__main__':
    app.run(debug=True)