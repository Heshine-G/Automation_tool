import os
import sys
import subprocess
from Bio import PDB
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue
from rdkit import Chem
from rdkit.Chem import AllChem


class NonStandardAminoAcidProcessor:
    def __init__(self, input_file):
        self.input_file = input_file
        self.parser = PDB.PDBParser(QUIET=True)
        self.io = PDBIO()
        self.standard_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL'
        }
        self.output_dir = "non_standard_residues"

    def create_output_directory(self):
        """Create output directory if it doesn't exist."""
        os.makedirs(self.output_dir, exist_ok=True)

    def is_amino_acid(self, residue):
        """Check if the residue has the basic amino acid backbone atoms (N, CA, C)."""
        backbone_atoms = {'N', 'CA', 'C'}
        residue_atoms = {atom.get_name() for atom in residue}
        return backbone_atoms.issubset(residue_atoms)

    def extract_non_standard_residues(self):
        """Extract non-standard residues and save them in separate directories."""
        structure = self.parser.get_structure('protein', self.input_file)
        non_standard_dirs = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() not in self.standard_residues and self.is_amino_acid(residue):
                        residue_name = f"{residue.get_resname()}_{residue.get_id()[1]}"
                        residue_dir = os.path.join(self.output_dir, residue_name)
                        os.makedirs(residue_dir, exist_ok=True)  # Create unique folder per residue

                        new_structure = Structure.Structure('non_standard')
                        new_model = Model.Model(0)
                        new_chain = Chain.Chain('A')
                        new_chain.add(residue.copy())
                        new_model.add(new_chain)
                        new_structure.add(new_model)

                        output_file = os.path.join(residue_dir, "residue.pdb")  # Standardized filename
                        self.io.set_structure(new_structure)
                        self.io.save(output_file)
                        non_standard_dirs.append(residue_dir)

        return non_standard_dirs

    def add_terminal_groups_with_pymol(self, residue_folder):
        """Run capping.py on a given residue folder."""
        try:
            script_path = os.path.join(os.path.dirname(__file__), "capping.py")
            print(f"Running: python {script_path} {residue_folder}")

            result = subprocess.run(
                [sys.executable, script_path, residue_folder],  # Use sys.executable for dynamic Python path
                capture_output=True,
                text=True,
                check=True
            )

            capped_file = os.path.join(residue_folder, "residue_capped.pdb")  # Corrected name
            return capped_file if os.path.exists(capped_file) else None
        except subprocess.CalledProcessError as e:
            print(f"Error running capping.py: {e}")
            return None

    def add_hydrogens(self, residue_folder):
        """Protonate the PDB structure using RDKit."""
        try:
            input_file = os.path.join(residue_folder, "residue_capped.pdb")  # Use capped file
            output_file = os.path.join(residue_folder, "residue_capped_protonated.pdb")

            mol = Chem.MolFromPDBFile(input_file, removeHs=False)
            if mol is None:
                raise ValueError(f"Error: Unable to parse {input_file} with RDKit.")

            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)

            Chem.MolToPDBFile(mol, output_file)
            return output_file if os.path.exists(output_file) else None
        except Exception as e:
            print(f"Error adding hydrogens: {str(e)}")
            return None

    def add_partial_charges(self, residue_folder):
        """Add Gasteiger partial charges using RDKit and convert to MOL2 format."""
        try:
            input_file = os.path.join(residue_folder, "residue_capped_protonated.pdb")  # Corrected filename
            sdf_file = os.path.join(residue_folder, "residue_capped_protonated_charged.sdf")  # Intermediate SDF file
            mol2_file = os.path.join(residue_folder, "residue_capped_protonated_charged.mol2")  # Final MOL2 output

            # 🚨 Check if input file exists
            if not os.path.exists(input_file):
                print(f"ERROR: {input_file} not found! Skipping charge addition.")
                return None

            mol = Chem.MolFromPDBFile(input_file, removeHs=False)
            if mol is None:
                print(f"ERROR: RDKit cannot read {input_file}. Skipping.")
                return None

            mol = Chem.AddHs(mol)
            AllChem.ComputeGasteigerCharges(mol)

            # Save in SDF format (intermediate step)
            writer = Chem.SDWriter(sdf_file)
            writer.write(mol)
            writer.close()

            # Convert SDF to MOL2 using Open Babel
            subprocess.run(["obabel", sdf_file, "-O", mol2_file], check=True)

            return mol2_file if os.path.exists(mol2_file) else None
        except Exception as e:
            print(f"Error adding partial charges: {str(e)}")
            return None

    def process(self):
        """Main processing function."""
        try:
            self.create_output_directory()

            non_standard_dirs = self.extract_non_standard_residues()
            if not non_standard_dirs:
                print("No non-standard amino acids found.")
                return

            print(f"Found {len(non_standard_dirs)} non-standard amino acids.")

            for residue_folder in non_standard_dirs:
                print(f"Processing {residue_folder}")

                capped_file = self.add_terminal_groups_with_pymol(residue_folder)
                hydrogenated_file = self.add_hydrogens(residue_folder)
                charged_file = self.add_partial_charges(residue_folder)

                print(f"Processed: {charged_file}")

        except Exception as e:
            print(f"Error during processing: {str(e)}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_pdb_file>")
        sys.exit(1)

    processor = NonStandardAminoAcidProcessor(sys.argv[1])
    processor.process()


if __name__ == "__main__":
    main()
