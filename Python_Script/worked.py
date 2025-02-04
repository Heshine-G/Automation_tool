import os
import sys
import subprocess
from Bio import PDB
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue
from rdkit import Chem
from rdkit.Chem import AllChem

class NonStandardAminoAcidProcessor:
    def __init__(self, input_file, is_peptide=False):
        self.input_file = input_file
        self.is_peptide = is_peptide
        self.parser = PDB.PDBParser(QUIET=True)
        self.io = PDBIO()
        self.standard_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL'
        }
        self.output_dir = "non_standard_residues"

    def create_output_directory(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def is_amino_acid(self, residue):
        """Check if the residue has the basic amino acid backbone atoms (N, CA, C)"""
        backbone_atoms = ['N', 'CA', 'C']
        residue_atoms = {atom.get_name() for atom in residue}

        return set(backbone_atoms).issubset(residue_atoms)

    def extract_non_standard_residues(self):
        """Extract non-standard residues from the input structure."""
        structure = self.parser.get_structure('protein', self.input_file)
        non_standard_files = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() not in self.standard_residues and self.is_amino_acid(residue):
                        new_structure = Structure.Structure('non_standard')
                        new_model = Model.Model(0)
                        new_chain = Chain.Chain('A')
                        new_chain.add(residue.copy())
                        new_model.add(new_chain)
                        new_structure.add(new_model)

                        output_file = f"{self.output_dir}/{residue.get_resname()}_{residue.get_id()[1]}.pdb"
                        self.io.set_structure(new_structure)
                        self.io.save(output_file)
                        non_standard_files.append(output_file)

        return non_standard_files

    def add_hydrogens(self, pdb_file):
        """Protonate the PDB structure using RDKit."""
        try:
            mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
            if mol is None:
                print(f"Error: Unable to parse {pdb_file} with RDKit.")
                return pdb_file

            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)

            output_file = pdb_file.replace('.pdb', '_protonated.pdb')
            Chem.MolToPDBFile(mol, output_file)
            return output_file
        except Exception as e:
            print(f"Error adding hydrogens: {str(e)}")
            return pdb_file

    def add_terminal_groups_with_pymol(self, input_file):
        """Call PyMOL script"""
        try:
            output_file = f"{os.path.splitext(input_file)[0]}_capped.pdb"
            pymol_command = f"pymol -c add_capping.pml {input_file} {output_file}"
            subprocess.run(pymol_command, shell=True, check=True)
            return output_file
        except Exception as e:
            print(f"Error in running PyMOL: {str(e)}")
            return input_file

        def add_partial_charges(self, pdb_file):
        """Add Gasteiger partial charges using RDKit and save as .mol2 using Open Babel."""
        try:
            mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
            if mol is None:
                print(f"Error: Unable to parse {pdb_file} for charge assignment.")
                return pdb_file

            # Ensure hydrogens are present
            mol = Chem.AddHs(mol)

            # Compute Gasteiger charges
            AllChem.ComputeGasteigerCharges(mol)

            # Save to temporary SDF format as RDKit does not support .mol2
            sdf_file = pdb_file.replace(".pdb", "_charged.sdf")
            Chem.SDWriter(sdf_file).write(mol)

            # Convert SDF to MOL2 using Open Babel
            mol2_file = pdb_file.replace(".pdb", "_charged.mol2")
            os.system(f"obabel {sdf_file} -O {mol2_file}")

            return mol2_file

        except Exception as e:
            print(f"Error adding partial charges: {str(e)}")
            return pdb_file

    def process(self):
        try:
            self.create_output_directory()

            if self.is_peptide:
                print("Input is a non-standard peptide. Skipping extraction step...")
                hydrogenated_file = self.add_hydrogens(self.input_file)
                capped_file = self.add_terminal_groups_with_pymol(hydrogenated_file)
                charged_file = self.add_partial_charges(capped_file)
                print(f"Processing complete. Check output PDB file: {charged_file}")
                return

            non_standard_files = self.extract_non_standard_residues()

            if not non_standard_files:
                print("No non-standard amino acids found in the input structure.")
                return

            print(f"Found {len(non_standard_files)} non-standard amino acids.")
            first_file = non_standard_files[0]
            print(f"Processing {first_file}")

            hydrogenated_file = self.add_hydrogens(first_file)
            capped_file = self.add_terminal_groups_with_pymol(hydrogenated_file)
            charged_file = self.add_partial_charges(capped_file)

            print(f"Processing complete. Check output PDB file: {charged_file}")

        except Exception as e:
            print(f"Error during processing: {str(e)}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_pdb_file> [--peptide]")
        sys.exit(1)

    input_file = sys.argv[1]
    is_peptide = '--peptide' in sys.argv

    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found.")
        sys.exit(1)

    processor = NonStandardAminoAcidProcessor(input_file, is_peptide)
    processor.process()


if __name__ == "__main__":
    main()
