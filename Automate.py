import os
import sys
from Bio import PDB
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd, editor

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
        """Check if the residue has the basic amino acid backbone atoms (N, CA, C)
        and ensure their correct order: N (head), CA (middle), C (tail)."""
        backbone_atoms = ['N', 'CA', 'C']
        residue_atoms = {atom.get_name() for atom in residue}

        # Ensure required backbone atoms exist
        if not set(backbone_atoms).issubset(residue_atoms):
            return False

        # Verify order of backbone atoms
        atom_names = [atom.get_name() for atom in residue.get_atoms()]
        return all(atom in atom_names for atom in backbone_atoms) and \
            atom_names.index('N') < atom_names.index('CA') < atom_names.index('C')

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
        """Protonate the PDB structure."""
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

    def create_pml_script(self, input_file, output_file):
        """Creates a PyMOL script to cap the terminals and save the output file."""
        pml_script = f"""
# add_capping.pml

set retain_order, 1
set pdb_retain_ids, 1

# Load the input PDB file
load {input_file}

# Select the terminal N atom for ACE capping
select terminal_N, name N 
editor.attach_amino_acid("terminal_N", "ace")

# Select the terminal C atom (not part of ACE) for NME capping
select terminal_C, name C and not resn ACE and name C
editor.attach_amino_acid("terminal_C", "nme")

# Save the final capped structure
save {output_file}

quit
"""
        pml_filename = "add_capping.pml"
        with open(pml_filename, "w") as file:
            file.write(pml_script)

        return pml_filename

    def add_terminal_groups_with_pymol(self, pdb_file):
        """Run PyMOL command to cap terminal residues."""
        try:
            output_pdb = pdb_file.replace('.pdb', '_capped.pdb')

            # Generate the .pml script
            pml_script = self.create_pml_script(pdb_file, output_pdb)

            # Run PyMOL from command line to execute the script
            pymol_command = f"pymol -c {pml_script}"
            os.system(pymol_command)

            # Check if the output file was saved
            if not os.path.exists(output_pdb):
                print(f"Error: Terminal capping failed. Check your PyMOL command.")
                return pdb_file

            return output_pdb
        except Exception as e:
            print(f"Error in adding terminal groups: {str(e)}")
            return pdb_file

    def reorder_pdb(self, pdb_file):
        """Reorder the PDB to place ACE and NME in the correct positions."""
        with open(pdb_file, 'r') as file:
            lines = file.readlines()

        ace_lines = [line for line in lines if "ACE" in line]
        nme_lines = [line for line in lines if "NME" in line]
        other_lines = [line for line in lines if "ACE" not in line and "NME" not in line]

        reordered_lines = other_lines + ace_lines + nme_lines

        output_file = pdb_file.replace('.pdb', '_fixed.pdb')
        with open(output_file, 'w') as file:
            file.writelines(reordered_lines)

        return output_file

    def process(self):
        try:
            self.create_output_directory()

            if self.is_peptide:
                print("Input is a non-standard peptide. Skipping extraction step...")
                hydrogenated_file = self.add_hydrogens(self.input_file)
                capped_file = self.add_terminal_groups_with_pymol(hydrogenated_file)
                print(f"Processing complete. Final structure saved as {capped_file}")
                return

            # Handle non-peptide case (for non-standard residues extraction, protonation, and capping)
            non_standard_files = self.extract_non_standard_residues()

            if not non_standard_files:
                print("No non-standard amino acids found in the input structure.")
                return

            print(f"Found {len(non_standard_files)} non-standard amino acids.")
            first_file = non_standard_files[0]
            print(f"Processing {first_file}")

            hydrogenated_file = self.add_hydrogens(first_file)
            capped_file = self.add_terminal_groups_with_pymol(hydrogenated_file)
            final_file = self.reorder_pdb(capped_file)
            print(f"Processing complete. Final structure saved as {final_file}")

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
