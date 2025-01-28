import os
import sys
from Bio import PDB
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue, Atom
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


class NonStandardAminoAcidProcessor:
    def __init__(self, input_file, is_peptide=False):
        """
        - input_file (str): Input PDB file.
        - is_peptide (bool): Flag to indicate if the input is a non-standard peptide.
        """
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
        """Check if the residue has the basic amino acid backbone atoms"""
        backbone_atoms = {'N', 'CA', 'C'}
        residue_atoms = {atom.get_name() for atom in residue}
        return backbone_atoms.issubset(residue_atoms)

    def verify_atom_positions(self, residue):
        """Verify if the residue contains backbone atoms N, CA, and C,
        with N as the head and C as the tail."""
        try:
            # Ensure required atoms exist
            required_atoms = {'N', 'CA', 'C'}
            residue_atoms = {atom.get_name() for atom in residue}
            if not required_atoms.issubset(residue_atoms):
                return False

            # Verify the order: N is head and C is tail
            atom_list = list(residue.get_atoms())
            atom_names = [atom.get_name() for atom in atom_list]

            return atom_names.index('N') < atom_names.index('CA') < atom_names.index('C')
        except KeyError:
            return False

    def extract_non_standard_residues(self):
        """Extract non-standard amino acids and save them as separate PDB files"""
        structure = self.parser.get_structure('protein', self.input_file)
        non_standard_files = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_name = residue.get_resname()
                    if (res_name not in self.standard_residues and
                            self.is_amino_acid(residue) and
                            self.verify_atom_positions(residue)):
                        # Create a new structure for this residue
                        new_structure = Structure.Structure('non_standard')
                        new_model = Model.Model(0)
                        new_chain = Chain.Chain('A')
                        new_chain.add(residue.copy())
                        new_model.add(new_chain)
                        new_structure.add(new_model)

                        # Save the residue
                        output_file = f"{self.output_dir}/{res_name}_{residue.get_id()[1]}.pdb"
                        self.io.set_structure(new_structure)
                        self.io.save(output_file)
                        non_standard_files.append(output_file)

        return non_standard_files

    def add_hydrogens(self, pdb_file):
        """Add hydrogens to the structure using RDKit"""
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

    def add_terminal_groups(self, pdb_file): # Add specific names to the atoms for each terminal. like ACE - Cx, Cy, Oz.
        # Generate just once as a file and add it whenever needed, do not have to call if existing already
        """
        Add acetyl group to N-terminal and N-methyl group to C-terminal using RDKit.
        """
        try:
            # Parse the PDB structure with Biopython
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('modified', pdb_file)
            model = structure[0]

            # Identify N-terminal and C-terminal residues
            chains = list(model.get_chains())
            if not chains:
                print("No chains found in the structure.")
                return pdb_file

            # Assuming modification on the first chain
            chain = chains[0]
            residues = list(chain.get_residues())
            if not residues:
                print("No residues found in the chain.")
                return pdb_file

            n_terminal = residues[0]
            c_terminal = residues[-1]

            # Check for backbone atoms
            try:
                n_atom = n_terminal['N']
                ca_atom = n_terminal['CA']
                c_atom = c_terminal['C']
            except KeyError as e:
                print(f"Missing backbone atom: {e}")
                return pdb_file

            # Generate Acetyl group (CH3-CO-) using RDKit
            acetyl_smiles = "CC(=O)"
            acetyl_mol = Chem.MolFromSmiles(acetyl_smiles)
            acetyl_mol = Chem.AddHs(acetyl_mol)
            AllChem.EmbedMolecule(acetyl_mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(acetyl_mol)
            acetyl_conf = acetyl_mol.GetConformer()

            # Generate N-Methyl group (-CH3) using RDKit
            n_methyl_smiles = "CN"
            n_methyl_mol = Chem.MolFromSmiles(n_methyl_smiles)
            n_methyl_mol = Chem.AddHs(n_methyl_mol)
            AllChem.EmbedMolecule(n_methyl_mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(n_methyl_mol)
            n_methyl_conf = n_methyl_mol.GetConformer()

            # Create new residues for acetyl and N-methyl groups
            acetyl_residue = Residue.Residue((' ', n_terminal.get_id()[1], ' '), 'ACE', n_terminal.get_id())
            n_methyl_residue = Residue.Residue((' ', c_terminal.get_id()[1], ' '), 'NME', c_terminal.get_id())

            # Define Atom names and types for acetyl group
            acetyl_atoms = ['C', 'CH3', 'O'] #check atoms in both ace and nme
            acetyl_atom_types = ['C', 'C', 'O'] #Assign atom names

            # Add atoms from acetyl group to the acetyl residue
            for i, atom_name in enumerate(acetyl_atoms):
                atom = acetyl_mol.GetAtomWithIdx(i)
                pos = acetyl_conf.GetAtomPosition(i)
                new_atom = Atom.Atom(atom_name, np.array([pos.x, pos.y, pos.z]), 0.0, 1.0, ' ', atom_name, i)
                acetyl_residue.add(new_atom)
#input and out
            # Define Atom names and types for N-methyl group
            n_methyl_atoms = ['C', 'N', 'CH3']
            n_methyl_atom_types = ['C', 'N', 'C']

            # Add atoms from N-methyl group to the N-methyl residue
            for i, atom_name in enumerate(n_methyl_atoms):
                atom = n_methyl_mol.GetAtomWithIdx(i)
                pos = n_methyl_conf.GetAtomPosition(i)
                new_atom = Atom.Atom(atom_name, np.array([pos.x, pos.y, pos.z]), 0.0, 1.0, ' ', atom_name, i)
                n_methyl_residue.add(new_atom)

            # Attach Acetyl group to N-terminus
            # Translate acetyl group to align with N-terminal N atom
            acetyl_translation = self.calculate_translation_vector(n_atom.get_coord(), acetyl_conf.GetAtomPosition(0))
            for atom in acetyl_residue:
                atom.set_coord(atom.get_coord() + acetyl_translation)

            # Attach N-Methyl group to C-terminus
            # Translate N-methyl group to align with C-terminal C atom
            n_methyl_translation = self.calculate_translation_vector(c_atom.get_coord(),
                                                                     n_methyl_conf.GetAtomPosition(1))
            for atom in n_methyl_residue:
                atom.set_coord(atom.get_coord() + n_methyl_translation)

            # Add the new residues to the chain
            chain.add(acetyl_residue)
            chain.add(n_methyl_residue)

            # Save the modified structure to a new PDB file
            output_file = pdb_file.replace('.pdb', '_modified.pdb')
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_file)

            print(f"Terminal groups added successfully. Modified structure saved as {output_file}")
            return output_file

        except Exception as e:
            print(f"Error in adding terminal groups: {str(e)}")
            return pdb_file  # Return original file if modification fails


    def process(self):
        """Main processing function"""
        try:
            self.create_output_directory()

            if self.is_peptide:
                print("Input is a non-standard peptide. Skipping extraction step...")
                hydrogenated_file = self.add_hydrogens(self.input_file)
                final_file = self.add_terminal_groups(hydrogenated_file)
                print(f"Processing complete. Final structure saved as {final_file}")
                return

            non_standard_files = self.extract_non_standard_residues()

            if not non_standard_files:
                print("No non-standard amino acids found in the input structure.")
                return

            print(f"Found {len(non_standard_files)} non-standard amino acids.")

            # Process the first non-standard residue
            first_file = non_standard_files[0]
            print(f"Processing {first_file}")

            hydrogenated_file = self.add_hydrogens(first_file)
            final_file = self.add_terminal_groups(hydrogenated_file)
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
