import subprocess
import os

def run_antechamber(input_mol2):

   # Runs antechamber to generate .ac and .mol2 files from an input .mol2 file.

    base_name = os.path.splitext(input_mol2)[0]
    residue_name = base_name.upper()

    # Define output file names
    ac_output = f"{base_name}.ac"
    mol2_output = f"{base_name}.mol2"

    # First antechamber command (.ac file)
    cmd1 = f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo ac -o {ac_output} -c bcc -at amber"

    # Second antechamber command (.mol2 file)
    cmd2 = f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo mol2 -o {mol2_output} -c bcc -at amber"

    try:
        # Run first command
        subprocess.run(cmd1, shell=True, check=True)
        print(f"Generated {ac_output} successfully.")

        # Run second command
        subprocess.run(cmd2, shell=True, check=True)
        print(f"Generated {mol2_output} successfully.")

    except subprocess.CalledProcessError as e:
        print(f"Error running antechamber: {e}")


run_antechamber("MVA.mol2")
