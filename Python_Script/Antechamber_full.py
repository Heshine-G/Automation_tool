import subprocess
import os
from rdkit import Chem


def run_antechamber(input_mol2):
    """Runs antechamber, generates .mc file, and marks capping atoms as OMIT_NAME."""

    base_name = os.path.splitext(input_mol2)[0]
    residue_name = base_name.upper()

    # Define output file names
    ac_output = f"{base_name}.ac"  # Change from .ac to .mc
    mol2_output = f"{base_name}.mol2"
    lib_output = f"{base_name}.lib"
    prepin_output = f"{base_name}.prepin"
    mc_output = f"{base_name}.mc"
    frcmod_output = f"{base_name}.frcmod"
    gaff_frcmod_output = f"{residue_name}_gaff.frcmod"
    ff14SB_frcmod_output = f"{residue_name}_ff14SB.frcmod"

    # First antechamber command (.ac file)
    cmd1 = f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo ac -o {ac_output} -c bcc -at amber"

    # Second antechamber command (.mol2 file)
    cmd2 = f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo mol2 -o {mol2_output} -c bcc -at amber"

    try:
        # Run first command
        subprocess.run(cmd1, shell=True, check=True)
        print(f"Generated {mc_output} successfully.")

        # Run second command
        subprocess.run(cmd2, shell=True, check=True)
        print(f"Generated {mol2_output} successfully.")

    except subprocess.CalledProcessError as e:
        print(f"Error running antechamber: {e}")
        return

    # Create Leap input script
    leap_script = f"""
    source leaprc.gaff
    {residue_name} = loadmol2 {mol2_output}
    edit {residue_name}
    desc {residue_name}
    set {residue_name} head {residue_name}.1.N
    set {residue_name} tail {residue_name}.1.C
    saveoff {residue_name} {lib_output}
    quit
    """

    leap_input_file = "leap.in"
    with open(leap_input_file, "w") as f:
        f.write(leap_script)

    # Run tleap
    try:
        subprocess.run("tleap -f leap.in", shell=True, check=True)
        print(f"Generated {lib_output} successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running Leap: {e}")
        return

    # Run prepgen command
    prepgen_cmd = f"prepgen -i {ac_output} -o {prepin_output} -m {mc_output} -rn {residue_name}"
    try:
        subprocess.run(prepgen_cmd, shell=True, check=True)
        print(f"Generated {prepin_output} successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running prepgen: {e}")
        return

    # Run parmchk2 commands
    parmchk_cmd1 = f"parmchk2 -i {prepin_output} -f prepi -o {frcmod_output} -a Y -p $AMBERHOME/dat/leap/parm/parm10.dat"
    parmchk_cmd2 = f"parmchk2 -i {ac_output} -f ac -o {gaff_frcmod_output}"
    parmchk_cmd3 = f"parmchk2 -i {ac_output} -f ac -o {ff14SB_frcmod_output} -a Y -p $AMBERHOME/dat/leap/parm/parm10.dat"

    try:
        subprocess.run(parmchk_cmd1, shell=True, check=True)
        print(f"Generated {frcmod_output} successfully.")
        subprocess.run(parmchk_cmd2, shell=True, check=True)
        print(f"Generated {gaff_frcmod_output} successfully.")
        subprocess.run(parmchk_cmd3, shell=True, check=True)
        print(f"Generated {ff14SB_frcmod_output} successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running parmchk2: {e}")
        return


run_antechamber("MVA.mol2")
