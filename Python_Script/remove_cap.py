import re

# Input and output file names
input_pdb = "residue_capped.pdb"
output_mc = "residue_capped.mc"

# Open the input file and process it
with open(input_pdb, "r") as pdb_file, open(output_mc, "w") as mc_file:
    for line in pdb_file:
        columns = line.split()
        if len(columns) > 12 and columns[3] in ["ACE", "NME"]:
            mc_file.write(f"OMIT_NAME {columns[11]} {columns[1]}\n")

print(f"Extraction complete. File saved as {output_mc}")
