set retain_order, 1
set pdb_retain_ids, 1

# Load PDB file
load r"MLU_10.pdb", prot

# Select and modify the structure
select pk1, name N
editor.attach_amino_acid("pk1", "ace")

select pk1, name C and not resn ace and name C
editor.attach_amino_acid("pk1", "nme")

# Save the modified structure
cmd.save(r"MLU_10_capped.pdb", "prot")
