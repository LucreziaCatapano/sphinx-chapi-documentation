# Read the molecular coordinates and map for 2pwt
imol = mc.read_coordinates("2pwt.pdb")
imol_map = mc.read_ccp4_map("2pwt.ccp4", False)

# Set the refinement map
mc.set_imol_refinement_map(imol_map)

# Add hydrogen atoms to the molecule
add_H_status = mc.add_hydrogen_atoms(imol)

# If hydrogen atoms were successfully added, refine the residues
if add_H_status == 1:
    
    # Import the CIF dictionary for refinement
    mc.import_cif_dictionary("/Applications/ccp4-9/lib/data/monomers/l/LHA.cif", imol) 

    # Refine the residues in the specified range
    mc.refine_residues_using_atom_cid(imol, "//B/24-44", "ALL", 10000)

    # Delete the hydrogen atoms after refinement
    mc.delete_hydrogen_atoms(imol)

    # Write the refined coordinates to a new PDB file
    mc.write_coordinates(imol, "2pwt_coot_refined.pdb")   