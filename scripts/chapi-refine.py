# Read the molecular coordinates and map for 2pwt
imol = chapi.read_coordinates("2pwt.pdb")
imol_map = chapi.read_ccp4_map("2pwt.ccp4", False)

# Set the refinement map
chapi.set_imol_refinement_map(imol_map)

# Add hydrogen atoms to the molecule
add_H_status = chapi.add_hydrogen_atoms(imol)

# If hydrogen atoms were successfully added, refine the residues
if add_H_status == 1:
    
    # Import the CIF dictionary for refinement
    chapi.import_cif_dictionary("/Applications/ccp4-9/lib/data/monomers/l/LHA.cif", imol) 

    # Refine the residues in the specified range
    chapi.refine_residues_using_atom_cid(imol, "//B/24-44", "ALL", 10000)

    # Delete the hydrogen atoms after refinement
    chapi.delete_hydrogen_atoms(imol)

    # Write the refined coordinates to a new PDB file
    chapi.write_coordinates(imol, "2pwt_coot_refined.pdb")   