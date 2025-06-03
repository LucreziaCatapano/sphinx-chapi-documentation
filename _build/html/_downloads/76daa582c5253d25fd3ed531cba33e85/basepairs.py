import pandas as pd
import numpy as np
import chapi


desired_atoms = {"O6", "N4", "N1", "N3", "N2", "O2"}

def check_d(pdb_path: str, residue_numbers: list[int], distance_threshold: float = 3.1):
    
    mc = chapi.molecules_container_t(False)
    mc.set_use_gemmi(False)
    
    # Read molecular coordinates and maps
    imol = mc.read_coordinates(pdb_path)
    distances_list = []

    for current_residue_res_no in residue_numbers:
        current_residue_cid = f"//B/{current_residue_res_no}"

        neighbours = mc.get_residues_near_residue(imol, 
                                                  current_residue_cid, 
                                                  distance_threshold)
        for n in neighbours:
            name = mc.get_residue_name(imol, n.chain_id, n.res_no, n.ins_code)
            if name == "HOH":
                continue  # Ignore water molecules
            delta_res_no = n.res_no - current_residue_res_no
            if np.abs(delta_res_no) == 1:
                continue  # Ignore adjacent residues
            neighbour_cid = f"//{n.chain_id}/{n.res_no}"
            distances = mc.get_distances_between_atoms_of_residues(imol, 
                                                                   current_residue_cid, 
                                                                   neighbour_cid, 
                                                                   distance_threshold)
            
            for d in distances:
                atom1_name = d.atom_1.atom_name.strip()  
                atom2_name = d.atom_2.atom_name.strip() 

                if atom1_name in desired_atoms and atom2_name in desired_atoms:
                    distance_data = {
                        "Chain ID1": d.atom_1.chain_id,
                        "Residue_no1": d.atom_1.res_no,
                        "Atom1": atom1_name,
                        "Chain ID2": d.atom_2.chain_id,
                        "Residue_no2": d.atom_2.res_no,
                        "Atom2": atom2_name,
                        "Distance": d.distance
                    }
                    distances_list.append(distance_data)

    # Return the distances as a pandas DataFrame
    return pd.DataFrame(distances_list)

residue_numbers = [i for i in range(24,45)]
df1 = check_d("2pwt/2pwt.pdb", residue_numbers)
df2 = check_d("2pwt_coot_refined.pdb", residue_numbers)

df1 = df1.rename(columns={"Distance": "Distance_original"})
df2 = df2.rename(columns={"Distance": "Distance_coot"})

# Merge the two DataFrames on the specified columns
df_merged = df1.merge(df2, on=["Chain ID1", "Chain ID2", "Residue_no1", "Residue_no2", "Atom1", "Atom2"])
df_merged.to_csv("distance_original_coot.csv", index=False)

# Calculate the differences and save to a new CSV
df_merged["Difference_coot"] = df_merged["Distance_coot"] - df_merged["Distance_original"]
df_merged.to_csv("basepairs_distance_differences.csv", index=False)