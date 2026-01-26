import coot_headless_api

chapi = coot_headless_api.molecules_container_t(True)
chapi.set_use_gemmi(False)

# read coordinates and map
imol = chapi.read_pdb("rnase.pdb")
imol_mtz = chapi.read_mtz("rnase-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", False, False)

# set the parameters for waters addition (the default values are given as arguments)
chapi.set_add_waters_water_to_protein_distance_lim_min(2.4)
chapi.set_add_waters_water_to_protein_distance_lim_max(3.4)
chapi.set_add_waters_variance_limit(0.1)
chapi.set_add_waters_sigma_cutoff(1.75)

# add waters
chapi.add_waters(imol, imol_mtz)