# read coordinates and map
imol = mc.read_pdb('tutorial-modern.pdb')
imol_mtz = mc.read_mtz("rnasa-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", False, False)

# set the parameters for waters addition (the default values are given as arguments)
mc.set_add_waters_water_to_protein_distance_lim_min(2.4)
mc.set_add_waters_water_to_protein_distance_lim_max(3.4)
mc.set_add_waters_variance_limit(0.1)
mc.set_add_waters_sigma_cutoff(1.75)

# add waters
mc.add_waters(imol, imol_mtz)
