# Read the coordinates and map for 2pwt
imol = mc.read_pdb("2pwt.pdb")
imol_map = mc.read_ccp4_map("2pwt.ccp4", False)

# Associate mtz file with map
mc.associate_data_mtz_file_with_map(imol_map,"2pwt.mtz" , "FP", "SIGFP", "FREE")

# Refine using servalcat
mc.servalcat_refine_xray(imol, imol_map, "ref1")
