import coot_headless_api

chapi = coot_headless_api.molecules_container_t(True)
chapi.set_use_gemmi(False)

# read coordinates and map
imol = chapi.read_pdb("rnase.pdb")
imol_mtz = chapi.read_mtz("rnase-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", False, False)

# delete water "outliers" - e.g., those with a distance to the protein less than 2.5 
# or more than 3.5

min_dist = 2.5
max_dist = 3.5
median_temperature_factor = chapi.get_median_temperature_factor(imol)
b_factor_limit = 2.0 * median_temperature_factor
outlier_map_rmsd_level = 1.0
ignore_part_occ_contact_flag = False
ignore_zero_occ_flag = False
water_outliers = chapi.find_water_baddies(imol,
                                        imol_mtz,
                                        b_factor_limit,
                                        outlier_map_rmsd_level,
                                        min_dist,
                                        max_dist,
                                        ignore_part_occ_contact_flag,
                                        ignore_zero_occ_flag)
for res in water_outliers:
    cid = "//" + res.chain_id + "/" + str(res.res_no)
    print("Deleting water", cid)
    chapi.delete_atom_using_cid(imol, cid)