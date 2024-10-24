class molecules_container_t():


    def simple_mesh_to_pythonic_mesh(mesh: int , mesh_mode: int):
        pass

    def get_pythonic_bonds_mesh(imol: int, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int):
        pass

    def get_pythonic_map_mesh(imol: int, x: float, y: float, z: float, radius: float, contour_level: float):
        pass

    def get_pythonic_molecular_representation_mesh(imol: int, atom_selection: str, colour_sheme: str, style: str):
        pass

    def get_pythonic_gaussian_surface_mesh(imol: int, sigma: float, contour_level: float, box_radius: float, grid_scale: float, fft_b_factor: float):
        pass

    def get_pythonic_simple_molecule(imol: int, cid: str, include_hydrogen_atoms_flag: bool):
        """ make a "proper" simple molecule python class one day. """
        """ 4: aromaticity flag (boolean) """
        pass


    def set_make_backups(state: bool):
        """ Allow the user to disable/enable backups ("""
        pass

    def get_make_backups():
        pass

    def file_name_to_string(file_name: str):
        pass

    def get_number_of_molecules():
        pass

    def create_empty_molecules(n_empty: int):
        """ The adds a number of empty molecules to the internal vector of molecules Note that his is not like """
        pass

    def set_imol_refinement_map(i: int):
        """ set the map used for refinement and fitting """
        pass

    def set_map_weight(w: float):
        """ set the map weight """
        pass

    def get_map_weight():
        pass

    def atom_cid_to_atom_spec(imol: int, cid: str):
        """ Convert atom cid string to a coot atom specifier. The test for these failing is """
        pass

    def residue_cid_to_residue_spec(imol: int, cid: str):
        """ Convert residue cid string to a coot residue specifier. """
        pass

    def set_show_timings(s: bool):
        """ this set the show_timings flag. Various (not all) functions in this class can calculate how long they took to run. Setting this will write the time to taken (in milliseconds) to stdout. The default is """
        pass

    def get_geom():
        pass

    def get_header_info(imol: int):
        """ get header info. """
        pass


    def read_coordinates(file_name: str):
        """ read a coordinates file (mmcif or PDB) """
        pass

    def read_pdb(file_name: str):
        """ read a PDB file (or mmcif coordinates file, despite the name). It does the same job as """
        pass

    def print_secondary_structure_info(imol: int):
        """ print the secondary structure information to standard out """
        pass

    def replace_molecule_by_model_from_file(imol: int, pdb_file_name: str):
        """ read a PDB file (or mmcif coordinates file, despite the name) to replace the current molecule. This will only work if the molecules is already a model molecule """
        pass

    def split_multi_model_molecule(imol: int):
        """ split an NMR model into multiple models - all in MODEL 1. """
        pass

    def make_ensemble(model_molecule_list: str):
        """ make a multi-model molecule given the input molecules """
        pass

    def molecule_to_PDB_string(imol: int):
        pass

    def molecule_to_mmCIF_string(imol: int):
        pass

    def get_active_atom(x: float, y: float, z: float, displayed_model_molecules_list: str):
        """ get the active atom given the screen centre"""
        pass

    def import_cif_dictionary(cif_file_name: str, imol_enc: int):
        """ IMOL_ENC_ANY = -999999 """
        pass

    def get_cif_file_name(comp_id: str, imol_enc: int):
        pass

    def get_cif_restraints_as_string(comp_id: str, imol_enc: int):
        pass

    def copy_dictionary(monomer_name: str, imol_current: int, imol_new: int):
        """ copy the dictionary that is specific for imol_current so that it can be used with imol_new """
        pass

    def get_monomer(monomer_name: str):
        """ get a monomer """
        pass

    def get_monomer_from_dictionary(comp_id: str, imol: int, idealised_flag: bool):
        """ get a monomer for a particular molecule - use -999999 (IMOL_ENC_ANY) if no molecule-specific dictionary is needed. """
        pass

    def get_monomer_and_position_at(comp_id: str, imol: int, x: float, y: float, z: float):
        """ get monomer and place it at the given position for a particular molecule - use -999999 if no molecule-specific dictionary is needed """
        pass

    def get_groups_for_monomers(residue_names: list):
        pass

    def get_group_for_monomer(residue_name: str):
        pass

    def get_hb_type(compound_id: str, imol_enc: int, atom_name: str):
        pass

    def get_gphl_chem_comp_info(compound_id: str, imol_enc: int):
        pass

    def write_png(compound_id: str, imol: int, file_name: str):
        """ write a PNG for the given compound_id. imol can be IMOL_ENC_ANY Currently this function does nothing (drawing is done with the not-allowed cairo) """
        pass

    def write_coordinates(imol: int, file_name: str):
        """ write the coordinate to the give file name """
        pass

    def set_draw_missing_residue_loops(state: bool):
        """ By default missing loops are drawn. This function allows missing loops to not be drawn. Sometimes that can clarify the representation. This is a lightweight function that sets a flag that is used by subsequent calls to """
        pass

    def get_bonds_mesh(imol: int, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int):
        """ get the bonds mesh."""
        pass

    def get_bonds_mesh_instanced(imol: int, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int):
        """ The arguments are as above:"""
        pass

    def get_bonds_mesh_for_selection_instanced(imol: int, atom_selection_cid: str, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int):
        """ As above, but only return the bonds for the atom selection. Typically one would call this with a wider bond_with than one would use for standards atoms (all molecule)"""
        pass

    def get_goodsell_style_mesh_instanced(imol: int, colour_wheel_rotation_step: float, saturation: float, goodselliness: float):
        pass

    def export_map_molecule_as_gltf(imol: int, pos_x: float, pos_y: float, pos_z: float, radius: float, contour_level: float, file_name: str):
        """ export map molecule as glTF """
        pass

    def export_model_molecule_as_gltf(imol: int, selection_cid: str, mode: str, against_a_dark_background: bool, bonds_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int, draw_hydrogen_atoms_flag: bool, draw_missing_residue_loops: bool, file_name: str):
        """ export model molecule as glTF - This API will change - we want to specify surfaces and ribbons too. """
        pass

    def export_molecular_represenation_as_gltf(imol: int, atom_selection_cid: str, colour_scheme: str, style: str, file_name: str):
        pass

    def get_colour_table(imol: int, against_a_dark_background: bool):
        """ return the colur table (for testing) """
        pass

    def set_colour_wheel_rotation_base(imol: int, r: float):
        """ set the colour wheel rotation base for the specified molecule (in degrees) """
        pass

    def set_base_colour_for_bonds(imol: int, r: float, g: float, b: float):
        """ set the base colour - to be used as a base for colour wheel rotation """
        pass

    def add_to_non_drawn_bonds(imol: int, atom_selection_cid: str):
        """ add a atom selection cid for atoms and bonds not to be drawn """
        pass

    def clear_non_drawn_bonds(imol: int):
        """ clear the set of non-drawn atoms (so that they can be displayed again) """
        pass

    def print_non_drawn_bonds(imol: int):
        pass

    def set_user_defined_bond_colours(imol: int, colour_map: list):
        """ user-defined colour-index to colour """
        pass

    def set_user_defined_atom_colour_by_selection(imol: int, indexed_residues_cids: list, colour_applies_to_non_carbon_atoms_also: bool):
        """ set the user-defined residue selections (CIDs) to colour index """
        pass

    def add_colour_rule(imol: int, selection_cid: str, colour: str):
        """ Add a colour rule for M2T representations. """
        pass

    def add_colour_rules_multi(imol: int, selections_and_colours_combo_string: str):
        """ add multiple colour rules, combined like the following "//A/1^#cc0000|//A/2^#cb0002|//A/3^#c00007" i.e. "|" is the separator for each rule and "^" is the separator for the selection string and the colour string """
        pass

    def delete_colour_rules(imol: int):
        """ delete the colour rules for the given molecule """
        pass

    def get_colour_rules(imol: int):
        """ get the colour rules """
        pass

    def print_colour_rules(imol: int):
        """ print the colour rules """
        pass

    def set_use_bespoke_carbon_atom_colour(imol: int, state: bool):
        """ use bespoke carbon atom colour """
        pass

    def set_bespoke_carbon_atom_colour(imol: int, col: int ):
        """ set bespoke carbon atom colour """
        pass

    def M2T_updateFloatParameter(imol: int, param_name: str, value: float):
        """ Update float parameter for MoleculesToTriangles molecular mesh. """
        pass

    def M2T_updateIntParameter(imol: int, param_name: str, value: int):
        """ Update int parameter for MoleculesToTriangles molecular mesh. """
        pass

    def get_molecular_representation_mesh(imol: int, cid: str, colour_scheme: str, style: str):
        """ get ribbon and surface representation """
        pass

    def get_gaussian_surface(imol: int, sigma: float, contour_level: float, box_radius: float, grid_scale: float, b_factor: float):
        """ b_factor = 100.0 (use 0.0 for no FFT-B-factor smoothing)"""
        pass

    def get_chemical_features_mesh(imol: int, cid: str):
        """ get chemical feaatures for the specified residue """
        pass

    def get_number_of_atoms(imol: int):
        pass

    def get_molecule_diameter(imol: int):
        pass

    def get_number_of_hydrogen_atoms(imol: int):
        pass

    def get_chains_in_model(imol: int):
        pass

    def get_ncs_related_chains(imol: int):
        """ Get the chains that are related by NCS or molecular symmetry: """
        pass

    def get_single_letter_codes_for_chain(imol: int, chain_id: str):
        pass

    def get_residue_names_with_no_dictionary(imol: int):
        pass

    def get_residue_name(imol: int, chain_id: str, res_no: int, ins_code: str):
        """ get residue name """
        pass

    def residues_with_missing_atoms(imol: int):
        pass

    def missing_atoms_info_raw(imol: int):
        """ Ths function is not const because missing_atoms() takes a non-const pointer to the geometry. """
        pass

    def get_residues_near_residue(imol: int, residue_cid: str, dist: float):
        pass

    def SSM_superpose(imol_ref: int, chain_id_ref: str, imol_mov: int, chain_id_mov: str):
        """ The specified chain of the moving molecule is superposed onto the chain in the reference molecule (if possible). There is some alignment screen output that would be better added to the return value. """
        pass

    def add_lsq_superpose_match(chain_id_ref: str, res_no_ref_start: int, res_no_ref_end: int, chain_id_mov: str, res_no_mov_start: int, res_no_mov_end: int, match_type: int):
        """ superpose using LSQ - setup the matches @params """
        pass

    def clear_lsq_matches():
        """ clear any existing lsq matchers """
        pass

    def lsq_superpose(imol_ref: int, imol_mov: int):
        """ apply the superposition using LSQ """
        pass

    def get_lsq_matrix(imol_ref: int, imol_mov: int):
        """ return the transformation matrix in a simple class - dont apply it to the coordinates """
        pass

    def get_lsq_matrix_internal(imol_ref: int, imol_mov: int):
        """ make this private """
        pass

    def get_symmetry(imol: int, symmetry_search_radius: float, centre_x: float, centre_y: float, centre_z: float):
        """ symmetry now comes in a simple container that also includes the cell """
        pass

    def get_cell(imol: int):
        """ Check that """
        pass

    def get_map_molecule_centre(imol: int):
        """ Get the middle of the "molecule blob" in cryo-EM reconstruction maps """
        pass

    def undo(imol: int):
        """ undo """
        pass

    def redo(imol: int):
        """ redo """
        pass


    def set_map_sampling_rate(msr: float):
        """ set the map sampling rate (default is 1.8). Higher numbers mean smoother maps, but they take longer to generate, longer to transfer, longer to parse and longer to draw """
        pass

    def read_mtz(file_name: str, f: str, phi: str, weight: str, use_weight: bool, is_a_difference_map: bool):
        """ Read the given mtz file. """
        pass

    def replace_map_by_mtz_from_file(imol: int, file_name: str, f: str, phi: str, weight: str, use_weight: bool):
        """ replace map """
        pass

    def auto_read_mtz(file_name: str):
        """ Read the given mtz file. """
        pass

    def read_ccp4_map(file_name: str, is_a_difference_map: bool):
        pass

    def write_map(imol: int, file_name: str):
        """ write a map. This function was be renamed from """
        pass

    def get_map_mean(imol: int):
        pass

    def get_map_rmsd_approx(imol_map: int):
        pass

    def get_map_histogram(imol: int, n_bins: int, zoom_factor: float):
        pass

    def get_suggested_initial_contour_level(imol: int):
        pass

    def is_EM_map(imol: int):
        pass

    def sharpen_blur_map(imol_map: int, b_factor: float, in_place_flag: bool):
        """ create a new map that is blurred/sharpened """
        pass

    def sharpen_blur_map_with_resample(imol_map: int, b_factor: float, resample_factor: float, in_place_flag: bool):
        """ create a new map that is blurred/sharpened and resampling. Note that resampling can be slow, a resample_factor of 1.5 is about the limit of the trade of of prettiness for speed. """
        pass

    def mask_map_by_atom_selection(imol_coords: int, imol_map: int, cid: str, atom_radius: float, invert_flag: bool):
        """ the """
        pass

    def flip_hand(imol_map: int):
        """ generate a new map which is the hand-flipped version of the input map. """
        pass

    def make_masked_maps_split_by_chain(imol: int, imol_map: int):
        """ Make a vector of maps that are split by chain-id of the input imol """
        pass

    def set_map_colour(imol: int, r: float, g: float, b: float):
        """ set the map colour. The next time a map mesh is requested, it will have this colour. This does not affect the colour of the difference maps. """
        pass

    def set_map_is_contoured_with_thread_pool(state: bool):
        """ set the state of the mode of the threading in map contouring """
        pass

    def get_map_contours_mesh(imol: int, position_x: float, position_y: float, position_z: float, radius: float, contour_level: float):
        """ This function is not """
        pass

    def get_map_contours_mesh_using_other_map_for_colours(imol_ref: int, imol_map_for_colouring: int, position_x: float, position_y: float, position_z: float, radius: float, contour_level: float, other_map_for_colouring_min_value: float, other_map_for_colouring_max_value: float, invert_colour_ramp: bool):
        """ get the mesh for the map contours using another map for colouring """
        pass

    def set_map_colour_saturation(imol: int, s: float):
        """ set the map saturation """
        pass

    def get_latest_sfcalc_stats():
        pass

    def get_r_factor_stats():
        pass

    def r_factor_stats_as_string(rfs: int ):
        pass

    def average_map(imol_maps: str, list):
        pass

    def regen_map(imol_map: int, imol_maps: str, scales: list):
        pass



    def testing_start_long_term_job(n_seconds: int):
        """ if """
        pass

    def testing_stop_long_term_job():
        """ stop the long-term job runnning (testing function) """
        pass

    def testing_interrogate_long_term_job():
        """ get the stats for the long-term job (testing function) """
        pass

    def get_contouring_time():
        """ get the time for conntouring in milliseconds """
        pass

    def set_max_number_of_threads(n_threads: int):
        """ set the maximum number of threads for both the thread pool and the vector of threads """
        pass

    def set_max_number_of_threads_in_thread_pool(n_threads: int):
        """ deprecated name for the above function """
        pass

    def test_the_threading(n_threads: int):
        """ get the time to run a test function in milliseconds """
        pass

    def test_launching_threads(n_threads_per_batch: int, n_batches: int):
        pass

    def test_thread_pool_threads(n_threads: int):
        pass

    def get_geometry():
        pass

    def make_mesh_for_map_contours_for_blender(imol: int, x: float, y: float, z: float, level: float, radius: float):
        pass

    def make_mesh_for_bonds_for_blender(imol: int, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int):
        pass

    def make_mesh_for_molecular_representation_for_blender(imol: int, cid: str, colour_scheme: str, style: str):
        pass

    def make_mesh_for_gaussian_surface_for_blender(imol: int, sigma: float, contour_level: float, box_radius: float, grid_scale: float, b_factor: float):
        pass

    def make_mesh_for_goodsell_style_for_blender(imol: int, colour_wheel_rotation_step: float, saturation: float, goodselliness: float):
        pass

    def get_colour_table_for_blender(imol: int):
        pass

    def get_vertices_for_blender(imol: int):
        pass

    def get_triangles_for_blender(imol: int):
        pass

    def get_molecule_name(imol: int):
        pass

    def set_molecule_name(imol: int, new_name: str):
        """ set the molecule name """
        pass

    def display_molecule_names_table():
        """ debugging function: display the table of molecule and names """
        pass

    def is_valid_model_molecule(imol: int):
        pass

    def is_valid_map_molecule(imol_map: int):
        pass

    def is_a_difference_map(imol_map: int):
        pass

    def new_molecule(name: str):
        """ create an empty molecule """
        pass

    def close_molecule(imol: int):
        """ close the molecule (and delete dynamically allocated memory) """
        pass

    def end_delete_closed_molecules():
        pass

    def pop_back():
        """ delete the most recent/last molecule in the molecule vector """
        pass

    def clear():
        """ delete all molecules """
        pass

    def get_eigenvalues(imol: int, chain_id: str, res_no: int, ins_code: str):
        pass

    def test_origin_cube():
        pass

    def fill_rotamer_probability_tables():
        """ fill the rotamer probability tables (currently not ARG and LYS) """
        pass

    def accept_rotamer_probability_tables_compressed_data(data_stream: str):
        """ the caller has access to a compressed file that contains the rotamer probabilities. libcootapi will fill the rotamer probabilities tables from this compressed data stream. (placeholder only) """
        pass

    def contains_unsaved_models():
        pass

    def save_unsaved_model_changes():
        """ Save the unsaved model - this function has not yet been written! """
        pass

    def geometry_init_standard():
        """ read the stardard list of residues """
        pass

    def non_standard_residue_types_in_model(imol: int):
        pass

    def auto_fit_rotamer(imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str, imol_map: int):
        """ auto-fit rotamer """
        pass

    def change_to_next_rotamer(imol: int, residue_cid: str, alt_conf: str):
        """ change to the next rotamer (rotamer cycling is implicit if needed)"""
        pass

    def change_to_previous_rotamer(imol: int, residue_cid: str, alt_conf: str):
        """ change to the next rotamer (rotamer cycling is implicit if needed)"""
        pass

    def change_to_first_rotamer(imol: int, residue_cid: str, alt_conf: str):
        """ change to the first (0th) rotamer """
        pass

    def delete_using_cid(imol: int, cid: str, scope: str):
        """ where scope is one of the strings: ["ATOM","WATER","RESIDUE","CHAIN","MOLECULE", "LITERAL"] """
        pass

    def delete_atom(imol: int, chain_id: str, res_no: int, ins_code: str, atom_name: str, alt_conf: str):
        """ delete atom """
        pass

    def delete_atom_using_cid(imol: int, cid: str):
        """ delete atom using atom cid """
        pass

    def delete_residue(imol: int, chain_id: str, res_no: int, ins_code: str):
        """ delete residue """
        pass

    def delete_residue_using_cid(imol: int, cid: str):
        """ delete residue using cid """
        pass

    def delete_residue_atoms_with_alt_conf(imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str):
        """ delete residue atoms using alt_conf """
        pass

    def delete_residue_atoms_using_cid(imol: int, cid: str):
        """ delete residue atoms using cid """
        pass

    def delete_side_chain(imol: int, chain_id: str, res_no: int, ins_code: str):
        """ delete side chain """
        pass

    def delete_side_chain_using_cid(imol: int, cid: str):
        """ delete side chain """
        pass

    def delete_chain_using_cid(imol: int, cid: str):
        """ delete chain. """
        pass

    def delete_literal_using_cid(imol: int, cid: str):
        """ delete the atoms specified in the CID selection """
        pass

    def add_terminal_residue_directly(imol: int, chain_id: str, res_no: int, ins_code: str):
        """ add a residue onto the end of the chain by fitting to density """
        pass

    def add_terminal_residue_directly_using_cid(imol: int, cid: str):
        """ the cid is for an atom. This used to return a pair, but I removed it so that I could compile the binding. """
        pass

    def add_terminal_residue_directly_using_bucca_ml_growing_using_cid(imol: int, cid: str):
        """ the cid is for an atom. buccaneer building """
        pass

    def add_terminal_residue_directly_using_bucca_ml_growing(imol: int, spec: int ):
        """ buccaneer building, called by the above """
        pass

    def set_add_waters_water_to_protein_distance_lim_min(d: float):
        """ parameter for """
        pass

    def set_add_waters_water_to_protein_distance_lim_max(d: float):
        """ parameter for """
        pass

    def set_add_waters_variance_limit(d: float):
        """ parameter for """
        pass

    def set_add_waters_sigma_cutoff(d: float):
        """ parameter for """
        pass

    def add_waters(imol_model: int, imol_map: int):
        """ add waters, updating imol_model (of course) """
        pass

    def add_hydrogen_atoms(imol_model: int):
        """ add hydrogen atoms, updating imol_model (of course) """
        pass

    def delete_hydrogen_atoms(imol_model: int):
        """ delete hydrogen atoms, updating imol_model (of course) """
        pass

    def add_alternative_conformation(imol_model: int, cid: str):
        """ add an alternative conformation for the specified residue """
        pass

    def fill_partial_residue(imol: int, chain_id: str, res_no: int, ins_code: str):
        """ fill the specified residue """
        pass

    def fill_partial_residue_using_cid(imol: int, cid: str):
        """ fill the specified residue """
        pass

    def fill_partial_residues(imol: int):
        """ fill all the the partially-filled residues in the molecule """
        pass

    def flip_peptide(imol: int, atom_spec: int , alt_conf: str):
        """ flip peptide """
        pass

    def flip_peptide_using_cid(imol: int, atom_cid: str, alt_conf: str):
        """ flip peptide using an atom CID """
        pass

    def eigen_flip_ligand(imol: int, chain_id: str, res_no: int, ins_code: str):
        """ eigen-flip ligand """
        pass

    def eigen_flip_ligand_using_cid(imol: int, residue_cid: str):
        """ eigen-flip ligand using CID """
        pass

    def mutate(imol: int, cid: str, new_residue_type: str):
        """ mutate residue """
        pass

    def side_chain_180(imol: int, atom_cid: str):
        """ rotate last chi angle of the side chain by 180 degrees """
        pass

    def jed_flip(imol: int, atom_cid: str, invert_selection: bool):
        """ JED-Flip the ligand (or residue) at the specified atom. """
        pass

    def move_molecule_to_new_centre(imol: int, x: float, y: float, z: float):
        """ move the molecule to the given centre """
        pass

    def multiply_residue_temperature_factors(imol: int, cid: str, factor: float):
        """ Interactive B-factor refinement (fun). "factor" might typically be say 0.9 or 1.1 """
        pass

    def get_molecule_centre(imol: int):
        """ get molecule centre """
        pass

    def copy_fragment_using_cid(imol: int, multi_cid: str):
        """ copy a fragment given the multi_cid selection string. """
        pass

    def copy_fragment_for_refinement_using_cid(imol: int, multi_cid: str):
        """ copy a fragment - use this in preference to """
        pass

    def copy_fragment_using_residue_range(imol: int, chain_id: str, res_no_start: int, res_no_end: int):
        """ copy a residue-range fragment """
        pass

    def apply_transformation_to_atom_selection(imol: int, atoms_selection_cid: str, n_atoms: int, m00: float, m01: float, m02: float, m10: float, m11: float, m12: float, m20: float, m21: float, m22: float, c0: float, c1: float, c2: float, t0: float, t1: float, t2: float):
        """ apply transformation to atom selection in the given molecule. """
        pass

    def new_positions_for_residue_atoms(imol: int, residue_cid: str, moved_atoms: list):
        """ update the positions of the atoms in the residue """
        pass

    def new_positions_for_atoms_in_residues(imol: int, moved_residues: list):
        """ update the positions of the atoms in the residues """
        pass

    def merge_molecules(imol: int, list_of_other_molecules: str):
        pass

    def cis_trans_convert(imol: int, atom_cid: str):
        """ Convert a cis peptide to a trans or vice versa. """
        pass

    def replace_fragment(imol_base: int, imol_reference: int, atom_selection: str):
        """ replace a fragment"""
        pass

    def rigid_body_fit(imol: int, multi_cid: str, imol_map: int):
        """ `multi_cids" is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||/B/56-66" """
        pass

    def change_chain_id(imol: int, from_chain_id: str, to_chain_id: str, use_resno_range: bool, start_resno: int, end_resno: int):
        """ change the chain id """
        pass

    def refine_residues_using_atom_cid(imol: int, cid: str, mode: str, n_cycles: int):
        """ refine the residues """
        pass

    def refine_residues(imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str, mode: str, n_cycles: int):
        """ refine the residues """
        pass

    def refine_residue_range(imol: int, chain_id: str, res_no_start: int, res_no_end: int, n_cycles: int):
        """ refine residue range """
        pass

    def minimize_energy(imol: int, atom_selection_cid: str, n_cycles: int, do_rama_plot_restraints: bool, rama_plot_weight: float, do_torsion_restraints: bool, torsion_weight: float, refinement_is_quiet: bool):
        pass

    def fix_atom_selection_during_refinement(imol: int, atom_selection_cid: str):
        """ fix atoms during refinement. Does nothing at the moment. """
        pass

    def add_target_position_restraint(imol: int, atom_cid: str, pos_x: float, pos_y: float, pos_z: float):
        """ add or update (if it has a pull restraint already) """
        pass

    def clear_target_position_restraint(imol: int, atom_cid: str):
        """ clear target_position restraint """
        pass

    def turn_off_when_close_target_position_restraint(imol: int):
        """ clear target_position restraint if it is (or they are) close to their target position """
        pass

    def set_use_rama_plot_restraints(state: bool):
        """ turn on or off rama restraints """
        pass

    def get_use_rama_plot_restraints():
        """ get the state of the rama plot restraints usage in refinement. """
        pass

    def set_rama_plot_restraints_weight(f: float):
        """ set the Ramachandran plot restraints weight """
        pass

    def get_rama_plot_restraints_weight():
        """ get the Ramachandran plot restraints weight """
        pass

    def set_use_torsion_restraints(state: bool):
        """ turn on or off torsion restraints """
        pass

    def get_use_torsion_restraints():
        """ get the state of the rama plot restraints usage in refinement. """
        pass

    def set_torsion_restraints_weight(f: float):
        """ set the Ramachandran plot restraints weight """
        pass

    def get_torsion_restraints_weight():
        """ get the Ramachandran plot restraints weight """
        pass

    def init_refinement_of_molecule_as_fragment_based_on_reference(imol_frag: int, imol_ref: int, imol_map: int):
        """ initialise the refinement of (all of) molecule """
        pass

    def refine(imol: int, n_cycles: int):
        """ Run some cycles of refinement and return a mesh. That way we can see the molecule animate as it refines """
        pass

    def add_target_position_restraint_and_refine(imol: int, atom_cid: str, pos_x: float, pos_y: float, pos_z: float, n_cycles: int):
        """ Create a new position for the given atom and create a new bonds mesh based on that. This is currently "heavyweight" as the bonds mesh is calculated from scratch (it is not (yet) merely a distortion of an internally-stored mesh). """
        pass

    def clear_target_position_restraints(imol: int):
        """ clear any and all drag-atom target position restraints """
        pass

    def clear_refinement(imol: int):
        """ call this after molecule refinement has finished (say when the molecule molecule is accepted into the original molecule) """
        pass

    def set_refinement_is_verbose(state: bool):
        """ for debugging the refinement - write out some diagnositics - some might be useful. API change 20240226 - this function now takes a boolean argument """
        pass

    def set_refinement_geman_mcclure_alpha(a: float):
        """ set the refinement Geman-McClure alpha """
        pass

    def get_geman_mcclure_alpha():
        """ set the refinement Geman-McClure alpha """
        pass

    def generate_self_restraints(imol: int, local_dist_max: float):
        """ generate GM self restraints for the whole molecule """
        pass

    def generate_chain_self_restraints(imol: int, local_dist_max: float, chain_id: str):
        """ generate GM self restraints for the given chain """
        pass

    def generate_local_self_restraints(imol: int, local_dist_max: float, residue_cids: str):
        """ generate GM self restraints for the given residues. `residue_cids" is a "||"-separated list of residues, e.g. "//A/12||//A/14||/B/56" """
        pass

    def add_parallel_plane_restraint(imol: int, residue_cid_1: str, residue_cid_2: str):
        """ generate parallel plane restraints (for RNA and DNA) """
        pass

    def get_extra_restraints_mesh(imol: int, mode: int):
        """ get the mesh for extra restraints currently mode is unused. """
        pass

    def read_extra_restraints(imol: int, file_name: str):
        """ read extra restraints (e.g. from ProSMART) """
        pass

    def clear_extra_restraints(imol: int):
        """ clear the extra restraints """
        pass

    def get_rotamer_dodecs(imol: int):
        """ get the rotamer dodecs for the model, not const because it regenerates the bonds. """
        pass

    def get_rotamer_dodecs_instanced(imol: int):
        """ get the rotamer dodecs for the model, not const because it regenerates the bonds. """
        pass

    def get_ramachandran_validation_markup_mesh(imol: int):
        """ 20221126-PE: the function was renamed from """
        pass

    def ramachandran_validation(imol: int):
        """ get the data for Ramachandran validation, which importantly contains probability information """
        pass

    def contact_dots_for_ligand(imol: int, cid: str, smoothness_factor: int):
        """ Recently (20230202) the smoothness factor has been added as an extra argument """
        pass

    def all_molecule_contact_dots(imol: int, smoothness_factor: int):
        """ Recently (20230202) the smoothness factor has been added as an extra argument """
        pass

    def get_simple_molecule(imol: int, residue_cid: str, draw_hydrogen_atoms_flag: bool):
        pass

    def make_exportable_environment_bond_box(imol: int):
        pass

    def get_h_bonds(imol: int, cid_str: str, mcdonald_and_thornton_mode: bool):
        pass

    def get_mesh_for_ligand_validation_vs_dictionary(imol: int, ligand_cid: str):
        """ get the mesh for ligand validation vs dictionary, coloured by badness. greater then 3 standard deviations is fully red. Less than 0.5 standard deviations is fully green. """
        pass

    def match_ligand_torsions(imol_ligand: int, imol_ref: int, chain_id_ref: str, resno_ref: int):
        """ match ligand torsions - return the success status """
        pass

    def match_ligand_position(imol_ligand: int, imol_ref: int, chain_id_ref: str, resno_ref: int):
        """ match ligand positions - return the success status """
        pass

    def match_ligand_torsions_and_position(imol_ligand: int, imol_ref: int, chain_id_ref: str, resno_ref: int):
        """ match ligand torsions and positions """
        pass

    def match_ligand_torsions_and_position_using_cid(imol_ligand: int, imol_ref: int, cid: str):
        """ match ligand torsions and positions, different api """
        pass

    def get_overlap_dots(imol: int):
        """ not const because it can dynamically add dictionaries """
        pass

    def get_overlap_dots_for_ligand(imol: int, cid_ligand: str):
        """ not const because it can dynamically add dictionaries """
        pass

    def get_overlaps(imol: int):
        """ not const because it can dynamically add dictionaries """
        pass

    def get_overlaps_for_ligand(imol: int, cid_ligand: str):
        """ not const because it can dynamically add dictionaries """
        pass

    def density_fit_analysis(imol_model: int, imol_map: int):
        """ density fit validation information """
        pass

    def density_correlation_analysis(imol_model: int, imol_map: int):
        """ density correlation validation information """
        pass

    def rotamer_analysis(imol_model: int):
        """ rotamer validation information """
        pass

    def ramachandran_analysis(imol_model: int):
        """ ramachandran validation information (formatted for a graph, not 3d) """
        pass

    def ramachandran_analysis_for_chain(imol_model: int, chain_id: str):
        """ ramachandran validation information (formatted for a graph, not 3d) for a given chain in a given molecule 20230127-PE This function does not exist yet."""
        pass

    def peptide_omega_analysis(imol_model: int):
        """ peptide omega validation information """
        pass

    def get_median_temperature_factor(imol: int):
        """ get the median temperature factor for the model """
        pass

    def get_interesting_places(imol: int, mode: str):
        """ get interesting places (does not work yet) """
        pass

    def difference_map_peaks(imol_map: int, imol_protein: int, n_rmsd: float):
        """ get difference map peaks """
        pass

    def pepflips_using_difference_map(imol_coords: int, imol_difference_map: int, n_sigma: float):
        """ get pepflips based on the difference map """
        pass

    def unmodelled_blobs(imol_model: int, imol_map: int):
        """ unmodelled blobs """
        pass

    def find_water_baddies(imol_model: int, imol_map: int, b_factor_lim: float, outlier_sigma_level: float, min_dist: float, max_dist: float, ignore_part_occ_contact_flag: bool, ignore_zero_occ_flag: bool):
        """ typical values for """
        pass

    def fourier_shell_correlation(imol_map_1: int, imol_map_2: int):
        """ Calculate the MMRRCC for the residues in the chain Multi Masked Residue Range Corellation Coefficient calculate the MMRRCC for the residues in the chain Multi Masked Residue Range Corellation Coefficient Fourier Shell Correlation (FSC) between maps """
        pass

    def calculate_new_rail_points():
        """ calling this adds to the rail_points history. Make this pairs when we add model scoring. """
        pass

    def rail_points_total():
        """ the total rail points """
        pass

    def associate_data_mtz_file_with_map(imol: int, data_mtz_file_name: str, f_col: str, sigf_col: str, free_r_col: str):
        """ call this before calling """
        pass

    def connect_updating_maps(imol_model: int, imol_with_data_info_attached: int, imol_map_2fofc: int, imol_map_fofc: int):
        """ reset the rail_points (calls reset_the_rail_points()), updates the maps (using internal/clipper SFC) so, update your contour lines meshes after calling this function. """
        pass

    def sfcalc_genmap(imol_model: int, imol_map_with_data_attached: int, imol_updating_difference_map: int):
        """ sfcalc and re-generate maps. This is a low-level function - generally one would use the updating maps method rather than this """
        pass

    def sfcalc_genmaps_using_bulk_solvent(imol_model: int, imol_2fofc_map: int, imol_updating_difference_map: int, imol_map_with_data_attached: int):
        """ Call this function after connecting maps for updating maps to set the initial R-factor and store the initial map flatness."""
        pass

    def shift_field_b_factor_refinement(imol: int, imol_with_data_attached: int):
        """ shift_field B-factor refinement. This function presumes that the Fobs,sigFobs and RFree data have been filled in the """
        pass

    def get_density_at_position(imol_map: int, x: float, y: float, z: float):
        """ get density at position """
        pass

    def get_diff_diff_map_peaks(imol_diff_map: int, screen_centre_x: float, screen_centre_y: float, screen_centre_z: float):
        pass

    def get_data_set_file_name(imol: int):
        """ the stored data set file name """
        pass

    def go_to_blob(x1: float, y1: float, z1: float, x2: float, y2: float, z2: float, contour_level: float):
        """ 20221022-PE: in future, this function should/will be provided with a list of displayed maps and their contour levels - but for now, it uses (only) imol_refinement_map. Use a string to pass the map information (map index and contour level), something like "0 0.45:1 1.2:2 0.1" """
        pass

    def fit_ligand_right_here(imol_protein: int, imol_map: int, imol_ligand: int, x: float, y: float, z: float, n_rmsd: float, use_conformers: bool, n_conformers: int):
        """ For trivial (i.e non-flexible) ligands you should instead use the jiggle-fit algorithm, which takes a fraction of a second. (That is the algorithm used for "Add Other Solvent Molecules" in Coot.)"""
        pass

    def fit_ligand(imol_protein: int, imol_map: int, imol_ligand: int, n_rmsd: float, use_conformers: bool, n_conformers: int):
        """ Fit ligand """
        pass

    def fit_ligand_multi_ligand(imol_protein: int, imol_map: int, multi_ligand_molecule_number_list: str, n_rmsd: float, use_conformers: bool, n_conformers: int):
        """ Fit ligands (place-holder) """
        pass

    def fit_to_map_by_random_jiggle(imol: int, res_spec: int , n_trials: int, translation_scale_factor: float):
        """ "Jiggle-Fit Ligand" if n_trials is 0, then a sensible default value will be used. if translation_scale_factor is negative then a sensible default value will be used. """
        pass

    def fit_to_map_by_random_jiggle_using_cid(imol: int, cid: str, n_trials: int, translation_scale_factor: float):
        """ "Jiggle-Fit Ligand" with a different interface - one that can use any atom selection (instead of just a ligand). As above, if n_trials is 0, then a sensible default value will be used. if translation_scale_factor is negative then a sensible default value will be used. """
        pass

    def fit_to_map_by_random_jiggle_with_blur_using_cid(imol: int, imol_map: int, cid: str, b_factor: float, n_trials: int, translation_scale_factor: float):
        """ Jiggle fit an atom selection, typically a whole molecule or a chain As above, if n_trials is 0, then a sensible default value will be used. if translation_scale_factor is negative then a sensible default value will be used. """
        pass

    def get_svg_for_residue_type(imol: int, comp_id: str, use_rdkit_svg: bool, dark_background_flag: bool):
        """ This function is not const because it caches the svgs if it can."""
        pass

    def add_compound(imol: int, tlc: str, imol_dict: int, imol_map: int, x: float, y: float, z: float):
        """ This function is for adding compounds/molecules like buffer agents and precipitants or anions and cations. """
        pass

    def get_non_standard_residues_in_molecule(imol: int):
        pass

    def get_map_section_texture(imol: int, section_id: int, axis: int, data_value_for_bottom: float, data_value_for_top: float):
        """ The new arguments, """
        pass

    def get_number_of_map_sections(imol_map: int, axis_id: int):
        pass

    def make_mesh_from_gltf_file(file_name: str):
        pass

    def get_octahemisphere(n_divisions: int):
        """ @params """
        pass

    def pae_png(pae_file_name: str):
        pass


    def set_updating_maps_need_an_update(imol: int):
        pass

    def update_updating_maps(imol: int):
        """ update the updating maps without generating a mesh """
        pass

    def adjust_refinement_residue_name(resname: str):
        pass


    # def find_serial_number_for_insert(seqnum_new: int, ins_code_for_new: str, chain_p: mmdb::Chain *):
    #     pass

    def read_standard_residues():
        """ read the standard protein, RNA, and DNA dictionaries. """
        pass

    def install_model(m: int ):
        pass

    def overlap_ligands_internal(imol_ligand: int, imol_ref: int, chain_id_ref: str, resno_ref: int, apply_rtop_flag: bool):
        pass

    # def superpose_with_atom_selection(asc_ref: atom_selection_container_t, asc_mov: atom_selection_container_t, imol_mov: int, moving_mol_name: str, reference_mol_name: str, move_copy_of_imol2_flag: bool):
    #     pass

    def valid_labels(mtz_file_name: str, f_col: str, phi_col: str, weight_col: str, use_weights: int):
        pass

    def init():
        pass

    def debug():
        pass

    def thread_for_refinement_loop_threaded():
        pass

    def get_restraints_lock(calling_function_name: str):
        pass

    def release_restraints_lock(calling_function_name: str):
        pass

    def all_atom_pulls_off():
        pass

    def atom_pull_off(spec: int ):
        pass

    def atom_pulls_off(specs: list):
        pass

    def molecules_container_t(verbose: bool):
        """ the one and only constructor """
        pass

    def molecules_container_t():
        pass

    def set_use_gemmi(state: bool):
        """ Set the state of using gemmi for coordinates parsing. The default is true. """
        pass

    def get_use_gemmi():
        """ get the state of using GEMMI for coordinates parsing """
        pass

