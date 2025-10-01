class lsq_results_t:
        pass

class molecules_container_t:

    def package_version(self) -> str:
        """ Get the package version

        :return: the package version, e.g. "1.1.11" - if this is a not yet a release version the version will end in a "+", such as "1.1.11+" """
        return 'a-string'

    def set_use_gemmi(self, state: bool) -> None:
        """ Set the state of using GEMMI for coordinates parsing

        :param state:  is True to mean that it is enabled. The default is True. """

    def get_use_gemmi(self) -> bool:
        """ Get the state of using GEMMI for coordinates parsing. """
        return True

    def set_make_backups(self, state: bool) -> None:
        """ Allow the user to disable/enable backups

        :param state:  is True to mean that it is enabled. The default is True. """

    def get_make_backups(self) -> bool:
        """ Get the state of the backups

        :return: the backup-enabled state """
        return True

    def file_name_to_string(self, file_name: str) -> str:
        """ File name to string

        :param file_name:  is the name of the given file

        :return: the string of the contents of the given file-name. """
        return 'a-string'

    def get_number_of_molecules(self) -> int:
        """ Get the number of molecules (either map or model)

        :return: the number of molecules """
        return 0

    def create_empty_molecules(self, n_empty: int) -> None:
        """ Add a number of empty molecules to the internal vector/list of molecules

        Note this is not like STL `reserve` as it will increase the molecule index of the next added molecule by `n_empty`.

        :param n_empty:  the number of empty molecules to create """

    def set_imol_refinement_map(self, i: int) -> None:
        """ Set the map used for refinement and fitting

        :param i:  the map molecule index used for refinement and fitting """

    def set_map_weight(self, w: float) -> None:
        """ Set the map weight

        :param w:  the map weight to be used for refinement, e.g. 50.0 """

    def get_map_weight(self) -> float:
        """ Get the map weight

        :return: the map weight """
        return 0.0

    def scale_map(self, imol_map: int, scale_factor: float) -> None:
        """ Scale map

        :param imol:  is the model molecule index 

        :param scale_factor:  is the scale factor """

    def atom_cid_to_atom_spec(self, imol: int, cid: str):
        """ Convert atom cid string to a coot atom specifier

        :param imol:  is the model molecule index 

        :param cid:  is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)

        :return: the atom spec., `spec.empty()` is true on failure. """
        pass

    def residue_cid_to_residue_spec(self, imol: int, cid: str):
        """ Convert residue cid string to a coot residue specifier

        :param imol:  is the model molecule index 

        :param cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)

        :return: the residues spec., `spec.empty()` is true on failure. """
        pass

    def set_show_timings(self, s: bool) -> None:
        """ Set the show_timings flag

        Various (not all) functions in this class can calculate how long they took to run. Setting this will write the time to taken (in milliseconds) to stdout.

        :param s:  is True to mean that it is enabled. The default is True. """

    def get_header_info(self, imol: int):
        """ Get header info

        (the header info is sparce at the moment)

        :param imol:  is the model molecule index

        :return: an object with header info. """
        pass

    def get_imol_enc_any(self) -> int:
        """ Get imol_enc_any (enc: encoded)

        imol_enc_any refers to the molecule number for dictionary that can be used with any molecule

        :return: the value of `imol_enc_any` """
        return 0


    def read_coordinates(self, file_name: str) -> int:
        """ Read a coordinates file (mmcif or PDB)

        :param file_name:  is the name of the coordinates file

        :return: the new molecule index on success and -1 on failure """
        return 0

    def read_pdb(self, file_name: str) -> int:
        """ Read a PDB file (or mmcif coordinates file, despite the name)

        It does the same job as `read_coordinates` but has (perhaps) a more familiar name

        :param file_name:  is the name of the coordinates file

        :return: the new molecule index on success and -1 on failure """
        return 0

    def read_small_molecule_cif(self, file_name: str) -> int:
        """ Read a small molecule CIF file

        :param file_name:  is the cif file-name

        :return: the new molecule index on success and -1 on failure """
        return 0

    def print_secondary_structure_info(self, imol: int) -> None:
        """ Print the secondary structure information

        :param imol:  is the model molecule index """

    def replace_molecule_by_model_from_file(self, imol: int, pdb_file_name: str) -> None:
        """ Read a PDB file (or mmcif coordinates file, despite the name) to replace the current molecule

        This will only work if the molecules is already a model molecule

        :param imol:  is the model molecule index 

        :param file_name:  is the name of the coordinates file """

    def split_multi_model_molecule(self, imol: int):
        """ Split an NMR model into multiple models

        :param imol:  is the model molecule index

        :return: the vector/list of new molecule indices """
        pass

    def make_ensemble(self, model_molecules_list: str) -> int:
        """ Make a multi-model molecule given the input molecules

        :param model_molecules_list:  is a colon-separated list of molecules, (e.g. "2:3:4")

        :return: the new molecule index, -1 if no models were found in the model molecules list """
        return 0

    def molecule_to_PDB_string(self, imol: int) -> str:
        """ Get the molecule as a PDB string

        :param imol:  is the model molecule index

        :return: the model molecule imol as a string. Return emtpy string on error """
        return 'a-string'

    def molecule_to_mmCIF_string(self, imol: int) -> str:
        """ Get the molecule as an mmCIF string

        :param imol:  is the model molecule index

        :return: the model molecule imol as a string. Return emtpy string on error """
        return 'a-string'

    def get_active_atom(self, x: float, y: float, z: float, displayed_model_molecules_list: str):
        """ Get the active atom given the screen centre

        :param x:  is x position of the centre of the screen 

        :param y:  is y position of the centre of the screen 

        :param z:  is z position of the centre of the screen 

        :param displayed_model_molecules_list:  is a colon-separated list of molecules, (e.g. "2:3:4")

        :return: the molecule index and the atom cid. On failure (no molecules with atoms in them) then return -1 and a blank string. """
        pass

    def import_cif_dictionary(self, cif_file_name: str, imol_enc: int) -> int:
        """ Import a dictionary cif

        :param imol_enc:  is used to specify to which molecule this dictionary should apply. Use IMOL_ENC_ANY to mean "it applies to all molecules." IMOL_ENC_ANY = -999999

        :return: 1 on success and 0 on failure """
        return 0

    def get_cif_file_name(self, comp_id: str, imol_enc: int) -> str:
        """ Get the cif file name

        :param comp_id:  is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine 

        :param imol_enc:  is the molecule index for the residue type/compound_id

        :return: the dictionary read for the given residue type, return an empty string on failure """
        return 'a-string'

    def get_cif_restraints_as_string(self, comp_id: str, imol_enc: int) -> str:
        """ Get the cif restraints as a string

        :param comp_id:  is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine 

        :param imol_enc:  is the molecule index for the residue type/compound_id

        :return: a string that is the contents of a dictionary cif file """
        return 'a-string'

    def copy_dictionary(self, monomer_name: str, imol_current: int, imol_new: int) -> bool:
        """ Copy the dictionary that is specific for `imol_current` so that it can be used with a new molecule

        :param monomer_name:  is the 3 letter code of the monomer in the dictionary, e.g. "ALA" for alanine 

        :param imol_current:  is the model molecule index with the dictionary to be copied from 

        :param imol_new:  is the model molecule index the dictionary will be copied into """
        return True

    def get_monomer(self, monomer_name: str) -> int:
        """ Get a monomer

        :param monomer_name:  is the 3-letter code of the monomer in the dictionary, e.g. "ALA" for alanine

        :return: the new molecule index on success and -1 on failure """
        return 0

    def get_monomer_from_dictionary(self, comp_id: str, imol: int, idealised_flag: bool) -> int:
        """ Get a monomer for a particular molecule

        :param comp_id:  is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine 

        :param imol:  is the model molecule index, use -999999 (IMOL_ENC_ANY) if no molecule-specific dictionary is needed 

        :param idealised_flag:  means that the coordinates have been minimised with a molecular modelling minimisation algo, usually the value is True

        :return: the new molecule index on success and -1 on failure """
        return 0

    def get_monomer_and_position_at(self, comp_id: str, imol: int, x: float, y: float, z: float) -> int:
        """ Get monomer and place it at the given position for a particular molecule

        :param comp_id:  is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine 

        :param imol:  is the model molecule index, use -999999 (IMOL_ENC_ANY) if no molecule-specific dictionary is needed 

        :param x:  is the x value of the target position 

        :param y:  is the y value of the target position 

        :param z:  is the z value of the target position

        :return: the new molecule index on success and -1 on failure """
        return 0

    def dictionary_atom_name_map(self, comp_id_1: str, imol_1: int, comp_id_2: str, imol_2: int):
        """ Match atom between 2 dictionaries

        :param comp_id_1:  is the 3-letter code for the residue/ligand in the first model, e.g. "ALA" for alanine 

        :param imol_1:  is the model molecule index of the first model 

        :param comp_id_2:  is the 3-letter code for the residue/ligand in the second model, e.g. "ALA" for alanine 

        :param imol_2:  is the model molecule index of the second model

        :return: the atom name match on superposing the atoms of the given dictionaries """
        pass

    def get_types_in_molecule(self, imol: int):
        """ get types """
        pass

    def get_groups_for_monomers(self, residue_names: list):
        """ Get the groups for a vector/list of monomers

        e.g. "NON POLYMER", "PEPTIDE", etc

        :param residue_names:  is a list of residue names, e.g. ["ALA", "TRP"]

        :return: the group for the given list of residue names """
        pass

    def get_group_for_monomer(self, residue_name: str) -> str:
        """ Get the group for a particular monomer

        e.g. "NON POLYMER", "PEPTIDE", etc 

        :param residue_name:  is the is the 3-letter code for the residue, e.g. "ALA" for alanine

        :return: the group for the given residue name """
        return 'a-string'

    def get_hb_type(self, compound_id: str, imol_enc: int, atom_name: str) -> str:
        """ Get the hydrogen bond type of a particular atom in a given residue type

        :param compound_id:  is the 3-letter code for the residue/ligand in the first model, e.g. "TYR" for tyrosine 

        :param imol_enc:  is the molecule index for the residue type/compound_id 

        :param atom_name:  is the name of the atom, e.g. "OH"

        :return: the hb_type for the given atom. On failure return an empty string. Valid types are: "HB_UNASSIGNED" ,"HB_NEITHER", "HB_DONOR", "HB_ACCEPTOR", "HB_BOTH", "HB_HYDROGEN". """
        return 'a-string'

    def get_gphl_chem_comp_info(self, compound_id: str, imol_enc: int):
        """ Get the GPhL extra restraint information (from the input cif file)

        :param compound_id:  is the 3-letter code for the residue/ligand in the first model, e.g. "TYR" for tyrosine 

        :param imol_enc:  is the molecule index for the residue type/compound_id

        :return: a vector/list of string pairs that were part of a gphl_chem_comp_info. return an empty vector/list on failure to find any such info. """
        pass

    def get_acedrg_atom_types(self, compound_id: str, imol_enc: int):
        """ Get a list of atom names and their associated AceDRG atom types

        :param compound_id:  is the 3-letter code for the residue/ligand in the first model, e.g. "TYR" for tyrosine 

        :param imol_enc:  is the molecule index for the residue type/compound_id

        :return: a list of atom names and their associated AceDRG atom types, return an empty list on failure (e.g. when atoms types are not in the dictionary) """
        pass

    def get_acedrg_atom_types_for_ligand(self, imol: int, residue_cid: str):
        """ Get AceDRG atom types for ligand bonds

        :param imol:  is the model molecule index 

        :param residue_cid:  is the atom selection CID e.g "//A/15" (residue 15 of chain A)

        :return: a `coot::acedrg_types_for_residue_t` - which contains a vector/list of bond descriptions. """
        pass

    def set_occupancy(self, imol: int, cid: str, occ_new: float) -> None:
        """ Set the occupancy for the given atom selection

        :param imol:  is the model molecule index 

        :param cid:  is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A) 

        :param occ_new:  is the new occupancy """

    def write_png(self, compound_id: str, imol: int, file_name: str) -> None:
        """ Write a PNG for the given compound_id.

        Currently this function does nothing (drawing is done with the not-allowed cairo)

        :param compound_id:  is the 3-letter code for the residue/ligand in the first model, e.g. "TYR" for tyrosine 

        :param imol:  is the model molecule index """

    def write_coordinates(self, imol: int, file_name: str) -> int:
        """ Write the coordinates to the given file name

        :param imol:  is the model molecule index 

        :param file_name:  is the name of the new model

        :return: 1 on success and 0 on failure """
        return 0

    def set_draw_missing_residue_loops(self, state: bool) -> None:
        """ Set the state for drawing missing loops

        By default missing loops are drawn. This function allows missing loops to not be drawn. Sometimes that can clarify the representation. This is a lightweight function that sets a flag that is used by subsequent calls to 

        :param state:  is True to mean it is enabled """

    def get_bonds_mesh(self, imol: int, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int):
        """ Get the bonds mesh

        :param mode:  is "COLOUR-BY-CHAIN-AND-DICTIONARY", "CA+LIGANDS" or "VDW-BALLS" 

        :param against_a_dark_background:  allows the bond colours to be relevant for the background. When the background is dark, the colours should (as a rule) be bright and pastelly. When the background is light/white, the colour are darker and more saturated. 

        :param smoothness_factor:  controls the number of triangles used to make the bond cylinders and spheres for the atoms - it rises in powers of 4. 1 is the smallest smoothness_factor, 2 looks nice and 3 is best. 

        :param bond_width:  is the bond width in Angstroms. 0.12 is a reasonable default value. 

        :param atom_radius_to_bond_width_ratio:  allows the representation of "ball and stick". To do so use a value between 1.5 and 3.0. The ratio for "liquorice" representation is 1.0.

        :return: a `simple_mesh_t` """
        pass

    def get_bonds_mesh_instanced(self, imol: int, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, show_atoms_as_aniso_flag: bool, show_aniso_atoms_as_ortep_flag: bool, draw_hydrogen_atoms_flag: bool, smoothness_factor: int):
        """ Get the instanced bonds mesh.

        :param mode:  is "COLOUR-BY-CHAIN-AND-DICTIONARY" - more modes to follow 

        :param against_a_dark_background:  allows the bond colours to be relevant for the background. When the background is dark, the colours should (as a rule) be bright and pastelly. When the background is light/white, the colour are darker and more saturated. 

        :param bond_width:  is the bond width in Angstroms. 0.12 is a reasonable default value. 

        :param atom_radius_to_bond_width_ratio:  allows the representation of "ball and stick". To do so use a value between 1.5 and 3.0. The ratio for "liquorice" representation is 1.0. 

        :param show_atoms_as_aniso_flag:  if true, if possible, show the atoms with thermal ellipsoids. 

        :param show_aniso_atoms_as_ortep_flag:  if true, show any anisotropic atoms with ortep style. 

        :param draw_hydrogen_atoms_flag:  if true, bonds to hydrogen atoms should be added. 

        :param smoothness_factor:  controls the number of triangles used to make the bond cylinders and spheres for the atoms - it rises in powers of 4. 1 is the smallest smoothness_factor, 2 looks nice and 3 is best. Instancing may mean that smoothness factor 3 should be used by default. 

        :return: a `instanced_mesh_t` """
        pass

    def get_bonds_mesh_for_selection_instanced(self, imol: int, atom_selection_cid: str, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, show_atoms_as_aniso_flag: bool, show_aniso_atoms_as_ortep_flag: bool, draw_hydrogen_atoms_flag: bool, smoothness_factor: int):
        """ As `get_bonds_mesh_instanced` above, but only return the bonds for the atom selection. Typically one would call this with a wider bond width than one would use for standards atoms (all molecule)

        :param atom_selection_cid:  e.g. "//A/15" (all the atoms in residue 15 of chain A) 

        :return: a `instanced_mesh_t` """
        pass

    def get_goodsell_style_mesh_instanced(self, imol: int, colour_wheel_rotation_step: float, saturation: float, goodselliness: float):
        """ Get the Goodsell style mesh

        :param colour_wheel_rotation_step':  the amount, in degrees, that the colour wheel advances between different chains. 97 degrees is a reasonable starting value 

        :param saturation':  a number between 0 and 1, where 0 is grey and 1 is "lego-like" colour scheme. 0.5 is a nice middle value 

        :param goodselliness:  is the degree to which the C-atoms are desaturated, 0.3 is a reasonable value

        :return: a `instanced_mesh_t` """
        pass

    def export_map_molecule_as_gltf(self, imol: int, pos_x: float, pos_y: float, pos_z: float, radius: float, contour_level: float, file_name: str) -> None:
        """ Export map molecule as glTF file

        glTF files can be imported into Blender or other 3D graphics applications

        :param imol:  is the model molecule index 

        :param x:  x of the center of the screen 

        :param y:  y of the center of the screen 

        :param z:  z of the center of the screen 

        :param radius:  e.g. 12.0 for X-ray map and 100.0 for Cryo-EM map 

        :param contour_level:  e.g. 1.5 sd for X-ray, 5.0 sd for cryo-EM 

        :param file_name:  extension should be .glb """

    def export_model_molecule_as_gltf(self, imol: int, selection_cid: str, mode: str, against_a_dark_background: bool, bonds_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int, draw_hydrogen_atoms_flag: bool, draw_missing_residue_loops: bool, file_name: str) -> None:
        """ Export model molecule as glTF file

        glTF files can be imported into Blender or other 3D graphics applications

        Same parameters as the `get_bonds_mesh` function. `draw_hydrogen_atoms_flag` and `draw_missing_residue_loops` are typically False. This API will change - we want to specify surfaces and ribbons too. """

    def export_molecular_representation_as_gltf(self, imol: int, atom_selection_cid: str, colour_scheme: str, style: str, secondary_structure_usage_flag: int, file_name: str) -> None:
        """ Export molecular representation as glTF file

        glTF files can be imported into Blender

        :param imol:  is the model molecule index 

        :param atom_selection_cid:  : e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param colour_scheme:  is one of "colorRampChainsScheme", "colorBySecondaryScheme", "Chain" 

        :param style:  "Ribbon" or "MolecularSurface" 

        :param secondary_structure_usage_flag:  0 (USE_HEADER), 1 (DONT_USE) or 2 (CALC_SECONDARY_STRUCTURE) 

        :param file_name:  of the glTF (the file will be compressed, so choose ".glb" as the extension) """

    def export_chemical_features_as_gltf(self, imol: int, cid: str, file_name: str) -> None:
        """ export chemical features for the specified residue """

    def set_gltf_pbr_roughness_factor(self, imol: int, roughness_factor: float) -> None:
        """ set the gltf PBR roughness factor

        :param imol:  is the model molecule index 

        :param roughness_factor:  is the factor for the roughness (0.0 to 1.0) """

    def set_gltf_pbr_metalicity_factor(self, imol: int, metalicity: float) -> None:
        """ set the gltf PBR metalicity factor

        :param imol:  is the model molecule index 

        :param metalicity:  is the factor for the roughness (0.0 to 1.0) """

    def get_colour_table(self, imol: int, against_a_dark_background: bool):
        """ Get colour table (for testing)

        :param imol:  is the model molecule index 

        :param against_a_dark_background:  allows the bond colours to be relevant for the background.

        :return: the colour table """
        pass

    def set_colour_wheel_rotation_base(self, imol: int, r: float) -> None:
        """ Set the colour wheel rotation base for the specified molecule (in degrees)

        :param imol:  is the model molecule index 

        :param r:  is the rotation angle in degrees """

    def set_base_colour_for_bonds(self, imol: int, r: float, g: float, b: float) -> None:
        """ Set the base colour - to be used as a base for colour wheel rotation

        RGB colour codes, e.g. green is r:0, g: 255, b:0

        :param imol:  is the model molecule index """

    def add_to_non_drawn_bonds(self, imol: int, atom_selection_cid: str) -> None:
        """ Add an atom selection cid for atoms and bonds not to be drawn

        :param imol:  is the model molecule index 

        :param atom_selection_cid:  e.g "//A/15" (all the atoms in residue 15 of chain A) """

    def clear_non_drawn_bonds(self, imol: int) -> None:
        """ Clear the set of non-drawn atoms (so that they can be displayed again)

        :param imol:  is the model molecule index """

    def print_non_drawn_bonds(self, imol: int) -> None:
        """ Print non-drawn bonds

        :param imol:  is the model molecule index """

    def set_user_defined_bond_colours(self, imol: int, colour_map: dict) -> None:
        """ User-defined colour-index to colour. """

    def set_user_defined_atom_colour_by_selection(self, imol: int, indexed_residues_cids: list, colour_applies_to_non_carbon_atoms_also: bool) -> None:
        """ Set the user-defined residue selections (CIDs) to colour index

        :param imol:  is the model molecule index """

    def add_colour_rule(self, imol: int, selection_cid: str, colour: str) -> None:
        """ Add a colour rule for M2T representations. """

    def add_colour_rules_multi(self, imol: int, selections_and_colours_combo_string: str) -> None:
        """ Add multiple colour rules

        :param selections_and_colours_combo_string:  e.g. "//A/1^#cc0000|//A/2^#cb0002|//A/3^#c00007", where "|" is the separator for each rule and "^" is the separator for the selection string and the colour string """

    def delete_colour_rules(self, imol: int) -> None:
        """ Delete the colour rules for the given molecule

        :param imol:  is the model molecule index """

    def get_colour_rules(self, imol: int):
        """ Get the colour rules

        :param imol:  is the model molecule index """
        pass

    def print_colour_rules(self, imol: int) -> None:
        """ Print the colour rules

        :param imol:  is the model molecule index """

    def set_use_bespoke_carbon_atom_colour(self, imol: int, state: bool) -> None:
        """ Use bespoke carbon atom colour

        :param imol:  is the model molecule index """

    def set_bespoke_carbon_atom_colour(self, imol: int, col: str) -> None:
        """ Set bespoke carbon atom colour

        :param imol:  is the model molecule index """

    def M2T_updateFloatParameter(self, imol: int, param_name: str, value: float) -> None:
        """ Update float parameter for MoleculesToTriangles molecular mesh. """

    def M2T_updateIntParameter(self, imol: int, param_name: str, value: int) -> None:
        """ Update int parameter for MoleculesToTriangles molecular mesh. """

    def get_molecular_representation_mesh(self, imol: int, cid: str, colour_scheme: str, style: str, secondary_structure_usage_flag: int):
        """ Get ribbon and molecular surface representation

        :param imol:  is the model molecule index 

        :param cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param colour_scheme:  should be one of "colorRampChainsScheme", "colorBySecondaryScheme", "Chain" 

        :param style:  "Ribbon" or "MolecularSurface" 

        :param secondary_structure_usage_flag:  0 (USE_HEADER), 1 (DONT_USE) or 2 (CALC_SECONDARY_STRUCTURE).

        :return: a `simple_mesh_t` """
        pass

    def get_gaussian_surface(self, imol: int, sigma: float, contour_level: float, box_radius: float, grid_scale: float, b_factor: float):
        """ Get a Gaussian surface representation

        :param imol:  is the model molecule index 

        :param sigma:  default 4.4 

        :param contour_level:  default 4.0 

        :param box_radius:  default 5.0 

        :param grid_scale:  default 0.7 

        :param b_factor:  default 100.0 (use 0.0 for no FFT-B-factor smoothing)

        :return: a `simple_mesh_t` composed of a number of Gaussian surfaces (one for each chain) """
        pass

    def get_chemical_features_mesh(self, imol: int, cid: str):
        """ Get chemical features for the specified residue

        :param imol:  is the model molecule index 

        :param cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)

        :return: a `simple_mesh_t` """
        pass

    def get_atom_using_cid(self, imol: int, cid: str):
        """ get an (mmdb-style) atom

        If more than one atom is selected by the selection cid, then the first atom is returned.

        Don't use this in emscript.

        :param imol:  is the model molecule index 

        :param cid:  is the coordinate-id for the atom. 

        :return: either the specified atom or nullopt (None) if not found """
        pass

    def get_residue_using_cid(self, imol: int, cid: str):
        """ get an (mmdb-style) residue

        If more than one residue is selected by the selection cid, then the first residue is returned.

        Don't use this in emscript.

        :param imol:  is the model molecule index 

        :param cid:  is the coordinate-id for the residue 

        :return: either the specified residue or nullopt (None) if not found """
        pass

    def residue_is_nucleic_acid(self, imol: int, cid: str) -> bool:
        """ Residue is nucleic acid?

        Every residue in the selection is checked

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A)

        :return: a bool """
        return True

    def get_residue_CA_position(self, imol: int, cid: str):
        """ Get the residue CA position

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A)

        :return: a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values """
        pass

    def get_residue_average_position(self, imol: int, cid: str):
        """ Get the average residue position

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A)

        :return: a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values """
        pass

    def get_residue_sidechain_average_position(self, imol: int, cid: str):
        """ Get the average residue side-chain position

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A)

        :return: a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values """
        pass

    def get_number_of_atoms(self, imol: int) -> int:
        """ Get number of atoms

        :param imol:  is the model molecule index

        :return: the number of atoms in the specified model, or 0 on error """
        return 0

    def get_molecule_diameter(self, imol: int) -> float:
        """ Get molecule diameter

        :param imol:  is the model molecule index

        :return: an estimate of the diameter of the model molecule (-1 on failure) """
        return 0.0

    def get_number_of_hydrogen_atoms(self, imol: int) -> int:
        """ Get number of hydrogen atoms

        :param imol:  is the model molecule index

        :return: the number of hydrogen atoms in the specified model, or -1 on error """
        return 0

    def get_chains_in_model(self, imol: int):
        """ Get the chain IDs in the given molecule

        :param imol:  is the model molecule index

        :return: vector/list of chain-ids for the given molecule """
        pass

    def get_ncs_related_chains(self, imol: int):
        """ Get the chains that are related by NCS or molecular symmetry

        :param imol:  is the model molecule index

        :return: a vector/list of vector/list of chain ids, e.g. [[A,C], [B,D]] (for hemoglobin). """
        pass

    def get_single_letter_codes_for_chain(self, imol: int, chain_id: str):
        """ Get the single letter codes for the residues in the specified chain

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A

        :return: vector/list of single letter codes - in a pair with the given residue spec """
        pass

    def get_residue_names_with_no_dictionary(self, imol: int):
        """ Get a list of residues that don't have a dictionary

        :param imol:  is the model molecule index

        :return: a list of residue names (compound_ids) for which there is no dictionary in the geometry store """
        pass

    def get_residue_name(self, imol: int, chain_id: str, res_no: int, ins_code: str) -> str:
        """ Get residue name

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A"

        :return: the residue name, return a blank string on residue not found. """
        return 'a-string'

    def get_SMILES_for_residue_type(self, residue_name: str, imol_enc: int) -> str:
        """ Get the SMILES string for the give residue type

        :param residue:  3 letter-code/name of the compound-id 

        :param imol_enc:  is the molecule index for the residue type/compound_id

        :return: the SMILES string if the residue type can be found in the dictionary store or the empty string on a failure. """
        return 'a-string'

    def residues_with_missing_atoms(self, imol: int):
        """ Get residues with missing atoms

        :param imol:  is the model molecule index

        :return: an object that has information about residues without dictionaries and residues with missing atom in the the specified molecule """
        pass

    def get_missing_residue_ranges(self, imol: int):
        """ Get missing residue ranges

        :param imol:  is the model molecule index 

        :return: missing residue ranges """
        pass

    def get_residues_near_residue(self, imol: int, residue_cid: str, dist: float):
        """ Get a list of residues specs that have atoms within distance of the atoms of the specified residue

        :param imol:  is the model molecule index 

        :param residue_cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param dist:  is the distance in Angstrom

        :return: a list of residue specs """
        pass

    def get_distances_between_atoms_of_residues(self, imol: int, cid_res_1: str, cid_res_2: str, dist_max: float):
        """ Get atom distances

        :param imol:  is the model molecule index 

        :param cid_res_1:  is the first atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A) 

        :param cid_res_2:  is the second atom selection CID e.g "//A/17/NH" (atom NH in residue 17 of chain A) 

        :param dist:  is the distance in Angstrom """
        pass

    def SSM_superpose(self, imol_ref: int, chain_id_ref: str, imol_mov: int, chain_id_mov: str):
        """ Superposition (using SSM)

        The specified chain of the moving molecule is superposed onto the chain in the reference molecule (if possible).

        :param imol_ref:  the reference model molecule index 

        :param chain_id_ref:  the chain ID for the reference chain 

        :param imol_mov:  the moving model molecule index 

        :param chain_id_mov:  the chain ID for the moving chain """
        pass

    def add_lsq_superpose_match(self, chain_id_ref: str, res_no_ref_start: int, res_no_ref_end: int, chain_id_mov: str, res_no_mov_start: int, res_no_mov_end: int, match_type: int) -> None:
        """ Superpose using LSQ - setup the matches

        :param chain_id_ref:  the chain ID for the reference chain 

        :param res_no_ref_start:  the starting residue number in the reference chain 

        :param res_no_ref_end:  the ending residue number in the reference chain 

        :param chain_id_mov:  the chain ID for the moving chain 

        :param res_no_mov_start:  the starting residue number in the moving chain 

        :param res_no_mov_end:  the ending residue number in the moving chain 

        :param match_type:  0: all, 1: main, 2: CAs, 3: N, CA, C, 4: N, CA, CB, C """

    def add_lsq_superpose_atom_match(self, chain_id_ref: str, res_no_ref: int, atom_name_ref: str, chain_id_mov: str, res_no_mov: int, atom_name_mov: str) -> None:
        """ Superpose using LSQ for a scpecific atom - setup the matches

        :param chain_id_ref:  the chain ID for the reference chain 

        :param res_no_ref:  the residue number in the reference chain 

        :param atom_name_ref:  the name of the reference atom 

        :param chain_id_mov:  the chain ID for the moving chain 

        :param res_no_mov:  the residue number in the moving chain 

        :param atom_name_mov:  the name of the moving atom """

    def clear_lsq_matches(self) -> None:
        """ Clear any existing lsq matchers. """

    def lsq_superpose(self, imol_ref: int, imol_mov: int) -> bool:
        """ Apply the superposition using LSQ

        :param imol_ref:  the reference model molecule index 

        :param imol_mov:  the moving model molecule index 

        :return: the success status, i.e. whether or not there were enough atoms to superpose """
        return True

    def transform_map_using_lsq_matrix(self, imol_map: int, lsq_matrix: lsq_results_t, x: float, y: float, z: float, radius: float) -> int:
        """ Transform a map and create a new map

        :param imol_map:  map molecule index 

        :param lsq_matrix:  is an object of type lsq_results_t, is the object returned by 

        :param x:  is the point in the map about which the map is transformed 

        :param y:  is the point in the map about which the map is transformed 

        :param z:  is the point in the map about which the map is transformed 

        :param radius:  the radius of the transformed map, typically between 10 and 100 A

        :return: the molecule index of the new map, -1 for failure """
        return 0

    def get_lsq_matrix(self, imol_ref: int, imol_mov: int, summary_to_screen: bool):
        """ Get LSQ matrix

        don't apply it to the coordinates

        :param imol_ref:  the reference model molecule index 

        :param imol_mov:  the moving model molecule index 

        :param summary_to_screen:  if True, write a summary of the statistics to the output

        :return: the transformation matrix in a simple class """
        pass

    def get_symmetry(self, imol: int, symmetry_search_radius: float, centre_x: float, centre_y: float, centre_z: float):
        """ Get symmetry

        now comes in a simple container that also includes the cell """
        pass

    def get_cell(self, imol: int):
        """ Get the cell

        Check that `is_set` is true before use.

        :param imol:  is the model molecule index

        :return: a `cell_t` """
        pass

    def get_map_molecule_centre(self, imol: int):
        """ Get the middle of the "molecule blob" in cryo-EM reconstruction maps

        :param imol:  is the map molecule index

        :return: a """
        pass

    def undo(self, imol: int) -> int:
        """ Undo

        :param imol:  is the model molecule index

        :return: 1 on successful undo, return 0 on failure """
        return 0

    def redo(self, imol: int) -> int:
        """ Redo

        :param imol:  is the model molecule index

        :return: 1 on successful redo, return 0 on failure """
        return 0

    def get_torsion(self, imol: int, cid: str, atom_names: list):
        """ Get the torsion of the specified atom in the specified residue

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID, e.g. //A/15 (residue 15 in chain A) 

        :param atom_names:  is a list of atom names, e.g. [" CA ", " CB ", " CG ", " CD1"]

        :return: a pair, the first of which is a succes status (1 success, 0 failure), the second is the torsion in degrees """
        pass

    def set_temperature_factors_using_cid(self, imol: int, cid: str, temp_fact: float) -> None:
        """ Change the B factors

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID, e.g. //A/15 (residue 15 in chain A) 

        :param temp_fact:  is the isotropic ADP/temperature factor, e.g., 22 """


    def get_map_sampling_rate(self) -> float:
        """ Get map sampling rate

        :return: the map sampling rate, the default is 1.8 """
        return 0.0

    def set_map_sampling_rate(self, msr: float) -> None:
        """ Set the map sampling rate

        Higher numbers mean smoother maps, but they take longer to generate, longer to transfer, longer to parse and longer to draw

        :param msr:  is the map sampling rate to set, the default is 1.8 """

    def read_mtz(self, file_name: str, f: str, phi: str, weight: str, use_weight: bool, is_a_difference_map: bool) -> int:
        """ Read the given mtz file

        :param file_name:  is the name of the MTZ file 

        :param f:  F column, "FWT" 

        :param phi:  phi column, "PHWT" 

        :param weight:  weight column, "W" 

        :param use_weight:  flag for weights usage, False 

        :param is_a_difference_map:  False

        :return: the new molecule number or -1 on failure """
        return 0

    def replace_map_by_mtz_from_file(self, imol: int, file_name: str, f: str, phi: str, weight: str, use_weight: bool) -> int:
        """ Replace map by mtz

        :param imol:  is the map molecule index 

        :param f:  F column, "FWT" 

        :param phi:  phi column, "PHWT" 

        :param weight:  weight column, "W" 

        :param use_weight:  flag for weights usage, False """
        return 0

    def auto_read_mtz(self, file_name: str):
        """ Auto read the given mtz file

        :param file_name:  is the name of the MTZ file

        :return: a vector of the maps created from reading the file """
        pass

    def read_ccp4_map(self, file_name: str, is_a_difference_map: bool) -> int:
        """ Read a CCP4 (or MRC) map

        There is currently a size limit of 1000 pixels per edge.

        :param file_name:  is the name of the map file 

        :param is_a_difference_map:  False

        :return: the new molecule number or -1 on failure """
        return 0

    def write_map(self, imol: int, file_name: str) -> int:
        """ Write a map

        :param imol:  is the map molecule index 

        :param file_name:  is the name of the new map file

        :return: 1 on a successful write, return 0 on failure. """
        return 0

    def get_map_mean(self, imol: int) -> float:
        """ Get map mean

        :param imol:  is the map molecule index

        :return: the mean of the map or -1 if imol is not a valid map molecule index """
        return 0.0

    def get_map_rmsd_approx(self, imol_map: int) -> float:
        """ Get map rmsd approx

        the function is approximate because the epsilon factor is not taken into account

        :param imol_map:  is the map molecule index

        :return: the map rmsd. -1 is returned if imol_map is not a valid map molecule index. """
        return 0.0

    def get_map_histogram(self, imol: int, n_bins: int, zoom_factor: float):
        """ Get map histogram

        :param imol:  is the map molecule index 

        :param n_bins:  is the number of bins - 200 is a reasonable default. 

        :param zoom_factor:  (reduces the range by the given factor) centred around the median (typically 1.0 but usefully can vary until ~20.0).

        :return: the map histogram """
        pass

    def get_suggested_initial_contour_level(self, imol: int) -> float:
        """ 

        :param imol:  is the map molecule index

        :return: the suggested initial contour level. Return -1 on not-a-map """
        return 0.0

    def is_EM_map(self, imol: int) -> bool:
        """ Check if a map is an EM map or not

        :param imol:  is the map molecule index

        :return: the "EM" status of this molecule. Return false on not-a-map. """
        return True

    def sharpen_blur_map(self, imol_map: int, b_factor: float, in_place_flag: bool) -> int:
        """ Create a new map that is blurred/sharpened

        :param imol:  is the map molecule index 

        :param b_factor:  e.g. 100.0 

        :param in_place_flag:  True if you want to replace the current map, False if you want to create a new map

        :return: the molecule index of the new map or -1 on failure or if `in_place_flag` was true. """
        return 0

    def sharpen_blur_map_with_resample(self, imol_map: int, b_factor: float, resample_factor: float, in_place_flag: bool) -> int:
        """ Create a new map that is blurred/sharpened and resampling

        Note that resampling can be slow, a resample_factor of 1.5 is about the limit of the trade of of prettiness for speed.

        :param imol_map:  is the map molecule index 

        :param b_factor:  e.g. 100.0 

        :param resample_factor:  e.g. 1.4 

        :param in_place_flag:  True if you want to replace the current map, False if you want to create a new map

        :return: the molecule index of the new map or -1 on failure or if `in_place_flag` was true. """
        return 0

    def mask_map_by_atom_selection(self, imol_coords: int, imol_map: int, cid: str, atom_radius: float, invert_flag: bool) -> int:
        """ Mask map by atom selection

        (note the argument order is reversed compared to the coot api).

        :param imol_coords:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param atom_radius:  is the atom radius. Use a negative number to mean "default". 

        :param invert_flag:  changes the parts of the map that are masked, so to highlight the density for a ligand one would pass the cid for the ligand and invert_flag as true, so that the parts of the map that are not the ligand are suppressed.

        :return: the index of the new map or -1 on failure """
        return 0

    def partition_map_by_chain(self, imol_map: int, imol_model: int):
        """ Partition the input map

        Each voxel in the map is assigned to the chain to which it is nearest. Unlike masking, the generated maps are not restricted to be "close" to the atoms in the atom selection.

        e.g. maskChains for ChimeraX - JiangLab

        :param imol_map:  is the map molecule index 

        :param imol_model:  is the model molecule index

        :return: a vector/list of the molecules indices of the newly created maps """
        pass

    def make_mask(self, imol_map_ref: int, imol_model: int, atom_selection_cid: str, radius: float) -> int:
        """ Make a masked map

        :param imol_map_ref:  is the map molecule index 

        :param imol_model:  is the model molecule index 

        :param atom_selection_cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param radius:  is the radius of the map, e.g. 12.0 for X-ray map and 100.0 for Cryo-EM map

        :return: the index of the newly created mask. Return -1 on failure. """
        return 0

    def flip_hand(self, imol_map: int) -> int:
        """ Generate a new map which is the hand-flipped version of the input map

        :param imol_map:  is the map molecule index

        :return: the molecule index of the new map, or -1 on failure. """
        return 0

    def make_masked_maps_split_by_chain(self, imol: int, imol_map: int):
        """ Make a vector/list of maps that are split by chain-id of the input imol

        :param imol:  is the model molecule index 

        :param imol_map:  is the map molecule index

        :return: a vector/list of the map molecule indices. """
        pass

    def dedust_map(self, imol: int) -> int:
        """ dedust map

        :param imol_map:  the map molecule index

        :return: the map molecule index of the dedusted map or -1 on failure """
        return 0

    def set_map_colour(self, imol: int, r: float, g: float, b: float) -> None:
        """ Set the map colour

        The next time a map mesh is requested, it will have this colour. This does not apply to/affect the colour of the difference maps.

        RGB colour codes, e.g. green is r:0, g: 255, b:0

        :param imol:  is the model molecule index """

    def set_map_is_contoured_with_thread_pool(self, state: bool) -> None:
        """ Set the state of the mode of the threading in map contouring

        :param state:  is True to mean that it is enabled. The default is True. """

    def get_map_contours_mesh(self, imol: int, position_x: float, position_y: float, position_z: float, radius: float, contour_level: float):
        """ Get the mesh for the map contours

        This function is not const because the internal state of a `coot_molecule_t` is changed.

        :param imol:  is the model molecule index 

        :param position_x:  is the x coordinate of the target position 

        :param position_y:  is the y coordinate of the target position 

        :param position_z:  is the z coordinate of the target position 

        :param radius:  is the radius of the map, e.g. 12.0 for X-ray map and 100.0 for Cryo-EM map 

        :param contour_level:  e.g. 1.5 sd for X-ray, 5.0 sd for cryo-EM

        :return: a `simple_mesh_t` for the map contours of the specified map """
        pass

    def get_map_contours_mesh_using_other_map_for_colours(self, imol_ref: int, imol_map_for_colouring: int, position_x: float, position_y: float, position_z: float, radius: float, contour_level: float, other_map_for_colouring_min_value: float, other_map_for_colouring_max_value: float, invert_colour_ramp: bool):
        """ Get the mesh for the map contours using another map for colouring

        :param imol_ref:  is the reference map index 

        :param imol_map_for_coloring:  is the target map index 

        :param position_x:  is the x coordinate of the target position 

        :param position_y:  is the y coordinate of the target position 

        :param position_z:  is the z coordinate of the target position 

        :param radius:  is the radius of the map, e.g. 12.0 for X-ray map and 100.0 for Cryo-EM map 

        :param contour_level:  e.g. 1.5 sd for X-ray, 5.0 sd for cryo-EM 

        :param other_map_for_colouring_min_value:  e.g. -1.0 in the case of correlation map 

        :param other_map_for_colouring_max_value:  e.g. 1.0 in the case of correlation map 

        :param invert_colour_ramp:  e.g. red to blue rather than blue to red

        :return: a `simple_mesh_t` for the map contours of the specified map """
        pass

    def set_map_colour_saturation(self, imol: int, s: float) -> None:
        """ Set the map saturation

        :param s:  is the map saturation, e.g. a number between 0 and 1, where 0 is grey and 1 is "lego-like" colour scheme. 0.5 is a nice middle value """


    def get_map_vertices_histogram(self, imol: int, imol_map_for_sampling: int, position_x: float, position_y: float, position_z: float, radius: float, contour_level: float, n_bins: int):
        """ Get map vertices histogram

        Note not const because 

        :param imol:  is the map molecule index 

        :param n_bins:  is the number of bins - 40 is a reasonable default.

        :return: the map vertices histogram """
        pass

    def get_latest_sfcalc_stats(self):
        """ Get the latest sfcalc stats

        :return: a sfcalc_genmap_stats_t object """
        pass

    def get_r_factor_stats(self):
        """ Get the R-factors

        :return: an object """
        pass

    def r_factor_stats_as_string(self, rfs: str) -> str:
        """ Get the R factor stats as a string

        :param rfs:  is the name of the string

        :return: a string with the R-factor stats """
        return 'a-string'

    def average_map(self, imol_maps: str, scales: list) -> int:
        """ Get the average map

        This function does no normalization of the scales, presuming that they are pre-normalized.

        :param imol_maps:  is a colon-separated list of map indices e.g. "2:3:4" 

        :param scales:  is the list of weights corresponding to the list of maps. The number of scales factors should match the number of maps

        :return: the index of the new map, or -1 on failure. """
        return 0

    def regen_map(self, imol_map: int, imol_maps: str, scales: list) -> bool:
        """ This function does no normalisation of the scales, presuming that they are pre-normalized.

        The number of maps in imol_maps should match the size of the scale vector. If not, nothing will happen and False will be returned

        :param imol_map:  is the map molecule index 

        :param imol_maps:  is a colon-separated list of map indices e.g. "2:3:4" 

        :param scales:  is the list of weights corresponding to the list of maps

        :return: the success status """
        return True



    def testing_start_long_term_job(self, n_seconds: int) -> None:
        """ Testing function

        start a long-term job.

        :param n_seconds:  is the number of seconds, if is 0, then run forever (or until interrupted) """

    def testing_stop_long_term_job(self) -> None:
        """ Testing function

        stop the long-term job runnning """

    def testing_interrogate_long_term_job(self):
        """ Testing function

        get the stats for the long-term job """
        pass

    def get_contouring_time(self):
        """ Testing function

        get the time for contouring in milliseconds """
        pass

    def set_max_number_of_threads(self, n_threads: int) -> None:
        """ Testing function

        set the maximum number of threads for both the thread pool and the vector of threads

        :param n_threads:  is the number of threads """

    def set_max_number_of_threads_in_thread_pool(self, n_threads: int) -> None:
        """ Testing function

        Deprecated name for the "set_max_number_of_threads()" function """

    def test_the_threading(self, n_threads: int):
        """ Testing function

        get the time to run a test function in milliseconds

        :param n_threads:  is the number of threads """
        pass

    def test_launching_threads(self, n_threads_per_batch: int, n_batches: int):
        """ Testing function

        :param n_threads_per_batch:  is the number of threads per batch 

        :param n_batches:  is the number batches

        :return: the time per batch in microseconds """
        pass

    def test_thread_pool_threads(self, n_threads: int):
        """ Testing function

        :param n_threads:  is the number of threads

        :return: time in microseconds """
        pass

    def mmcif_tests(self, last_test_only: bool) -> int:
        """ Testing function

        a test for mmdb/gemmi/mmcif functionality

        :param last_test_only:  is True to mean that only that last test should be run. The default is False. This is useful to set to True while a test is being developed.

        :return: the success status: 1 means that all the tests passed. """
        return 0

    def get_molecule_name(self, imol: int) -> str:
        """ Get the molecule name

        :param imol:  is the model molecule index

        :return: the name of the molecule """
        return 'a-string'

    def set_molecule_name(self, imol: int, new_name: str) -> None:
        """ Set the molecule name

        :param imol:  is the model or map molecule index 

        :param new_name:  is the new name of the model or map """

    def display_molecule_names_table(self) -> None:
        """ Debugging function: display the table of molecule and names. """

    def is_valid_model_molecule(self, imol: int) -> bool:
        """ Check if the model index is valid

        e.g. if the molecule is a map you will have an invalid model

        :param imol:  is the model molecule index

        :return: True or False """
        return True

    def is_valid_map_molecule(self, imol_map: int) -> bool:
        """ Check if the map index is valid

        e.g. if the map is a model you will have an invalid map

        :param imol_map:  is the map molecule index

        :return: True or False """
        return True

    def is_a_difference_map(self, imol_map: int) -> bool:
        """ Check if it the map is a difference map

        :param imol_map:  is the map molecule index

        :return: True or False """
        return True

    def new_molecule(self, name: str) -> int:
        """ Create an empty molecule

        :return: the index of the new molecule """
        return 0

    def close_molecule(self, imol: int) -> int:
        """ Close the molecule (and delete dynamically allocated memory)

        :param imol:  is the model molecule index

        :return: 1 on successful closure and 0 on failure to close """
        return 0

    def end_delete_closed_molecules(self) -> None:
        """ Delete the most recent/last closed molecule in the molecule vector, until the first non-closed molecule is found (working from the end) """

    def pop_back(self) -> None:
        """ Delete the most recent/last molecule in the molecule vector. """

    def clear(self) -> None:
        """ Delete all molecules. """

    def get_eigenvalues(self, imol: int, chain_id: str, res_no: int, ins_code: str):
        """ Get the eigenvalues of the specified residue

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A"

        :return: the eigenvalues of the atoms in the specified residue """
        pass

    def test_origin_cube(self):
        """ Get a simple test mesh

        :return: the mesh of a unit solid cube at the origin """
        pass

    def fill_rotamer_probability_tables(self) -> None:
        """ Fill the rotamer probability tables (currently not ARG and LYS) """

    def accept_rotamer_probability_tables_compressed_data(self, data_stream: str) -> None:
        """ Access to a compressed file that contains the rotamer probabilities

        libcootapi will fill the rotamer probabilities tables from this compressed data stream (placeholder only) """

    def contains_unsaved_models(self) -> bool:
        """ Check if there are unsaved changes for this model

        e.g. as yet not written to disk

        :return: a flag of unsaved models state - e.g. if any of them are unsaved, then this returns True. """
        return True

    def save_unsaved_model_changes(self) -> None:
        """ Save the unsaved model - this function has not yet been written! """

    def geometry_init_standard(self) -> None:
        """ Read the standard list of residues. """

    def non_standard_residue_types_in_model(self, imol: int):
        """ Get a list of non-standard residues in the given molecule

        :param imol:  is the model molecule index

        :return: a vector/list of non-standard residues """
        pass

    def auto_fit_rotamer(self, imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str, imol_map: int) -> int:
        """ Auto-fit rotamer

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A" 

        :param alt_conf:  is the alternate conformation, e.g. "A" or "B" 

        :param imol_map:  is the map molecule index

        :return: 1 on successful modification, return 0 on failure """
        return 0

    def change_to_next_rotamer(self, imol: int, residue_cid: str, alt_conf: str):
        """ Change to the next rotamer (rotamer cycling is implicit if needed)

        :param imol:  is the model molecule index 

        :param residue_cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param alt_conf:  is the alternate conformation, e.g. "A" or "B"

        :return: the change information. """
        pass

    def change_to_previous_rotamer(self, imol: int, residue_cid: str, alt_conf: str):
        """ Change to the previous rotamer (rotamer cycling is implicit if needed)

        :param imol:  is the model molecule index 

        :param residue_cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param alt_conf:  is the alternate conformation, e.g. "A" or "B"

        :return: the change information. """
        pass

    def change_to_first_rotamer(self, imol: int, residue_cid: str, alt_conf: str):
        """ Change to the first (0th) rotamer

        :param imol:  is the model molecule index 

        :param residue_cid:  is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param alt_conf:  is the alternate conformation, e.g. "A" or "B"

        :return: the change information. """
        pass

    def delete_using_cid(self, imol: int, cid: str, scope: str) -> tuple:
        """ Delete item

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (all the atoms in residue 15 of chain A) 

        :param scope:  is one of the strings: ["ATOM", "WATER", "RESIDUE"," CHAIN"," MOLECULE", "LITERAL"]

        :return: 1 on successful deletion, return 0 on failure """

    def delete_atom(self, imol: int, chain_id: str, res_no: int, ins_code: str, atom_name: str, alt_conf: str) -> tuple:
        """ Delete atom

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A" 

        :param atom_name:  is the name of the atom, e.g. "OH" 

        :param alt_conf:  is the alternate conformation, e.g. "A" or "B"

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_atom_using_cid(self, imol: int, cid: str) -> tuple:
        """ Delete atom using cid

        :param imol:  is the model molecule index 

        :param cid:  is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_residue(self, imol: int, chain_id: str, res_no: int, ins_code: str) -> tuple:
        """ Delete residue

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A"

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_residue_using_cid(self, imol: int, cid: str) -> tuple:
        """ Delete residue using cid

        :param imol:  is the model molecule index 

        :param cid:  is the residue selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_residue_atoms_with_alt_conf(self, imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str) -> tuple:
        """ Delete residue atoms using alt_conf

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A" 

        :param alt_conf:  is the alternate conformation, e.g. "A" or "B"

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_residue_atoms_using_cid(self, imol: int, cid: str) -> tuple:
        """ Delete residue atoms using cid

        This is the same as `delete_atom_using_cid`. It will be deleted in the future

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_side_chain(self, imol: int, chain_id: str, res_no: int, ins_code: str) -> tuple:
        """ Delete side chain

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A"

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_side_chain_using_cid(self, imol: int, cid: str) -> tuple:
        """ Delete side chain using cid

        :param imol:  is the model molecule index 

        :param cid:  is the residue selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_chain_using_cid(self, imol: int, cid: str) -> tuple:
        """ Delete chain using chain cid

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A" (chain A), "//*" (all chains)

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_literal_using_cid(self, imol: int, cid: str) -> tuple:
        """ Delete the atoms specified in the cid selection

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15/OH" (atom OH in residue 15)

        :return: 1 on successful deletion, return 0 on failure to delete. """

    def delete_all_carbohydrate(self, imol: int) -> bool:
        """ delete all carbohydrate

        :param imol:  is the model molecule index

        :return: true on successful deletion, return false on no deletion. """
        return True

    def change_alt_locs(self, imol: int, cid: str, change_mode: str) -> int:
        """ Change alternate conformation

        Note that this function only deals with (swaps) alt confs "A" and "B" - any alt-conf other than that is ignored.

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 in chain A) 

        :param change_mode:  is either "residue", "main-chain", "side-chain" or a comma-separated atom-name pairs (e.g "N,CA") - you can (of course) specify just one atom, e.g.: "N".

        :return: the success status (1 is done, 0 means failed to do) """
        return 0

    def add_terminal_residue_directly(self, imol: int, chain_id: str, res_no: int, ins_code: str):
        """ Add a residue onto the end of the chain by fitting to density

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A"

        :return: first: 1 on success, second is failure message """
        pass

    def add_terminal_residue_directly_using_cid(self, imol: int, cid: str) -> int:
        """ Add a residue onto the end of the chain by fitting to density using cid

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15/OH" (atom OH in residue 15) 

        :return: success status (1 for good, 0 for not done) """
        return 0

    def add_terminal_residue_directly_using_bucca_ml_growing_using_cid(self, imol: int, cid: str) -> int:
        """ Add a residue onto the end of the chain by fitting to density using Buccaneer building and cid

        This function has been removed - is is now a noop.

        :param imol:  is the model molecule index 

        :param cid:  is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15) """
        return 0

    def add_terminal_residue_directly_using_bucca_ml_growing(self, imol: int, spec: str) -> int:
        """ Add a residue onto the end of the chain by fitting to density using Buccaneer building

        :param imol:  is the model molecule index 

        :param spec:  is the residue specifier, residue_spec_t("A", 10, "") """
        return 0

    def set_add_waters_water_to_protein_distance_lim_min(self, d: float) -> None:
        """ Parameter for `add_waters` function

        :param d:  is the min distance, default 2.4 """

    def set_add_waters_water_to_protein_distance_lim_max(self, d: float) -> None:
        """ Parameter for `add_waters` function

        :param d:  is the max distance, default 3.4 """

    def set_add_waters_variance_limit(self, d: float) -> None:
        """ Parameter for `add_waters` function

        :param d:  is the variance limit, default is 0.1 """

    def set_add_waters_sigma_cutoff(self, d: float) -> None:
        """ Parameter for `add_waters` function

        :param d:  is the sigma cutoff, default is 1.75 """

    def add_waters(self, imol_model: int, imol_map: int) -> int:
        """ Add waters

        :param imol:  is the model molecule index 

        :param imol_map:  is the map molecule index

        :return: the number of waters added on a success, -1 on failure. """
        return 0

    def flood(self, imol_model: int, imol_map: int, n_rmsd: float) -> int:
        """ Flood with dummy atoms

        :param imol:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param n_rmsd:  e.g., 4.0

        :return: the number of waters added on a success, -1 on failure. """
        return 0

    def add_hydrogen_atoms(self, imol_model: int) -> int:
        """ Add hydrogen atoms

        :param imol_model:  is the model molecule index

        :return: 1 on success, 0 on failure. """
        return 0

    def delete_hydrogen_atoms(self, imol_model: int) -> int:
        """ Delete hydrogen atoms

        :param imol_model:  is the model molecule index

        :return: 1 on a successful deletion, 0 on failure. """
        return 0

    def add_alternative_conformation(self, imol_model: int, cid: str) -> int:
        """ Add an alternative conformation for the specified residue

        :param imol_model:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 in chain A)

        :return: 1 on a successful addition, 0 on failure. """
        return 0

    def fill_partial_residue(self, imol: int, chain_id: str, res_no: int, ins_code: str) -> int:
        """ Fill the specified residue

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A"

        :return: 1 on a successful fill, 0 on failure. """
        return 0

    def fill_partial_residue_using_cid(self, imol: int, cid: str) -> int:
        """ Fill the specified residue using cid

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 in chain A)

        :return: 1 on a successful fill, 0 on failure. """
        return 0

    def fill_partial_residues(self, imol: int) -> int:
        """ Fill all the the partially-filled residues in the molecule

        :param imol:  is the model molecule index

        :return: 1 on a successful fill, 0 on failure. """
        return 0

    def add_named_glyco_tree(self, imol_model: int, imol_map: int, glycosylation_name: str, asn_chain_id: str, asn_res_no: int) -> None:
        """ Add N-linked glycosylation

        :param imol_model:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param glycosylation_name:  is the type of glycosylation, one of: "NAG-NAG-BMA" or "high-mannose" or "hybrid" or "mammalian-biantennary" or "plant-biantennary" 

        :param asn_chain_id:  is the chain-id of the ASN to which the carbohydrate is to be added 

        :param asn_res_no:  is the residue number of the ASN to which the carbohydrate is to be added """

    def flip_peptide_using_cid(self, imol: int, atom_cid: str, alt_conf: str) -> int:
        """ Flip peptide using cid

        :param imol:  is the model molecule index 

        :param atom_cid:  is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A) 

        :param alt_conf:  is the alternate conformation, e.g. "A" or "B"

        :return: 1 on a successful flip """
        return 0

    def eigen_flip_ligand(self, imol: int, chain_id: str, res_no: int, ins_code: str) -> None:
        """ Eigen-flip the specified ligand

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A" """

    def eigen_flip_ligand_using_cid(self, imol: int, residue_cid: str) -> None:
        """ Eigen-flip ligand using cid

        :param imol:  is the model molecule index 

        :param residue_cid:  is the residue selection CID e.g "//A/15" (residue 15 of chain A) """

    def mutate(self, imol: int, cid: str, new_residue_type: str) -> int:
        """ Mutate residue

        :param imol:  is the model molecule index 

        :param residue_cid:  is the residue selection CID e.g "//A/15" (residue 15 of chain A) 

        :param new_residue_type:  is the 3-letter code of the new residue, e.g. "TYR" for tyrosine

        :return: 1 on a successful move, 0 on failure. """
        return 0

    def side_chain_180(self, imol: int, atom_cid: str) -> int:
        """ Rotate last chi angle of the side chain by 180 degrees

        :param imol:  is the model molecule index 

        :param atom_cid:  is the atom selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A)

        :return: 1 on a successful move, 0 on failure. """
        return 0

    def jed_flip(self, imol: int, atom_cid: str, invert_selection: bool) -> str:
        """ JED-Flip the ligand (or residue) at the specified atom

        :param imol:  is the model molecule index 

        :param atom_cid:  is the residue selection CID e.g "//A/15" (residue 15 of chain A) 

        :param invert_selection:  is True if you want to use the larger fragment

        :return: a non-blank message if there is a problem """
        return 'a-string'

    def move_molecule_to_new_centre(self, imol: int, x: float, y: float, z: float) -> int:
        """ Move the molecule to the given centre

        :param imol:  is the model molecule index 

        :param x:  is the x coordinate of the new centre of the screen 

        :param y:  is the y coordinate of the new centre of the screen 

        :param z:  is the z coordinate of the new centre of the screen

        :return: 1 on a successful move, 0 on failure. """
        return 0

    def multiply_residue_temperature_factors(self, imol: int, cid: str, factor: float) -> None:
        """ Interactive B-factor refinement

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A) 

        :param factor:  might typically be 0.9 or 1.1 """

    def get_molecule_centre(self, imol: int):
        """ Get molecule centre

        :param imol:  is the model molecule index

        :return: the molecule centre """
        pass

    def get_radius_of_gyration(self, imol: int):
        """ Get Radius of Gyration

        :param imol:  is the model molecule index

        :return: the molecule centre. If the number is less than zero, there was a problem finding the molecule or atoms. """
        pass

    def copy_molecule(self, imol: int) -> int:
        """ Copy the molecule

        :param imol:  the specified molecule

        :return: the new molecule number """
        return 0

    def copy_fragment_using_cid(self, imol: int, multi_cid: str) -> int:
        """ Copy a fragment given the multi_cid selection string

        :param imol:  is the model molecule index 

        :param multi_cids:  is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||/B/56-66"

        :return: the new molecule number (or -1 on no atoms selected) """
        return 0

    def copy_fragment_for_refinement_using_cid(self, imol: int, multi_cid: str) -> int:
        """ Copy a fragment given the multi_cid selection string for refinement

        Use this in preference to `copy_fragment_using_cid` when copying a molecule fragment to make a molten zone for refinement. That is because this version quietly also copies the residues near the residues of the selection, so that those residues can be used for links and non-bonded contact restraints.

        :param imol:  is the model molecule index 

        :param multi_cids:  is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||//B/56-66"

        :return: the new molecule number (or -1 on no atoms selected) """
        return 0

    def copy_fragment_using_residue_range(self, imol: int, chain_id: str, res_no_start: int, res_no_end: int) -> int:
        """ Copy a residue-range fragment

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" 

        :param res_no_start:  the starting residue number 

        :param res_no_ref_end:  the ending residue number

        :return: the new molecule number (or -1 on no atoms selected) """
        return 0

    def apply_transformation_to_atom_selection(self, imol: int, atoms_selection_cid: str, n_atoms: int, m00: float, m01: float, m02: float, m10: float, m11: float, m12: float, m20: float, m21: float, m22: float, c0: float, c1: float, c2: float, t0: float, t1: float, t2: float) -> int:
        """ Apply transformation to atom selection in the given molecule

        :return: the number of atoms moved. """
        return 0

    def new_positions_for_residue_atoms(self, imol: int, residue_cid: str, moved_atoms: list) -> int:
        """ Update the positions of the atoms in the residue

        :param imol:  is the model molecule index 

        :param residue_cid:  is the residue selection CID e.g "//A/15" (residue 15 of chain A) 

        :param moved_atoms:  is a list of the atoms moved in the specified residue, e.g. moved_atom_t(" CA ", 1, 2, 3) """
        return 0

    def new_positions_for_atoms_in_residues(self, imol: int, moved_residues: list) -> int:
        """ Update the positions of the atoms in the residues

        :param imol:  is the model molecule index 

        :param moved_residue:  is a list of the residues with moved atoms, e.g. moved_residue_t("A", 10, "") """
        return 0

    def merge_molecules(self, imol: int, list_of_other_molecules: str):
        """ Merge molecules

        :param imol:  is the model molecule index 

        :param list_of_other_molecules:  is a colon-separated list of molecules, e.g. "2:3:4"

        :return: the first is a flag set to 1 if a merge occurred (and 0 if it did not) the second is a vector of merge results, i.e. if you merged a ligand, what is the new residue spec of the ligand, and if you merged a (polymer) chain, what is the new chain-id of that chain. """
        pass

    def cis_trans_convert(self, imol: int, atom_cid: str) -> int:
        """ Convert a cis peptide to a trans or vice versa

        :param imol:  is the model molecule index 

        :param atom_cid:  is the atom selection CID e.g "//A/15/OH" (atom OH residue 15 of chain A)

        :return: 1 on a successful conversion. """
        return 0

    def replace_residue(self, imol: int, residue_cid: str, new_residue_type: str, imol_enc: int) -> None:
        """ Replace a residue

        This has a different meaning of "replace" to 

        Change the type of a residue (for example, "TYR" to "CYS") The algorithm will superpose the mainchain CA, C and N and try to set matching torsion to the angle that they were in the reference structure.

        :param imol:  is the model molecule index 

        :param residue_cid:  is the residue selection CID e.g "//A/15" (residue 15 of chain A) 

        :param new_residue_type:  is the 3-letter code of the new residue, e.g "CYS" 

        :param imol_enc:  is the molecule index for the residue type/compound_id """

    def replace_fragment(self, imol_base: int, imol_reference: int, atom_selection: str) -> int:
        """ Replace a fragment

        :param imol_base:  is the base model index 

        :param imol_reference:  is the reference model index 

        :param atom_selection:  is the selection CID e.g "//A/15-17" (residue 15, 16 and 17 of chain A)

        :return: the success status """
        return 0

    def rigid_body_fit(self, imol: int, multi_cid: str, imol_map: int) -> int:
        """ Rigid-body fitting

        :param imol:  is the model molecule index 

        :param multi_cids:  is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||/B/56-66" 

        :param imol_map:  is the map molecule index """
        return 0

    def rotate_around_bond(self, imol: int, residue_cid: str, atom_name_1: str, atom_name_2: str, atom_name_3: str, atom_name_4: str, torsion_angle: float) -> int:
        """ Rotate atoms around torsion

        the bond is presumed to be between atom-2 and atom-3. Atom-1 and atom-4 are used to define the absolute torsion angle.

        :param imol:  is the model molecule index 

        :param residue_cid:  is the residue selection CID e.g "//A/15" (residue 15 of chain A) 

        :param atom_name_1:  e.g. " CA " 

        :param atom_name_2:  e.g. " CB " 

        :param atom_name_3:  e.g. " CG " 

        :param atom_name_4:  e.g. " CD1" 

        :param torsion_angle:  e.g. 12.3 degrees

        :return: status 1 if successful, 0 if not. """
        return 0

    def change_chain_id(self, imol: int, from_chain_id: str, to_chain_id: str, use_resno_range: bool, start_resno: int, end_resno: int):
        """ Change the chain ID

        :param imol:  is the model molecule index 

        :param from_chain_id:  e.g. "A" 

        :param to_chain_id:  e.g. "C" 

        :param use_resno_range:  use residue number range, typically True 

        :param start_resno:  the starting residue number of the range 

        :param end_resno:  the ending residue number of the range

        :return: -1 on a conflict, 1 on good, 0 on did nothing, return also an information/error message """
        pass

    def split_residue_using_map(self, imol: int, residue_cid: str, imol_diff_map: int) -> int:
        """ Split a residue into alt-confs

        do nothing if the residue already has alt-confs.

        :param imol:  the modified model 

        :param residue_cid:  the modified residue, the residue selection CID e.g "//A/15" (residue 15 of chain A) 

        :param imol_diff_map:  is the difference map that is used to determine the residue split

        :return: split success status """
        return 0

    def associate_sequence(self, imol: int, name_or_chain_id: str, sequence: str) -> None:
        """ Associate a sequence with a molecule

        :param imol:  is the model molecule index 

        :param name_or_chain_id:  e.g. "A" 

        :param sequence:  is the model sequence """

    def assign_sequence(self, imol_model: int, imol_map: int) -> None:
        """ Assign a sequence to a molecule

        Often one might copy out a fragment from a more complete molecule (and then copy it back after the sequence has been added). This runs `backrub_rotamer()` on the newly assigned residues

        :param imol:  is the model molecule index 

        :param imol_map:  is the map molecule index """

    def get_sequence_info(self, imol: int):
        """ Get the sequence information

        :param imol:  is the molecule index 

        :return: the sequence information """
        pass

    def get_mutation_info(self, imol: int):
        """ Get mutation information

        The reference sequece is that which has been provided using the 

        :param imol:  is the model molecule index 

        :return: the mismatches/mutations as insertions, deletions or mutations """
        pass

    def refine_residues_using_atom_cid(self, imol: int, cid: str, mode: str, n_cycles: int) -> int:
        """ Refine the residues using cid

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A) 

        :param mode:  is the mode of real space refinement e.g. "SINGLE", "TRIPLE", "QUINTUPLE", "HEPTUPLE", "SPHERE", "BIG_SPHERE", "CHAIN", "ALL" 

        :param n_cycles:  is the number of refinement cycles

        :return: a value of 1 if the refinement was performed and 0 if it was not. """
        return 0

    def refine_residues(self, imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str, mode: str, n_cycles: int) -> int:
        """ Refine the residues

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no:  is the residue number, e.g. 12 

        :param ins_code:  is the insertion code, e.g. "A" 

        :param alt_conf:  is the alternate conformation, e.g. "A" or "B" 

        :param mode:  is the mode of real space refinement e.g. "SINGLE", "TRIPLE", "QUINTUPLE", "HEPTUPLE", "SPHERE", "BIG_SPHERE", "CHAIN", "ALL" 

        :param n_cycles:  is the number of refinement cycles

        :return: a value of 1 if the refinement was performed and 0 if it was not. """
        return 0

    def refine_residue_range(self, imol: int, chain_id: str, res_no_start: int, res_no_end: int, n_cycles: int) -> int:
        """ Refine residue range

        :param imol:  is the model molecule index 

        :param chain_id:  e.g. "A" for chain A 

        :param res_no_start:  the starting residue number 

        :param res_no_ref_end:  the ending residue number 

        :param n_cycles:  is the number of refinement cycles

        :return: a value of 1 if the refinement was performed and 0 if it was not. """
        return 0

    def minimize_energy(self, imol: int, atom_selection_cid: str, n_cycles: int, do_rama_plot_restraints: bool, rama_plot_weight: float, do_torsion_restraints: bool, torsion_weight: float, refinement_is_quiet: bool):
        """ Minimise/optimise the geometry of the specified residue(s)

        The use of "energy" should not be taken literally here

        :param imol:  is the model molecule index 

        :param atom_selection_cid:  is the selection CID e.g. "//A/15" (residue 15 of chain A) 

        :param n_cycles:  is the number of refinement cycles. If you pass n_cycles = 100 (or some such) then you can get the mesh for the partially optimized ligand/residues 

        :param do_rama_plot_restraints:  is the flag for the usage of Ramachandran plot restraints 

        :param rama_plot_weight:  is the flag to set the Ramachandran plot restraints weight 

        :param do_torsion_restraints:  is the flag for the usage of torsion restraints 

        :param torsion_weight:  is the flag to set the torsion restraints weight 

        :param refinement_is_quiet:  is used to reduce the amount of diagnostic text written to the output

        :return: the success status 1 if the minimization was performed and 0 if it was not. """
        pass

    def minimize(self, imol: int, atom_selection_cid: str, n_cycles: int, do_rama_plot_restraints: bool, rama_plot_weight: float, do_torsion_restraints: bool, torsion_weight: float, refinement_is_quiet: bool) -> float:
        """ 

        :param imol:  is the model molecule index 

        :param atom_selection_cid:  is the selection CID e.g. "//A/15" (residue 15 of chain A) 

        :param n_cycles:  is the number of refinement cycles. If you pass n_cycles = 100 (or some such) then you can get the mesh for the partially optimized ligand/residues 

        :param do_rama_plot_restraints:  is the flag for the usage of Ramachandran plot restraints 

        :param rama_plot_weight:  is the flag to set the Ramachandran plot restraints weight 

        :param do_torsion_restraints:  is the flag for the usage of torsion restraints 

        :param torsion_weight:  is the flag to set the torsion restraints weight 

        :param refinement_is_quiet:  is used to reduce the amount of diagnostic text written to the output

        :return: the function value at termination """
        return 0.0

    def fix_atom_selection_during_refinement(self, imol: int, atom_selection_cid: str) -> None:
        """ Fix atoms during refinement

        Does nothing at the moment

        :param imol:  is the model molecule index 

        :param atom_selection_cid:  is the selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A) """

    def add_target_position_restraint(self, imol: int, atom_cid: str, pos_x: float, pos_y: float, pos_z: float) -> None:
        """ Add or update restraint (if it has a pull restraint already)

        :param imol:  is the model molecule index 

        :param atom_cid:  is the selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A) 

        :param pos_x:  is the x coordinate of the target position of the specified atom 

        :param pos_y:  is the y coordinate of the target position of the specified atom 

        :param pos_z:  is the z coordinate of the target position of the specified atom """

    def clear_target_position_restraint(self, imol: int, atom_cid: str) -> None:
        """ Clear target_position restraint

        :param imol:  is the model molecule index 

        :param atom_cid:  is the selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A) """

    def turn_off_when_close_target_position_restraint(self, imol: int) -> None:
        """ Clear target_position restraint if it is (or they are) close to their target position

        :param imol:  is the model molecule index """

    def set_logging_level(self, level: str) -> None:
        """ Control the logging

        :param level:  is the logging level, level is either "LOW" or "HIGH" or "DEBUGGING" """

    def set_logging_file(self, file_name: str) -> None:
        """ make the logging output go to a file

        :param file_name:  the looging file name """

    def set_use_rama_plot_restraints(self, state: bool) -> None:
        """ Turn on or off ramachandran restraints

        :param state:  is True to mean that it is enabled """

    def get_use_rama_plot_restraints(self) -> bool:
        """ Get the state of the rama plot restraints usage in refinement

        :return: the state """
        return True

    def set_rama_plot_restraints_weight(self, f: float) -> None:
        """ Set the Ramachandran plot restraints weight

        :param f:  is the weight to set, default 1.0 """

    def get_rama_plot_restraints_weight(self) -> float:
        """ Get the Ramachandran plot restraints weight

        :return: the Ramachandran plot restraints weight """
        return 0.0

    def set_use_torsion_restraints(self, state: bool) -> None:
        """ Turn on or off torsion restraints

        :param state:  is True to mean that it is enabled """

    def get_use_torsion_restraints(self) -> bool:
        """ Get the state of the rama plot restraints usage in refinement

        :return: the state """
        return True

    def set_torsion_restraints_weight(self, f: float) -> None:
        """ Set the torsiont restraints weight

        :param f:  is the weight to set, default value is 1.0 """

    def get_torsion_restraints_weight(self) -> float:
        """ Get the torsion restraints weight

        :return: the torsion restraints weight """
        return 0.0

    def init_refinement_of_molecule_as_fragment_based_on_reference(self, imol_frag: int, imol_ref: int, imol_map: int) -> None:
        """ Initialise the refinement of (all of) molecule imol_frag

        :param imol_frag:  is the model molecule index of the fragment 

        :param imol_ref:  is the model molecule index of the reference 

        :param imol_map:  is the map molecule index """

    def refine(self, imol: int, n_cycles: int):
        """ Run some cycles of refinement and return a mesh

        That way we can see the molecule animate as it refines

        :param imol:  is the model molecule index 

        :param n_cycles:  is the number of refinement cycles

        :return: a pair: the first of which is the status of the refinement: GSL_CONTINUE, GSL_SUCCESS, GSL_ENOPROG (no progress). i.e. don't call this function again unless the status is GSL_CONTINUE (-2); The second is a """
        pass

    def add_target_position_restraint_and_refine(self, imol: int, atom_cid: str, pos_x: float, pos_y: float, pos_z: float, n_cycles: int):
        """ Create a new position for the given atom and create a new bonds mesh based on that

        This is currently "heavyweight" as the bonds mesh is calculated from scratch (it is not (yet) merely a distortion of an internally-stored mesh).

        :param imol:  is the model molecule index 

        :param atom_cid:  is the selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A) 

        :param pos_x:  is the x coordinate of the target position of the specified atom 

        :param pos_y:  is the y coordinate of the target position of the specified atom 

        :param pos_z:  is the z coordinate of the target position of the specified atom 

        :param n_cycles:  specifies the number of refinement cyles to run after the target position of the atom has been applied. If n_cycles is -1 then, no cycles are done and the mesh is bonds merely calculated.

        :return: a `instanced_mesh_t` """
        pass

    def clear_target_position_restraints(self, imol: int) -> None:
        """ Clear any and all drag-atom target position restraints

        :param imol:  is the model molecule index """

    def clear_refinement(self, imol: int) -> None:
        """ Call this after molecule refinement has finished (say when the molecule molecule is accepted into the original molecule)

        :param imol:  is the model molecule index """

    def set_refinement_is_verbose(self, state: bool) -> None:
        """ For debugging the refinement - write out some diagnositics - some might be useful.

        API change 20240226 - this function now takes a boolean argument

        :param state:  is True to mean that it is enabled """

    def set_refinement_geman_mcclure_alpha(self, a: float) -> None:
        """ Set the refinement Geman-McClure alpha

        :param a:  is the Geman-McClure alpha, e.g. 0.01 """

    def get_geman_mcclure_alpha(self) -> float:
        """ Get the refinement Geman-McClure alpha

        :return: the Geman-McClure alpha """
        return 0.0

    def generate_self_restraints(self, imol: int, local_dist_max: float) -> int:
        """ Generate GM self restraints for the whole molecule

        :param imol:  is the model molecule index 

        :param local_dist_max:  is the maximum distance, e.g. 4.6 """
        return 0

    def generate_chain_self_restraints(self, imol: int, local_dist_max: float, chain_id: str) -> None:
        """ Generate GM self restraints for the given chain

        :param imol:  is the model molecule index 

        :param local_dist_max:  is the maximum distance, e.g. 4.6 

        :param chain_id:  e.g. "A" for chain A """

    def generate_local_self_restraints(self, imol: int, local_dist_max: float, residue_cids: str) -> None:
        """ Generate GM self restraints for the given residues

        :param imol:  is the model molecule index 

        :param local_dist_max:  is the maximum distance, e.g. 4.6 

        :param residue_cids:  is a "||"-separated list of residues, e.g. "//A/12||//A/14||//B/56" """

    def add_parallel_plane_restraint(self, imol: int, residue_cid_1: str, residue_cid_2: str) -> None:
        """ Generate parallel plane restraints (for RNA and DNA)

        :param imol:  is the model molecule index 

        :param residue_cid_1:  is the selection CID e.g "//A/15" (residue 15 of chain A) 

        :param residue_cid_2:  is the selection CID e.g "//A/17" (residue 17 of chain A) """

    def get_extra_restraints_mesh(self, imol: int, mode: int):
        """ Get the mesh for extra restraints

        :param imol:  is the model molecule index 

        :param mode:  is currently unused """
        pass

    def read_extra_restraints(self, imol: int, file_name: str) -> int:
        """ Read extra restraints (e.g. from ProSMART)

        :param imol:  is the model molecule index """
        return 0

    def clear_extra_restraints(self, imol: int) -> None:
        """ Clear the extra restraints

        :param imol:  is the model molecule index """

    def servalcat_refine_xray(self, imol: int, imol_map: int, output_prefix: str) -> int:
        """ External refinement using servalcat, using data that has already been associated.

        :param imol:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param output_prefix:  is the prefix of the output filename, e.g. "ref-1"

        :return: the imol of the refined model. """
        return 0

    def servalcat_refine_xray_with_keywords(self, imol: int, imol_map: int, output_prefix: str, key_value_pairs: dict) -> int:
        """ Use servalcat keywords

        :param imol:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param output_prefix:  is the prefix of the output filename, e.g. "ref-1" 

        :param key_value_pairs:  is a dictionary of key-value pairs for the servalcat keywords, e.g. resolution: 2.05

        :return: the imol of the refined model. """
        return 0

    def get_rotamer_dodecs(self, imol: int):
        """ Get the rotamer dodecs for the model

        :param imol:  is the model molecule index

        :return: a `simple_mesh_t` """
        pass

    def get_rotamer_dodecs_instanced(self, imol: int):
        """ Get the rotamer dodecs for the model instanced

        :param imol:  is the model molecule index

        :return: an `instanced_mesh_t` """
        pass

    def get_ramachandran_validation_markup_mesh(self, imol: int):
        """ Get the Ramachandran validation markup mesh

        20221126-PE: the function was renamed from ramachandran_validation_markup_mesh().

        :param imol:  is the model molecule index

        :return: a `simple_mesh_t` """
        pass

    def ramachandran_validation(self, imol: int):
        """ Get the data for Ramachandran validation, which importantly contains probability information

        :param imol:  is the model molecule index

        :return: a vector/list of `phi_psi_prob_t` """
        pass

    def contact_dots_for_ligand(self, imol: int, cid: str, smoothness_factor: int):
        """ Contact dots for ligand

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A) 

        :param smoothness_factor:  is 1, 2 or 3 (3 is the most smooth). Recently added (20230202)

        :return: the instanced mesh for the specified ligand """
        pass

    def all_molecule_contact_dots(self, imol: int, smoothness_factor: int):
        """ Contact dots for the whole molecule/model

        :param imol:  is the model molecule index 

        :param smoothness_factor:  is 1, 2 or 3 (3 is the most smooth). Recently added (20230202)

        :return: the instanced mesh for the specified molecule. """
        pass

    def get_simple_molecule(self, imol: int, residue_cid: str, draw_hydrogen_atoms_flag: bool):
        """ Get a simple molecule

        :param imol:  is the model molecule index 

        :param residue_cid:  is the selection CID e.g "//A/15" (residue 15 of chain A) 

        :param draw_hydrogen_atoms_flag:  is the flag for drawing H atoms

        :return: a simple::molecule_t for the specified residue. """
        pass

    def make_exportable_environment_bond_box(self, imol: int, max_dist: float):
        """ 

        :param imol:  is the model molecule index 

        :param spec:  is the residue specifier, e.g. residue_spec_t("A", 10, "") 

        :param max_dist:  specifies the maximum distance of the interaction, typically 3.8

        :return: a vector of lines for non-bonded contacts and hydrogen bonds """
        pass

    def get_h_bonds(self, imol: int, cid_str: str, mcdonald_and_thornton_mode: bool):
        """ Get hydrogen bonds

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A) 

        :param mcdonald_and_thornton_mode:  turns on the McDonald & Thornton algorithm - using explicit hydrogen atoms

        :return: a vector of hydrogen bonds around the specified residue (typically a ligand) """
        pass

    def get_mesh_for_ligand_validation_vs_dictionary(self, imol: int, ligand_cid: str):
        """ Get the mesh for ligand validation vs dictionary, coloured by badness

        Greater then 3 standard deviations is fully red. Less than 0.5 standard deviations is fully green

        :param imol:  is the model molecule index 

        :param ligand_cid:  is the ligand selection CID e.g "//A/15" (ligand 15 of chain A) """
        pass

    def get_ligand_validation_vs_dictionary(self, imol: int, ligand_cid: str, include_non_bonded_contacts: bool):
        """ Ligand validation

        :param imol:  is the model molecule index 

        :param ligand_cid:  is the ligand selection CID e.g "//A/15" (ligand 15 of chain A) 

        :param include_non_bonded_contacts:  is the flag to include non bonded contacts

        :return: a vector/list of interesting geometry - one for each chain involved """
        pass

    def get_validation_vs_dictionary_for_selection(self, imol: int, selection_cid: str, include_non_bonded_contacts: bool):
        """ General fragment distortion analysis

        :param imol:  is the model molecule index 

        :param selection_cid:  is the selection CID e.g "//A/15-23" 

        :param include_non_bonded_contacts:  is the flag to include non bonded contacts

        :return: a vector/list of interesting geometry - one for each chain involved """
        pass

    def get_ligand_distortion(self, imol: int, ligand_cid: str, include_non_bonded_contacts: bool):
        """ Get ligand distortion

        a more simple interface to the above

        :param imol:  is the model molecule index 

        :param selection_cid:  is the selection CID e.g "//A/15-23" 

        :param include_non_bonded_contacts:  is the flag to include non bonded contacts

        :return: a pair: the first is the status (1 for OK, 0 for failed to determine the distortion) """
        pass

    def match_ligand_torsions(self, imol_ligand: int, imol_ref: int, chain_id_ref: str, resno_ref: int) -> bool:
        """ Match ligand torsions

        :param imol_ligand:  is the ligand molecule index 

        :param imol_ref:  is the reference model molecule index 

        :param chain_id_ref:  is the reference chain, e.g. "A" 

        :param resno_ref:  is the reference residue number, e.g. 12

        :return: the success status """
        return True

    def match_ligand_position(self, imol_ligand: int, imol_ref: int, chain_id_ref: str, resno_ref: int) -> bool:
        """ Match ligand positions

        i.e. do a least-squares superposition of the atoms that match in the graphs of the two specified ligands - typically one would use this function after matching ligand torsions.

        :param imol_ligand:  is the ligand molecule index 

        :param imol_ref:  is the reference model molecule index 

        :param chain_id_ref:  is the reference chain, e.g. "A" 

        :param resno_ref:  is the reference residue number, e.g. 12

        :return: the success status """
        return True

    def match_ligand_torsions_and_position(self, imol_ligand: int, imol_ref: int, chain_id_ref: str, resno_ref: int) -> bool:
        """ Match ligand torsions and positions

        :param imol_ligand:  is the ligand molecule index 

        :param imol_ref:  is the reference model molecule index 

        :param chain_id_ref:  is the reference chain, e.g. "A" 

        :param resno_ref:  is the reference residue number, e.g. 12

        :return: the success status. """
        return True

    def match_ligand_torsions_and_position_using_cid(self, imol_ligand: int, imol_ref: int, cid: str) -> bool:
        """ Match ligand torsions and positions using cid

        :param imol_ligand:  is the ligand molecule index 

        :param imol_ref:  is the reference model molecule index 

        :param cid:  is the selection CID e.g "//A/15" (residue 15 of chain A) """
        return True

    def get_overlap_dots(self, imol: int):
        """ 

        :param imol:  is the model molecule index """
        pass

    def get_overlap_dots_for_ligand(self, imol: int, cid_ligand: str):
        """ This function not const because it can dynamically add dictionaries

        :param imol:  is the model molecule index 

        :param cid_ligand:  is the ligand selection CID e.g "//A/15" (ligand 15 of chain A) """
        pass

    def get_overlaps(self, imol: int):
        """ Get Atom Overlaps. 

        :param imol:  is the model molecule index 

        :return: a vector of atom overlap objects """
        pass

    def get_atom_overlap_score(self, imol: int) -> float:
        """ Get the atom overlap score

        :param imol:  the model molecule index 

        :return: the overlap score - a negative number indicates failure """
        return 0.0

    def get_overlaps_for_ligand(self, imol: int, cid_ligand: str):
        """ Gat Atom Overlaps for a ligand or residue. 

        :param imol:  is the model molecule index 

        :param cid_ligand:  is the ligand selection CID e.g "//A/15" (ligand 15 of chain A) 

        :return: a vector of atom overlap objects """
        pass

    def get_atom_differences(self, imol1: int, imol2: int):
        """ Get the atom differences between two molecules typically after refinement

        :param imol1:  is the first model molecule index 

        :param imol2:  is the second model molecule index

        :return: a vector/list of `positioned_atom_spec_t` """
        pass

    def density_fit_analysis(self, imol_model: int, imol_map: int):
        """ Density fit validation information.

        This function returns the sum of the densiy of the atoms in the residue

        :param imol_model:  is the model molecule index 

        :param imol_map:  is the map molecule index

        :return: an object `validation_information_t` """
        pass

    def get_sum_density_for_atoms_in_residue(self, imol: int, cid: str, atom_names: list, imol_map: int):
        """ 

        :return: the sum of the density of the given atoms in the specified CID return -1001 on failure to find the residue or any atoms in the residue or if imol_map is not a map """
        pass

    def get_number_of_atoms_in_residue(self, imol: int, residue_cid: str) -> int:
        """ get the number of atoms in a given residue

        :param imol:  is the model molecule index 

        :param residue_cid:  is the selection CID e.g "//A/15" (residue 15 of chain A) 

        :return: the number of atoms in the residue, or -1 on failure """
        return 0

    def density_correlation_analysis(self, imol_model: int, imol_map: int):
        """ Get the density correlation validation information

        :param imol_model:  is the model molecule index 

        :param imol_map:  is the map molecule index

        :return: an object `validation_information_t` """
        pass

    def rotamer_analysis(self, imol_model: int):
        """ Get the rotamer validation information

        :param imol_model:  is the model molecule index

        :return: an object `validation_information_t` """
        pass

    def ramachandran_analysis(self, imol_model: int):
        """ Get the ramachandran validation information (formatted for a graph, not 3D)

        :param imol_model:  is the model molecule index

        :return: an object `validation_information_t` """
        pass

    def ramachandran_analysis_for_chain(self, imol_model: int, chain_id: str):
        """ Get the ramachandran validation information (formatted for a graph, not 3D) for a given chain in a given molecule

        This function does not exist yet (20230127-PE)

        :param imol_model:  is the model molecule index 

        :param chain_id:  e.g. "A"

        :return: an object `validation_information_t` """
        pass

    def peptide_omega_analysis(self, imol_model: int):
        """ Peptide omega validation information

        :param imol_model:  is the model molecule index

        :return: an object `validation_information_t` """
        pass

    def get_median_temperature_factor(self, imol: int) -> float:
        """ Get the median temperature factor for the model

        :param imol:  is the model molecule index

        :return: a negative number on failure """
        return 0.0

    def get_temperature_factor_of_atom(self, imol: int, atom_cid: str) -> float:
        """ Get the atom temperature factor

        :param imol:  is the model molecule index 

        :param atom_cid:  is the selection cid for the atom

        :return: a negative number on failure, otherwise the temperature factor """
        return 0.0

    def get_interesting_places(self, imol: int, mode: str):
        """ Get interesting places

        This function does not work yet

        :return: a vector/list of `validation_information_t` """
        pass

    def difference_map_peaks(self, imol_map: int, imol_protein: int, n_rmsd: float):
        """ Get difference map peaks

        :param imol_map:  is the map molecule index 

        :param imol_protein:  is the model molecule index 

        :param n_rmsd:  number of sd, e.g. 4.8

        :return: a vector/list of `validation_information_t` """
        pass

    def pepflips_using_difference_map(self, imol_coords: int, imol_difference_map: int, n_sigma: float):
        """ Get pepflips based on the difference map

        :param imol_coords:  is the model molecule index 

        :param imol_difference_map:  is the difference map molecule index 

        :param n_sigma:  number of sd, e.g. 4.8

        :return: a vector/list of `validation_information_t` """
        pass

    def unmodelled_blobs(self, imol_model: int, imol_map: int, rmsd_cut_off: float):
        """ Unmodelled blobs

        :param imol_model:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param rmsd_cut_off:  is the low map limit for cluster generation 1.4 is a reasonable value.

        :return: a vector/list of `validation_information_t` """
        pass

    def find_water_baddies(self, imol_model: int, imol_map: int, b_factor_lim: float, outlier_sigma_level: float, min_dist: float, max_dist: float, ignore_part_occ_contact_flag: bool, ignore_zero_occ_flag: bool):
        """ Check waters, using implicit logical OR

        :param imol_model:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param b_factor_lim:  typical value is 60.0 

        :param outlier_sigma_level:  typical value is 0.8 

        :param min_dist:  typical value is 2.3 

        :param max_dist:  typical value is 3.5 

        :param ignore_part_occ_contact_flag:  typical value is False 

        :param ignore_zero_occ_flag:  typical value is False

        :return: a vector/list of atom specifiers """
        pass

    def get_HOLE(self, imol: int, start_pos_x: float, start_pos_y: float, start_pos_z: float, end_pos_x: float, end_pos_y: float, end_pos_z: float):
        """ Get HOLE

        HOLE is a program for the analysis of the pore dimesions of ion channels. See Smart et al., 1996.

        :return: a list of spheres on the surface of the pore """
        pass

    def mmrrcc(self, imol: int, chain_id: str, imol_map: int):
        """         Sphinx-Doc-Placeholder"""
        pass


    def fourier_shell_correlation(self, imol_map_1: int, imol_map_2: int):
        """ Fourier Shell Correlation (FSC) between maps

        :param imol_map_1:  is the first map molecule index 

        :param imol_map_2:  is the second map molecule index

        :return: a vector/list or pairs of graph points (resolution, correlation). The resolution is in inverse Angstroms squared. An empty list is returned on failure """
        pass

    def make_power_scaled_map(self, imol_ref: int, imol_map_for_scaling: int) -> int:
        """ Make a FSC-scaled map

        :param imol_ref:  is the reference map molecule index 

        :param imol_map_for_scaling:  is the second map molecule index

        :return: the molecule index of the new map """
        return 0

    def get_q_score(self, imol_model: int, imol_map: int):
        """ Get the Q Score (Pintilie et al.)

        :param imol_model:  is the model molecule index 

        :param imol_map:  is the map molecule index

        :return: a `validation_information_t` object """
        pass

    def get_q_score_for_cid(self, imol_model: int, cid: str, imol_map: int):
        """ Get the Pintile et al. Q Score for a particular residue (typically a ligand)

        :param cid:  If the `cid` matches more than one residue the score will be returned for all of the residues covered in the `cid`. Typically, of course the `cid` will be something like "//A/301".

        :return: a `validation_information_t` object """
        pass

    def get_mean_and_variance_of_density_for_non_water_atoms(self, imol_coords: int, imol_map: int):
        """ get mean and variance of map at non-waters

        :param imol_model:  is the model molecule index 

        :param imol_map:  is the map molecule index

        :return: the mean and variance or a negative number on failure """
        pass

    def get_spherical_variance(self, imol_map: int, imol_model: int, atom_cid: str, mean_density_other_atoms: float) -> float:
        """ Get spherical variance - typically for water atoms

        :param imol_model:  is the model molecule index 

        :param imol_map:  is the map molecule index

        :return: the variance or a negative number on failure """
        return 0.0

    def calculate_new_rail_points(self) -> int:
        """ Calling this adds to the rail_points history. Make this pairs when we add model scoring

        :return: the new rail points (since last modification) """
        return 0

    def rail_points_total(self) -> int:
        """ The total rail points

        :return: the sum of all rail points accumulated since the maps were connected. """
        return 0

    def associate_data_mtz_file_with_map(self, imol: int, data_mtz_file_name: str, f_col: str, sigf_col: str, free_r_col: str) -> None:
        """ Associate a data mtz file with a molecule

        This function is called before calling "connect_updating_maps()"

        :param imol:  is the map molecule index 

        :param data_mtz_file_name:  is the name of the mtz file 

        :param f_col:  is the F column, e.g. "FOBS" 

        :param sigf_col:  e.g. "SIGFOBS" 

        :param free_r_col:  e.g. "RFREE" """

    def connect_updating_maps(self, imol_model: int, imol_with_data_info_attached: int, imol_map_2fofc: int, imol_map_fofc: int) -> int:
        """ Connect updating maps

        Reset the rail_points (calls "reset_the_rail_points()"), updates the maps (using internal/clipper SFC). Update your contour lines meshes after calling this function.

        :param imol_model:  is the model molecule index 

        :param imol_with_data_info_attached:  is the map index with the data have been attached by the previous function (associate_data_mtz_file_with_map) 

        :param imol_map_2fofc:  is the map molecule index of the 2FO-FC map 

        :param imol_map_fofc:  is the map molecule index of the FO-FC map

        :return: 1 if the connection was successful """
        return 0

    def sfcalc_genmap(self, imol_model: int, imol_map_with_data_attached: int, imol_updating_difference_map: int) -> None:
        """ Calculate SF and re-generate maps

        This is a low-level function - generally one would use the updating maps method rather than this

        :param imol_model:  is the model molecule index 

        :param imol_map_with_data_attached:  is the map index with the data have been attached by the previous function (associate_data_mtz_file_with_map) 

        :param imol_updating_difference_map:  is the index of the difference map that you want to update when the model updates """

    def sfcalc_genmaps_using_bulk_solvent(self, imol_model: int, imol_2fofc_map: int, imol_updating_difference_map: int, imol_map_with_data_attached: int):
        """ Calculate SF and re-generate maps using bulk solvent

        This is a low-level function. Call this function after connecting maps for updating maps to set the initial R-factor and store the initial map flatness.

        :param imol_model:  is the model molecule index 

        :param imol_2fofc_map:  is the map molecule index of the 2FO-FC map 

        :param imol_updating_difference_map:  is the index of the difference map that you want to update when the model updates 

        :param imol_map_with_data_attached:  is the map index with the data have been attached by the previous function (associate_data_mtz_file_with_map)

        :return: a class of interesting statistics. On failure to calculate SFs and generate the maps the returned r_factor in the returned stats will be set to -1. """
        pass

    def shift_field_b_factor_refinement(self, imol: int, imol_with_data_attached: int) -> bool:
        """ Shift-field B-factor refinement

        This function presumes that the Fobs, sigFobs and RFree data have been filled in the `imol_map_with_data_attached` molecule

        :param imol:  is the model molecule index 

        :param imol_map_with_data_attached:  is the map index with the data have been attached by the previous function (associate_data_mtz_file_with_map)

        :return: success status """
        return True

    def get_density_at_position(self, imol_map: int, x: float, y: float, z: float) -> float:
        """ Get density at position 

        :param imol_map:  is the map molecule index 

        :param x:  is the x coordinate of the target position 

        :param y:  is the y coordinate of the target position 

        :param z:  is the z coordinate of the target position

        :return: the density value """
        return 0.0

    def get_diff_diff_map_peaks(self, imol_diff_map: int, screen_centre_x: float, screen_centre_y: float, screen_centre_z: float):
        """ 

        :param imol_diff_map:  is the map molecule index of the difference map 

        :param screen_centre_x:  is the position x of the center of the screen 

        :param screen_centre_y:  is the position y of the center of the screen 

        :param screen_centre_z:  is the position z of the center of the screen

        :return: a vector of the position where the difference map has been flattened. The associated float value is the amount that the map has been flattened. """
        pass

    def get_data_set_file_name(self, imol: int) -> str:
        """ Get the stored data set file name

        :param imol:  is the model molecule index """
        return 'a-string'

    def go_to_blob(self, x1: float, y1: float, z1: float, x2: float, y2: float, z2: float, contour_level: float):
        """ Given a point on the front clipping plane (x1, y1, z1) and a point on the back clipping plane (x2, y2, z2) this function searches imol_refinement_map (if set) to find a the centre of a blob above the contour level. Blobs at the "front" are selected in preference to blobs at the back. If no blob is found, then the first of the pair is false. If it is found, then the second is the centre of the blob.

        In future, this function should/will be provided with a list of displayed maps and their contour levels - but for now, it uses (only) imol_refinement_map.

        Use a string to pass the map information (map index and contour level), something like "0 0.45:1 1.2:2 0.1" 

        :param x1:  is the x point of the front clipping plane 

        :param y1:  is the y point of the front clipping plane 

        :param z1:  is the z point of the front clipping plane 

        :param x2:  is the x point of the back clipping plane 

        :param y2:  is the y point of the back clipping plane 

        :param z2:  is the z point of the back clipping plane """
        pass

    def fit_ligand_right_here(self, imol_protein: int, imol_map: int, imol_ligand: int, x: float, y: float, z: float, n_rmsd: float, use_conformers: bool, n_conformers: int):
        """ Fit the ligand at specified position. You can expect this to take about 20 seconds. For trivial (i.e non-flexible) ligands you should instead use the jiggle-fit algorithm, which takes a fraction of a second. (That is the algorithm used for "Add Other Solvent Molecules" in Coot.)

        :param imol_protein:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param imol_ligand:  is the ligand molecule index 

        :param x:  is the x position of the blob 

        :param y:  is the y position of the blob 

        :param z:  is the z position of the blob 

        :param n_rmsd:  number of sd, e.g. 4.8 

        :param use_conformers:  is True for flexible ligands 

        :param n_conformers:  set the number of conformers

        :return: a vector/list of indices of molecules for the best fitting ligands to this blob. """
        pass

    def fit_ligand(self, imol_protein: int, imol_map: int, imol_ligand: int, n_rmsd: float, use_conformers: bool, n_conformers: int):
        """ Fit ligand

        :param imol_protein:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param imol_ligand:  is the ligand molecule index 

        :param n_rmsd:  the number of sd used as a cut-off for the map level when finding clusters, e.g. 1.2 

        :param use_conformers:  is True for flexible ligands 

        :param n_conformers:  set the number of conformers

        :return: a vector/list of interesting information about the fitted ligands """
        pass

    def fit_ligand_multi_ligand(self, imol_protein: int, imol_map: int, multi_ligand_molecule_number_list: str, n_rmsd: float, use_conformers: bool, n_conformers: int):
        """ Fit multiple ligands (place-holder)

        :param imol_protein:  is the model molecule index 

        :param imol_map:  is the map molecule index 

        :param multi_ligand_molecule_number_list:  is a colon-separated list of molecules, e.g. "2:3:4" 

        :param n_rmsd:  is number of sd, e.g. 4.8 

        :param use_conformers:  is True for flexible ligands 

        :param n_conformers:  set the number of conformers

        :return: an empty vector (at the moment) """
        pass

    def fit_to_map_by_random_jiggle(self, imol: int, res_spec: str, n_trials: int, translation_scale_factor: float) -> float:
        """ Jiggle-Fit Ligand

        :param imol:  is the model molecule index 

        :param res_spec:  is the residue specifier, e.g. residue_spec_t("A", 10, "") 

        :param n_trials:  is the number of trials, if n_trials is 0, then a sensible default value will be used. 

        :param translation_scale_factor:  if is negative then a sensible default value will be used.

        :return: a value less than -99.9 on failure to fit. """
        return 0.0

    def fit_to_map_by_random_jiggle_using_cid(self, imol: int, cid: str, n_trials: int, translation_scale_factor: float) -> float:
        """ Jiggle-Fit Ligand using cid

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID, e.g "//A/15" (ligand 15 of chain A) 

        :param n_trials:  is the number of trials, if n_trials is 0, then a sensible default value will be used. 

        :param translation_scale_factor:  if is negative then a sensible default value will be used.

        :return: a value less than -99.9 on failure to fit. """
        return 0.0

    def fit_to_map_by_random_jiggle_with_blur_using_cid(self, imol: int, imol_map: int, cid: str, b_factor: float, n_trials: int, translation_scale_factor: float) -> float:
        """ Jiggle-Fit an atom selection, typically a whole molecule or a chain

        :param imol:  is the model molecule index 

        :param cid:  is the selection CID, e.g. "//A" (chain A) 

        :param b_factor:  e.g. 100.0 

        :param n_trials:  is the number of trials, if n_trials is 0, then a sensible default value will be used. 

        :param translation_scale_factor:  if is negative then a sensible default value will be used.

        :return: a value less than -99.9 on failure to fit. """
        return 0.0

    def add_compound(self, imol: int, tlc: str, imol_dict: int, imol_map: int, x: float, y: float, z: float) -> int:
        """ This function is for adding compounds/molecules like buffer agents and precipitants or anions and cations. e.g. those ligands that can be positioned without need for internal torsion angle manipulation.

        :param imol:  is the model molecule index 

        :param tlc:  is the 3-letter-code/compound-id 

        :param imol_dict:  is the molecule to which the ligand is attached (if any). Typically this will be IMOL_ENC_ANY (-999999). 

        :param imol_map:  is the map molecule index 

        :param x:  is the x position 

        :param y:  is the y position 

        :param z:  is the z position

        :return: the success status, 1 for good, 0 for not good. """
        return 0

    def get_svg_for_residue_type(self, imol: int, comp_id: str, use_rdkit_svg: bool, background_type: str) -> str:
        """ Get svg for residue type

        It won't work unless the dictionary for that ligand has been imported. The native output renderings are not very good at the moment. (The RDKit renderings are pretty good).

        :param imol:  is the model molecule index, except for unusual cases, it will be IMOL_ENC_ANY (-999999) 

        :param comp_id:  is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine 

        :param use_rdkit_svg:  is the flag for using the rdkit svg renderer 

        :param background_type:  is one of:

        This function is not const because it caches the svgs.

        :return: the string for the SVG representation. """
        return 'a-string'

    def get_svg_for_2d_ligand_environment_view(self, imol: int, residue_cid: str, add_key: bool) -> str:
        """ Get SVG for 2d ligand environment view (FLEV)

        The caller should make sure that the dictionary for the ligand has been loaded - this function won't do that. It will add hydrogen atoms if needed.

        From time to time (depending on the ligand) this function will fail to produce a result.

        Not const because get_monomer_restraints_at_least_minimal() is called. Hmm.

        :param imol:  is the model molecule index 

        :param residue_cid:  is the cid for the residue 

        :param add_key:  should a key be added to the figure? 

        :return: an svg string of the representation. On failure, return an empty string. """
        return 'a-string'

    def get_non_standard_residues_in_molecule(self, imol: int):
        """ Get non-standard residues in a model

        :param imol:  is the model molecule index

        :return: a vector/list of residue specifiers - the residue name is encoded in the `string_user_data` data item of the residue specifier """
        pass

    def try_read_dictionaries_for_new_residue_types(self, imol: int) -> bool:
        """ Try to read the dictionaries for any residue type in imol that as yet does not have a dictionary

        :param imol:  is the model molecule index 

        :return: true if there were no dictionary for new types that couldn't be read. """
        return True

    def get_dictionary_conformers(self, comp_id: str, imol_enc: int, remove_internal_clash_conformers: bool):
        """ Get the conformers that can be generated by variation around rotatable bonds as described in the dictionary.

        Torsions that are marked as "const" are excluded from the variation, as are pyranose ring torsions and torsions that rotate hydrogen atoms

        :param comp_id:  is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine 

        :param imol_enc:  is the molecule index for the residue type/compound_id 

        :param remove_internal_clash_conformers:  is the flag for removing internal clash

        :return: a vector/list of indices of the new molecules """
        pass

    def get_map_section_texture(self, imol: int, section_id: int, axis: int, data_value_for_bottom: float, data_value_for_top: float):
        """ 

        :param imol:  is the map molecule index 

        :param section_id:  e.g. 2 

        :param axis:  e.g. 0 for X-axis, 1 for Y-axis, 2 for Z-axis`data_value_for_bottom`, `data_value_for_top` should be pre-calculated (don't calculate them for every call to this function).

        :return: a texture_as_floats_t object for the given section. On failure, the image_data vector is empty. """
        pass

    def get_number_of_map_sections(self, imol_map: int, axis_id: int) -> int:
        """ 

        :param imol_map:  is the map molecule index 

        :param axis_id:  is 0 for X-axis, 1 for Y-axis, 2 for Z-axis

        :return: the number of sections in the map along the given axis, -1 on failure. """
        return 0

    def make_mesh_from_gltf_file(self, file_name: str):
        """ 

        :param file_name:  is the glTF file name

        :return: a `simple_mesh_t` from the given file. """
        pass

    def get_octahemisphere(self, n_divisions: int):
        """ Get octahemisphere

        :param n_divisions:  is a number divisible by 2, at least 4 (typically 16)

        :return: a unit-vector end-cap octohemisphere mesh """
        pass

    def get_max_number_of_simple_mesh_vertices(self) -> int:
        """         Sphinx-Doc-Placeholder"""
        return 0

    def set_max_number_of_simple_mesh_vertices(self, n: int) -> None:
        """         Sphinx-Doc-Placeholder"""

    def pae_png(self, pae_file_name: str) -> str:
        """ Predicted alignment error (AlphaFold) 

        :return: a string of a png """
        return 'a-string'

    def make_mesh_for_map_contours_for_blender(self, imol: int, x: float, y: float, z: float, level: float, radius: float) -> None:
        """ Function for Blender interface. """

    def make_mesh_for_bonds_for_blender(self, imol: int, mode: str, against_a_dark_background: bool, bond_width: float, atom_radius_to_bond_width_ratio: float, smoothness_factor: int) -> None:
        """ Function for Blender interface. """

    def make_mesh_for_molecular_representation_for_blender(self, imol: int, cid: str, colour_scheme: str, style: str, secondary_structure_usage_flag: int) -> None:
        """ Function for Blender interface

        Make an (internal) mesh

        This function doesn't return a value, instead it stores a `blender_mesh_t` blender_mesh in this model. One then (shortly later) uses get_triangles_for_blender(imol) (etc) to import this mesh into blender.

        @modifies internal state to fill the internal `blender_mesh` object """

    def make_mesh_for_gaussian_surface_for_blender(self, imol: int, sigma: float, contour_level: float, box_radius: float, grid_scale: float, b_factor: float) -> None:
        """ Function for Blender interface. """

    def make_mesh_for_goodsell_style_for_blender(self, imol: int, colour_wheel_rotation_step: float, saturation: float, goodselliness: float) -> None:
        """ blender Function for Blender interface """

    def get_colour_table_for_blender(self, imol: int):
        """ Function for Blender interface. """
        pass

    def get_vertices_for_blender(self, imol: int):
        """ Function for Blender interface. """
        pass

    def get_triangles_for_blender(self, imol: int):
        """ Function for Blender interface. """
        pass

    def test_function(self, s: str) -> None:
        """         Sphinx-Doc-Placeholder"""





























    def adjust_refinement_residue_name(self, resname: str) -> str:
        """ 

        :param resname:  is the 3 letter code for the residue, e.g. "ALA" for alanine

        :return: the state of having found restraints """
        return 'a-string'

