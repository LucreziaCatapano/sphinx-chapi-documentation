# --- The Virtual Trackball ---

def vt_surface(mode: int) -> None:
        """ 

        :param mode:  1 for "Flat", 2 for "Spherical Surface" """

def vt_surface_status() -> int:
        """ 

        :return: the status, mode=1 for "Flat", mode=2 for "Spherical Surface" """
        return 0

# --- NCS ---

def set_draw_ncs_ghosts(imol: int, istate: int) -> None:
        """ set drawing state of NCS ghosts for molecule number imol """

def draw_ncs_ghosts_state(imol: int) -> int:
        """ return the drawing state of NCS ghosts for molecule number imol. Return -1 on imol is a bad molecule or no ghosts. """
        return 0

def set_ncs_ghost_bond_thickness(imol: int, f: float) -> None:
        """ set bond thickness of NCS ghosts for molecule number imol """

def ncs_update_ghosts(imol: int) -> None:
        """ update ghosts for molecule number imol """

def make_dynamically_transformed_ncs_maps(imol_model: int, imol_map: int, overwrite_maps_of_same_name_flag: int) -> int:
        """ make NCS map """
        return 0

def make_ncs_ghosts_maybe(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def add_ncs_matrix(imol: int, this_chain_id: str, target_chain_id: str, m11: float, m12: float, m13: float, m21: float, m22: float, m23: float, m31: float, m32: float, m33: float, t1: float, t2: float, t3: float) -> None:
        """ Add NCS matrix. """

def clear_ncs_ghost_matrices(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def add_strict_ncs_matrix(imol: int, this_chain_id: str, target_chain_id: str, m11: float, m12: float, m13: float, m21: float, m22: float, m23: float, m31: float, m32: float, m33: float, t1: float, t2: float, t3: float) -> int:
        """ add an NCS matrix for strict NCS molecule representation for CNS strict NCS usage: expand like normal symmetry does """
        return 0

def add_strict_ncs_from_mtrix_from_self_file(imol: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def show_strict_ncs_state(imol: int) -> int:
        """ return the state of NCS ghost molecules for molecule number imol """
        return 0

def set_show_strict_ncs(imol: int, state: int) -> None:
        """ set display state of NCS ghost molecules for molecule number imol """

def set_ncs_homology_level(flev: float) -> None:
        """ At what level of homology should we say that we can't see homology for NCS calculation? (default 0.8) """

def copy_chain(imol: int, from_chain: str, to_chain: str) -> None:
        """ Copy single NCS chain. """

def copy_from_ncs_master_to_others(imol: int, chain_id: str) -> None:
        """ Copy chain from master to all related NCS chains. """

def copy_residue_range_from_ncs_master_to_others(imol: int, master_chain_id: str, residue_range_start: int, residue_range_end: int) -> None:
        """ Copy residue range to all related NCS chains. If the target residues do not exist in the peer chains, then create them. """

def ncs_control_change_ncs_master_to_chain(imol: int, ichain: int) -> None:
        """ change the NCS master chain (by number) """

def ncs_control_change_ncs_master_to_chain_id(imol: int, chain_id: str) -> None:
        """ change the NCS master chain (by chain_id) """

def ncs_control_display_chain(imol: int, ichain: int, state: int) -> None:
        """ display the NCS master chain """

def set_ncs_matrix_type(flag: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_ncs_matrix_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

# --- Startup Functions ---

def set_prefer_python() -> None:
        """ tell coot that you prefer to run python scripts if/when there is an option to do so. """

def prefer_python() -> int:
        """ the python-prefered mode. This is available so that the scripting functions know whether on not to put themselves onto in as menu items.

        If you consider using this, consider in preference use_gui_qm == 2, which is used elsewhere to stop python functions adding to the gui, when guile-gtk functions have alread done so. We should clean up this (rather obscure) interface at some stage.

        return 1 for python is prefered, 0 for not. """
        return 0

# --- File System Functions ---

def make_directory_maybe(dir: str) -> int:
        """ make a directory dir (if it doesn't exist) and return error code If it can be created, create the directory dir, return the success status like mkdir: mkdir

        :return: zero on success, or -1 if an error occurred. If dir already exists as a directory, return 0 of course. """
        return 0

def set_show_paths_in_display_manager(i: int) -> None:
        """ Show Paths in Display Manager? Some people don't like to see the full path names in the display manager here is the way to turn them off, with an argument of 1. """

def show_paths_in_display_manager_state() -> int:
        """ return the internal state What is the internal flag?

        :return: 1 for "yes, display paths" , 0 for not """
        return 0

def add_coordinates_glob_extension(ext: str) -> None:
        """ add an extension to be treated as coordinate files """

def add_data_glob_extension(ext: str) -> None:
        """ add an extension to be treated as data (reflection) files """

def add_dictionary_glob_extension(ext: str) -> None:
        """ add an extension to be treated as geometry dictionary files """

def add_map_glob_extension(ext: str) -> None:
        """ add an extension to be treated as geometry map files """

def remove_coordinates_glob_extension(ext: str) -> None:
        """ remove an extension to be treated as coordinate files """

def remove_data_glob_extension(ext: str) -> None:
        """ remove an extension to be treated as data (reflection) files """

def remove_dictionary_glob_extension(ext: str) -> None:
        """ remove an extension to be treated as geometry dictionary files """

def remove_map_glob_extension(ext: str) -> None:
        """ remove an extension to be treated as geometry map files """

def set_sticky_sort_by_date() -> None:
        """ sort files in the file selection by date? some people like to have their files sorted by date by default """

def unset_sticky_sort_by_date() -> None:
        """ do not sort files in the file selection by date? removes the sorting of files by date """

def set_filter_fileselection_filenames(istate: int) -> None:
        """ on opening a file selection dialog, pre-filter the files. set to 1 to pre-filter, [0 (off, non-pre-filtering) is the default """

def filter_fileselection_filenames_state() -> int:
        """ , return the state of the above variable """
        return 0

def file_type_coords(file_name: str):
        """ is the given file name suitable to be read as coordinates? """
        pass

def open_coords_dialog() -> None:
        """ display the open coordinates dialog """

def set_file_chooser_selector(istate: int) -> None:
        """ this flag set chooser as default for windows, otherwise use selector 0 is selector 1 is chooser """

def file_chooser_selector_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_file_chooser_overwrite(istate: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def file_chooser_overwrite_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def export_map_gui(export_map_fragment: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Widget Utilities ---

def set_main_window_title(s: str) -> None:
        """ set the main window title. function added for Lothar Esser """

# --- MTZ and data handling utilities ---

def manage_column_selector(filename: str) -> None:
        """ given a filename, try to read it as a data file We try as .phs and .cif files first """

# --- Molecule Info Functions ---

def chain_n_residues(chain_id: str, imol: int) -> int:
        """ the number of residues in chain chain_id and molecule number imol 

        :return: the number of residues """
        return 0

def molecule_centre_internal(imol: int, iaxis: int) -> float:
        """ internal function for molecule centre 

        :return: status, less than -9999 is for failure (eg. bad imol); """
        return 0.0

def seqnum_from_serial_number(imol: int, chain_id: str, serial_num: int) -> int:
        """ a residue seqnum (normal residue number) from a residue serial number 

        :return: < -9999 on failure """
        return 0

def insertion_code_from_serial_number(imol: int, chain_id: str, serial_num: int):
        """ the insertion code of the residue. 

        :return: NULL (scheme False) on failure. """
        pass

def n_models(imol: int) -> int:
        """ the chain_id (string) of the ichain-th chain molecule number imol 

        :return: the chain-id

        useful for NMR or other such multi-model molecules.

        return the number of models or -1 if there was a problem with the given molecule. """
        return 0

def n_chains(imol: int) -> int:
        """ number of chains in molecule number imol 

        :return: the number of chains """
        return 0

def is_solvent_chain_p(imol: int, chain_id: str) -> int:
        """ is this a solvent chain? [Raw function] This is a raw interface function, you should generally not use this, but instead use (is-solvent-chain? imol chain-id)

        :return: -1 on error, 0 for no, 1 for is "a solvent chain". We wouldn't want to be doing rotamer searches and the like on such a chain."""
        return 0

def is_protein_chain_p(imol: int, chain_id: str) -> int:
        """ is this a protein chain? [Raw function] This is a raw interface function, you should generally not use this, but instead use (is-protein-chain? imol chain-id)

        :return: -1 on error, 0 for no, 1 for is "a protein chain". We wouldn't want to be doing rotamer searches and the like on such a chain."""
        return 0

def is_nucleotide_chain_p(imol: int, chain_id: str) -> int:
        """ is this a nucleic acid chain? [Raw function] This is a raw interface function, you should generally not use this, but instead use (is-nucleicacid-chain? imol chain-id)

        :return: -1 on error, 0 for no, 1 for is "a nucleicacid chain". We wouldn't want to be doing rotamer searches and the like on such a chain."""
        return 0

def n_residues(imol: int) -> int:
        """ return the number of residues in the molecule, return -1 if this is a map or closed. """
        return 0

def n_atoms(imol: int) -> int:
        """ return the ATOMs of residues in the molecule, return -1 if this is a map or closed. HETATMs are not counted. """
        return 0

def sort_chains(imol: int) -> None:
        """ return a list of the remarks of hte molecule number imol sort the chain ids of the imol-th molecule in lexographical order """

def sort_residues(imol: int) -> None:
        """ sort the residues of the imol-th molecule """

def remarks_dialog(imol: int) -> None:
        """ a gui dialog showing remarks header info (for a model molecule). """

def print_header_secondary_structure_info(imol: int) -> None:
        """ simply print secondary structure info to the terminal/console. In future, this could/should return the info. """

def add_header_secondary_structure_info(imol: int) -> None:
        """ add secondary structure info to the internal representation of the model """

def write_header_secondary_structure_info(imol: int, file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def copy_molecule(imol: int) -> int:
        """ copy molecule imol 

        :return: the new molecule number. Return -1 on failure to copy molecule (out of range, or molecule is closed) """
        return 0

def add_ligand_delete_residue_copy_molecule(imol_ligand_new: int, chain_id_ligand_new: str, resno_ligand_new: int, imol_current: int, chain_id_ligand_current: str, resno_ligand_current: int) -> int:
        """ Copy a molecule with addition of a ligand and a deletion of current ligand. This function is used when adding a new (modified) ligand to a structure. It creates a new molecule that is a copy of the current molecule except that the new ligand is added and the current ligand/residue is deleted. """
        return 0

def exchange_chain_ids_for_seg_ids(imol: int) -> int:
        """ Experimental interface for Ribosome People. Ribosome People have many chains in their pdb file, they prefer segids to chainids (chainids are only 1 character). But coot uses the concept of chain ids and not seg-ids. mmdb allow us to use more than one char in the chainid, so after we read in a pdb, let's replace the chain ids with the segids. Will that help? """
        return 0

def show_remarks_browswer() -> None:
        """ show the remarks browser """

# --- Library and Utility Functions ---

def git_revision_count() -> int:
        """ return the git revision count for for this build. """
        return 0

def svn_revision() -> int:
        """ an alias to """
        return 0

def molecule_name(imol: int):
        """ return the name of molecule number imol 

        :return: 0 if not a valid name ( -> False in scheme) e.g. "/a/b/c.pdb" for "d/e/f.mtz FWT PHWT" """
        pass

def set_molecule_name(imol: int, new_name: str) -> None:
        """ set the molecule name of the imol-th molecule """

def coot_checked_exit(retval: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def coot_real_exit(retval: int) -> None:
        """ exit from coot, give return value retval back to invoking process. """

def coot_no_state_real_exit(retval: int) -> None:
        """ exit without writing a state file """

def coot_clear_backup_or_real_exit(retval: int) -> None:
        """ exit coot doing clear-backup maybe """

def coot_save_state_and_exit(retval: int, save_state_flag: int) -> None:
        """ exit coot, write a state file """

def first_coords_imol() -> int:
        """ What is the molecule number of first coordinates molecule? return -1 when there is none. """
        return 0

def first_small_coords_imol() -> int:
        """ molecule number of first small (<400 atoms) molecule. return -1 on no such molecule """
        return 0

def first_unsaved_coords_imol() -> int:
        """ What is the molecule number of first unsaved coordinates molecule? return -1 when there is none. """
        return 0

def mmcif_sfs_to_mtz(cif_file_name: str, mtz_file_name: str) -> int:
        """ convert the structure factors in cif_file_name to an mtz file. Return 1 on success. Return 0 on a file without Rfree, return -1 on complete failure to write a file. """
        return 0

# --- Graphics Utility Functions ---

def set_do_anti_aliasing(state: int) -> None:
        """ set the bond lines to be antialiased """

def do_anti_aliasing_state() -> int:
        """ return the flag for antialiasing the bond lines """
        return 0

def set_do_GL_lighting(state: int) -> None:
        """ turn the GL lighting on (state = 1) or off (state = 0) slows down the display of simple lines """

def do_GL_lighting_state() -> int:
        """ return the flag for GL lighting """
        return 0

def use_graphics_interface_state():
        """ shall we start up the Gtk and the graphics window? if passed the command line argument 

        An interface function for Ralf. """
        pass

def set_use_dark_mode(state: int) -> None:
        """ set the GUI dark mode state """

def python_at_prompt_at_startup_state():
        """ is the python interpreter at the prompt? 

        :return: 1 for yes, 0 for no. """
        pass

def reset_view() -> int:
        """ "Reset" the view return 1 if we moved, else return 0.

        centre on last-read molecule with zoom 100. If we are there, then go to the previous molecule, if we are there, then go to the origin. """
        return 0

def set_view_rotation_scale_factor(f: float) -> None:
        """ set the view rotation scale factor Useful/necessary for high resolution displayed, where, without this factor the view doesn't rotate enough """

def get_number_of_molecules() -> int:
        """ return the number of molecules (coordinates molecules and map molecules combined) that are currently in coot 

        :return: the number of molecules (closed molecules are not counted) """
        return 0

def graphics_n_molecules() -> int:
        """ As above, return the number of molecules (coordinates molecules and map molecules combined) that are currently in coot. This is the old name for the function.

        :return: the number of molecules (closed molecules are not counted) """
        return 0

def molecule_has_hydrogens_raw(imol: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def own_molecule_number(imol: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def toggle_idle_spin_function() -> None:
        """ Spin spin spin (or not) """

def toggle_idle_rock_function() -> None:
        """ Rock (not roll) (self-timed) """

def get_idle_function_rock_target_angle():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_rocking_factors(width_scale: float, frequency_scale: float) -> None:
        """ Settings for the inevitable discontents who dislike the default rocking rates (defaults 1 and 1) """

def set_idle_function_rotate_angle(f: float) -> None:
        """ how far should we rotate when (auto) spinning? Fast computer? set this to 0.1 """

def idle_function_rotate_angle() -> float:
        """ what is the idle function rotation angle? """
        return 0.0

def make_updating_model_molecule(filename: str) -> int:
        """ make a model molecule from the give file name. If the file updates, then the model will be updated. """
        return 0

def show_calculate_updating_maps_pythonic_gui() -> None:
        """ show the updating maps gui this function is called from callbacks.c and calls a python gui function """

def allow_duplicate_sequence_numbers() -> None:
        """ enable reading PDB/pdbx files with duplicate sequence numbers """

def set_convert_to_v2_atom_names(state: int) -> None:
        """ shall we convert nucleotides to match the old dictionary names? Usually (after 2006 or so) we do not want to do this (given current Coot architecture). Coot should handle the residue synonyms transparently.

        default off (0). """

def assign_hetatms(imol: int) -> int:
        """ some programs produce PDB files with ATOMs where there should be HETATMs. This is a function to assign HETATMs as per the PDB definition. """
        return 0

def hetify_residue(imol: int, chain_id: str, resno: int, ins_code: str) -> int:
        """ if this is not a standard group, then turn the atoms to HETATMs. Return 1 on atoms changes, 0 on not. Return -1 if residue not found. """
        return 0

def residue_has_hetatms(imol: int, chain_id: str, resno: int, ins_code: str) -> int:
        """ residue has HETATMs? return 1 if all atoms of the specified residue are HETATMs, else, return 0. If residue not found, return -1. """
        return 0

def het_group_n_atoms(comp_id: str) -> int:
        """ return the number of non-hydrogen atoms in the given het-group (comp-id). Return -1 on comp-id not found in dictionary. """
        return 0

def replace_fragment(imol_target: int, imol_fragment: int, atom_selection: str) -> int:
        """ replace the parts of molecule number imol that are duplicated in molecule number imol_frag """
        return 0

def copy_residue_range(imol_target: int, chain_id_target: str, imol_reference: int, chain_id_reference: str, resno_range_start: int, resno_range_end: int) -> int:
        """ copy the given residue range from the reference chain to the target chain resno_range_start and resno_range_end are inclusive. """
        return 0

def clear_and_update_model_molecule_from_file(molecule_number: int, file_name: str) -> int:
        """ replace the given residues from the reference molecule to the target molecule replace pdb. Fail if molecule_number is not a valid model molecule. Return -1 on failure. Else return molecule_number """
        return 0

def screendump_image(filename: str) -> None:
        """ dump the current screen image to a file. Format ppm You can use this, in conjunction with spinning and view moving functions to make movies """

def check_for_dark_blue_density() -> None:
        """ give a warning dialog if density it too dark (blue) """

def set_draw_solid_density_surface(imol: int, state: int) -> None:
        """ sets the density map of the given molecule to be drawn as a (transparent) solid surface. """

def set_draw_map_standard_lines(imol: int, state: int) -> None:
        """ toggle for standard lines representation of map. This turns off/on standard lines representation of map. transparent surface is another representation type.

        If you want to just turn off a map, don't use this, use """

def set_solid_density_surface_opacity(imol: int, opacity: float) -> None:
        """ set the opacity of density surface representation of the given map. 0.0 is totally transparent, 1.0 is completely opaque and (because the objects are no longer depth sorted) considerably faster to render. 0.3 is a reasonable number. """

def get_solid_density_surface_opacity(imol: int) -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def set_flat_shading_for_solid_density_surface(state: int) -> None:
        """ set the flag to do flat shading rather than smooth shading for solid density surface. Default is 1 (on. """

# --- Interface Preferences ---

def set_scroll_by_wheel_mouse(istate: int) -> None:
        """ Some people (like Phil Evans) don't want to scroll their map with the mouse-wheel. To turn off mouse wheel recontouring call this with istate value of 0 """

def scroll_by_wheel_mouse_state() -> int:
        """ return the internal state of the scroll-wheel map contouring """
        return 0

def set_auto_recontour_map(state: int) -> None:
        """ turn off (0) or on (1) auto recontouring (on screen centre change) (default it on) """

def get_auto_recontour_map() -> int:
        """ return the auto-recontour state """
        return 0

def set_default_initial_contour_level_for_map(n_sigma: float) -> None:
        """ set the default inital contour for 2FoFc-style map in sigma """

def set_default_initial_contour_level_for_difference_map(n_sigma: float) -> None:
        """ set the default inital contour for FoFc-style map in sigma """

def print_view_matrix() -> None:
        """ print the view matrix to the console, useful for molscript, perhaps """

def get_view_matrix_element(row: int, col: int) -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def get_view_quaternion_internal(element: int) -> float:
        """ internal function to get an element of the view quaternion. The whole quaternion is returned by the scheme function view-quaternion """
        return 0.0

def set_view_quaternion(i: float, j: float, k: float, l: float) -> None:
        """ Set the view quaternion. """

def apply_ncs_to_view_orientation(imol: int, current_chain: str, next_ncs_chain: str) -> None:
        """ Given that we are in chain current_chain, apply the NCS operator that maps current_chain on to next_ncs_chain, so that the relative view is preserved. For NCS skipping. """

def apply_ncs_to_view_orientation_and_screen_centre(imol: int, current_chain: str, next_ncs_chain: str, forward_flag: int) -> None:
        """ as above, but shift the screen centre also. """

def set_show_fps(t: int) -> None:
        """ set show frame-per-second flag """

def set_fps_flag(t: int) -> None:
        """ the old name for """

def get_fps_flag() -> int:
        """ set the state of show frames-per-second flag """
        return 0

def set_show_origin_marker(istate: int) -> None:
        """ set a flag: is the origin marker to be shown? 1 for yes, 0 for no. """

def show_origin_marker_state() -> int:
        """ return the origin marker shown? state """
        return 0

def hide_main_toolbar() -> None:
        """ hide the horizontal main toolbar in the GTK2 version """

def show_main_toolbar() -> None:
        """ show the horizontal main toolbar in the GTK2 version (the toolbar is shown by default) """

def suck_model_fit_dialog() -> int:
        """ reparent the Model/Fit/Refine dialog so that it becomes part of the main window, next to the GL graphics context """
        return 0

def suck_model_fit_dialog_bl() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_model_fit_refine_dialog_stays_on_top(istate: int) -> None:
        """ model-fit-refine dialog stays on top """

def model_fit_refine_dialog_stays_on_top_state() -> int:
        """ return the state model-fit-refine dialog stays on top """
        return 0

def set_accept_reject_dialog_docked(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def accept_reject_dialog_docked_state() -> int:
        """ the accept/reject dialog docked state """
        return 0

def set_accept_reject_dialog_docked_show(state: int) -> None:
        """ set the accept/reject dialog docked show state """

def accept_reject_dialog_docked_show_state() -> int:
        """ what is the accept/reject dialog docked show state? """
        return 0

def set_main_toolbar_style(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def main_toolbar_style_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

# --- Mouse Buttons ---

def quanta_buttons() -> None:
        """ quanta-like buttons Note, when you have set these, there is no way to turn them of again (other than restarting). """

def quanta_like_zoom() -> None:
        """ quanta-like zoom buttons Note, when you have set these, there is no way to turn them of again (other than restarting). """

def set_control_key_for_rotate(state: int) -> None:
        """ Alternate mode for rotation. Prefered by some, including Dirk Kostrewa. I don't think this mode works properly yet """

def control_key_for_rotate_state() -> int:
        """ return the control key rotate state """
        return 0

def blob_under_pointer_to_screen_centre() -> int:
        """ Put the blob under the cursor to the screen centre. Check only positive blobs. Useful function if bound to a key. The refinement map must be set. (We can't check all maps because they are not (or may not be) on the same scale).

        :return: 1 if successfully found a blob and moved there. return 0 if no move. """
        return 0

# --- Cursor Function ---

def normal_cursor() -> None:
        """ normal cursor """

def fleur_cursor() -> None:
        """ fleur cursor """

def pick_cursor_maybe() -> None:
        """ pick cursor maybe """

def rotate_cursor() -> None:
        """ rotate cursor """

def set_pick_cursor_index(icursor_index: int) -> None:
        """ let the user have a different pick cursor sometimes (the default) GDK_CROSSHAIR is hard to see, let the user set their own """

# --- Model/Fit/Refine Functions ---

def show_select_map_frame() -> None:
        """ display the Display Manager dialog """

def set_model_fit_refine_rotate_translate_zone_label(txt: str) -> None:
        """ Allow the changing of Model/Fit/Refine button label from "Rotate/Translate Zone". """

def set_model_fit_refine_place_atom_at_pointer_label(txt: str) -> None:
        """ Allow the changing of Model/Fit/Refine button label from "Place Atom at Pointer". """

def set_refinement_move_atoms_with_zero_occupancy(state: int) -> None:
        """ shall atoms with zero occupancy be moved when refining? (default 1, yes) """

def refinement_move_atoms_with_zero_occupancy_state() -> int:
        """ return the state of "shall atoms with zero occupancy be moved
when refining?" """
        return 0

# --- Backup Functions ---

def make_backup(imol: int) -> None:
        """ make backup for molecule number imol """

def turn_off_backup(imol: int) -> None:
        """ turn off backups for molecule number imol """

def turn_on_backup(imol: int) -> None:
        """ turn on backups for molecule number imol """

def backup_state(imol: int) -> int:
        """ return the backup state for molecule number imol return 0 for backups off, 1 for backups on, -1 for unknown """
        return 0

def apply_undo() -> int:
        """ apply undo - the "Undo" button callback """
        return 0

def apply_redo() -> int:
        """ apply redo - the "Redo" button callback """
        return 0

def set_have_unsaved_changes(imol: int) -> None:
        """ set the molecule number imol to be marked as having unsaved changes """

def have_unsaved_changes_p(imol: int) -> int:
        """ does molecule number imol have unsaved changes? 

        :return: -1 on bad imol, 0 on no unsaved changes, 1 on has unsaved changes """
        return 0

def set_undo_molecule(imol: int) -> None:
        """ set the molecule to which undo operations are done to molecule number imol """

def show_set_undo_molecule_chooser() -> None:
        """ show the Undo Molecule chooser - i.e. choose the molecule to which the "Undo" button applies. """

def set_unpathed_backup_file_names(state: int) -> None:
        """ set the state for adding paths to backup file names by default directories names are added into the filename for backup (with / to _ mapping). call this with state=1 to turn off directory names """

def unpathed_backup_file_names_state() -> int:
        """ return the state for adding paths to backup file names """
        return 0

def set_decoloned_backup_file_names(state: int) -> None:
        """ set the state for adding paths to backup file names by default directories names are added into the filename for backup (with / to _ mapping). call this with state=1 to turn off directory names """

def decoloned_backup_file_names_state() -> int:
        """ return the state for adding paths to backup file names """
        return 0

def backup_compress_files_state() -> int:
        """ return the state for compression of backup files """
        return 0

def set_backup_compress_files(state: int) -> None:
        """ set if backup files will be compressed or not using gzip """

# --- Recover Session Function ---

def recover_session() -> None:
        """ recover session After a crash, we provide this convenient interface to restore the session. It runs through all the molecules with models and looks at the coot backup directory looking for related backup files that are more recent that the read file. (Not very good, because you need to remember which files you read in before the crash - should be improved.) """

# --- Map Functions ---

def calc_phases_generic(mtz_file_name: str) -> None:
        """ fire up a GUI, which asks us which model molecule we want to calc phases from. On "OK" button there, we call """

def map_from_mtz_by_refmac_calc_phases(mtz_file_name: str, f_col: str, sigf_col: str, imol_coords: int) -> int:
        """ Calculate SFs (using refmac optionally) from an MTZ file and generate a map. Get F and SIGF automatically (first of their type) from the mtz file. 

        :return: the new molecule number, -1 on a problem. """
        return 0

def map_from_mtz_by_calc_phases(mtz_file_name: str, f_col: str, sigf_col: str, imol_coords: int) -> int:
        """ Calculate SFs from an MTZ file and generate a map. 

        :return: the new molecule number. """
        return 0

def calculate_maps_and_stats_py(imol_model: int, imol_map_with_data_attached: int, imol_map_2fofc: int, imol_map_fofc: int):
        """ Calculate structure factors and make a 2FoFC map and a Fo-Fc map updating the given molecule numbers for those maps - if thase molecule ids are not valid maps, them generate new maps (return the model number information in the returned object) """
        pass

def sfcalc_genmap(imol_model: int, imol_map_with_data_attached: int, imol_updating_difference_map: int) -> None:
        """ Calculate structure factors from the model and update the given difference map accordingly. """

def set_auto_updating_sfcalc_genmap(imol_model: int, imol_map_with_data_attached: int, imol_updating_difference_map: int) -> None:
        """ As above, calculate structure factors from the model and update the given difference map accordingly - but difference map gets updated automatically on modification of the imol_model molecule. """

def set_auto_updating_sfcalc_genmaps(imol_model: int, imol_map_with_data_attached: int, imol_updating_2fofc_map: int, imol_updating_difference_map: int) -> None:
        """ As above, calculate structure factors from the model and update the given difference map accordingly - but the 2fofc and difference map get updated automatically on modification of the imol_model molecule. """

def set_scroll_wheel_map(imap: int) -> None:
        """ set the map that is moved by changing the scroll wheel and """

def set_scrollable_map(imol: int) -> None:
        """ return the molecule number to which the mouse scroll wheel is attached set the map that has its contour level changed by the scrolling the mouse wheel to molecule number imol (same as """

def scroll_wheel_map() -> int:
        """ the contouring of which map is altered when the scroll wheel changes? """
        return 0

def save_previous_map_colour(imol: int) -> None:
        """ save previous colour map for molecule number imol """

def restore_previous_map_colour(imol: int) -> None:
        """ restore previous colour map for molecule number imol """

def set_active_map_drag_flag(t: int) -> None:
        """ set the state of immediate map upate on map drag. By default, it is on (t=1). On slower computers it might be better to set t=0. """

def get_active_map_drag_flag():
        """ return the state of the dragged map flag """
        pass

def set_last_map_colour(f1: float, f2: float, f3: float) -> None:
        """ set the colour of the last (highest molecule number) map """

def set_map_colour(imol: int, red: float, green: float, blue: float) -> None:
        """ set the colour of the imolth map """

def set_map_hexcolour(imol: int, hex_colour: str) -> None:
        """ set the colour of the imolth map using a (7-character) hex colour """

def set_contour_level_absolute(imol_map: int, level: float) -> None:
        """ set the contour level, direct control """

def set_contour_level_in_sigma(imol_map: int, level: float) -> None:
        """ set the contour level, direct control in r.m.s.d. (if you like that sort of thing) """

def get_contour_level_absolute(imol: int) -> float:
        """ get the contour level """
        return 0.0

def get_contour_level_in_sigma(imol: int) -> float:
        """ get the contour level in rmd above 0. """
        return 0.0

def set_last_map_sigma_step(f: float) -> None:
        """ set the sigma step of the last map to f sigma """

def set_contour_by_sigma_step_by_mol(imol: int, f: float, state: int) -> None:
        """ set the contour level step set the contour level step of molecule number imol to f and variable state (setting state to 0 turns off contouring by sigma level) """

def data_resolution(imol: int) -> float:
        """ return the resolution of the data for molecule number imol. Return negative number on error, otherwise resolution in A (eg. 2.0) """
        return 0.0

def model_resolution(imol: int) -> float:
        """ return the resolution set in the header of the model/coordinates file. If this number is not available, return a number less than 0. """
        return 0.0

def export_map(imol: int, filename: str) -> int:
        """ export (write to disk) the map of molecule number imol to filename. Return 0 on failure, 1 on success. """
        return 0

def export_map_fragment(imol: int, x: float, y: float, z: float, radius: float, filename: str) -> int:
        """ export a fragment of the map about (x,y,z) """
        return 0

def export_map_fragment_with_text_radius(imol: int, radius_text: str, filename: str) -> None:
        """ convenience function, called from callbacks.c """

def export_map_fragment_with_origin_shift(imol: int, x: float, y: float, z: float, radius: float, filename: str) -> int:
        """ export a fragment of the map about (x,y,z) """
        return 0

def export_map_fragment_to_plain_file(imol: int, x: float, y: float, z: float, radius: float, filename: str) -> int:
        """ tmp interface for Hamish """
        return 0

def transform_map_raw(imol: int, r00: float, r01: float, r02: float, r10: float, r11: float, r12: float, r20: float, r21: float, r22: float, t0: float, t1: float, t2: float, pt0: float, pt1: float, pt2: float, box_half_size: float, ref_space_group: str, cell_a: float, cell_b: float, cell_c: float, alpha: float, beta: float, gamma: float) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def difference_map(imol1: int, imol2: int, map_scale: float) -> int:
        """ make a difference map, taking map_scale * imap2 from imap1, on the grid of imap1. Return the new molecule number. Return -1 on failure. """
        return 0

def set_map_has_symmetry(imol: int, state: int) -> None:
        """ by default, maps that are P1 and have 90 degree angles are considered as maps without symmetry (i.e. EM maps). In some cases though P1 maps do/should have symmetry - and this is the means by you can tell Coot that. 

        :param imol:  is the moleculle number to be acted on 

        :param state:  the desired state, a value of 1 turns on map symmetry """

def reinterp_map(map_no: int, reference_map_no: int) -> int:
        """ make a new map (a copy of map_no) that is in the cell, spacegroup and gridding of the map in reference_map_no. Return the new map molecule number - return -1 on failure """
        return 0

def smooth_map(map_no: int, sampling_multiplier: float) -> int:
        """ make a new map (a copy of map_no) that is in the cell, spacegroup and a multiple of the sampling of the input map (a sampling factor of more than 1 makes the output maps smoother) """
        return 0

# --- Density Increment ---

def get_text_for_iso_level_increment_entry(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def get_text_for_diff_map_iso_level_increment_entry(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_iso_level_increment(val: float) -> None:
        """ set the contour scroll step (in absolute e/A3) for 2Fo-Fc-style maps to val The is only activated when scrolling by sigma is turned off """

def get_iso_level_increment() -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def set_iso_level_increment_from_text(text: str, imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_diff_map_iso_level_increment(val: float) -> None:
        """ set the contour scroll step for difference map (in absolute e/A3) to val The is only activated when scrolling by sigma is turned off """

def get_diff_map_iso_level_increment() -> float:
        """ return difference maps iso-map level increment """
        return 0.0

def set_diff_map_iso_level_increment_from_text(text: str, imol: int) -> None:
        """ set the difference maps iso-map level increment """

def set_map_sampling_rate_text(text: str) -> None:
        """ sampling rate find the molecule for which the single map dialog applies and set the contour level and redraw """

def set_map_sampling_rate(r: float) -> None:
        """ set the map sampling rate (default 1.5) Set to something like 2.0 or 2.5 for more finely sampled maps. Useful for baton-building low resolution maps. """

def get_text_for_map_sampling_rate_text():
        """ Sphinx-Doc-Placeholder"""
        pass

def get_map_sampling_rate() -> float:
        """ return the map sampling rate """
        return 0.0

def change_contour_level(is_increment: int) -> None:
        """ change the contour level of the current map by a step if is_increment=1 the contour level is increased. If is_increment=0 the map contour level is decreased. """

def set_last_map_contour_level(level: float) -> None:
        """ set the contour level of the map with the highest molecule number to level """

def set_last_map_contour_level_by_sigma(n_sigma: float) -> None:
        """ set the contour level of the map with the highest molecule number to n_sigma sigma """

def set_stop_scroll_diff_map(i: int) -> None:
        """ create a lower limit to the "Fo-Fc-style" map contour level changing (default 1 on) """

def set_stop_scroll_iso_map(i: int) -> None:
        """ create a lower limit to the "2Fo-Fc-style" map contour level changing (default 1 on) """

def set_stop_scroll_iso_map_level(f: float) -> None:
        """ set the actual map level changing limit (default 0.0) """

def set_stop_scroll_diff_map_level(f: float) -> None:
        """ set the actual difference map level changing limit (default 0.0) """

def set_residue_density_fit_scale_factor(f: float) -> None:
        """ set the scale factor for the Residue Density fit analysis """

# --- Density Functions ---

def set_map_line_width(w: int) -> None:
        """ draw the lines of the chickenwire density in width w """

def map_line_width_state() -> int:
        """ return the width in which density contours are drawn """
        return 0

def make_and_draw_map(mtz_file_name: str, f_col: str, phi_col: str, weight: str, use_weights: int, is_diff_map: int) -> int:
        """ make a map from an mtz file (simple interface) given mtz file mtz_file_name and F column f_col and phases column phi_col and optional weight column weight_col (pass use_weights=0 if weights are not to be used). Also mark the map as a difference map (is_diff_map=1) or not (is_diff_map=0) because they are handled differently inside coot.

        :return: -1 on error, else return the new molecule number """
        return 0

def read_mtz(mtz_file_name: str, f_col: str, phi_col: str, weight: str, use_weights: int, is_diff_map: int) -> int:
        """ the function is a synonym of the above function - which now has an archaic-style name """
        return 0

def make_and_draw_map_with_refmac_params(mtz_file_name: str, a: str, b: str, weight: str, use_weights: int, is_diff_map: int, have_refmac_params: int, fobs_col: str, sigfobs_col: str, r_free_col: str, sensible_f_free_col: int) -> int:
        """ as the above function, execpt set refmac parameters too pass along the refmac column labels for storage (not used in the creation of the map)

        :return: -1 on error, else return imol """
        return 0

def make_and_draw_map_with_reso_with_refmac_params(mtz_file_name: str, a: str, b: str, weight: str, use_weights: int, is_diff_map: int, have_refmac_params: int, fobs_col: str, sigfobs_col: str, r_free_col: str, sensible_f_free_col: int, is_anomalous: int, use_reso_limits: int, low_reso_limit: float, high_reso_lim: float) -> int:
        """ as the above function, except set expert options too. """
        return 0

def make_updating_map(mtz_file_name: str, f_col: str, phi_col: str, weight: str, use_weights: int, is_diff_map: int) -> int:
        """ make a map molecule from the give file name. If the file updates, then the map will be updated. """
        return 0

def stop_updating_molecule(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def mtz_file_has_phases_p(mtz_file_name: str) -> int:
        """ does the mtz file have phases? """
        return 0

def is_mtz_file_p(filename: str) -> int:
        """ is the given filename an mtz file? """
        return 0

def cns_file_has_phases_p(cns_file_name: str) -> int:
        """ does the given file have cns phases? """
        return 0

def wrapped_auto_read_make_and_draw_maps(filename: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_auto_read_do_difference_map_too(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def auto_read_do_difference_map_too_state() -> int:
        """ return the flag to do a difference map (too) on auto-read MTZ 

        :return: 0 means no, 1 means yes. """
        return 0

def set_auto_read_column_labels(fwt: str, phwt: str, is_for_diff_map_flag: int) -> None:
        """ set the expected MTZ columns for Auto-reading MTZ file. Not every program uses the default refmac labels ("FWT"/"PHWT") for its MTZ file. Here we can tell coot to expect other labels so that coot can "Auto-open" such MTZ files.

        e.g. (set-auto-read-column-labels "2FOFCWT" "PH2FOFCWT" 0) """

def get_text_for_density_size_widget():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_density_size_from_widget(text: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_text_for_density_size_em_widget():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_density_size_em_from_widget(text: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_map_radius(f: float) -> None:
        """ set the extent of the box/radius of electron density contours for x-ray maps """

def set_map_radius_em(radius: float) -> None:
        """ set the extent of the box/radius of electron density contours for EM map """

def set_density_size(f: float) -> None:
        """ another (old) way of setting the radius of the map """

def set_map_radius_slider_max(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_display_intro_string(str: str) -> None:
        """ Give me this nice message str when I start coot. """

def get_map_radius() -> float:
        """ return the extent of the box/radius of electron density contours """
        return 0.0

def set_esoteric_depth_cue(istate: int) -> None:
        """ not everone likes coot's esoteric depth cueing system Pass an argument istate=1 to turn it off

        (this function is currently disabled). """

def esoteric_depth_cue_state() -> int:
        """ native depth cueing system return the state of the esoteric depth cueing flag """
        return 0

def set_swap_difference_map_colours(i: int) -> None:
        """ not everone likes coot's default difference map colouring. Pass an argument i=1 to swap the difference map colouring so that red is positive and green is negative. """

def swap_difference_map_colours_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_map_is_difference_map(imol: int, bool_flag: int) -> int:
        """ post-hoc set the map of molecule number imol to be a difference map 

        :return: success status, 0 -> failure (imol does not have a map) """
        return 0

def map_is_difference_map(imol: int) -> int:
        """ map is difference map? """
        return 0

def another_level() -> int:
        """ Add another contour level for the last added map. Currently, the map must have been generated from an MTZ file. 

        :return: the molecule number of the new molecule or -1 on failure """
        return 0

def another_level_from_map_molecule_number(imap: int) -> int:
        """ Add another contour level for the given map. Currently, the map must have been generated from an MTZ file. 

        :return: the molecule number of the new molecule or -1 on failure """
        return 0

def residue_density_fit_scale_factor() -> float:
        """ return the scale factor for the Residue Density fit analysis """
        return 0.0

def density_at_point(imol_map: int, x: float, y: float, z: float) -> float:
        """ return the density at the given point for the given map. Return 0 for bad imol """
        return 0.0

# --- Parameters from map ---

def mtz_hklin_for_map(imol_map: int):
        """ return the mtz file that was use to generate the map return 0 when there is no mtz file associated with that map (it was generated from a CCP4 map file say). """
        pass

def mtz_fp_for_map(imol_map: int):
        """ return the FP column in the file that was use to generate the map return 0 when there is no mtz file associated with that map (it was generated from a CCP4 map file say).

        Caller should dispose of returned pointer. """
        pass

def mtz_phi_for_map(imol_map: int):
        """ return the phases column in mtz file that was use to generate the map return 0 when there is no mtz file associated with that map (it was generated from a CCP4 map file say). Caller should dispose of returned pointer. """
        pass

def mtz_weight_for_map(imol_map: int):
        """ return the weight column in the mtz file that was use to generate the map return 0 when there is no mtz file associated with that map (it was generated from a CCP4 map file say) or no weights were used. Caller should dispose of returned pointer. """
        pass

def mtz_use_weight_for_map(imol_map: int):
        """ return flag for whether weights were used that was use to generate the map return 0 when no weights were used or there is no mtz file associated with that map. """
        pass

# --- PDB Functions ---

def write_pdb_file(imol: int, file_name: str) -> int:
        """ write molecule number imol as a PDB to file file_name """
        return 0

def write_cif_file(imol: int, file_name: str) -> int:
        """ write molecule number imol as a mmCIF to file file_name """
        return 0

def write_residue_range_to_pdb_file(imol: int, chainid: str, resno_start: int, resno_end: int, filename: str) -> int:
        """ write molecule number imol's residue range as a PDB to file file_name """
        return 0

def write_chain_to_pdb_file(imol: int, chainid: str, filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def quick_save() -> int:
        """ save all modified coordinates molecules to the default names and save the state too. """
        return 0

def get_write_conect_record_state() -> int:
        """ return the state of the write_conect_records_flag. """
        return 0

def set_write_conect_record_state(state: int) -> None:
        """ set the flag to write (or not) conect records to the PDB file. """

# --- Info Dialog ---

def info_dialog(txt: str) -> None:
        """ create a dialog with information create a dialog with information string txt. User has to click to dismiss it, but it is not modal (nothing in coot is modal). """

def info_dialog_and_text(txt: str) -> None:
        """ create a dialog with information and print to console as info_dialog but print to console as well. """

def info_dialog_with_markup(txt: str) -> None:
        """ as above, create a dialog with information This dialog is left-justified and can use markup such as angled bracketted tt or i """

# --- Refmac Functions ---

def set_refmac_counter(imol: int, refmac_count: int) -> None:
        """ set counter for runs of refmac so that this can be used to construct a unique filename for new output """

def swap_map_colours(imol1: int, imol2: int) -> None:
        """ swap the colours of maps swap the colour of maps imol1 and imol2. Useful to some after running refmac, so that the map to be build into is always the same colour """

def set_keep_map_colour_after_refmac(istate: int) -> None:
        """ flag to enable above call this with istate=1 """

def keep_map_colour_after_refmac_state() -> int:
        """ the keep-map-colour-after-refmac internal state 

        :return: 1 for "yes", 0 for "no" """
        return 0

# --- Symmetry Functions ---

def get_text_for_symmetry_size_widget():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_symmetry_size_from_widget(text: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_symmetry_size(f: float) -> None:
        """ set the size of the displayed symmetry """

def get_symmetry_bonds_colour(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def get_show_symmetry():
        """ is symmetry master display control on? """
        pass

def set_show_symmetry_master(state: int) -> None:
        """ set display symmetry, master controller """

def set_show_symmetry_molecule(mol_no: int, state: int) -> None:
        """ set display symmetry for molecule number mol_no pass with state=0 for off, state=1 for on """

def symmetry_as_calphas(mol_no: int, state: int) -> None:
        """ display symmetry as CAs? pass with state=0 for off, state=1 for on """

def get_symmetry_as_calphas_state(imol: int):
        """ what is state of display CAs for molecule number mol_no? return state=0 for off, state=1 for on """
        pass

def set_symmetry_molecule_rotate_colour_map(imol: int, state: int) -> None:
        """ set the colour map rotation (i.e. the hue) for the symmetry atoms of molecule number imol """

def symmetry_molecule_rotate_colour_map_state(imol: int) -> int:
        """ should there be colour map rotation (i.e. the hue) change for the symmetry atoms of molecule number imol? return state=0 for off, state=1 for on """
        return 0

def set_symmetry_colour_by_symop(imol: int, state: int) -> None:
        """ set symmetry colour by symop mode """

def set_symmetry_whole_chain(imol: int, state: int) -> None:
        """ set symmetry colour for the chain """

def set_symmetry_atom_labels_expanded(state: int) -> None:
        """ set use expanded symmetry atom labels """

def has_unit_cell_state(imol: int) -> int:
        """ molecule number imol has a unit cell? 

        :return: 1 on "yes, it has a cell", 0 for "no" """
        return 0

def add_symmetry_on_to_preferences_and_apply() -> None:
        """ Sphinx-Doc-Placeholder"""

def undo_symmetry_view() -> int:
        """ Undo symmetry view. Translate back to main molecule from this symmetry position. """
        return 0

def first_molecule_with_symmetry_displayed() -> int:
        """ return the molecule number. 

        :return: -1 if there is no molecule with symmetry displayed. """
        return 0

def save_symmetry_coords(imol: int, filename: str, symop_no: int, shift_a: int, shift_b: int, shift_c: int, pre_shift_to_origin_na: int, pre_shift_to_origin_nb: int, pre_shift_to_origin_nc: int) -> None:
        """ save the symmetry coordinates of molecule number imol to filename Allow a shift of the coordinates to the origin before symmetry expansion is apllied (this is how symmetry works in Coot internals). """

def new_molecule_by_symmetry(imol: int, name: str, m11: float, m12: float, m13: float, m21: float, m22: float, m23: float, m31: float, m32: float, m33: float, tx: float, ty: float, tz: float, pre_shift_to_origin_na: int, pre_shift_to_origin_nb: int, pre_shift_to_origin_nc: int) -> int:
        """ create a new molecule (molecule number is the return value) from imol. The rotation/translation matrix components are given in 

        Allow a shift of the coordinates to the origin before symmetry expansion is aplied.

        Pass "" as the name-in and a name will be constructed for you.

        Return -1 on failure. """
        return 0

def new_molecule_by_symmetry_with_atom_selection(imol: int, name: str, mmdb_atom_selection_string: str, m11: float, m12: float, m13: float, m21: float, m22: float, m23: float, m31: float, m32: float, m33: float, tx: float, ty: float, tz: float, pre_shift_to_origin_na: int, pre_shift_to_origin_nb: int, pre_shift_to_origin_nc: int) -> int:
        """ create a new molecule (molecule number is the return value) from imol, but only for atom that match the mmdb_atom_selection_string. The rotation/translation matrix components are given in 

        Allow a shift of the coordinates to the origin before symmetry expansion is aplied.

        Pass "" as the name-in and a name will be constructed for you.

        Return -1 on failure. """
        return 0

def new_molecule_by_symop(imol: int, symop_string: str, pre_shift_to_origin_na: int, pre_shift_to_origin_nb: int, pre_shift_to_origin_nc: int) -> int:
        """ create a new molecule (molecule number is the return value) from imol. """
        return 0

def n_symops(imol: int) -> int:
        """ return the number of symmetry operators for the given molecule return -1 on no-symmetry for molecule or inappropriate imol number """
        return 0

def move_reference_chain_to_symm_chain_position() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def setup_save_symmetry_coords() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_space_group(imol: int, spg: str):
        """ set the space group for a coordinates molecule for shelx FA pdb files, there is no space group. So allow the user to set it. This can be initted with a HM symbol or a symm list for clipper.

        This will only work on model molecules.

        :return: the success status of the setting (1 good, 0 fail). """
        pass

def set_unit_cell_and_space_group(imol: int, a: float, b: float, c: float, alpha: float, beta: float, gamma: float, space_group: str) -> int:
        """ set the unit cell for a given model molecule Angles in degress, cell lengths in Angstroms.

        :return: the success status of the setting (1 good, 0 fail). """
        return 0

def set_unit_cell_and_space_group_using_molecule(imol: int, imol_from: int) -> int:
        """ set the unit cell for a given model molecule using the cell of moecule imol_from This will only work on model molecules. 

        :return: the success status of the setting (1 good, 0 fail). """
        return 0

def set_symmetry_shift_search_size(shift: int) -> None:
        """ set the cell shift search size for symmetry searching. When the coordinates for one (or some) symmetry operator are missing (which happens sometimes, but rarely), try changing setting this to 2 (default is 1). It slows symmetry searching, which is why it is not set to 2 by default. """

# --- History Functions ---

def print_all_history_in_scheme() -> None:
        """ print the history in scheme format """

def print_all_history_in_python() -> None:
        """ print the history in python format """

def set_console_display_commands_state(istate: int) -> None:
        """ set a flag to show the text command equivalent of gui commands in the console as they happen. 1 for on, 0 for off. """

def set_console_display_commands_hilights(bold_flag: int, colour_flag: int, colour_index: int) -> None:
        """ set a flag to show the text command equivalent of gui commands in the console as they happen in bold and colours. colour_flag: pass 1 for on, 0 for off.

        colour_index 0 to 7 inclusive for various different colourings. """

# --- State Functions ---

def save_state() -> None:
        """ save the current state to the default filename """

def save_state_file(filename: str) -> None:
        """ save the current state to file filename """

def save_state_file_py(filename: str) -> None:
        """ save the current state to file filename """

def set_save_state_file_name(filename: str) -> None:
        """ set the default state file name (default 0-coot.state.scm) """

def save_state_file_name_raw():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_run_state_file_status(istat: int) -> None:
        """ set run state file status 0: never run it 1: ask to run it 2: run it, no questions """

def run_state_file() -> None:
        """ run the state file (reading from default filenname) """

def run_state_file_py() -> None:
        """ Sphinx-Doc-Placeholder"""

def run_state_file_maybe() -> None:
        """ run the state file depending on the state variables """

# --- Clipping Functions ---

def increase_clipping_front() -> None:
        """ increase the amount of clipping, that is (independent of projection matrix) """

def increase_clipping_back() -> None:
        """ increase the amount of clipping, that is (independent of projection matrix) """

def decrease_clipping_front() -> None:
        """ decrease the amount of clipping, that is (independent of projection matrix) """

def decrease_clipping_back() -> None:
        """ decrease the amount of clipping, that is (independent of projection matrix) """

def set_clipping_back(v: float) -> None:
        """ set clipping plane back - this goes in differnent directions for orthographics vs perspective """

def set_clipping_front(v: float) -> None:
        """ set clipping plane front - this goes in differnent directions for orthographics vs perspective """

def get_clipping_plane_front() -> float:
        """ get clipping plane front """
        return 0.0

def get_clipping_plane_back() -> float:
        """ get clipping plane back """
        return 0.0

# --- Unit Cell interface ---

def get_show_unit_cell(imol: int):
        """ return the stage of show unit cell for molecule number imol """
        pass

def set_show_unit_cells_all(istate: int) -> None:
        """ set the state of show unit cell for all molecules 1 for displayed 0 for undisplayed """

def set_show_unit_cell(imol: int, istate: int) -> None:
        """ set the state of show unit cell for the particular molecule number imol 1 for displayed 0 for undisplayed """

def set_unit_cell_colour(red: float, green: float, blue: float) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Colour ---

def set_symmetry_colour_merge(v: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_colour_map_rotation_on_read_pdb(f: float) -> None:
        """ set the hue change step on reading a new molecule """

def set_colour_map_rotation_on_read_pdb_flag(i: int) -> None:
        """ shall the hue change step be used? 

        :param i:  0 for no, 1 for yes """

def set_colour_map_rotation_on_read_pdb_c_only_flag(i: int) -> None:
        """ shall the colour map rotation apply only to C atoms? 

        :param i:  0 for no, 1 for yes """

def set_colour_by_chain(imol: int) -> None:
        """ colour molecule number imol by chain type """

def set_colour_by_ncs_chain(imol: int, goodsell_mode: int) -> None:
        """ colour molecule number imol by chain type """

def set_colour_by_chain_goodsell_mode(imol: int) -> None:
        """ colour molecule number imol by chain type, goodsell-like colour scheme """

def set_goodsell_chain_colour_wheel_step(s: float) -> None:
        """ set the goodsell chain colour colour wheel step (default 0.22) """

def set_colour_by_molecule(imol: int) -> None:
        """ colour molecule number imol by molecule """

def get_colour_map_rotation_on_read_pdb_c_only_flag() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_symmetry_colour(r: float, g: float, b: float) -> None:
        """ set the symmetry colour base """

# --- Map colour ---

def set_colour_map_rotation_for_map(f: float) -> None:
        """ set the colour map rotation (hue change) for maps default: for maps is 14 degrees. """

def set_molecule_bonds_colour_map_rotation(imol: int, theta: float) -> None:
        """ set the colour map rotation for molecule number imol theta is in degrees """

def get_molecule_bonds_colour_map_rotation(imol: int) -> float:
        """ Get the colour map rotation for molecule number imol. """
        return 0.0

# --- Anisotropic Atoms Interface ---

def get_limit_aniso() -> float:
        """ get the aniso radius limit """
        return 0.0

def get_show_limit_aniso():
        """ get show the aniso limit """
        pass

def get_show_aniso():
        """ return show-aniso-atoms state - FIXME- per molecule """
        pass

def set_limit_aniso(state: int) -> None:
        """ set the aniso atom limit """

def set_show_aniso(state: int) -> None:
        """ does nothing """

def set_show_aniso_atoms(imol: int, state: int) -> None:
        """ set show aniso atoms """

def set_show_aniso_atoms_as_ortep(imol: int, state: int) -> None:
        """ set show aniso atoms as ortep """

def set_aniso_limit_size_from_widget(text: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_text_for_aniso_limit_radius_entry():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_aniso_probability(f: float) -> None:
        """ DELETE-ME """

def get_aniso_probability() -> float:
        """ DELETE-ME """
        return 0.0

# --- Display Functions ---

def set_graphics_window_size(x_size: int, y_size: int) -> None:
        """ set the window size """

def set_graphics_window_size_internal(x_size: int, y_size: int, as_widget_flag: int) -> None:
        """ set the window size as gtk_widget (flag=1) or gtk_window (flag=0) """

def set_graphics_window_position(x_pos: int, y_pos: int) -> None:
        """ set the graphics window position """

def store_graphics_window_position(x_pos: int, y_pos: int) -> None:
        """ store the graphics window position """

def graphics_window_size_and_position_to_preferences() -> None:
        """ store the graphics window position and size to zenops-graphics-window-size-and-postion.scm in the preferences directory. """

def graphics_draw() -> None:
        """ draw a frame """

def zalman_stereo_mode() -> None:
        """ try to turn on Zalman stereo mode """

def hardware_stereo_mode() -> None:
        """ try to turn on stereo mode """

def set_stereo_style(mode: int) -> None:
        """ set the stereo mode (the relative view of the eyes) 0 is 2010-mode 1 is modern mode """

def stereo_mode_state() -> int:
        """ what is the stero state? 

        :return: 1 for in hardware stereo, 2 for side by side stereo, else return 0. """
        return 0

def mono_mode() -> None:
        """ try to turn on mono mode """

def side_by_side_stereo_mode(use_wall_eye_mode: int) -> None:
        """ turn on side bye side stereo mode """

def set_dti_stereo_mode(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_hardware_stereo_angle_factor(f: float) -> None:
        """ how much should the eyes be separated in stereo mode? 

        :param f:  the angular difference (in multiples of 4.5 degrees) """

def hardware_stereo_angle_factor_state() -> float:
        """ return the hardware stereo angle factor """
        return 0.0

def set_model_display_radius(state: int, radius: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_model_fit_refine_dialog_position(x_pos: int, y_pos: int) -> None:
        """ set position of Model/Fit/Refine dialog """

def set_display_control_dialog_position(x_pos: int, y_pos: int) -> None:
        """ set position of Display Control dialog """

def set_go_to_atom_window_position(x_pos: int, y_pos: int) -> None:
        """ set position of Go To Atom dialog """

def set_delete_dialog_position(x_pos: int, y_pos: int) -> None:
        """ set position of Delete dialog """

def set_rotate_translate_dialog_position(x_pos: int, y_pos: int) -> None:
        """ set position of the Rotate/Translate Residue Range dialog """

def set_accept_reject_dialog_position(x_pos: int, y_pos: int) -> None:
        """ set position of the Accept/Reject dialog """

def set_ramachandran_plot_dialog_position(x_pos: int, y_pos: int) -> None:
        """ set position of the Ramachadran Plot dialog """

def set_edit_chi_angles_dialog_position(x_pos: int, y_pos: int) -> None:
        """ set edit chi angles dialog position """

def set_rotamer_selection_dialog_position(x_pos: int, y_pos: int) -> None:
        """ set rotamer selection dialog position """

# --- Smooth Scrolling ---

def set_smooth_scroll_flag(v: int) -> None:
        """ set smooth scrolling 

        :param v:  use v=1 to turn on smooth scrolling, v=0 for off (default on). """

def get_smooth_scroll() -> int:
        """ return the smooth scrolling state """
        return 0

def set_smooth_scroll_steps_str(t: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_smooth_scroll_steps(i: int) -> None:
        """ set the number of steps in the smooth scroll Set more steps (e.g. 50) for more smoothness (default 10). """

def get_text_for_smooth_scroll_steps():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_smooth_scroll_limit_str(t: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_smooth_scroll_limit(lim: float) -> None:
        """ do not scroll for distances greater this limit """

def get_text_for_smooth_scroll_limit():
        """ Sphinx-Doc-Placeholder"""
        pass

# --- Font Parameters ---

def set_font_size(i: int) -> None:
        """ set the font size 

        :param i:  1 (small) 2 (medium, default) 3 (large) """

def get_font_size() -> int:
        """ return the font size 

        :return: 1 (small) 2 (medium, default) 3 (large) """
        return 0

def set_font_colour(red: float, green: float, blue: float) -> None:
        """ set the colour of the atom label font - the arguments are in the range 0->1 """

def set_use_stroke_characters(state: int) -> None:
        """ set use stroke characters """

# --- Rotation Centre ---

def set_rotation_centre_size_from_widget(text: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_text_for_rotation_centre_cube_size():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_rotation_centre_size(f: float) -> None:
        """ set rotoation centre marker size """

def set_rotation_centre_cross_hairs_colour(r: float, g: float, b: float, alpha: float) -> None:
        """ set rotation centre colour This is the colour for a dark background - if the background colour is not dark, then the cross-hair colour becomes the inverse colour """

def recentre_on_read_pdb():
        """ return the recentre-on-pdb state """
        pass

def set_rotation_centre(x: float, y: float, z: float) -> None:
        """ set the rotation centre """

def set_rotation_centre_internal(x: float, y: float, z: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def rotation_centre_position(axis: int) -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def go_to_ligand() -> None:
        """ centre on the ligand of the "active molecule", if we are already there, centre on the next hetgroup (etc) """

def set_go_to_ligand_n_atoms_limit(n_atom_min: int) -> None:
        """ go to the ligand that has more than n_atom_min atoms """

def set_reorienting_next_residue_mode(state: int) -> None:
        """ rotate the view so that the next main-chain atoms are oriented in the same direction as the previous - hence side-chain always seems to be "up" - set this mode to 1 for reorientation-mode - and 0 for off (standard translation) """

# --- Orthogonal Axes ---

def set_draw_axes(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Atom Selection Utilities ---

def atom_index(imol: int, chain_id: str, iresno: int, atom_id: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def atom_index_full(imol: int, chain_id: str, iresno: int, inscode: str, atom_id: str, altconf: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def atom_index_first_atom_in_residue(imol: int, chain_id: str, iresno: int, ins_code: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def atom_index_first_atom_in_residue_with_altconf(imol: int, chain_id: str, iresno: int, ins_code: str, alt_conf: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def min_resno_in_chain(imol: int, chain_id: str) -> int:
        """ return the minimum residue number for imol chain chain_id """
        return 0

def max_resno_in_chain(imol: int, chain_id: str) -> int:
        """ return the maximum residue number for imol chain chain_id """
        return 0

def median_temperature_factor(imol: int) -> float:
        """ return the median temperature factor for imol """
        return 0.0

def average_temperature_factor(imol: int) -> float:
        """ return the average temperature factor for the atoms in imol """
        return 0.0

def standard_deviation_temperature_factor(imol: int) -> float:
        """ return the standard deviation of the atom temperature factors for imol """
        return 0.0

def clear_pending_picks() -> None:
        """ clear pending picks (stop coot thinking that the user is about to pick an atom). """

def centre_of_mass_string(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def centre_of_mass_string_py(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_default_temperature_factor_for_new_atoms(new_b: float) -> None:
        """ set the default temperature factor for newly created atoms (initial default 20) """

def default_new_atoms_b_factor() -> float:
        """ return the default temperature factor for newly created atoms """
        return 0.0

def set_reset_b_factor_moved_atoms(state: int) -> None:
        """ reset temperature factor for all moved atoms to the default for new atoms (usually 30) """

def get_reset_b_factor_moved_atoms_state() -> int:
        """ return the state if temperature factors shoudl be reset for moved atoms """
        return 0

def set_atom_attribute(imol: int, chain_id: str, resno: int, ins_code: str, atom_name: str, alt_conf: str, attribute_name: str, val: float) -> int:
        """ set a numberical attibute to the atom with the given specifier. Attributes can be "x", "y","z", "B", "occ" and the attribute val is a floating point number """
        return 0

def set_atom_string_attribute(imol: int, chain_id: str, resno: int, ins_code: str, atom_name: str, alt_conf: str, attribute_name: str, attribute_value: str) -> int:
        """ set a string attibute to the atom with the given specifier. Attributes can be "atom-name", "alt-conf", "element" or "segid". """
        return 0

def set_residue_name(imol: int, chain_id: str, res_no: int, ins_code: str, new_residue_name: str) -> None:
        """ set lots of atom attributes at once by-passing the rebonding and redrawing of the above 2 functions set the residue name of the specified residue """

# --- Skeletonization Interface ---

def skel_greer_on() -> None:
        """ Sphinx-Doc-Placeholder"""

def skel_greer_off() -> None:
        """ Sphinx-Doc-Placeholder"""

def skeletonize_map(imol: int, prune_flag: int) -> int:
        """ skeletonize molecule number imol the prune_flag should almost always be 0.

        NOTE:: The arguments to have been reversed for coot 0.8.3 and later (now the molecule number comes first). """
        return 0

def unskeletonize_map(imol: int) -> int:
        """ undisplay the skeleton on molecule number imol """
        return 0

def set_initial_map_for_skeletonize() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_max_skeleton_search_depth(v: int) -> None:
        """ set the skeleton search depth, used in baton building For high resolution maps, you need to search deeper down the skeleton tree. This limit needs to be increased to 20 or so for high res maps (it is 10 by default) """

def get_text_for_skeletonization_level_entry():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_skeletonization_level_from_widget(txt: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_text_for_skeleton_box_size_entry():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_skeleton_box_size_from_widget(txt: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_skeleton_box_size(f: float) -> None:
        """ the box size (in Angstroms) for which the skeleton is displayed """

# --- Save Coordinates ---

def save_coordinates(imol: int, filename: str) -> int:
        """ save coordinates of molecule number imol in filename 

        :return: status 1 is good (success), 0 is fail. """
        return 0

def set_save_coordinates_in_original_directory(i: int) -> None:
        """ set save coordinates in the starting directory """

def save_molecule_number_from_option_menu() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_save_molecule_number(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Read Phases File Functions ---

def read_phs_and_coords_and_make_map(pdb_filename: str) -> None:
        """ read phs file use coords to get cell and symm to make map uses pending data to make the map. """

def read_phs_and_make_map_using_cell_symm_from_previous_mol(phs_filename: str) -> int:
        """ read a phs file, the cell and symm information is from previously read (most recently read) coordinates file For use with phs data filename provided on the command line """
        return 0

def read_phs_and_make_map_using_cell_symm_from_mol(phs_filename: str, imol: int) -> int:
        """ read phs file and use a previously read molecule to provide the cell and symmetry information 

        :return: the new molecule number, return -1 if problem creating the map (e.g. not phs data, file not found etc). """
        return 0

def read_phs_and_make_map_using_cell_symm_from_mol_using_implicit_phs_filename(imol: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_phs_and_make_map_using_cell_symm(phs_file_name: str, hm_spacegroup: str, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> int:
        """ read phs file use coords to use cell and symm to make map in degrees """
        return 0

def read_phs_and_make_map_with_reso_limits(imol: int, phs_file_name: str, reso_lim_low: float, reso_lim_high: float) -> int:
        """ read a phs file and use the cell and symm in molecule number imol and use the resolution limits reso_lim_high (in Angstroems). 

        :param imol:  is the molecule number of the reference (coordinates) molecule from which the cell and symmetry can be obtained.

        :param phs_file_name:  is the name of the phs data file.

        :param reso_lim_high:  is the high resolution limit in Angstroems.

        :param reso_lim_low:  the low resoluion limit (currently ignored). """
        return 0

def graphics_store_phs_filename(phs_filename: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def possible_cell_symm_for_phs_file():
        """ Sphinx-Doc-Placeholder"""
        pass

def get_text_for_phs_cell_chooser(imol: int, field: str):
        """ Sphinx-Doc-Placeholder"""
        pass

# --- Graphics Move ---

def undo_last_move() -> None:
        """ undo last move """

def translate_molecule_by(imol: int, x: float, y: float, z: float) -> None:
        """ translate molecule number imol by (x,y,z) in Angstroms """

def transform_molecule_by(imol: int, m11: float, m12: float, m13: float, m21: float, m22: float, m23: float, m31: float, m32: float, m33: float, x: float, y: float, z: float) -> None:
        """ transform molecule number imol by the given rotation matrix, then translate by (x,y,z) in Angstroms """

def transform_zone(imol: int, chain_id: str, resno_start: int, resno_end: int, ins_code: str, m11: float, m12: float, m13: float, m21: float, m22: float, m23: float, m31: float, m32: float, m33: float, x: float, y: float, z: float) -> None:
        """ transform fragment of molecule number imol by the given rotation matrix, then translate by (x,y,z) in Angstroms """

# --- Go To Atom Widget Functions ---

def post_go_to_atom_window() -> None:
        """ Post the Go To Atom Window. """

def go_to_atom_molecule_number() -> int:
        """ the go-to-atom molecule number """
        return 0

def go_to_atom_chain_id():
        """ the go-to-atom chain-id """
        pass

def go_to_atom_atom_name():
        """ the go-to-atom atom name """
        pass

def go_to_atom_residue_number() -> int:
        """ the go-to-atom residue number """
        return 0

def go_to_atom_ins_code():
        """ the go-to-atom insertion code """
        pass

def go_to_atom_alt_conf():
        """ the go-to-atom alt conf """
        pass

def set_go_to_atom_chain_residue_atom_name(t1_chain_id: str, iresno: int, t3_atom_name: str) -> int:
        """ set the go to atom specification It seems important for swig that the char * arguments are const char *, not const gchar * (or else we get wrong type of argument error on (say) "A"

        :return: the success status of the go to. 0 for fail, 1 for success. """
        return 0

def set_go_to_atom_chain_residue_atom_name_full(chain_id: str, resno: int, ins_code: str, atom_name: str, alt_conf: str) -> int:
        """ set the go to (full) atom specification It seems important for swig that the char * arguments are const char *, not const gchar * (or else we get wrong type of argument error on (say) "A"

        :return: the success status of the go to. 0 for fail, 1 for success. """
        return 0

def set_go_to_atom_chain_residue_atom_name_no_redraw(t1: str, iresno: int, t3: str, make_the_move_flag: int) -> int:
        """ set go to atom but don't redraw """
        return 0

def set_go_to_atom_chain_residue_atom_name_strings(t1: str, t2: str, txt: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def update_go_to_atom_from_current_position() -> None:
        """ update the Go To Atom widget entries to atom closest to screen centre. """

def update_go_to_atom_residue_list(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def atom_spec_to_atom_index(mol: int, chain: str, resno: int, atom_name: str) -> int:
        """ what is the atom index of the given atom? """
        return 0

def full_atom_spec_to_atom_index(imol: int, chain: str, resno: int, inscode: str, atom_name: str, altloc: str) -> int:
        """ what is the atom index of the given atom? """
        return 0

def update_go_to_atom_window_on_changed_mol(imol: int) -> None:
        """ update the Go To Atom window """

def update_go_to_atom_window_on_new_mol() -> None:
        """ update the Go To Atom window. This updates the option menu for the molecules. """

def update_go_to_atom_window_on_other_molecule_chosen(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_go_to_atom_molecule(imol: int) -> None:
        """ set the molecule for the Go To Atom For dynarama callback sake. The widget/class knows which molecule that it was generated from, so in order to go to the molecule from dynarama, we first need to the the molecule - because 

        Also used in scripting, where go-to-atom-chain-residue-atom-name does not mention the molecule number.

        20090914-PE set-go-to-atom-molecule can be used in a script and it should change the go-to-atom-molecule in the Go To Atom dialog (if it is being displayed). This does mean, of course that using the ramachandran plot to centre on atoms will change the Go To Atom dialog. Maybe that is surprising (maybe not). """

def unset_go_to_atom_widget() -> None:
        """ Sphinx-Doc-Placeholder"""

def autobuild_ca_off() -> None:
        """ Sphinx-Doc-Placeholder"""

def test_fragment() -> None:
        """ Sphinx-Doc-Placeholder"""

def do_skeleton_prune() -> None:
        """ Sphinx-Doc-Placeholder"""

def test_function(i: int, j: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def glyco_tree_test() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Map and Molecule Control ---

def post_display_control_window() -> None:
        """ display the Display Constrol window """

def add_map_display_control_widgets() -> None:
        """ Sphinx-Doc-Placeholder"""

def add_mol_display_control_widgets() -> None:
        """ Sphinx-Doc-Placeholder"""

def add_map_and_mol_display_control_widgets() -> None:
        """ Sphinx-Doc-Placeholder"""

def reset_graphics_display_control_window() -> None:
        """ Sphinx-Doc-Placeholder"""

def close_graphics_display_control_window() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_map_displayed(imol: int, state: int) -> None:
        """ make the map displayed/undisplayed, 0 for off, 1 for on """

def set_mol_displayed(imol: int, state: int) -> None:
        """ make the coordinates molecule displayed/undisplayed, 0 for off, 1 for on """

def set_display_only_model_mol(imol: int) -> None:
        """ from all the model molecules, display only imol This stops flashing/delayed animations with many molecules """

def set_mol_active(imol: int, state: int) -> None:
        """ make the coordinates molecule active/inactve (clickable), 0 for off, 1 for on """

def mol_is_displayed(imol: int) -> int:
        """ return the display state of molecule number imol 

        :return: 1 for on, 0 for off """
        return 0

def mol_is_active(imol: int) -> int:
        """ return the active state of molecule number imol 

        :return: 1 for on, 0 for off """
        return 0

def map_is_displayed(imol: int) -> int:
        """ return the display state of molecule number imol 

        :return: 1 for on, 0 for off """
        return 0

def set_all_maps_displayed(on_or_off: int) -> None:
        """ if on_or_off is 0 turn off all maps displayed, for other values of on_or_off turn on all maps """

def set_all_models_displayed_and_active(on_or_off: int) -> None:
        """ if on_or_off is 0 turn off all models displayed and active, for other values of on_or_off turn on all models. """

def set_only_last_model_molecule_displayed() -> None:
        """ only display the last model molecule """

def display_only_active() -> None:
        """ display only the active mol and the refinement map """

def show_spacegroup(imol: int):
        """ return the spacegroup of molecule number imol . Deprecated. 

        :return: "No Spacegroup" when the spacegroup of a molecule has not been set. """
        pass

# --- Align and Mutate ---

def align_and_mutate(imol: int, chain_id: str, fasta_maybe: str, renumber_residues_flag: int) -> None:
        """ align and mutate the given chain to the given sequence """

def set_alignment_gap_and_space_penalty(wgap: float, wspace: float) -> None:
        """ set the penalty for affine gap and space when aligning, defaults -3.0 and -0.4 """

# --- Renumber Residue Range ---

def renumber_residue_range(imol: int, chain_id: str, start_res: int, last_res: int, offset: int) -> int:
        """ renumber the given residue range by offset residues """
        return 0

def change_residue_number(imol: int, chain_id: str, current_resno: int, current_inscode: str, new_resno: int, new_inscode: str) -> int:
        """ change chain id, residue number or insertion code for given residue """
        return 0

# --- Change Chain ID ---

def change_chain_id(imol: int, from_chain_id: str, to_chain_id: str, use_res_range_flag: int, from_resno: int, to_resno: int) -> None:
        """ change the chain id of the specified residue """

# --- Scripting Interface ---

def probe_available_p() -> int:
        """ Can we run probe (was the executable variable set properly?) (predicate). 

        :return: 1 for yes, 2 for no """
        return 0

def probe_available_p_py() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def post_scripting_window() -> None:
        """ do nothing - compatibility function """

def post_scheme_scripting_window() -> None:
        """ pop-up a scripting window for scheming """

def run_command_line_scripts() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_guile_gui_loaded_flag() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_python_gui_loaded_flag() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_found_coot_gui() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_found_coot_python_gui() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Monomer ---

def get_monomer_for_molecule_by_index(dict_idx: int, imol_enc: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def run_script(filename: str) -> None:
        """ run script file """

def run_guile_script(filename: str) -> None:
        """ guile run script file """

def run_python_script(filename: str) -> None:
        """ run python script file """

def import_python_module(module_name: str, use_namespace: int) -> int:
        """ import python module """
        return 0

# --- Regularization and Refinement ---

def do_regularize(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def do_refine(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def add_planar_peptide_restraints() -> None:
        """ add a restraint on peptides to make them planar This adds a 5 atom restraint that includes both CA atoms of the peptide. Use this rather than editting the mon_lib_list.cif file. """

def remove_planar_peptide_restraints() -> None:
        """ remove restraints on peptides to make them planar. """

def make_tight_planar_peptide_restraints() -> None:
        """ make the planar peptide restraints tight Useful when refining models with cryo-EM maps """

def planar_peptide_restraints_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_use_trans_peptide_restraints(on_off_state: int) -> None:
        """ add a restraint on peptides to keep trans peptides trans i.e. omega in trans-peptides is restraints to 180 degrees. """

def add_omega_torsion_restriants() -> None:
        """ add restraints on the omega angle of the peptides (that is the torsion round the peptide bond). Omega angles that are closer to 0 than to 180 will be refined as cis peptides (and of course if omega is greater than 90 then the peptide will be refined as a trans peptide (this is the normal case). """

def remove_omega_torsion_restriants() -> None:
        """ remove omega restraints on CIS and TRANS linked residues. """

def set_refine_hydrogen_bonds(state: int) -> None:
        """ add or remove auto H-bond restraints """

def set_refinement_immediate_replacement(istate: int) -> None:
        """ set immediate replacement mode for refinement and regularization. You need this (call with istate=1) if you are scripting refinement/regularization """

def refinement_immediate_replacement_state() -> int:
        """ query the state of the immediate replacement mode """
        return 0

def set_refine_use_noughties_physics(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_refine_use_noughties_physics_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_residue_selection_flash_frames_number(i: int) -> None:
        """ set the number of frames for which the selected residue range flashes On fast computers, this can be set to higher than the default for more aesthetic appeal. """

def c_accept_moving_atoms() -> None:
        """ accept the new positions of the regularized or refined residues If you are scripting refinement and/or regularization, this is the function that you need to call after refine-zone or regularize-zone. """

def accept_regularizement() -> None:
        """ a hideous alias for the above """

def clear_up_moving_atoms() -> None:
        """ Sphinx-Doc-Placeholder"""

def clear_moving_atoms_object() -> None:
        """ Sphinx-Doc-Placeholder"""

def stop_refinement_internal() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_refinement_use_soft_mode_nbc_restraints(flag: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def shiftfield_b_factor_refinement(imol: int) -> None:
        """ shiftfield B-factor refinement """

def shiftfield_xyz_factor_refinement(imol: int) -> None:
        """ shiftfield xyz refinement """

def set_refine_with_torsion_restraints(istate: int) -> None:
        """ turn on (or off) torsion restraints Pass with istate=1 for on, istate=0 for off. """

def refine_with_torsion_restraints_state() -> int:
        """ return the state of above """
        return 0

def set_matrix(f: float) -> None:
        """ set the relative weight of the geometric terms to the map terms The default is 60.

        The higher the number the more weight that is given to the map terms but the resulting chi squared values are higher). This will be needed for maps generated from data not on (or close to) the absolute scale or maps that have been scaled (for example so that the sigma level has been scaled to 1.0). """

def matrix_state() -> float:
        """ return the relative weight of the geometric terms to the map terms. """
        return 0.0

def get_map_weight() -> float:
        """ return the relative weight of the geometric terms to the map terms. A more sensible name for the """
        return 0.0

def estimate_map_weight(imol_map: int) -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def set_refine_auto_range_step(i: int) -> None:
        """ change the +/- step for autoranging (default is 1) Auto-ranging alow you to select a range from one button press, this allows you to set the number of residues either side of the clicked residue that becomes the selected zone """

def set_refine_max_residues(n: int) -> None:
        """ set the heuristic fencepost for the maximum number of residues in the refinement/regularization residue range Default is 20 """

def refine_zone_atom_index_define(imol: int, ind1: int, ind2: int) -> None:
        """ refine a zone based on atom indexing """

def refine_zone(imol: int, chain_id: str, resno1: int, resno2: int, altconf: str) -> None:
        """ refine a zone presumes that imol_Refinement_Map has been set """

def repeat_refine_zone() -> None:
        """ repeat the previous (user-selected) refine zone """

def refine_auto_range(imol: int, chain_id: str, resno1: int, altconf: str) -> None:
        """ refine a zone using auto-range presumes that imol_Refinement_Map has been set """

def regularize_zone(imol: int, chain_id: str, resno1: int, resno2: int, altconf: str) -> int:
        """ regularize a zone 

        :return: a status, whether the regularisation was done or not. 0 for no, 1 for yes. """
        return 0

def set_dragged_refinement_steps_per_frame(v: int) -> None:
        """ set the number of refinement steps applied to the intermediate atoms each frame of graphics. smaller numbers make the movement of the intermediate atoms slower, smoother, more elegant.

        Default: 80. """

def dragged_refinement_steps_per_frame() -> int:
        """ return the number of steps per frame in dragged refinement """
        return 0

def set_refinement_refine_per_frame(istate: int) -> None:
        """ allow refinement of intermediate atoms after dragging, before displaying (default: 0, off). An attempt to do something like xfit does, at the request of Frank von Delft.

        Pass with istate=1 to enable this option. """

def refinement_refine_per_frame_state() -> int:
        """ query the state of the above option """
        return 0

def set_refinement_drag_elasticity(e: float) -> None:
        """ Default 0.33

        Bigger numbers mean bigger movement of the other atoms. """

def set_refine_ramachandran_angles(state: int) -> None:
        """ turn on Ramachandran angles refinement in refinement and regularization name consistent with """

def set_refine_ramachandran_torsion_angles(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_refine_ramachandran_restraints_type(type: int) -> None:
        """ change the target function type """

def set_refine_ramachandran_restraints_weight(w: float) -> None:
        """ change the target function weight a big number means bad things """

def refine_ramachandran_restraints_weight() -> float:
        """ ramachandran restraints weight

        :return: weight as a float """
        return 0.0

def set_torsion_restraints_weight(w: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_refine_rotamers(state: int) -> None:
        """ set the state for using rotamer restraints "drive" mode 1 in on, 0 is off (off by default) """

def set_refinement_geman_mcclure_alpha_from_text(combobox_item_idx: int, t: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_refinement_lennard_jones_epsilon_from_text(combobox_item_idx: int, t: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_refinement_ramachandran_restraints_weight_from_text(combobox_item_idx: int, t: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_refinement_overall_weight_from_text(t: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_refinement_torsion_weight_from_text(combobox_item_index: int, t: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_refine_params_dialog_more_control_frame_is_active(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def refine_ramachandran_angles_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_numerical_gradients(istate: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_debug_refinement(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_fix_chiral_volumes_before_refinement(istate: int) -> None:
        """ correct the sign of chiral volumes before commencing refinement? Do we want to fix chiral volumes (by moving the chiral atom to the other side of the chiral plane if necessary). Default yes (1). Note: doesn't work currently. """

def check_chiral_volumes(imol: int) -> None:
        """ query the state of the above option """

def set_show_chiral_volume_errors_dialog(istate: int) -> None:
        """ For experienced Cooters who don't like Coot nannying about chiral volumes during refinement. """

def set_secondary_structure_restraints_type(itype: int) -> None:
        """ set the type of secondary structure restraints 0 no sec str restraints

        1 alpha helix restraints

        2 beta strand restraints """

def secondary_structure_restraints_type() -> int:
        """ return the secondary structure restraints type """
        return 0

def imol_refinement_map() -> int:
        """ the molecule number of the map used for refinement 

        :return: the map number, if it has been set or there is only one map, return -1 on no map set (ambiguous) or no maps. """
        return 0

def set_imol_refinement_map(imol: int) -> int:
        """ set the molecule number of the map to be used for refinement/fitting. 

        :return: imol on success, -1 on failure """
        return 0

def does_residue_exist_p(imol: int, chain_id: str, resno: int, inscode: str) -> int:
        """ Does the residue exist? (Raw function) 

        :return: 0 on not-exist, 1 on does exist. """
        return 0

def delete_restraints(comp_id: str) -> int:
        """ delete the restraints for the given comp_id (i.e. residue name) 

        :return: success status (0 is failed, 1 is success) """
        return 0

def add_extra_bond_restraint(imol: int, chain_id_1: str, res_no_1: int, ins_code_1: str, atom_name_1: str, alt_conf_1: str, chain_id_2: str, res_no_2: int, ins_code_2: str, atom_name_2: str, alt_conf_2: str, bond_dist: float, esd: float) -> int:
        """ add a user-define bond restraint this extra restraint is used when the given atoms are selected in refinement or regularization.

        :return: the index of the new restraint.

        :return: -1 when the atoms were not found and no extra bond restraint was stored. """
        return 0

def add_extra_geman_mcclure_restraint(imol: int, chain_id_1: str, res_no_1: int, ins_code_1: str, atom_name_1: str, alt_conf_1: str, chain_id_2: str, res_no_2: int, ins_code_2: str, atom_name_2: str, alt_conf_2: str, bond_dist: float, esd: float) -> int:
        """ add a user-define GM distance restraint this extra restraint is used when the given atoms are selected in refinement or regularization.

        :return: the index of the new restraint.

        :return: -1 when the atoms were not found and no extra bond restraint was stored. """
        return 0

def set_show_extra_distance_restraints(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def add_extra_angle_restraint(imol: int, chain_id_1: str, res_no_1: int, ins_code_1: str, atom_name_1: str, alt_conf_1: str, chain_id_2: str, res_no_2: int, ins_code_2: str, atom_name_2: str, alt_conf_2: str, chain_id_3: str, res_no_3: int, ins_code_3: str, atom_name_3: str, alt_conf_3: str, torsion_angle: float, esd: float) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def add_extra_torsion_restraint(imol: int, chain_id_1: str, res_no_1: int, ins_code_1: str, atom_name_1: str, alt_conf_1: str, chain_id_2: str, res_no_2: int, ins_code_2: str, atom_name_2: str, alt_conf_2: str, chain_id_3: str, res_no_3: int, ins_code_3: str, atom_name_3: str, alt_conf_3: str, chain_id_4: str, res_no_4: int, ins_code_4: str, atom_name_4: str, alt_conf_4: str, torsion_angle: float, esd: float, period: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def add_extra_start_pos_restraint(imol: int, chain_id_1: str, res_no_1: int, ins_code_1: str, atom_name_1: str, alt_conf_1: str, esd: float) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def add_extra_target_position_restraint(imol: int, chain_id: str, res_no: int, ins_code: str, atom_name: str, alt_conf: str, x: float, y: float, z: float, weight: float) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def delete_all_extra_restraints(imol: int) -> None:
        """ clear out all the extra/user-defined restraints for molecule number imol """

def delete_extra_restraints_for_residue(imol: int, chain_id: str, res_no: int, ins_code: str) -> None:
        """ clear out all the extra/user-defined restraints for this residue in molecule number imol """

def delete_extra_restraints_worse_than(imol: int, n_sigma: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def add_refmac_extra_restraints(imol: int, file_name: str) -> None:
        """ read in prosmart (typically) extra restraints """

def set_show_extra_restraints(imol: int, state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def extra_restraints_are_shown(imol: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_extra_restraints_prosmart_sigma_limits(imol: int, limit_high: float, limit_low: float) -> None:
        """ often we don't want to see all prosmart restraints, just the (big) violations """

def generate_local_self_restraints(imol: int, chain_id: str, local_dist_max: float) -> None:
        """ generate external distance local self restraints """

def generate_self_restraints(imol: int, local_dist_max: float) -> None:
        """ generate external distance all-molecule self restraints """

def write_interpolated_extra_restraints(imol_1: int, imol_2: int, n_steps: int, file_name_stub: str) -> None:
        """ proSMART interpolated restraints for model morphing """

def write_interpolated_models_and_extra_restraints(imol_1: int, imol_2: int, n_steps: int, file_name_stub: str, interpolation_mode: int) -> None:
        """ proSMART interpolated restraints for model morphing and write interpolated model interpolation_mode is currently dummy - in due course I will addd torion angle interpolation. """

def set_show_parallel_plane_restraints(imol: int, state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def parallel_plane_restraints_are_shown(imol: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def add_parallel_plane_restraint(imol: int, chain_id_1: str, re_no_1: int, ins_code_1: str, chain_id_2: str, re_no_2: int, ins_code_2: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_extra_restraints_representation_for_bonds_go_to_CA(imol: int, state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_use_only_extra_torsion_restraints_for_torsions(state: int) -> None:
        """ set use only extra torsion restraints for torsions """

def use_only_extra_torsion_restraints_for_torsions_state() -> int:
        """ return only-use-extra-torsion-restraints-for-torsions state """
        return 0

def clear_all_atom_pull_restraints() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_auto_clear_atom_pull_restraint(state: int) -> None:
        """ set auto-clear atom pull restraint """

def get_auto_clear_atom_pull_restraint_state() -> int:
        """ get auto-clear atom pull restraint state """
        return 0

def increase_proportional_editing_radius() -> None:
        """ iscrease the proportional editing radius """

def decrease_proportional_editing_radius() -> None:
        """ descrease the proportional editing radius """

# --- Simplex Refinement Interface ---

def fit_residue_range_to_map_by_simplex(res1: int, res2: int, altloc: str, chain_id: str, imol: int, imol_for_map: int) -> None:
        """ refine residue range using simplex optimization """

def score_residue_range_fit_to_map(res1: int, res2: int, altloc: str, chain_id: str, imol: int, imol_for_map: int) -> float:
        """ simply score the residue range fit to map """
        return 0.0

# --- Nomenclature Errors ---

def fix_nomenclature_errors(imol: int) -> int:
        """ fix nomenclature errors in molecule number imol 

        :return: the number of resides altered. """
        return 0

def set_nomenclature_errors_on_read(mode: str) -> None:
        """ set way nomenclature errors should be handled on reading coordinates. mode should be "auto-correct", "ignore", "prompt". The default is "prompt" """

# --- Atom Info  Interface ---

def output_atom_info_as_text(imol: int, chain_id: str, resno: int, ins_code: str, atname: str, altconf: str) -> None:
        """ output to the terminal the Atom Info for the give atom specs """

# --- Residue Info ---

def do_residue_info_dialog() -> None:
        """ Sphinx-Doc-Placeholder"""

def output_residue_info_dialog(imol: int, atom_index: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def residue_info_dialog(imol: int, chain_id: str, resno: int, ins_code: str) -> None:
        """ show residue info dialog for given residue """

def residue_info_dialog_is_displayed() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def output_residue_info_as_text(atom_index: int, imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def do_distance_define() -> None:
        """ Sphinx-Doc-Placeholder"""

def do_angle_define() -> None:
        """ Sphinx-Doc-Placeholder"""

def do_torsion_define() -> None:
        """ Sphinx-Doc-Placeholder"""

def residue_info_apply_all_checkbutton_toggled() -> None:
        """ Sphinx-Doc-Placeholder"""

def clear_residue_info_edit_list() -> None:
        """ Sphinx-Doc-Placeholder"""

def unset_residue_info_widget() -> None:
        """ Sphinx-Doc-Placeholder"""

def clear_measure_distances() -> None:
        """ Sphinx-Doc-Placeholder"""

def clear_last_measure_distance() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Edit Fuctions ---

def do_edit_copy_molecule() -> None:
        """ Sphinx-Doc-Placeholder"""

def do_edit_copy_fragment() -> None:
        """ Sphinx-Doc-Placeholder"""

def do_edit_replace_fragment() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Residue Environment Functions ---

def set_show_environment_distances(state: int) -> None:
        """ show environment distances. If state is 0, distances are turned off, otherwise distances are turned on. """

def set_show_environment_distances_bumps(state: int) -> None:
        """ show bumps environment distances. If state is 0, bump distances are turned off, otherwise bump distances are turned on. """

def set_show_environment_distances_h_bonds(state: int) -> None:
        """ show H-bond environment distances. If state is 0, bump distances are turned off, otherwise H-bond distances are turned on. """

def show_environment_distances_state() -> int:
        """ show the state of display of the environment distances """
        return 0

def set_environment_distances_distance_limits(min_dist: float, max_dist: float) -> None:
        """ min and max distances for the environment distances """

def set_show_environment_distances_as_solid(state: int) -> None:
        """ show the environment distances with solid modelling """

def set_environment_distances_label_atom(state: int) -> None:
        """ Label the atom on Environment Distances start/change. """

def label_neighbours() -> None:
        """ Label the atoms in the residues around the central residue. """

def label_atoms_in_residue() -> None:
        """ Label the atoms in the central residue. """

def set_show_local_b_factors(state: int) -> None:
        """ Label the atoms with their B-factors. """

def add_geometry_distance(imol_1: int, x_1: float, y_1: float, z_1: float, imol_2: int, x_2: float, y_2: float, z_2: float):
        """ Add a geometry distance between points in a given molecule. 

        :return: the distance between the points """
        pass

# --- Pointer Functions ---

def set_show_pointer_distances(istate: int) -> None:
        """ turn on (or off) the pointer distance by passing 1 (or 0). """

def show_pointer_distances_state() -> int:
        """ show the state of display of the pointer distances """
        return 0

# --- Zoom Functions ---

def scale_zoom(f: float) -> None:
        """ scale the view by f external (scripting) interface (with redraw) 

        :param f:  the smaller f, the bigger the zoom, typical value 1.3. Values outside the range 0.5 to 1.8 are filtered out """

def scale_zoom_internal(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def zoom_factor() -> float:
        """ return the current zoom factor i.e. get_zoom_factor() """
        return 0.0

def set_smooth_scroll_do_zoom(i: int) -> None:
        """ set smooth scroll with zoom 

        :param i:  0 means no, 1 means yes: (default 0) """

def smooth_scroll_do_zoom() -> int:
        """ return the state of the above system """
        return 0

def smooth_scroll_zoom_limit() -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def set_smooth_scroll_zoom_limit(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_zoom(f: float) -> None:
        """ set the zoom factor (absolute value) - maybe should be called set_zoom_factor() """

# --- CNS Data Functions ---

def handle_cns_data_file(filename: str, imol: int) -> int:
        """ read CNS data (currently only a placeholder) """
        return 0

def handle_cns_data_file_with_cell(filename: str, imol: int, a: float, b: float, c: float, alpha: float, beta: float, gamma: float, spg_info: str) -> int:
        """ read CNS data (currently only a placeholder) a, b,c are in Angstroems. alpha, beta, gamma are in degrees. spg is the space group info, either ;-delimited symmetry operators or the space group name """
        return 0

# --- mmCIF Functions ---

def auto_read_cif_data_with_phases(filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data_with_phases_sigmaa(filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data_with_phases_diff_sigmaa(filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data(filename: str, imol_coords: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data_2fofc_map(filename: str, imol_coords: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data_fofc_map(filename: str, imol_coords: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data_with_phases_fo_fc(filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data_with_phases_2fo_fc(filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data_with_phases_nfo_fc(filename: str, map_type: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_data_with_phases_fo_alpha_calc(filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def write_connectivity(monomer_name: str, filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def open_cif_dictionary_file_selector_dialog() -> None:
        """ open the cif dictionary file selector dialog """

def import_all_refmac_cifs() -> None:
        """ Sphinx-Doc-Placeholder"""

def read_small_molecule_cif(file_name: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_small_molecule_data_cif(file_name: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_small_molecule_data_cif_and_make_map_using_coords(file_name: str, imol_coords: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

# --- Validation Functions ---

def deviant_geometry(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def is_valid_model_molecule(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def is_valid_map_molecule(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def difference_map_peaks(imol: int, imol_coords: int, level: float, max_closeness: float, do_positive_level_flag: int, do_negative_level_flag: int, around_model_only_flag: int) -> None:
        """ generate a list of difference map peaks peaks within max_closeness (2.0 A typically) of a larger peak are not listed.

        the flag around_model_only_flag limits the peak list to those only within 4A of the selected model (useful for maps with molecular symmetry). """

def set_difference_map_peaks_max_closeness(m: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def difference_map_peaks_max_closeness() -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def clear_diff_map_peaks() -> None:
        """ Sphinx-Doc-Placeholder"""

def gln_asn_b_factor_outliers(imol: int) -> None:
        """ Make a gui for GLN adn ASN B-factor outiers, compairing the O and N temperatur factors difference to the distribution of temperature factors from the other atoms. """

def gln_asn_b_factor_outliers_py(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def clear_multi_residue_torsion_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_multi_residue_torsion_reverse_mode(mode: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def show_multi_residue_torsion_dialog() -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_multi_residue_torsion() -> None:
        """ Sphinx-Doc-Placeholder"""

def atom_overlap_score(imol: int) -> float:
        """ return the atom overlap score """
        return 0.0

def set_show_chiral_volume_outliers(imol: int, state: int) -> None:
        """ set the state of showing chiral volume outlier markers - of a model molecule that is, not the intermediate atoms (derived from restraints) set the state of showing non-bonded contact markers - of a model molecule that is, not the intermediate atoms (derived from restraints) """

# --- Ramachandran Plot Functions ---

def do_ramachandran_plot(imol: int) -> None:
        """ Ramachandran plot for molecule number imol. """

def set_kleywegt_plot_n_diffs(n_diffs: int) -> None:
        """ set the number of biggest difference arrows on the Kleywegt plot. """

def set_ramachandran_plot_contour_levels(level_prefered: float, level_allowed: float) -> None:
        """ set the contour levels for the ramachandran plot, default values are 0.02 (prefered) 0.002 (allowed) """

def set_ramachandran_plot_background_block_size(blocksize: float) -> None:
        """ set the ramachandran plot background block size. Smaller is smoother but slower. Should be divisible exactly into"""

def set_ramachandran_psi_axis_mode(mode: int) -> None:
        """ set the psi axis for the ramachandran plot. Default (0) from -180 to 180. Alternative (1) from -120 to 240. """

def ramachandran_psi_axis_mode() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_moving_atoms(phi: float, psi: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def accept_phi_psi_moving_atoms() -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_edit_phi_psi(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_dynamic_distances(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def destroy_edit_backbone_rama_plot() -> None:
        """ Sphinx-Doc-Placeholder"""

def ramachandran_plot_differences(imol1: int, imol2: int) -> None:
        """ 2 molecule ramachandran plot (NCS differences) a.k.a. A Kleywegt Plot. """

def ramachandran_plot_differences_by_chain(imol1: int, imol2: int, a_chain: str, b_chain: str) -> None:
        """ A chain-specific Kleywegt Plot. """

# --- Sequence View Interface ---

def sequence_view(imol: int) -> None:
        """ display the sequence view dialog for molecule number imol """

def do_sequence_view(imol: int) -> None:
        """ old name for the above function """

def update_sequence_view_current_position_highlight_from_active_atom() -> None:
        """ update the sequnce view current position highlight based on active atom """

def change_peptide_carbonyl_by(angle: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def change_peptide_peptide_by(angle: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def execute_setup_backbone_torsion_edit(imol: int, atom_index: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_backbone_torsion_edit(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_backbone_torsion_peptide_button_start_pos(ix: int, iy: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def change_peptide_peptide_by_current_button_pos(ix: int, iy: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_backbone_torsion_carbonyl_button_start_pos(ix: int, iy: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def change_peptide_carbonyl_by_current_button_pos(ix: int, iy: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Atom Labelling ---

def add_atom_label(imol: int, chain_id: str, iresno: int, atom_id: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def remove_atom_label(imol: int, chain_id: str, iresno: int, atom_id: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def remove_all_atom_labels() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_label_on_recentre_flag(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def centre_atom_label_status() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_brief_atom_labels(istat: int) -> None:
        """ use brief atom names for on-screen labels call with istat=1 to use brief labels, istat=0 for normal labels """

def brief_atom_labels_state() -> int:
        """ the brief atom label state """
        return 0

def set_seg_ids_in_atom_labels(istat: int) -> None:
        """ set if brief atom labels should have seg-ids also """

# --- Screen Rotation ---

def rotate_y_scene(nsteps: int, stepsize: float) -> None:
        """ rotate view round y axis stepsize degrees for nstep such steps """

def rotate_x_scene(nsteps: int, stepsize: float) -> None:
        """ rotate view round x axis stepsize degrees for nstep such steps """

def rotate_z_scene(nsteps: int, stepsize: float) -> None:
        """ rotate view round z axis stepsize degrees for nstep such steps """

def spin_zoom_trans(axis: int, nstep: int, stepsize: float, zoom_by: float, x_rel: float, y_rel: float, z_rel: float) -> None:
        """ Bells and whistles rotation. spin, zoom and translate.

        where axis is either x,y or z, stepsize is in degrees, zoom_by and x_rel etc are how much zoom, x,y,z should have changed by after nstep steps. """

# --- Screen Translation ---

def translate_scene_x(nsteps: int) -> None:
        """ translate rotation centre relative to screen axes for nsteps """

def translate_scene_y(nsteps: int) -> None:
        """ translate rotation centre relative to screen axes for nsteps """

def translate_scene_z(nsteps: int) -> None:
        """ translate rotation centre relative to screen axes for nsteps """

# --- Views Interface ---

def add_view_here(view_name: str) -> int:
        """ return the view number """
        return 0

def add_view_raw(rcx: float, rcy: float, rcz: float, quat1: float, quat2: float, quat3: float, quat4: float, zoom: float, view_name: str) -> int:
        """ return the view number """
        return 0

def play_views() -> None:
        """ Sphinx-Doc-Placeholder"""

def remove_this_view() -> None:
        """ Sphinx-Doc-Placeholder"""

def remove_named_view(view_name: str) -> int:
        """ the view with the given name """
        return 0

def remove_view(view_number: int) -> None:
        """ the given view number """

def go_to_first_view(snap_to_view_flag: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def go_to_view_number(view_number: int, snap_to_view_flag: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def add_spin_view(view_name: str, n_steps: int, degrees_total: float) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def add_view_description(view_number: int, description: str) -> None:
        """ Add a view description/annotation to the give view number. """

def add_action_view(view_name: str, action_function: str) -> int:
        """ add a view (not add to an existing view) that 

        :return: the view number for this (new) view. """
        return 0

def insert_action_view_after_view(view_number: int, view_name: str, action_function: str) -> int:
        """ add an action view after the view of the given view number 

        :return: the view number for this (new) view. """
        return 0

def n_views() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def save_views(view_file_name: str) -> None:
        """ save views to view_file_name """

def views_play_speed() -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def set_views_play_speed(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def clear_all_views() -> None:
        """ Clear the view list. """

# --- Movies Interface ---

def set_movie_file_name_prefix(file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_movie_frame_number(frame_number: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def movie_frame_number() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_make_movie_mode(make_movies_flag: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Background Colour ---

def set_background_colour(red: float, green: float, blue: float) -> None:
        """ set the background colour red, green and blue are numbers between 0.0 and 1.0 """

def redraw_background() -> None:
        """ re draw the background colour when switching between mono and stereo """

def background_is_black_p() -> int:
        """ is the background black (or nearly black)? 

        :return: 1 if the background is black (or nearly black), else return 0. """
        return 0

# --- Ligand Fitting Functions ---

def set_ligand_acceptable_fit_fraction(f: float) -> None:
        """ set the fraction of atoms which must be in positive density after a ligand fit """

def set_ligand_cluster_sigma_level(f: float) -> None:
        """ set the default sigma level that the map is searched to find potential ligand sites """

def set_ligand_flexible_ligand_n_samples(i: int) -> None:
        """ set the number of conformation samples big ligands require more samples. Default 10. """

def set_ligand_verbose_reporting(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_find_ligand_n_top_ligands(n: int) -> None:
        """ search the top n sites for ligands. Default 10. """

def set_find_ligand_do_real_space_refinement(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_find_ligand_multi_solutions_per_cluster(lim_1: float, lim_2: float) -> None:
        """ allow multiple ligand solutions per cluster. The first limit is the fraction of the top scored positions that go on to correlation scoring (closer to 1 means less and faster - default 0.7).

        The second limit is the fraction of the top correlation score that is considered interesting. Limits the number of solutions displayed to user. Default 0.9.

        There is currently no chi-angle set redundancy filtering - I suspect that there should be.

        Nino-mode. """

def set_find_ligand_mask_waters(istate: int) -> None:
        """ how shall we treat the waters during ligand fitting? pass with istate=1 for waters to mask the map in the same way that protein atoms do. """

def set_ligand_search_protein_molecule(imol: int) -> None:
        """ set the protein molecule for ligand searching """

def set_ligand_search_map_molecule(imol_map: int) -> None:
        """ set the map molecule for ligand searching """

def add_ligand_search_ligand_molecule(imol_ligand: int) -> None:
        """ add a rigid ligand molecule to the list of ligands to search for in ligand searching """

def add_ligand_search_wiggly_ligand_molecule(imol_ligand: int) -> None:
        """ add a flexible ligand molecule to the list of ligands to search for in ligand searching """

def set_find_ligand_here_cluster(state: int) -> None:
        """ Allow the user a scripting means to find ligand at the rotation centre. """

def execute_ligand_search() -> None:
        """ Sphinx-Doc-Placeholder"""

def add_ligand_clear_ligands() -> None:
        """ Sphinx-Doc-Placeholder"""

def ligand_expert() -> None:
        """ this sets the flag to have expert option ligand entries in the Ligand Searching dialog """

def do_find_ligands_dialog() -> None:
        """ display the find ligands dialog if maps, coords and ligands are available, that is. """

def match_ligand_atom_names(imol_ligand: int, chain_id_ligand: str, resno_ligand: int, ins_code_ligand: str, imol_reference: int, chain_id_reference: str, resno_reference: int, ins_code_reference: str) -> None:
        """ Overlap residue with "template"-based matching. Overlap the first residue in imol_ligand onto the residue specified by the reference parameters. Use graph matching, not atom names.

        :return: success status, False = failed to find residue in either imol_ligand or imo_ref. If success, return the RT operator.

        By using graph matching, make the names of the atoms of the given ligand/residue match those of the reference residue/ligand as closely as possible - where there would be an atom name clash, invent a new atom name. """

def match_ligand_atom_names_to_comp_id(imol_ligand: int, chain_id_ligand: str, resno_ligand: int, ins_code_ligand: str, comp_id_ref: str) -> None:
        """ Match ligand atom names to a reference ligand type (comp_id) By using graph matching, make the names of the atoms of the given ligand/residue match those of the reference ligand from the geometry store as closely as possible. Where there would be an atom name clash, invent a new atom name.

        This doesn't create a new dictionary for the selected ligand - and that's a big problem (see match_residue_and_dictionary). """

def exchange_ligand(imol_lig: int, chain_id_lig: str, resno_lig: int, ins_code_lig: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def flip_ligand(imol: int, chain_id: str, resno: int) -> None:
        """ flip the ligand (usually active residue) around its eigen vectors to the next flip number. Immediate replacement (like flip peptide). """

def jed_flip(imol: int, chain_id: str, res_no: int, ins_code: str, atom_name: str, alt_conf: str, invert_selection: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Water Fitting Functions ---

def show_create_find_waters_dialog() -> None:
        """ create a dialog for water fitting """

def renumber_waters(imol: int) -> None:
        """ Renumber the waters of molecule number imol with consecutive numbering. """

def execute_find_waters_real(imol_for_map: int, imol_for_protein: int, new_waters_mol_flag: int, rmsd_cut_off: float) -> None:
        """ find waters """

def find_waters(imol_for_map: int, imol_for_protein: int, new_waters_mol_flag: int, rmsd_cut_off: float, show_blobs_dialog: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def move_waters_to_around_protein(imol: int) -> int:
        """ move waters of molecule number imol so that they are around the protein. 

        :return: the number of moved waters. """
        return 0

def move_hetgroups_to_around_protein(imol: int) -> None:
        """ move all hetgroups (including waters) of molecule number imol so that they are around the protein. """

def max_water_distance(imol: int) -> float:
        """ return the maximum minimum distance of any water atom to any protein atom - used in validation of """
        return 0.0

def get_text_for_find_waters_sigma_cut_off():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_value_for_find_waters_sigma_cut_off(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_water_check_spherical_variance_limit(f: float) -> None:
        """ set the limit of interesting variance, above which waters are listed (otherwise ignored) default 0.12. """

def set_ligand_water_to_protein_distance_limits(f1: float, f2: float) -> None:
        """ set ligand to protein distance limits f1 is the minimum distance, f2 is the maximum distance """

def set_ligand_water_n_cycles(i: int) -> None:
        """ set the number of cycles of water searching """

def set_write_peaksearched_waters() -> None:
        """ Sphinx-Doc-Placeholder"""

def execute_find_blobs(imol_model: int, imol_for_map: int, cut_off: float, interactive_flag: int) -> None:
        """ find blobs """

def split_water(imol: int, chain_id: str, res_no: int, ins_code: str) -> None:
        """ split the given water and fit to map. If refinement map is not defined, don't do anything.

        If there is more than one atom in the specified resiue, don't do anything.

        If the given atom does not have an alt conf of "", don't do anything. """

# --- Bond Representation ---

def set_default_bond_thickness(t: int) -> None:
        """ set the default thickness for bonds (e.g. in ~/.coot) """

def set_bond_thickness(imol: int, t: float) -> None:
        """ set the thickness of the bonds in molecule number imol to t pixels """

def set_bond_thickness_intermediate_atoms(t: float) -> None:
        """ set the thickness of the bonds of the intermediate atoms to t pixels """

def set_use_variable_bond_thickness(state: int) -> None:
        """ allow lines that are further away to be thinner """

def set_bond_colour_rotation_for_molecule(imol: int, f: float) -> None:
        """ set bond colour for molecule """

def set_draw_stick_mode_atoms_default(state: int) -> None:
        """ set default for the drawing of atoms in stick mode (default is on (1)) """

def get_bond_colour_rotation_for_molecule(imol: int) -> float:
        """ get the bond colour for molecule. Return -1 on err (bad molecule number) """
        return 0.0

def set_unbonded_atom_star_size(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_default_representation_type(type: int) -> None:
        """ set the default representation type (default 1). """

def get_default_bond_thickness() -> int:
        """ get the default thickness for bonds """
        return 0

def set_draw_zero_occ_markers(status: int) -> None:
        """ set status of drawing zero occupancy markers. default status is 1. """

def set_draw_cis_peptide_markups(status: int) -> None:
        """ set status of drawing cis-peptide markups default status is 1. """

def set_draw_hydrogens(imol: int, istat: int) -> None:
        """ set the hydrogen drawing state. istat = 0 is hydrogens off, istat = 1: show hydrogens """

def draw_hydrogens_state(imol: int) -> int:
        """ the state of draw hydrogens for molecule number imol. return -1 on bad imol. """
        return 0

def set_draw_stick_mode_atoms(imol: int, state: int) -> None:
        """ draw little coloured balls on atoms turn off with state = 0

        turn on with state = 1 """

def set_draw_missing_residues_loops(state: int) -> None:
        """ set the state for drawing missing resiude loops For taking screenshots, we often don't want to see them. """

def graphics_to_ca_representation(imol: int) -> None:
        """ draw molecule number imol as CAs """

def graphics_to_colour_by_chain(imol: int) -> None:
        """ draw molecule number imol coloured by chain """

def graphics_to_ca_plus_ligands_representation(imol: int) -> None:
        """ draw molecule number imol as CA + ligands """

def graphics_to_ca_plus_ligands_and_sidechains_representation(imol: int) -> None:
        """ draw molecule number imol as CA + ligands + sidechains """

def graphics_to_bonds_no_waters_representation(imol: int) -> None:
        """ draw molecule number imol with no waters """

def graphics_to_bonds_representation(mol: int) -> None:
        """ draw molecule number imol with normal bonds """

def graphics_to_colour_by_molecule(imol: int) -> None:
        """ draw molecule with colour-by-molecule colours """

def graphics_to_ca_plus_ligands_sec_struct_representation(imol: int) -> None:
        """ draw molecule number imol with CA bonds in secondary structure representation and ligands """

def graphics_to_sec_struct_bonds_representation(imol: int) -> None:
        """ draw molecule number imol with bonds in secondary structure representation """

def graphics_to_rainbow_representation(imol: int) -> None:
        """ draw molecule number imol in Jones' Rainbow """

def graphics_to_b_factor_representation(imol: int) -> None:
        """ draw molecule number imol coloured by B-factor """

def graphics_to_b_factor_cas_representation(imol: int) -> None:
        """ draw molecule number imol coloured by B-factor, CA + ligands """

def graphics_to_occupancy_representation(imol: int) -> None:
        """ draw molecule number imol coloured by occupancy """

def graphics_to_user_defined_atom_colours_representation(imol: int) -> None:
        """ draw molecule number imol in CA+Ligands mode coloured by user-defined atom colours """

def graphics_to_user_defined_atom_colours_all_atoms_representation(imol: int) -> None:
        """ draw molecule number imol all atoms coloured by user-defined atom colours """

def get_graphics_molecule_bond_type(imol: int) -> int:
        """ what is the bond drawing state of molecule number imol """
        return 0

def set_b_factor_bonds_scale_factor(imol: int, f: float) -> int:
        """ scale the colours for colour by b factor representation """
        return 0

def change_model_molecule_representation_mode(up_or_down: int) -> None:
        """ change the representation of the model molecule closest to the centre of the screen """

def set_use_grey_carbons_for_molecule(imol: int, state: int) -> None:
        """ make the carbon atoms for molecule imol be grey """

def set_grey_carbon_colour(imol: int, r: float, g: float, b: float) -> None:
        """ set the colour for the carbon atoms can be not grey if you desire, r, g, b in the range 0 to 1. """

def set_draw_moving_atoms_restraints(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_draw_moving_atoms_restraints():
        """ Sphinx-Doc-Placeholder"""
        pass

def make_ball_and_stick(imol: int, atom_selection_str: str, bond_thickness: float, sphere_size: float, do_spheres_flag: int) -> int:
        """ make a ball and stick representation of imol given atom selection e.g. (make-ball-and-stick 0 "/1" 0.15 0.25 1) """
        return 0

def clear_ball_and_stick(imol: int) -> int:
        """ clear ball and stick representation of molecule number imol """
        return 0

def set_model_molecule_representation_style(imol: int, mode: int) -> None:
        """ set the model molecule representation stye 0 for ball-and-stick/licorice (default) and 1 for ball """

def set_show_molecular_representation(imol: int, mesh_index: int, state: int) -> None:
        """ set show a ribbon/mesh for a given molecule """

def set_show_additional_representation(imol: int, representation_number: int, on_off_flag: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_show_all_additional_representations(imol: int, on_off_flag: int) -> None:
        """ display/undisplay all the additional representations for the given molecule """

def all_additional_representations_off_except(imol: int, representation_number: int, ball_and_sticks_off_too_flag: int) -> None:
        """ removed from API brief undisplay all the additional representations for the given molecule, except the given representation number (if it is off, leave it off) """

def delete_additional_representation(imol: int, representation_number: int) -> None:
        """ removed from API brief delete a given additional representation """

def additional_representation_by_string(imol: int, atom_selection: str, representation_type: int, bonds_box_type: int, bond_width: float, draw_hydrogens_flag: int) -> int:
        """ removed from API brief return the index of the additional representation. Return -1 on error """
        return 0

def additional_representation_by_attributes(imol: int, chain_id: str, resno_start: int, resno_end: int, ins_code: str, representation_type: int, bonds_box_type: int, bond_width: float, draw_hydrogens_flag: int) -> int:
        """ return the index of the additional representation. 

        :return: -1 on error. """
        return 0

def set_flev_idle_ligand_interactions(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def toggle_flev_idle_ligand_interactions() -> None:
        """ Sphinx-Doc-Placeholder"""

def calculate_hydrogen_bonds(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_draw_hydrogen_bonds(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Dots Representation ---

def dots(imol: int, atom_selection_str: str, dots_object_name: str, dot_density: float, sphere_size_scale: float) -> int:
        """ display dotted surface return a generic objects handle (which can be used to remove later) """
        return 0

def set_dots_colour(imol: int, r: float, g: float, b: float) -> None:
        """ set the colour of the surface dots of the imol-th molecule to be the given single colour r,g,b are values between 0.0 and 1.0 """

def unset_dots_colour(imol: int) -> None:
        """ no longer set the dots of molecule imol to a single colour i.e. go back to element-based colours. """

def clear_dots(imol: int, dots_handle: int) -> None:
        """ clear dots in imol with dots_handle """

def clear_dots_by_name(imol: int, dots_object_name: str) -> None:
        """ clear the first dots object for imol with given name """

def n_dots_sets(imol: int) -> int:
        """ return the number of dots sets for molecule number imol """
        return 0

# --- Pep-flip Interface ---

def do_pepflip(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def pepflip(imol: int, chain_id: str, resno: int, inscode: str, altconf: str) -> None:
        """ pepflip the given residue """

def pepflip_intermediate_atoms() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def pepflip_intermediate_atoms_other_peptide() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

# --- Rigid Body Refinement Interface ---

def do_rigid_body_refine(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def rigid_body_refine_zone(imol: int, chain_id: str, reso_start: int, resno_end: int) -> None:
        """ setup rigid body refine zone where we set the atom selection holders according to the arguments and then call execute_rigid_body_refine() """

def rigid_body_refine_by_atom_selection(imol: int, atom_selection_string: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def execute_rigid_body_refine(auto_range_flag: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_rigid_body_fit_acceptable_fit_fraction(f: float) -> None:
        """ set rigid body fraction of atoms in positive density 

        :param f:  in the range 0.0 -> 1.0 (default 0.75) """

# --- Dynamic Map ---

def toggle_dynamic_map_display_size() -> None:
        """ Sphinx-Doc-Placeholder"""

def toggle_dynamic_map_sampling() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_dynamic_map_size_display_on() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_dynamic_map_size_display_off() -> None:
        """ Sphinx-Doc-Placeholder"""

def get_dynamic_map_size_display() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_dynamic_map_sampling_on() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_dynamic_map_sampling_off() -> None:
        """ Sphinx-Doc-Placeholder"""

def get_dynamic_map_sampling() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_dynamic_map_zoom_offset(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Add Terminal Residue Functions ---

def do_add_terminal_residue(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_add_terminal_residue_n_phi_psi_trials(n: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_add_terminal_residue_add_other_residue_flag(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_add_terminal_residue_do_rigid_body_refine(v: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_terminal_residue_do_rigid_body_refine(v: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_add_terminal_residue_debug_trials(debug_state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def add_terminal_residue_immediate_addition_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_add_terminal_residue_immediate_addition(i: int) -> None:
        """ set immediate addition of terminal residue call with i=1 for immediate addtion """

def add_terminal_residue(imol: int, chain_id: str, residue_number: int, residue_type: str, immediate_add: int) -> int:
        """ Add a terminal residue. residue type can be "auto" and immediate_add is recommended to be 1.

        :return: 0 on failure, 1 on success """
        return 0

def add_nucleotide(imol: int, chain_id: str, res_no: int) -> int:
        """ Add a terminal nucleotide. No fitting is done """
        return 0

def add_terminal_residue_using_phi_psi(imol: int, chain_id: str, res_no: int, residue_type: str, phi: float, psi: float) -> int:
        """ Add a terminal residue using given phi and psi angles. 

        :return: the success status, 0 on failure, 1 on success """
        return 0

def set_add_terminal_residue_default_residue_type(type: str) -> None:
        """ set the residue type of an added terminal residue. """

def set_add_terminal_residue_do_post_refine(istat: int) -> None:
        """ set a flag to run refine zone on terminal residues after an addition. """

def add_terminal_residue_do_post_refine_state() -> int:
        """ what is the value of the previous flag? """
        return 0

# --- Delete Residues ---

def delete_atom_by_atom_index(imol: int, index: int, do_delete_dialog: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def delete_residue_by_atom_index(imol: int, index: int, do_delete_dialog: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def delete_residue_hydrogens_by_atom_index(imol: int, index: int, do_delete_dialog: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def delete_residue_range(imol: int, chain_id: str, resno_start: int, end_resno: int) -> None:
        """ delete residue range """

def delete_residue(imol: int, chain_id: str, resno: int, inscode: str) -> None:
        """ delete residue """

def delete_residue_with_full_spec(imol: int, imodel: int, chain_id: str, resno: int, inscode: str, altloc: str) -> None:
        """ delete residue with altconf """

def delete_residue_hydrogens(imol: int, chain_id: str, resno: int, inscode: str, altloc: str) -> None:
        """ delete hydrogen atoms in residue """

def delete_atom(imol: int, chain_id: str, resno: int, ins_code: str, at_name: str, altloc: str) -> None:
        """ delete atom in residue """

def delete_residue_sidechain(imol: int, chain_id: str, resno: int, ins_code: str, do_delete_dialog: int) -> None:
        """ delete all atoms in residue that are not main chain or CB """

def delete_hydrogen_atoms(imol: int) -> int:
        """ delete all hydrogens in molecule, 

        :return: number of hydrogens deleted. """
        return 0

def delete_hydrogens(imol: int) -> int:
        """ delete all hydrogens in molecule, 

        :return: number of hydrogens deleted. """
        return 0

def delete_waters(imol: int) -> int:
        """ delete all waters in molecule, 

        :return: number of waters deleted. """
        return 0

def post_delete_item_dialog() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_delete_atom_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_delete_residue_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_delete_residue_zone_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_delete_residue_hydrogens_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_delete_water_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_delete_sidechain_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_delete_sidechain_range_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_delete_chain_mode() -> None:
        """ Sphinx-Doc-Placeholder"""

def delete_item_mode_is_atom_p():
        """ Sphinx-Doc-Placeholder"""
        pass

def delete_item_mode_is_residue_p():
        """ Sphinx-Doc-Placeholder"""
        pass

def delete_item_mode_is_water_p():
        """ Sphinx-Doc-Placeholder"""
        pass

def delete_item_mode_is_sidechain_p():
        """ Sphinx-Doc-Placeholder"""
        pass

def delete_item_mode_is_sidechain_range_p():
        """ Sphinx-Doc-Placeholder"""
        pass

def delete_item_mode_is_chain_p():
        """ Sphinx-Doc-Placeholder"""
        pass

def clear_pending_delete_item() -> None:
        """ Sphinx-Doc-Placeholder"""

def do_rot_trans_setup(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def rot_trans_reset_previous() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_rotate_translate_zone_rotates_about_zone_centre(istate: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_rot_trans_object_type(rt_type: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_rot_trans_object_type() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def do_cis_trans_conversion_setup(istate: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def cis_trans_convert(imol: int, chain_id: str, resno: int, altconf: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def cis_trans_convert_intermediate_atoms() -> int:
        """ cis-trans convert the active residue of the active atom in the inermediate atoms, and continue with the refinement """
        return 0

# --- Mainchain Building Functions ---

def do_db_main(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def db_mainchain(imol: int, chain_id: str, iresno_start: int, iresno_end: int, direction: str) -> int:
        """ CA -> mainchain conversion. direction is either "forwards" or "backwards"

        See also the function below.

        return the new molecule number """
        return 0

def db_mainchains_fragment(imol: int, chain_id: str, res_no: int) -> int:
        """ CA-Zone to Mainchain for a fragment based on the given residue. Both directions are built. This is the modern interface. """
        return 0

# --- Close Molecule Functions ---

def close_molecule(imol: int) -> None:
        """ close the molecule """

# --- Rotamer Functions ---

def set_rotamer_search_mode(mode: int) -> None:
        """ set the mode of rotamer search, options are (ROTAMERSEARCHAUTOMATIC), (ROTAMERSEARCHLOWRES) (aka. "backrub rotamers"), (ROTAMERSEARCHHIGHRES) (with rigid body fitting) """

def rotamer_search_mode_state() -> int:
        """ \ brief get the mode of rotamer search """
        return 0

def setup_rotamers(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def do_rotamers(atom_index: int, imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def show_rotamers_dialog(imol: int, chain_id: str, resno: int, ins_code: str, altconf: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_rotamer_lowest_probability(f: float) -> None:
        """ For Dunbrack rotamers, set the lowest probability to be considered. Set as a percentage i.e. 1.00 is quite low. For Richardson Rotamers, this has no effect. """

def set_rotamer_check_clashes(i: int) -> None:
        """ set a flag: 0 is off, 1 is on """

def auto_fit_best_rotamer(imol_coords: int, chain_id: str, resno: int, insertion_code: str, altloc: str, imol_map: int, clash_flag: int, lowest_probability: float) -> float:
        """ auto fit by rotamer search. return the score, for some not very good reason. clash_flag determines if we use clashes with other residues in the score for this rotamer (or not). It would be cool to call this from a script that went residue by residue along a (newly-built) chain (now available). """
        return 0.0

def auto_fit_rotamer_active_residue() -> float:
        """ auto-fit the rotamer for the active residue """
        return 0.0

def set_auto_fit_best_rotamer_clash_flag(i: int) -> None:
        """ set the clash flag for rotamer search And this functions for [pre-setting] the variables for auto_fit_best_rotamer called interactively (using a graphics_info_t function). 0 off, 1 on. """

def rotamer_score(imol: int, chain_id: str, res_no: int, insertion_code: str, alt_conf: str) -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def setup_auto_fit_rotamer(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def n_rotamers(imol: int, chain_id: str, resno: int, ins_code: str) -> int:
        """ return the number of rotamers for this residue - return -1 on no residue found. """
        return 0

def set_residue_to_rotamer_number(imol: int, chain_id: str, resno: int, ins_code: str, alt_conf: str, rotamer_number: int) -> int:
        """ set the residue specified to the rotamer number specifed. """
        return 0

def set_residue_to_rotamer_name(imol: int, chain_id: str, resno: int, ins_code: str, alt_conf: str, rotamer_name: str) -> int:
        """ set the residue specified to the rotamer name specified. (rotamer names are the Richardson rotamer names.)

        return value is 0 if atoms were not moved (e.g. because rotamer-name was not know) """
        return 0

def fill_partial_residues(imol: int) -> None:
        """ fill all the residues of molecule number imol that have missing atoms. To be used to remove the effects of chainsaw. """

def fill_partial_residue(imol: int, chain_id: str, resno: int, inscode: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def simple_fill_partial_residues(imol: int) -> None:
        """ Fill amino acid residues. do backrub rotamer search for residues, but don't do refinement """

# --- 180 Flip Side chain ---

def do_180_degree_side_chain_flip(imol: int, chain_id: str, resno: int, inscode: str, altconf: str) -> None:
        """ rotate 180 degrees around the last chi angle """

def setup_180_degree_flip(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def side_chain_flip_180_intermediate_atoms() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

# --- Mutate Functions ---

def setup_mutate(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_mutate_auto_fit(state: int) -> None:
        """ Mutate then fit to map. that we have a map define is checked first """

def do_mutation(type: str, is_stub_flag: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def mutate_active_residue() -> None:
        """ display a dialog that allows the choice of residue type to which to mutate """

def progressive_residues_in_chain_check(chain_id: str, imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def mutate(imol: int, chain_id: str, ires: int, inscode: str, target_res_type: str) -> int:
        """ mutate a given residue target_res_type is a three-letter-code.

        Return 1 on a good mutate. """
        return 0

def mutate_base(imol: int, chain_id: str, res_no: int, ins_code: str, res_type: str) -> int:
        """ mutate a base. return success status, 1 for a good mutate. """
        return 0

def nudge_residue_sequence(imol: int, chain_id: str, res_no_range_start: int, res_no_range_end: int, nudge_by: int, nudge_residue_numbers_also: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_mutate_auto_fit_do_post_refine(istate: int) -> None:
        """ Do you want Coot to automatically run a refinement after every mutate and autofit? 1 for yes, 0 for no. """

def mutate_auto_fit_do_post_refine_state() -> int:
        """ what is the value of the previous flag? """
        return 0

def set_rotamer_auto_fit_do_post_refine(istate: int) -> None:
        """ Do you want Coot to automatically run a refinement after every rotamer autofit? 1 for yes, 0 for no. """

def rotamer_auto_fit_do_post_refine_state() -> int:
        """ what is the value of the previous flag? """
        return 0

def mutate_single_residue_by_serial_number(ires_ser: int, chain_id: str, imol: int, target_res_type: str) -> int:
        """ an alternate interface to mutation of a singe residue. 

        :return: 1 on success, 0 on failure

        Note that the target_res_type is a char, not a string (or a char *). So from the scheme interface you'd use (for example) hash backslash A for ALA. """
        return 0

def mutate_single_residue_by_seqno(imol: int, chain_id: str, ires: int, inscode: str, target_res_type: str) -> int:
        """ ires is the seqnum of the residue (conventional) """
        return 0

def mutate_and_autofit_residue_range(imol: int, chain_id: str, start_res_no: int, stop_res_no: int, sequence: str) -> int:
        """ mutate and auto-fit (Move this and the above function into """
        return 0

def do_base_mutation(type: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_residue_type_chooser_stub_state(istat: int) -> None:
        """ set a flag saying that the residue chosen by mutate or auto-fit mutate should only be added as a stub (mainchain + CB) """

def handle_residue_type_chooser_entry_chose_type(entry_text: str, stub_mode: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Alternative Conformation ---

def alt_conf_split_type_number():
        """ Sphinx-Doc-Placeholder"""
        pass

def set_add_alt_conf_split_type_number(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def unset_add_alt_conf_dialog() -> None:
        """ Sphinx-Doc-Placeholder"""

def unset_add_alt_conf_define() -> None:
        """ Sphinx-Doc-Placeholder"""

def altconf() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_add_alt_conf_new_atoms_occupancy(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_add_alt_conf_new_atoms_occupancy() -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def set_show_alt_conf_intermediate_atoms(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def show_alt_conf_intermediate_atoms_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def zero_occupancy_residue_range(imol: int, chain_id: str, ires1: int, ires2: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def fill_occupancy_residue_range(imol: int, chain_id: str, ires1: int, ires2: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_occupancy_residue_range(imol: int, chain_id: str, ires1: int, ires2: int, occ: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_b_factor_residue_range(imol: int, chain_id: str, ires1: int, ires2: int, bval: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def reset_b_factor_residue_range(imol: int, chain_id: str, ires1: int, ires2: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Pointer Atom Functions ---

def place_atom_at_pointer() -> None:
        """ Sphinx-Doc-Placeholder"""

def place_atom_at_pointer_by_window() -> None:
        """ Sphinx-Doc-Placeholder"""

def place_typed_atom_at_pointer(type: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_pointer_atom_is_dummy(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def display_where_is_pointer() -> None:
        """ Sphinx-Doc-Placeholder"""

def create_pointer_atom_molecule_maybe() -> int:
        """ Return the current pointer atom molecule, create a pointer atom molecule if necessary (i.e. when the user has not set it). """
        return 0

def pointer_atom_molecule() -> int:
        """ Return the current pointer atom molecule. """
        return 0

def set_pointer_atom_molecule(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Baton Build Interface Functions ---

def set_baton_mode(i: int) -> None:
        """ toggle so that mouse movement moves the baton not rotates the view. """

def try_set_draw_baton(i: int) -> int:
        """ draw the baton or not """
        return 0

def accept_baton_position() -> None:
        """ accept the baton tip position - a prime candidate for a key binding """

def baton_tip_try_another() -> None:
        """ move the baton tip position - another prime candidate for a key binding """

def baton_tip_previous() -> None:
        """ move the baton tip to the previous position """

def shorten_baton() -> None:
        """ shorten the baton length """

def lengthen_baton() -> None:
        """ lengthen the baton """

def baton_build_delete_last_residue() -> None:
        """ delete the most recently build CA position """

def set_baton_build_params(istart_resno: int, chain_id: str, direction: str) -> None:
        """ set the parameters for the start of a new baton-built fragment. direction can either be "forwards" or "backwards" """

# --- Terminal OXT Atom ---

def add_OXT_to_residue(imol: int, chain_id: str, reso: int, insertion_code: str):
        """ Sphinx-Doc-Placeholder"""
        pass

# --- Crosshairs  Interface ---

def set_draw_crosshairs(i: int) -> None:
        """ draw the distance crosshairs, 0 for off, 1 for on. """

def draw_crosshairs_state():
        """ Sphinx-Doc-Placeholder"""
        pass

# --- Edit Chi Angles ---

def setup_edit_chi_angles(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def rotate_chi(am: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_find_hydrogen_torsions(state: int) -> None:
        """ show torsions that rotate hydrogens in the torsion angle manipulation dialog. Note that this may be needed if, in the dictionary cif file torsion which have as a 4th atom both a hydrogen and a heavier atom bonding to the 3rd atom, but list the 4th atom as a hydrogen (not a heavier atom). """

def set_graphics_edit_current_chi(ichi: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def unset_moving_atom_move_chis() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_moving_atom_move_chis() -> None:
        """ Sphinx-Doc-Placeholder"""

def edit_chi_angles(imol: int, chain_id: str, resno: int, ins_code: str, altconf: str) -> int:
        """ display the edit chi angles gui for the given residue return a status of 0 if it failed to fined the residue, return a value of 1 if it worked. """
        return 0

def set_show_chi_angle_bond(imode: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_edit_chi_angles_reverse_fragment_state(istate: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_torsion_general(state: int) -> None:
        """ beloved torsion general at last makes an entrance onto the Coot scene... """

def toggle_torsion_general_reverse() -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_residue_partial_alt_locs(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Backrubbing function ---

def backrub_rotamer(imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def backrub_rotamer_intermediate_atoms() -> int:
        """ apply rotamer backrub to the active atom of the intermediate atoms """
        return 0

# --- Masks ---

def mask_map_by_molecule(map_mol_no: int, coord_mol_no: int, invert_flag: int) -> int:
        """ generate a new map that has been masked by some coordinates (mask-map-by-molecule map-no mol-no invert?) creates and displays a masked map, cuts down density where the coordinates are (invert is 0). If invert? is 1, cut the density down where there are no atoms atoms. """
        return 0

def mask_map_by_atom_selection(map_mol_no: int, coords_mol_no: int, mmdb_atom_selection: str, invert_flag: int) -> int:
        """ mask map by atom selection """
        return 0

def make_masked_maps_split_by_chain(imol: int, imol_map: int) -> int:
        """ make chain masked maps needs to return a list of values """
        return 0

def set_map_mask_atom_radius(rad: float) -> None:
        """ set the atom radius for map masking """

def map_mask_atom_radius() -> float:
        """ get the atom radius for map masking """
        return 0.0

# --- check Waters Interface ---

def set_check_waters_b_factor_limit(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_check_waters_map_sigma_limit(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_check_waters_min_dist_limit(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_check_waters_max_dist_limit(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def delete_checked_waters_baddies(imol: int, b_factor_lim: float, map_sigma_lim: float, min_dist: float, max_dist: float, part_occ_contact_flag: int, zero_occ_flag: int, logical_operator_and_or_flag: int) -> None:
        """ Delete waters that are fail to meet the given criteria. """

def check_waters_by_difference_map(imol_waters: int, imol_diff_map: int, interactive_flag: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def check_waters_by_difference_map_sigma_level_state() -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def set_check_waters_by_difference_map_sigma_level(f: float) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Least-Squares matching ---

def clear_lsq_matches() -> None:
        """ Sphinx-Doc-Placeholder"""

def add_lsq_match(reference_resno_start: int, reference_resno_end: int, chain_id_reference: str, moving_resno_start: int, moving_resno_end: int, chain_id_moving: str, match_type: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def apply_lsq_matches_simple(imol_reference: int, imol_moving: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def setup_lsq_deviation(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_lsq_plane_define(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def unset_lsq_plane_dialog() -> None:
        """ Sphinx-Doc-Placeholder"""

def remove_last_lsq_plane_atom() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Trim ---

def trim_molecule_by_map(imol_coords: int, imol_map: int, map_level: float, delete_or_zero_occ_flag: int) -> None:
        """ cut off (delete or give zero occupancy) atoms in the given molecule if they are below the given map (absolute) level. """

def trim_molecule_by_b_factor(imol: int, limit: float, keep_higher: int) -> None:
        """ trim the molecule by the value in the B-factor column. If an atom in a residue has a "B-factor" above (or below, if keep_higher is true) limit, then the whole residue is deleted """

def pLDDT_to_b_factor(imol: int) -> None:
        """ convert the value in the B-factor column (typically pLDDT for AlphaFold models) to a temperature factor """

# --- External Ray-Tracing ---

def raster3d(rd3_filename: str) -> None:
        """ create a r3d file for the current view """

def povray(filename: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def renderman(rib_filename: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def make_image_raster3d(filename: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def make_image_povray(filename: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def make_image_raster3d_py(filename: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def make_image_povray_py(filename: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_raster3d_bond_thickness(f: float) -> None:
        """ set the bond thickness for the Raster3D representation """

def set_raster3d_atom_radius(f: float) -> None:
        """ set the atom radius for the Raster3D representation """

def set_raster3d_density_thickness(f: float) -> None:
        """ set the density line thickness for the Raster3D representation """

def set_renderer_show_atoms(istate: int) -> None:
        """ set the flag to show atoms for the Raster3D representation """

def set_raster3d_bone_thickness(f: float) -> None:
        """ set the bone (skeleton) thickness for the Raster3D representation """

def set_raster3d_shadows_enabled(state: int) -> None:
        """ turn off shadows for raster3d output - give argument 0 to turn off """

def set_raster3d_water_sphere(istate: int) -> None:
        """ set the flag to show waters as spheres for the Raster3D representation. 1 show as spheres, 0 the usual stars. """

def set_raster3d_font_size(size_in: str) -> None:
        """ set the font size (as a string) for raster3d """

def raster_screen_shot() -> None:
        """ run raster3d and display the resulting image. """

def raster_screen_shot_py() -> None:
        """ Sphinx-Doc-Placeholder"""

def citation_notice_off() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Superposition (SSM) ---

def superpose(imol1: int, imol2: int, move_imol2_flag: int) -> None:
        """ simple interface to superposition. Superpose all residues of imol2 onto imol1. imol1 is reference, we can either move imol2 or copy it to generate a new molecule depending on the vaule of move_imol2_flag (1 for copy 0 for move). """

def superpose_with_chain_selection(imol1: int, imol2: int, chain_imol1: str, chain_imol2: str, chain_used_flag_imol1: int, chain_used_flag_imol2: int, move_imol2_copy_flag: int) -> None:
        """ chain-based interface to superposition. Superpose the given chains of imol2 onto imol1. imol1 is reference, we can either move imol2 or copy it to generate a new molecule depending on the vaule of move_imol2_flag (1 for move 0 for copy). """

def superpose_with_atom_selection(imol1: int, imol2: int, mmdb_atom_sel_str_1: str, mmdb_atom_sel_str_2: str, move_imol2_copy_flag: int) -> int:
        """ detailed interface to superposition. Superpose the given atom selection (specified by the mmdb atom selection strings) of imol2 onto imol1. imol1 is reference, we can either move imol2 or copy it to generate a new molecule depending on the vaule of move_imol2_flag (1 for move 0 for copy).

        :return: the index of the superposed molecule - which could either be a new molecule (if move_imol2_flag was 1) or the imol2 or -1 (signifying failure to do the SMM superposition). """
        return 0

# --- Helices and Strands ---

def place_helix_here() -> int:
        """ add a helix Add a helix somewhere close to this point in the map, try to fit the orientation. Add to a molecule called "Helix", create it if needed. Create another moecule called "Reverse Helix" if the helix orientation isn't completely unequivocal.

        :return: the index of the new molecule. """
        return 0

def place_strand_here(n_residues: int, n_sample_strands: int) -> int:
        """ add a strands Add a strand close to this point in the map, try to fit the orientation. Add to a molecule called "Strand", create it if needed. n_residues is the estimated number of residues in the strand.

        n_sample_strands is the number of strands from the database tested to fit into this strand density. 8 is a suggested number. 20 for a more rigourous search, but it will be slower.

        :return: the index of the new molecule. """
        return 0

def set_place_helix_here_fudge_factor(ff: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def place_strand_here_dialog() -> None:
        """ show the strand placement gui. Choose the python version in there, if needed. Call scripting function, display it in place, don't return a widget. """

def find_helices() -> int:
        """ autobuild helices Find secondary structure in the current map. Add to a molecule called "Helices", create it if needed.

        :return: the index of the new molecule. """
        return 0

def find_strands() -> int:
        """ autobuild strands Find secondary structure in the current map. Add to a molecule called "Strands", create it if needed.

        :return: the index of the new molecule. """
        return 0

def find_secondary_structure(use_helix: int, helix_length: int, helix_target: int, use_strand: int, strand_length: int, strand_target: int) -> int:
        """ autobuild secondary structure Find secondary structure in the current map. Add to a molecule called "SecStruc", create it if needed.

        :return: the index of the new molecule. """
        return 0

def find_secondary_structure_local(use_helix: int, helix_length: int, helix_target: int, use_strand: int, strand_length: int, strand_target: int, radius: float) -> int:
        """ autobuild secondary structure Find secondary structure local to current view in the current map. Add to a molecule called "SecStruc", create it if needed.

        :return: the index of the new molecule. """
        return 0

# --- Nucleotides ---

def find_nucleic_acids_local(radius: float) -> int:
        """ autobuild nucleic acid chains Find secondary structure local to current view in the current map. Add to a molecule called "NuclAcid", create it if needed.

        :return: the index of the new molecule. """
        return 0

# --- New Molecule by Section Interface ---

def new_molecule_by_residue_type_selection(imol: int, residue_type: str) -> int:
        """ create a new molecule that consists of only the residue of type residue_type in molecule number imol 

        :return: the new molecule number, -1 means an error. """
        return 0

def new_molecule_by_atom_selection(imol: int, atom_selection: str) -> int:
        """ create a new molecule that consists of only the atoms specified by the mmdb atoms selection string in molecule number imol 

        :return: the new molecule number, -1 means an error. """
        return 0

def new_molecule_by_sphere_selection(imol: int, x: float, y: float, z: float, r: float, allow_partial_residues: int) -> int:
        """ create a new molecule that consists of only the atoms within the given radius (r) of the given position. 

        :return: the new molecule number, -1 means an error. """
        return 0

# --- RNA/DNA ---

def ideal_nucleic_acid(RNA_or_DNA: str, form: str, single_stranged_flag: int, sequence: str) -> int:
        """ create a molecule of idea nucleotides use the given sequence (single letter code)

        RNA_or_DNA is either "RNA" or "DNA"

        form is either "A" or "B"

        :return: the new molecule number or -1 if a problem """
        return 0

def watson_crick_pair(imol: int, chain_id: str, resno: int) -> int:
        """ Return a molecule that contains a residue that is the WC pair partner of the clicked/picked/selected residue. """
        return 0

def watson_crick_pair_for_residue_range(imol: int, chain_id: str, resno_start: int, resno_end: int) -> int:
        """ add base pairs for the given residue range, modify molecule imol by creating a new chain """
        return 0

def setup_base_pairing(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Sequence File (Assignment/Association) ---

def print_sequence_chain(imol: int, chain_id: str) -> None:
        """ Print the sequence to the console of the given molecule. """

def print_sequence_chain_general(imol: int, chain_id: str, pir_format: int, file_output: int, file_name: str) -> None:
        """ optionally write the sequence to the file for the given molecule, optionally in PIR format """

def assign_fasta_sequence(imol: int, chain_id_in: str, seq: str) -> None:
        """ Assign a FASTA sequence to a given chain in the molecule. """

def assign_pir_sequence(imol: int, chain_id_in: str, seq: str) -> None:
        """ Assign a PIR sequence to a given chain in the molecule. If the chain of the molecule already had a chain assigned to it, then this will overwrite that old assignment with the new one. """

def assign_sequence(imol_model: int, imol_map: int, chain_id: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def assign_sequence_from_file(imol: int, file: str) -> None:
        """ Assign a sequence to a given molecule from (whatever) sequence file by alignment. """

def assign_sequence_from_string(imol: int, chain_id_in: str, seq: str) -> None:
        """ Assign a sequence to a given molecule from a simple string. """

def delete_all_sequences_from_molecule(imol: int) -> None:
        """ Delete all the sequences from a given molecule. """

def delete_sequence_by_chain_id(imol: int, chain_id_in: str) -> None:
        """ Delete the sequence for a given chain_id from a given molecule. """

def associate_sequence_from_file(imol: int, file_name: str) -> None:
        """ Associate the sequence to the molecule - to be used later for sequence assignment (.c.f assign_pir_sequence) """

# --- Surface Interface ---

def do_surface(imol: int, istate: int) -> None:
        """ draw surface of molecule number imol if state = 1 draw the surface (normal representation goes away)

        if state = 0 don't draw surface """

def molecule_is_drawn_as_surface_int(imol: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def make_molecular_surface(imol: int, selection_string: str) -> None:
        """ make molecular surface for given atom selection

        per-chain functions can be added later """

def make_electrostatic_surface(imol: int, selection_string: str) -> None:
        """ make electrostatics surface for given atom selection

        per-chain functions can be added later """

def set_electrostatic_surface_charge_range(v: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_electrostatic_surface_charge_range() -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def set_transparent_electrostatic_surface(imol: int, opacity: float) -> None:
        """ simple on/off screendoor transparency at the moment, an opacity > 0.0 will turn on screendoor transparency (stippling). """

def get_electrostatic_surface_opacity(imol: int) -> float:
        """ return 1.0 for non transparent and 0.5 if screendoor transparency has been turned on. """
        return 0.0

# --- FFFearing ---

def fffear_search(imol_model: int, imol_map: int) -> int:
        """ fffear search model in molecule number imol_model in map number imol_map """
        return 0

def set_fffear_angular_resolution(f: float) -> None:
        """ set and return the fffear angular resolution in degrees """

def fffear_angular_resolution() -> float:
        """ return the fffear angular resolution in degrees """
        return 0.0

# --- Remote Control ---

def make_socket_listener_maybe() -> None:
        """ try to make socket listener """

def set_coot_listener_socket_state_internal(sock_state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_socket_string_waiting(s: str) -> None:
        """ feed the main thread a scheme script to evaluate """

def set_socket_python_string_waiting(s: str) -> None:
        """ feed the main thread a python script to evaluate """

def set_remote_control_port(port_number: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_remote_control_port_number() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_tip_of_the_day_flag(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Display Lists for Maps ---

def set_display_lists_for_maps(i: int) -> None:
        """ Should display lists be used for maps? It may speed things up if these are turned on (or off) - depends on graphics card and drivers. Pass 1 for on, 0 for off. """

def display_lists_for_maps_state() -> int:
        """ return the state of display_lists_for_maps. """
        return 0

def update_maps() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Browser Interface ---

def browser_url(url: str) -> None:
        """ try to open given url in Web browser """

def set_browser_interface(browser: str) -> None:
        """ set command to open the web browser, examples are "open" or "mozilla" """

def handle_online_coot_search_request(entry_text: str) -> None:
        """ the search interface find words, construct a url and open it. """

# --- Molprobity Interface ---

def handle_read_draw_probe_dots(dots_file: str) -> None:
        """ pass a filename that contains molprobity's probe output in XtalView format """

def handle_read_draw_probe_dots_unformatted(dots_file: str, imol: int, show_clash_gui_flag: int) -> None:
        """ pass a filename that contains molprobity's probe output in unformatted format """

def set_do_probe_dots_on_rotamers_and_chis(state: int) -> None:
        """ shall we run molprobity for on edit chi angles intermediate atoms? """

def do_probe_dots_on_rotamers_and_chis_state():
        """ return the state of if run molprobity for on edit chi angles intermediate atoms? """
        pass

def set_do_probe_dots_post_refine(state: int) -> None:
        """ shall we run molprobity after a refinement has happened? """

def do_probe_dots_post_refine_state():
        """ show the state of shall we run molprobity after a refinement has happened? """
        pass

def set_do_coot_probe_dots_during_refine(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_do_coot_probe_dots_during_refine():
        """ Sphinx-Doc-Placeholder"""
        pass

def unmangle_hydrogen_name(pdb_hydrogen_name: str):
        """ make an attempt to convert pdb hydrogen name to the name used in Coot (and the refmac dictionary, perhaps). """
        pass

def set_interactive_probe_dots_molprobity_radius(r: float) -> None:
        """ set the radius over which we can run interactive probe, bigger is better but slower. default is 6.0 """

def interactive_probe_dots_molprobity_radius() -> float:
        """ return the radius over which we can run interactive probe. """
        return 0.0

# --- Map Sharpening Interface ---

def sharpen(imol: int, b_factor: float) -> None:
        """ Sharpen map imol by b_factor (note (of course) that positive numbers blur the map). """

def sharpen_with_gompertz_scaling(imol: int, b_factor: float, try_gompertz: int, gompertz_factor: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_map_sharpening_scale_limit(f: float) -> None:
        """ set the limit of the b-factor map sharpening slider (default 30) """

# --- Marking Fixed Atom Interface ---

def setup_fixed_atom_pick(ipick: int, is_unpick: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def clear_all_fixed_atoms(imol: int) -> None:
        """ clear all fixed atoms """

def clear_fixed_atoms_all() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_debug_atom_picking(istate: int) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Partial Charges ---

def show_partial_charge_info(imol: int, chain_id: str, resno: int, ins_code: str) -> None:
        """ show the partial charges for the residue of the given specs (charges are read from the dictionary) """

# --- EM interface ---

def scale_cell(imol_map: int, fac_u: float, fac_v: float, fac_w: float) -> int:
        """ Scale the cell, for use with EM maps, where the cell needs to be adjusted. Use like: (scale-cell 2 1.012 1.012 1.012). Return error status, 1 means it worked, 0 means it did not work. """
        return 0

def segment_map(imol_map: int, low_level: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def segment_map_multi_scale(imol_map: int, low_level: float, b_factor_inc: float, n_rounds: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def map_histogram(imol_map: int) -> None:
        """ make a map histogram """

def set_ignore_pseudo_zeros_for_map_stats(state: int) -> None:
        """ ignore pseudo-zeros when calculationg maps stats (default 1 = true) """

# --- CCP4mg Interface ---

def set_add_ccp4i_projects_to_file_dialogs(state: int) -> None:
        """ allow the user to not add ccp4i directories to the file choosers use state=0 to turn it off """

def write_ccp4mg_picture_description(filename: str) -> None:
        """ write a ccp4mg picture description file """

# --- Dipoles ---

def delete_dipole(imol: int, dipole_number: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def make_and_draw_patterson(mtz_file_name: str, f_col: str, sigf_col: str) -> int:
        """ Make a patterson molecule. 

        :return: a new molecule number or -1 on failure """
        return 0

def make_and_draw_patterson_using_intensities(mtz_file_name: str, i_col: str, sigi_col: str) -> int:
        """ Make a patterson molecule. 

        :return: a new molecule number or -1 on failure """
        return 0

# --- Aux functions ---

def laplacian(imol: int) -> int:
        """ Create the "Laplacian" (-ve second derivative) of the given map. """
        return 0

# --- SMILES ---

def do_smiles_gui() -> None:
        """ display the SMILES string dialog """

def do_tw() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- PHENIX Support ---

def set_button_label_for_external_refinement(button_label: str) -> None:
        """ set the button label of the external Refinement program """

# --- Graphics Text ---

def place_text(text: str, x: float, y: float, z: float, size: int) -> int:
        """ Put text at x,y,z. 

        :return: a text handle"""
        return 0

def remove_text(text_handle: int) -> None:
        """ Remove "3d" text item. """

def edit_text(text_handle: int, new_text: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def text_index_near_position(x: float, y: float, z: float, r: float) -> int:
        """ return the closest text that is with r A of the given position. If no text item is close, then return -1 """
        return 0

# --- PISA Interaction ---

def pisa_interaction(imol_1: int, imol_2: int) -> int:
        """ return the molecule number of the interacting residues. Return -1 if no new model was created. Old, not very useful. """
        return 0

# --- Jiggle Fit ---

def fit_to_map_by_random_jiggle(imol: int, chain_id: str, resno: int, ins_code: str, n_trials: int, jiggle_scale_factor: float) -> float:
        """ jiggle fit to the current refinment map. return < -100 if not possible, else return the new best fit for this residue. """
        return 0.0

def fit_molecule_to_map_by_random_jiggle(imol: int, n_trials: int, jiggle_scale_factor: float) -> float:
        """ jiggle fit the molecule to the current refinment map. return < -100 if not possible, else return the new best fit for this molecule. """
        return 0.0

def fit_molecule_to_map_by_random_jiggle_and_blur(imol: int, n_trials: int, jiggle_scale_factor: float, map_blur_factor: float) -> float:
        """ jiggle fit the molecule to the current refinment map. return < -100 if not possible, else return the new best fit for this molecule - create a map that is blurred by the given factor for fitting """
        return 0.0

def fit_chain_to_map_by_random_jiggle(imol: int, chain_id: str, n_trials: int, jiggle_scale_factor: float) -> float:
        """ jiggle fit the chain to the current refinment map. return < -100 if not possible, else return the new best fit for this chain. """
        return 0.0

def fit_chain_to_map_by_random_jiggle_and_blur(imol: int, chain_id: str, n_trials: int, jiggle_scale_factor: float, map_blur_factor: float) -> float:
        """ jiggle fit the chain to the current refinment map Use a map that is blurred by the give factor for fitting. 

        :return: < -100 if not possible, else return the new best fit for this chain. """
        return 0.0

# --- SBase interface ---

def get_ccp4srs_monomer_and_dictionary(comp_id: str) -> int:
        """ return the new molecule number of the monomer. The monomer will have chainid "A" and residue number 1.

        Return -1 on failure to get monomer. """
        return 0

def get_sbase_monomer(comp_id: str) -> int:
        """ same as above but using old name for back-compatibility """
        return 0

# --- FLE-View ---

def fle_view(imol: int, chain_id: str, res_no: int, ins_code: str, dist_max: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def fle_view_internal(imol: int, chain_id: str, res_no: int, ins_code: str, imol_ligand_fragment: int, prodrg_output_flat_mol_file_name: str, prodrg_output_flat_pdb_file_name: str, prodrg_output_3d_pdb_file_name: str, prodrg_output_dict_cif_file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def fle_view_internal_to_png(imol: int, chain_id: str, res_no: int, ins_code: str, imol_ligand_fragment: int, prodrg_output_flat_mol_file_name: str, prodrg_output_flat_pdb_file_name: str, prodrg_output_3d_pdb_file_name: str, prodrg_output_dict_cif_file_name: str, output_to_png_file_flag: int, png_file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def fle_view_with_rdkit(imol: int, chain_id: str, res_no: int, ins_code: str, residues_near_radius: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def fle_view_with_rdkit_to_png(imol: int, chain_id: str, res_no: int, ins_code: str, residues_near_radius: float, png_file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def fle_view_with_rdkit_to_svg(imol: int, chain_id: str, res_no: int, ins_code: str, residues_near_radius: float, svg_file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def fle_view_with_rdkit_internal(imol: int, chain_id: str, res_no: int, ins_code: str, residues_near_radius: float, file_format: str, file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def fle_view_set_water_dist_max(dist_max: float) -> None:
        """ set the maximum considered distance to water default 3.25 A. """

def fle_view_set_h_bond_dist_max(h_bond_dist_max: float) -> None:
        """ set the maximum considered hydrogen bond distance default 3.9 A. """

def sprout_hydrogens(imol: int, chain_id: str, res_no: int, ins_code: str) -> int:
        """ Add hydrogens to specificied residue. 

        :return: success status."""
        return 0

# --- LSQ-improve ---

def lsq_improve(imol_ref: int, ref_selection: str, imol_moving: int, moving_selection: str, n_res: int, dist_crit: float) -> None:
        """ an slightly-modified implementation of the "lsq_improve" algorithm of Kleywegt and Jones (1997). Note that if a residue selection is specified in the residue selection(s), then the first residue of the given range must exist in the molecule (if not, then mmdb will not select any atoms from that molecule).

        Kleywegt and Jones set n_res to 4 and dist_crit to 6.0. """

# --- single-model view ---

def single_model_view_model_number(imol: int, imodel: int) -> None:
        """ put molecule number imol to display only model number imodel """

def single_model_view_this_model_number(imol: int) -> int:
        """ the current model number being displayed return 0 on non-multimodel-molecule. """
        return 0

def single_model_view_next_model_number(imol: int) -> int:
        """ change the representation to the next model number to be displayed return 0 on non-multimodel-molecule. """
        return 0

def single_model_view_prev_model_number(imol: int) -> int:
        """ change the representation to the previous model number to be displayed return 0 on non-multimodel-molecule. """
        return 0

# --- graphics 2D ligand view ---

def set_show_graphics_ligand_view(state: int) -> None:
        """ set the graphics ligand view state (default is 1 (on)). """

# --- Experimental ---

def fetch_and_superpose_alphafold_models_using_active_molecule() -> None:
        """ Sphinx-Doc-Placeholder"""

def start_ligand_builder_gui() -> None:
        """ display the ligand builder dialog """

def globularize(imol: int) -> None:
        """ globularize the molecule. This is not guaranteed to generate the correct biological entity, but will bring together molecules (chains/domains) that are dispersed throughout the unit cell. """

# --- Uncategorized ---

def try_load_python_extras_dir() -> None:
        """ Sphinx-Doc-Placeholder"""

def reverse_direction_of_fragment(imol: int, chain_id: str, resno: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_reverse_direction(i: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_axis_orientation_matrix(m11: float, m12: float, m13: float, m21: float, m22: float, m23: float, m31: float, m32: float, m33: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_axis_orientation_matrix_usage(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def optimal_B_kurtosis(imol: int) -> float:
        """ Sphinx-Doc-Placeholder"""
        return 0.0

def add_linked_residue(imol: int, chain_id: str, resno: int, ins_code: str, new_residue_comp_id: str, link_type: str, n_trials: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_add_linked_residue_do_fit_and_refine(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def toolbar_multi_refine_stop() -> None:
        """ Sphinx-Doc-Placeholder"""

def toolbar_multi_refine_continue() -> None:
        """ Sphinx-Doc-Placeholder"""

def toolbar_multi_refine_cancel() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_visible_toolbar_multi_refine_stop_button(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_visible_toolbar_multi_refine_continue_button(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_visible_toolbar_multi_refine_cancel_button(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def toolbar_multi_refine_button_set_sensitive(button_type: str, state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def load_tutorial_model_and_data() -> None:
        """ load tutorial model and data """

def run_update_self_maybe() -> None:
        """ Sphinx-Doc-Placeholder"""

def show_go_to_residue_keyboarding_mode_window() -> None:
        """ Sphinx-Doc-Placeholder"""

def handle_go_to_residue_keyboarding_mode(text: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def git_commit() -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def filtered_by_glob(pre_directory: str, data_type: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def string_member(search: str, list: list):
        """ Sphinx-Doc-Placeholder"""
        pass

def compare_strings(a: str, b: str) -> bool:
        """ Sphinx-Doc-Placeholder"""
        return True

def add_cablam_markup(imol: int, cablam_file_name: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def add_cablam_markup_py(imol: int, cablam_log_file_name: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_rotation_centre(pos: list) -> None:
        """ Sphinx-Doc-Placeholder"""

def goto_next_atom_maybe_py(chain_id: str, resno: int, ins_code: str, atom_name: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def goto_prev_atom_maybe_py(chain_id: str, resno: int, ins_code: str, atom_name: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_go_to_atom_from_spec(atom_spec: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_go_to_atom_from_res_spec(spec: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_go_to_atom_from_res_spec_py(residue_spec: object) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_go_to_atom_from_atom_spec_py(residue_spec: object) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def active_atom_spec():
        """ Sphinx-Doc-Placeholder"""
        pass

def active_atom_spec_py():
        """ Sphinx-Doc-Placeholder"""
        pass

def get_symmetry_py(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def clashes_with_symmetry(imol: int, chain_id: str, res_no: int, ins_code: str, clash_dist: float) -> int:
        """ return 1 if this residue clashes with the symmetry-related atoms of the same molecule. 0 means that it did not clash, -1 means that the residue or molecule could not be found or that there was no cell and symmetry. """
        return 0

def add_molecular_symmetry(imol: int, r_00: float, r_01: float, r_02: float, r_10: float, r_11: float, r_12: float, r_20: float, r_21: float, r_22: float, about_origin_x: float, about_origin_y: float, about_origin_z: float) -> None:
        """ Add molecular symmetry

        You will need to know how to expand your point group molecular symmetry to a set of 3x3 matrices. Call this function for every matrix. """

def add_molecular_symmetry_from_mtrix_from_file(imol: int, file_name: str) -> int:
        """ Add molecular symmetry.

        Often molecular symmetry is descibed using MTRIX card in a PDB file header. Use this function to extract and apply such molecular symmmetry """
        return 0

def add_molecular_symmetry_from_mtrix_from_self_file(imol: int) -> int:
        """ Add molecular symmetry

        This is a convenience function for the above - where you don't need to specify the PDB file name. """
        return 0

def regen_map_internal(imol_map: int, weighted_map_indices: list) -> None:
        """ We overwrite the imol_map and we also presume that the grid sampling of the contributing maps match. This makes it much faster to generate than an average map. """

def make_weighted_map_simple_internal(weighted_map_indices: list) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def colour_map_by_other_map(imol_map: int, imol_map_used_for_colouring: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def colour_map_by_other_map_py(imol_map: int, imol_map_used_for_colouring: int, table_bin_start: float, table_bin_size: float, colour_table_list: object) -> None:
        """ the colour_table should be a list of colours So, if the map has 4 entries covering the range from 0 to 1, then the table_bin_size would be 0.25 and the colour_table list would have 4 entries covering the range 0->0.25, 0.25->0.5, 0.5->0.75, 0.75->1.0 """

def export_molecule_as_x3d(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def export_molecule_as_obj(imol: int, file_name: str) -> bool:
        """ Sphinx-Doc-Placeholder"""
        return True

def export_molecule_as_gltf(imol: int, file_name: str) -> bool:
        """ Sphinx-Doc-Placeholder"""
        return True

def colour_map_by_other_map_turn_off(imol_map: int) -> None:
        """ turn of colour map by other map """

def add_density_map_cap() -> None:
        """ Add map caps. """

def recolour_mesh_by_map(imol_model: int, imol_map: int, scale: float, offset: float) -> None:
        """ colour meshes (e.g. Ribbon diagrams) by map scale might be 2 and offset 1 (for example) """

def import_bild(file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def servalcat_fofc(imol_model: int, imol_fofc_map: int, half_map_1: str, half_map_2: str, resolution: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def servalcat_refine(imol_model: int, half_map_1: str, half_map_2: str, mask_map: str, resolution: float) -> None:
        """ resolution in A. """

def run_acedrg_link_generation(acedrg_link_command: str) -> None:
        """ run acedrg link """

def add_toolbar_subprocess_button(button_label: str, subprocess_command: str, arg_list: object, on_completion_function: object, on_completion_args: object) -> None:
        """ run generic process - doesn't work at the moment - on_completion_args is wrongly interpretted. """

def compare_mtimes() -> bool:
        """ Sphinx-Doc-Placeholder"""
        return True

def parse_ccp4i_defs(filename: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def ccp4_project_directory(ccp4_project_name: str) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def add_to_history(ls: list) -> None:
        """ Sphinx-Doc-Placeholder"""

def add_to_history_simple(cmd: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def add_to_history_typed(command: str, args: list) -> None:
        """ Sphinx-Doc-Placeholder"""

def single_quote(s: str) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def pythonize_command_name(s: str) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def schemize_command_name(s: str) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def languagize_command(command_parts: list) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def add_to_database(command_strings: list) -> None:
        """ Sphinx-Doc-Placeholder"""

def merge_molecules_by_vector(add_molecules: list, imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def monomer_restraints_py(monomer_type: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def monomer_restraints_for_molecule_py(monomer_type: str, imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_monomer_restraints_py(monomer_type: str, restraints: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def show_restraints_editor(monomer_type: str) -> None:
        """ show restraints editor """

def show_restraints_editor_by_index(menu_item_index: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def write_restraints_cif_dictionary(monomer_type: str, file_name: str) -> None:
        """ write cif restraints for monomer """

def list_nomenclature_errors(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def list_nomenclature_errors_py(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def show_fix_nomenclature_errors_gui(imol: int, nomenclature_errors: list) -> None:
        """ Sphinx-Doc-Placeholder"""

def dipole_to_py(dp: list):
        """ Sphinx-Doc-Placeholder"""
        pass

def coot_has_guile():
        """ Sphinx-Doc-Placeholder"""
        pass

def coot_can_do_lidia_p() -> bool:
        """ Sphinx-Doc-Placeholder"""
        return True

def run_scheme_command(scheme_command: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def pyrun_simple_string(python_command: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def residue_spec_to_py(res: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def residue_spec_make_triple_py(residue_spec_py: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def residue_spec_from_py(residue_in: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def get_residue_by_type(imol: int, residue_type: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def get_residue_specs_in_mol(imol: int, residue_type: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def get_residue_specs_in_mol_py(imol: int, residue_type: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def resname_from_serial_number(imol: int, chain_id: str, serial_num: int) -> str:
        """ return the rename from a residue serial number 

        :return: blank ("") on failure. """
        return 'a-string'

def residue_name(imol: int, chain_id: str, resno: int, ins_code: str) -> str:
        """ return the residue name of the specified residue """
        return 'a-string'

def serial_number_from_residue_specs(imol: int, chain_id: str, res_no: int, ins_code: str) -> int:
        """ return the serial number of the specified residue 

        :return: -1 on failure to find the residue """
        return 0

def hydrogenate_region(radius: float) -> None:
        """ find the active residue, find the near residues (within radius) create a new molecule, run reduce on that, import hydrogens from the result and apply them to the molecule of the active residue. """

def add_hydrogens_from_file(imol: int, pdb_with_Hs_file_name: str) -> None:
        """ Add hydrogens to imol from the given pdb file. """

def add_hydrogen_atoms_to_residue(imol: int, chain_id: str, res_no: int, ins_code: str) -> None:
        """ add hydrogen atoms to the specified residue """

def add_hydrogen_atoms_to_residue_py(imol: int, residue_spec_py: object) -> None:
        """ add hydrogen atoms to the specified residue """

def atom_info_string_py(imol: int, chain_id: str, resno: int, ins_code: str, atname: str, altconf: str):
        """ output atom info in a python list for use in scripting: in this format [occ, temp_factor, element, x, y, z]. Return empty list if atom not found. *"""
        pass

def molecule_to_pdb_string_py(imol: int):
        """ Return the molecule as a PDB string. """
        pass

def residue_info_py(imol: int, chain_id: str, resno: int, ins_code: str):
        """ Return a list of atom info for each atom in the specified residue: output is like this: [ [[atom-name,alt-conf] [occ,temp_fact,element] [x,y,z]]] """
        pass

def residue_name_py(imol: int, chain_id: str, resno: int, ins_code: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def residue_centre_from_spec_py(imol: int, spec_py: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def chain_fragments_py(imol: int, screen_output_also: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_b_factor_residues_py(imol: int, residue_specs_b_value_tuple_list_py: object) -> None:
        """ Sphinx-Doc-Placeholder"""

def rigid_body_fit_with_residue_ranges(imol: int, ranges: list) -> int:
        """ return 0 on fail to refine (no sensible place to put atoms) and 1 on fitting happened. """
        return 0

def morph_fit_all(imol: int, transformation_averaging_radius: float) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def morph_fit_chain(imol: int, chain_id: str, transformation_averaging_radius: float) -> int:
        """ Morph the given chain. """
        return 0

def morph_fit_residues_py(imol: int, residue_specs: object, transformation_averaging_radius: float) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def morph_fit_residues(imol: int, residue_specs: list, transformation_averaging_radius: float) -> int:
        """ morph the given residues. """
        return 0

def morph_fit_by_secondary_structure_elements(imol: int, chain_id: str) -> int:
        """ morph transformation are based primarily on rigid body refinement of the secondary structure elements. """
        return 0

def check_waters_baddies(imol: int, b_factor_lim: float, map_sigma_lim: float, min_dist: float, max_dist: float, part_occ_contact_flag: int, zero_occ_flag: int, logical_operator_and_or_flag: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def find_blobs(imol_model: int, imol_map: int, cut_off_density_level: float):
        """ Sphinx-Doc-Placeholder"""
        pass

def find_blobs_py(imol_model: int, imol_map: int, cut_off_density_level: float):
        """ Sphinx-Doc-Placeholder"""
        pass

def b_factor_distribution_graph(imol: int) -> None:
        """ B-factor distribution histogram. """

def monomer_lib_3_letter_codes_matching(search_string: str, allow_minimal_descriptions_flag: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def mutate_residue_range(imol: int, chain_id: str, res_no_start: int, res_no_end: int, target_sequence: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def mutate_internal(ires: int, chain_id: str, imol: int, target_res_type: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def mutate_active_residue_to_single_letter_code(slc: str) -> None:
        """ mutate active residue to single letter code slc """

def show_keyboard_mutate_dialog() -> None:
        """ show keyboard mutate dialog """

def mutate_by_overlap(imol: int, chain_id: str, res_no: int, new_type: str) -> int:
        """ mutate by overlap """
        return 0

def overlap_ligands_internal(imol_ligand: int, imol_ref: int, chain_id_ref: str, resno_ref: int, apply_rtop_flag: bool):
        """ Sphinx-Doc-Placeholder"""
        pass

def do_smiles_to_simple_3d_overlay_frame() -> None:
        """ display the SMILES entry. This is the simple version - no dictionary is generated. """

def ligand_search_make_conformers_py():
        """ Sphinx-Doc-Placeholder"""
        pass

def ligand_search_make_conformers_internal():
        """ Sphinx-Doc-Placeholder"""
        pass

def add_animated_ligand_interaction(imol: int, lb: list) -> None:
        """ Sphinx-Doc-Placeholder"""

def cootaneer_internal(imol_map: int, imol_model: int, atom_spec: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def cootaneer_py(imol_map: int, imol_model: int, atom_in_fragment_atom_spec: object) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def is_interesting_dots_object_next_p(vs: list):
        """ Sphinx-Doc-Placeholder"""
        pass

def generic_string_vector_to_list_internal_py(v: list):
        """ Sphinx-Doc-Placeholder"""
        pass

def generic_int_vector_to_list_internal_py(v: list):
        """ Sphinx-Doc-Placeholder"""
        pass

def generic_list_to_string_vector_internal_py(l: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def rtop_to_python(rtop: list):
        """ Sphinx-Doc-Placeholder"""
        pass

def inverse_rtop_py(rtop_py: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def atom_spec_from_python_expression(expr: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def atom_spec_to_py(spec: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_display_control_button_state(imol: int, button_type: str, state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def fullscreen() -> None:
        """ Sphinx-Doc-Placeholder"""

def unfullscreen() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_use_trackpad(state: int) -> None:
        """ Use left-mouse for view rotation. """

def set_use_primary_mouse_button_for_view_rotation(state: int) -> None:
        """ this is an alias for the above """

def get_coords_for_accession_code(code: str) -> None:
        """ if possible, read in the new coords getting coords via web. (no return value because get-url-str does not return one). """

def coot_get_url_as_string_internal(url: str) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def stop_curl_download(file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_drug_mdl_via_wikipedia_and_drugbank(drugname: str) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def fetch_and_superpose_alphafold_models(imol: int) -> None:
        """ fetch and superpose AlphaFold models corresponding to model model must have Uniprot DBREF info in the header. """

def fetch_alphafold_model_for_uniprot_id(uniprot_id: str) -> int:
        """ return the model number """
        return 0

def fetch_emdb_map(emd_accession_code: str) -> None:
        """ Loads up map frmo emdb. """

def fetch_cod_entry(cod_entry_id: str) -> int:
        """ return the COD entry, return a molecule index """
        return 0

def orient_view(imol: int, central_residue_spec: str, neighbour_residue_spec: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def topological_equivalence_chiral_centres(residue_type: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def screendump_tga(file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_framebuffer_scale_factor(sf: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_use_perspective_projection(state: int) -> None:
        """ set use perspective mode """

def use_perspective_projection_state() -> int:
        """ query if perspective mode is being used """
        return 0

def set_perspective_fov(degrees: float) -> None:
        """ set the perspective fov. Default 20 degrees. """

def set_use_ambient_occlusion(state: int) -> None:
        """ set use ambient occlusion """

def use_ambient_occlusion_state() -> int:
        """ query use ambient occlusion """
        return 0

def set_use_depth_blur(state: int) -> None:
        """ set use depth blur """

def use_depth_blur_state() -> int:
        """ query use depth blur """
        return 0

def set_use_fog(state: int) -> None:
        """ set use fog """

def use_fog_state() -> int:
        """ query use fog """
        return 0

def set_use_outline(state: int) -> None:
        """ set use ourline """

def use_outline_state() -> int:
        """ query use outline """
        return 0

def set_map_shininess(imol: int, shininess: float) -> None:
        """ set the map shininess """

def set_map_specular_strength(imol: int, specular_strength: float) -> None:
        """ set the map specular strength """

def set_draw_normals(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def draw_normals_state() -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def set_draw_mesh(imol: int, mesh_index: int, state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def draw_mesh_state(imol: int, mesh_index: int) -> int:
        """ return -1 on unable to lookup mesh """
        return 0

def set_map_material_specular(imol: int, specular_strength: float, shininess: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_model_material_specular(imol: int, specular_strength: float, shininess: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_model_material_ambient(imol: int, r: float, g: float, b: float, alpha: float) -> None:
        """ set the ambient """

def set_model_material_diffuse(imol: int, r: float, g: float, b: float, alpha: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_model_goodselliness(pastelization_factor: float) -> None:
        """ set the goodselliness (pastelization_factor) 0.3 is about right, but "the right value" depends on the renderer, may be some personal choice. """

def set_map_fresnel_settings(imol: int, state: int, bias: float, scale: float, power: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def reload_map_shader() -> None:
        """ Sphinx-Doc-Placeholder"""

def reload_model_shader() -> None:
        """ Sphinx-Doc-Placeholder"""

def set_atom_radius_scale_factor(imol: int, scale_factor: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_use_fancy_lighting(state: int) -> None:
        """ set use fancy lighting (default 1 = true); """

def set_use_simple_lines_for_model_molecules(state: int) -> None:
        """ set use simple lines for model molecule """

def set_fresnel_colour(imol: int, red: float, green: float, blue: float, opacity: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_focus_blur_z_depth(z: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_focus_blur_strength(st: float) -> None:
        """ set focus blur strength """

def set_shadow_strength(s: float) -> None:
        """ set shadow stren """

def set_shadow_resolution(reso_multiplier: int) -> None:
        """ set the shadow resolution (1,2,3,4) """

def set_shadow_box_size(size: float) -> None:
        """ set shadow box size - default 66; """

def set_ssao_kernel_n_samples(n_samples: int) -> None:
        """ set SSAO kernel n samples """

def set_ssao_strength(strength: float) -> None:
        """ set SSAO strength """

def set_ssao_radius(radius: float) -> None:
        """ set SSAO strength """

def set_ssao_bias(bias: float) -> None:
        """ set SSAO bias """

def set_ssao_blur_size(blur_size: int) -> None:
        """ set SSAO blur size (0, 1, or 2) """

def set_shadow_softness(softness: int) -> None:
        """ set the shadow softness (1, 2 or 3) """

def set_shadow_texture_resolution_multiplier(m: int) -> None:
        """ set the shadow softness (1, 2 or 3) """

def set_effects_shader_output_type(type: int) -> None:
        """ adjust the effects shader output type (for debugging effects) """

def set_effects_shader_brightness(f: float) -> None:
        """ adjust the effects shader brightness """

def set_effects_shader_gamma(f: float) -> None:
        """ adjust the effects shader gamma """

def set_bond_smoothness_factor(fac: int) -> None:
        """ set bond smoothness (default 1 (not smooth)) """

def set_draw_gl_ramachandran_plot_during_refinement(state: int) -> None:
        """ set the draw state of the Ramachandran plot display during Real Space Refinement """

def set_fps_timing_scale_factor(f: float) -> None:
        """ set the FPS timing scale factor - default 0.0025 """

def set_draw_background_image(state: bool) -> None:
        """ draw background image """

def read_test_gltf_models() -> None:
        """ Sphinx-Doc-Placeholder"""

def load_gltf_model(gltf_file_name: str) -> None:
        """ load a gltf model """

def scale_model(model_index: int, scale_factor: float) -> None:
        """ load a gltf model """

def reset_framebuffers() -> None:
        """ reset the frame buffers """

# --- Extra Map Functions ---

def auto_read_make_and_draw_maps(filename: str):
        """ read MTZ file filename and from it try to make maps Useful for reading the output of refmac. The default labels (FWT/PHWT and DELFWT/PHDELFWT) can be changed using ...[something]

        :return: a list of molecule numbers for the new maps """
        pass

def auto_read_make_and_draw_maps_from_mtz(file_name: str):
        """ set the flag to do a difference map (too) on auto-read MTZ """
        pass

def auto_read_make_and_draw_maps_from_cns(file_name: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def valid_labels(mtz_file_name: str, f_col: str, phi_col: str, weight_col: str, use_weights_flag: bool) -> int:
        """ does the mtz file have the columms that we want it to have? """
        return 0

def sharpen_blur_map(imol_map: int, b_factor: float) -> int:
        """ make a sharpened or blurred map blurred maps are generated by using a positive value of b_factor.

        :return: the index of the map created by applying a b-factor to the given map. Return -1 on failure. """
        return 0

def sharpen_blur_map_with_resampling(imol_map: int, b_factor: float, resample_factor: float) -> int:
        """ make a sharpened or blurred map with resampling resampling factor might typically be 1.3

        blurred maps are generated by using a positive value of b_factor.

        :return: the index of the map created by applying a b-factor to the given map. Return -1 on failure. """
        return 0

def sharpen_blur_map_with_resampling_threaded_version(imol_map: int, b_factor: float, resample_factor: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def multi_sharpen_blur_map_py(imol_map: int, b_factors_list: object) -> None:
        """ make many sharpened or blurred maps blurred maps are generated by using a positive value of b_factor. """

def amplitude_vs_resolution_py(mol_map: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def flip_hand(imol_map: int) -> int:
        """ Flip the hand of the map. in case it was accidentally generated on the wrong one. 

        :return: the molecule number of the flipped map. """
        return 0

def analyse_map_point_density_change(map_number_list: list, imol_map_mask: int) -> int:
        """ test function for analysis of multiple map """
        return 0

def analyse_map_point_density_change_py(map_number_list: object, imol_map_mask: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def go_to_map_molecule_centre(imol_map: int) -> None:
        """ Go to the centre of the molecule - for Cryo-EM Molecules. """

def b_factor_from_map(imol_map: int) -> float:
        """ b-factor from map calculate structure factors and use the amplitudes to estimate the B-factor of the data using a wilson plot using a low resolution limit of 4.5A. 

        :return: -1 when given a bad map or there were no data beyond 4.5A """
        return 0.0

def map_colour_components_py(imol: int):
        """ return the colour triple of the imolth map (e.g.: (list 0.4 0.6 0.8). If invalid imol return scheme false.

        return the colour triple of the imolth map e.g.: [0.4, 0.6, 0.8]. If invalid imol return Py_False. """
        pass

def read_ccp4_map(filename: str, is_diff_map_flag: int) -> int:
        """ read a CCP4 map or a CNS map (despite the name). """
        return 0

def handle_read_ccp4_map(filename: str, is_diff_map_flag: int) -> int:
        """ same function as above - old name for the function. Deleted from the API at some stage """
        return 0

def handle_read_emdb_data(dir_name: str) -> int:
        """ this reads a EMDB bundle - I don't think they exist any more """
        return 0

def show_map_partition_by_chain_dialog() -> None:
        """ Sphinx-Doc-Placeholder"""

def map_partition_by_chain(imol_map: int, imol_model: int):
        """ Use the function for scriptng. """
        pass

def map_partition_by_chain_threaded(imol_map: int, imol_model: int) -> None:
        """ Use the function for use in the GUI (non-blocking, no results returned) """

def set_use_vertex_gradients_for_map_normals(imol: int, state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def use_vertex_gradients_for_map_normals_for_latest_map() -> None:
        """ the map should be displayed and not a difference map """

def set_use_vertex_gradients_for_map_normals_for_latest_map() -> None:
        """ alias for the above (more canonical naming) """

# --- Multi-Residue Torsion ---

def multi_residue_torsion_fit(imol: int, specs: list, n_trials: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def multi_residue_torsion_fit_py(imol: int, residues_specs_py: object, n_trials: int):
        """ fit residues (note: fit to the current-refinement map) """
        pass

# --- Add an Atom ---

def add_an_atom(element: str) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Nudge the B-factors ---

def nudge_the_temperature_factors_py(imol: int, residue_spec_py: object, amount: float) -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Merge Fragments ---

def merge_fragments(imol: int) -> int:
        """ each fragment is presumed to be in its own chain. """
        return 0

# --- Delete Items ---

def delete_chain(imol: int, chain_id: str) -> None:
        """ delete the chain """

def delete_sidechains_for_chain(imol: int, chain_id: str) -> None:
        """ delete the side chains in the chain """

# --- Execute Refmac ---

def execute_refmac_real(pdb_in_filename: str, pdb_out_filename: str, mtz_in_filename: str, mtz_out_filename: str, cif_lib_filename: str, fobs_col_name: str, sigfobs_col_name: str, r_free_col_name: str, have_sensible_free_r_flag: int, make_molecules_flag: int, refmac_count_string: str, swap_map_colours_post_refmac_flag: int, imol_refmac_map: int, diff_map_flag: int, phase_combine_flag: int, phib_string: str, fom_string: str, ccp4i_project_dir: str) -> None:
        """ if swap_map_colours_post_refmac_flag is not 1 thenn imol_refmac_map is ignored. """

def refmac_name(imol: int) -> str:
        """ the name for refmac 

        :return: a stub name used in the construction of filename for refmac output """
        return 'a-string'

# --- Dictionary Functions ---

def handle_cif_dictionary(filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def read_cif_dictionary(filename: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def handle_cif_dictionary_for_molecule(filename: str, imol_enc: int, new_molecule_from_dictionary_cif_checkbutton_state: int) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def dictionary_entries():
        """ dictionary entries """
        pass

def debug_dictionary() -> None:
        """ debug dictionary information """

def SMILES_for_comp_id(comp_id: str) -> str:
        """ this can throw an exception """
        return 'a-string'

def dictionaries_read_py():
        """ return a list of all the dictionaries read """
        pass

def cif_file_for_comp_id_py(comp_id: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def dictionary_entries_py():
        """ Sphinx-Doc-Placeholder"""
        pass

def SMILES_for_comp_id_py(comp_id: str):
        """ Sphinx-Doc-Placeholder"""
        pass

# --- Using S-expression molecules ---

def clear_and_update_molecule_py(molecule_number: int, molecule_expression: object) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def add_molecule_py(molecule_expression: object, name: str) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def active_residue_py():
        """ Return a list of [imol, chain-id, resno, ins-code, atom-name, alt-conf] for atom that is closest to the screen centre. If there are multiple models with the same coordinates at the screen centre, return the attributes of the atom in the highest number molecule number. return False if no active residue """
        pass

def closest_atom_simple_py():
        """ return the spec of the closest displayed atom Return a list of [imol, chain-id, resno, ins-code, atom-name, alt-conf, [x, y, z]] for atom that is closest to the screen centre in the given molecule (unlike active-residue, potential CA substition is not performed). If there is no atom, or if imol is not a valid model molecule, return False. """
        pass

def closest_atom_py(imol: int):
        """ return closest atom in imolth molecule Return a list of [imol, chain-id, resno, ins-code, atom-name, alt-conf, [x, y, z]] for atom that is closest to the screen centre in the given molecule (unlike active-residue, no account is taken of the displayed state of the molecule). If there is no atom, or if imol is not a valid model molecule, return False. """
        pass

def closest_atom_raw_py():
        """ return the specs of the closest atom to the centre of the screen Return a list of (list imol chain-id resno ins-code atom-name alt-conf (list x y z)) for atom that is closest to the screen for displayed molecules. If there is no atom, return scheme false. Don't choose the CA of the residue if there is a CA in the residue of the closest atom """
        pass

def residues_near_residue_py(imol: int, residue_in: object, radius: float):
        """ Sphinx-Doc-Placeholder"""
        pass

def residues_near_residues_py(imol: int, residues_in: object, radius: float):
        """ Sphinx-Doc-Placeholder"""
        pass

def residues_near_position_py(imol: int, pos_in: object, radius: float):
        """ Return residue specs for residues that have atoms that are closer than radius Angstroems to the given position. """
        pass

def label_closest_atoms_in_neighbour_residues_py(imol: int, residue_spec_py: object, radius: float) -> None:
        """ label the closest atoms in the residues that neighbour residue_spec """

def get_bonds_representation(imol: int):
        """ return a Python object for the bonds """
        pass

def set_new_non_drawn_bonds(imol: int, cid: str) -> None:
        """ replace any current non-drawn bonds with these - and regen bonds """

def add_to_non_drawn_bonds(imol: int, cid: str) -> None:
        """ add to non-drawn bonds - and regen bonds """

def clear_non_drawn_bonds(imol: int) -> None:
        """ clear the non-drawn bonds - force regen bonds to restore all """

def get_dictionary_radii():
        """ return a Python object for the radii of the atoms in the dictionary """
        pass

def get_environment_distances_representation_py(imol: int, residue_spec_py: object):
        """ return a Python object for the representation of bump and hydrogen bonds oft """
        pass

def get_intermediate_atoms_bonds_representation():
        """ return a Python object for the intermediate atoms bonds """
        pass

def get_continue_updating_refinement_atoms_state() -> int:
        """ return the continue-updating-refinement-atoms state """
        return 0

# --- status bar string functions ---

def atom_info_as_text_for_statusbar(atom_index: int, imol: int) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def atom_info_as_text_for_statusbar(atom_index: int, imol: int, sts: list) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

# --- Refinement with specs ---

def all_residues_with_serial_numbers_py(imol: int):
        """ a utility to return the specs of all the residues, each spec prefixed by the serial number """
        pass

def regularize_residues(imol: int, residues: list) -> None:
        """ regularize the given residues """

def mtz_file_name(imol: int) -> str:
        """ presumes that imol_Refinement_Map has been set """
        return 'a-string'

def refine_zone_with_full_residue_spec_py(imol: int, chain_id: str, resno1: int, inscode_1: str, resno2: int, inscode_2: str, altconf: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_draw_moving_atoms_rota_markup(state: int) -> None:
        """ set display of rotamer markup during interactive real space refinement """

def set_draw_moving_atoms_rama_markup(state: int) -> None:
        """ set display of ramachandran markup during interactive real space refinement """

def set_show_intermediate_atoms_rota_markup(state: int) -> None:
        """ the old names for the above functions: """

def set_show_intermediate_atoms_rama_markup(state: int) -> None:
        """ the old names for the above functions: """

def get_draw_moving_atoms_rota_markup_state() -> int:
        """ the geters for the rota markup """
        return 0

def get_draw_moving_atoms_rama_markup_state() -> int:
        """ the geters for the rama markup """
        return 0

def get_show_intermediate_atoms_rota_markup() -> int:
        """ the old names for the above functions: """
        return 0

def get_show_intermediate_atoms_rama_markup() -> int:
        """ the old names for the above functions: """
        return 0

def set_cryo_em_refinement(mode: bool) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_cryo_em_refinement() -> bool:
        """ Sphinx-Doc-Placeholder"""
        return True

def accept_moving_atoms_py():
        """ Sphinx-Doc-Placeholder"""
        pass

def register_post_intermediate_atoms_moved_hook(function_name: object) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_regenerate_bonds_needs_make_bonds_type_checked(state: bool) -> None:
        """ Sphinx-Doc-Placeholder"""

def get_regenerate_bonds_needs_make_bonds_type_checked_state() -> bool:
        """ Sphinx-Doc-Placeholder"""
        return True

# --- Water Chain Functions ---

def water_chain_from_shelx_ins_py(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def water_chain_py(imol: int):
        """ return the chain id of the water chain. Raw interface """
        pass

# --- Interface Utils ---

def add_status_bar_text(s: str) -> None:
        """ Put text s into the status bar. use this to put info for the user in the statusbar (less intrusive than popup). """

def set_logging_level(level: str) -> None:
        """ set the logging level

        :param level:  is either "LOW" or "HIGH" or "DEBUGGING" """

# --- Glyco Tools ---

def print_glyco_tree(imol: int, chain_id: str, resno: int, ins_code: str) -> None:
        """ print the glycosylation tree that contains the specified residue """

# --- Variance Map ---

def make_variance_map(map_molecule_number_vec: list) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def make_variance_map_py(map_molecule_number_list: object) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

# --- Spin Search Functions ---

def spin_search_by_atom_vectors(imol_map: int, imol: int, chain_id: str, resno: int, ins_code: str, direction_atoms_list: list, moving_atoms_list: list) -> None:
        """ Sphinx-Doc-Placeholder"""

def spin_search_py(imol_map: int, imol: int, chain_id: str, resno: int, ins_code: str, direction_atoms_list: object, moving_atoms_list: object) -> None:
        """ for the given residue, spin the atoms in moving_atom_list... around the bond defined by direction_atoms_list looking for the best fit to density of imom_map map of the first atom in moving_atom_list. Works (only) with atoms in altconf "" """

def spin_N_py(imol: int, residue_spec: object, angle: float) -> None:
        """ Spin N and CB (and the rest of the side chain if extant) Sometime on N-terminal addition, then N ends up pointing the wrong way. The allows us to (more or less) interchange the positions of the CB and the N. angle is in degrees. """

def CG_spin_search_py(imol_model: int, imol_map: int):
        """ Spin search the density based on possible positions of CG of a side-chain. """
        pass

# --- Sequence from Map ---

def sequence_from_map(imol: int, chain_id: str, resno_start: int, resno_end: int, imol_map: int) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def apply_sequence_to_fragment(imol: int, chain_id: str, resno_start: int, resno_end: int, imol_map: int, file_name_for_sequences: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def assign_sequence_to_active_fragment() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Rotamer Scoring ---

def score_rotamers(imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str, imol_map: int, clash_flag: int, lowest_probability: float):
        """ Sphinx-Doc-Placeholder"""
        pass

def score_rotamers_py(imol: int, chain_id: str, res_no: int, ins_code: str, alt_conf: str, imol_map: int, clash_flag: int, lowest_probability: float):
        """ Sphinx-Doc-Placeholder"""
        pass

# --- protein-db ---

def protein_db_loops(imol_coords: int, residue_specs: list, imol_map: int, nfrags: int, preserve_residue_names: bool):
        """ Cowtan's protein_db loops. """
        pass

def protein_db_loop_specs_to_atom_selection_string(specs: list) -> str:
        """ Sphinx-Doc-Placeholder"""
        return 'a-string'

def protein_db_loops_py(imol_coords: int, residues_specs: object, imol_map: int, nfrags: int, preserve_residue_names: bool):
        """ Sphinx-Doc-Placeholder"""
        pass

# --- Coot's Hole implementation ---

def hole(imol: int, start_x: float, start_y: float, start_z: float, end_x: float, end_y: float, end_z: float, colour_map_multiplier: float, colour_map_offset: float, n_runs: int, show_probe_radius_graph_flag: bool, export_surface_dots_file_name: str) -> None:
        """ starting point and end point, colour map multiplier and shall the probe radius graph be shown (dummy value currently). if export_dots_file_name string length is zero, then don't try to export the surface dots. """

# --- Coot's Gaussian Surface ---

def gaussian_surface(imol: int) -> int:
        """ The surface is separated into chains to make Generid Display Objects There is no colour, contour, grid or sigma or material control yet for Gaussian surfaces. there should be one day... """
        return 0

def set_gaussian_surface_sigma(s: float) -> None:
        """ set the sigma for gaussian surface (default 4.0) """

def set_gaussian_surface_contour_level(s: float) -> None:
        """ set the contour_level for gaussian surface (default 4.4) """

def set_gaussian_surface_box_radius(s: float) -> None:
        """ set the box_radius for gaussian surface (defautl 5) """

def set_gaussian_surface_grid_scale(s: float) -> None:
        """ set the grid_scale for gaussian surface (default 0.7) """

def set_gaussian_surface_fft_b_factor(f: float) -> None:
        """ set the fft B-factor for gaussian surface. Use 0 for no B-factor (default 100) """

def set_gaussian_surface_chain_colour_mode(mode: int) -> None:
        """ set the chain colour mode for Gaussian surfaces mode = 1 means each chain has its own colour mode = 2 means the chain colour is determined from NCS/molecular symmetry (so that, in this mode, chains with the same sequence have the same colour """

def show_gaussian_surface_overlay() -> None:
        """ Sphinx-Doc-Placeholder"""

def make_acedrg_dictionary_via_CCD_dictionary(imol: int, spec: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def make_link(imol: int, spec_1: str, spec_2: str, link_name: str, length: float) -> None:
        """ make a link between the specified atoms """

def make_link_py(imol: int, spec_1: object, spec_2: object, link_name: str, length: float) -> None:
        """ Sphinx-Doc-Placeholder"""

def link_info_py(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def show_acedrg_link_interface_overlay() -> None:
        """ Sphinx-Doc-Placeholder"""

# --- Drag and Drop Functions ---

def handle_drag_and_drop_string(uri: str) -> int:
        """ handle the string that get when a file or URL is dropped. """
        return 0

# --- Map Display Control ---

def undisplay_all_maps_except(imol_map: int) -> None:
        """ undisplay all maps except the given one """

# --- Map Contouring Functions ---

def map_contours(imol: int, contour_level: float):
        """ return a list of pairs of vertices for the lines """
        pass

def map_contours_as_triangles(imol: int, contour_level: float):
        """ return two lists: a list of vertices and a list of index-triples for connection """
        pass

def set_radial_map_colouring_enabled(imol: int, state: int) -> None:
        """ enable radial map colouring """

def set_radial_map_colouring_centre(imol: int, x: float, y: float, z: float) -> None:
        """ radial map colouring centre """

def set_radial_map_colouring_min_radius(imol: int, r: float) -> None:
        """ radial map colouring min """

def set_radial_map_colouring_max_radius(imol: int, r: float) -> None:
        """ radial map colouring max """

def set_radial_map_colouring_invert(imol: int, invert_state: int) -> None:
        """ radial map colouring inverted colour map """

def set_radial_map_colouring_saturation(imol: int, saturation: float) -> None:
        """ radial map colouring saturation saturation is a number between 0 and 1, typically 0.5 """

# --- Map to Model Correlation ---

def set_map_correlation_atom_radius(r: float) -> None:
        """ The atom radius is not passed as a parameter to correlation. """

def map_to_model_correlation_py(imol: int, residue_specs: object, neighb_residue_specs: object, atom_mask_mode: int, imol_map: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def map_to_model_correlation_stats_py(imol: int, residue_specs: object, neighb_residue_specs: object, atom_mask_mode: int, imol_map: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def map_to_model_correlation_stats_per_residue_range_py(imol: int, chain_id: str, imol_map: int, n_residue_per_residue_range: int, exclude_NOC_flag: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def map_to_model_correlation(imol: int, residue_specs: list, neigh_residue_specs: list, atom_mask_mode: int, imol_map: int) -> float:
        """ atom-mask-mode is as follows: """
        return 0.0

def map_to_model_correlation_stats(imol: int, residue_specs: list, neigh_residue_specs: list, atom_mask_mode: int, imol_map: int):
        """ map to model density correlation stats atom-mask-mode is as follows: """
        pass

def map_to_model_correlation_per_residue(imol: int, specs: list, atom_mask_mode: int, imol_map: int):
        """ map to model density correlation, reported per residue atom-mask-mode is as follows: """
        pass

def map_to_model_correlation_stats_per_residue(imol: int, residue_specs: list, atom_mask_mode: int, atom_radius_for_masking: float, imol_map: int):
        """ map to model density statistics, reported per residue """
        pass

def map_to_model_correlation_stats_per_residue_range(imol: int, chain_id: str, imol_map: int, n_residue_per_residue_range: int, exclude_NOC_flag: int):
        """ map to model density statistics, reported per residue, the middle residue of a range of residues 

        :return: the all-atom stats first and side chains stats second """
        pass

def map_to_model_correlation_per_residue_py(imol: int, residue_specs: object, atom_mask_mode: int, imol_map: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def qq_plot_map_and_model_py(imol: int, residue_specs_py: object, neigh_residue_specs_py: object, atom_mask_mode: int, imol_map: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def density_score_residue(imol: int, chain_id: str, res_no: int, ins_code: str, imol_map: int) -> float:
        """ simple density score for given residue (over-ridden by scripting function) """
        return 0.0

def map_mean_py(imol: int):
        """ return sigma for the given map. Return Python False if not a valid map molecule number. """
        pass

def map_sigma_py(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def map_statistics_py(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

# --- Get Sequence ---

def get_sequence_as_fasta_for_chain(imol: int, chain_id: str) -> str:
        """ get the sequence for chain_id in imol """
        return 'a-string'

def write_sequence(imol: int, file_name: str) -> None:
        """ write the sequence for imol as fasta """

def res_tracer(imol_map: int, pir_file_name: str) -> None:
        """ trace the given map and try to apply the sequence in the given pir file """

def register_interesting_positions_list_py(pos_list: object) -> None:
        """ Sphinx-Doc-Placeholder"""

def molecule_atom_overlaps_py(imol: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def prodrg_import_function(file_name: str, comp_id: str) -> None:
        """ import given mdl file into prodrg or other 3d generation program the function passed to lbg, so that it calls it when a new prodrg-in.mdl file has been made. We no longer have a timeout function waiting for prodrg-in.mdl to be updated/written. """

def sbase_import_function(comp_id: str) -> None:
        """ import molecule from CCP4 SRS (or SBase, as it used to be called). the function passed to lbg, so that it calls it when a new SBase comp_id is required. We no longer have a timeout function waiting for prodrg-in.mdl to be updated/written. """

def align_to_closest_chain(target_seq: str, match_fraction: float):
        """ align sequence to closest chain (compare across all chains in all molecules). Typically match_fraction is 0.95 or so.

        Return the molecule number and chain id if successful, return -1 as the molecule number if not. """
        pass

def resolve_clashing_sidechains_by_deletion(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def resolve_clashing_sidechains_by_rebuilding(imol: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def simple_text_dialog(dialog_title: str, text: str, geom_x: int, geom_y: int) -> None:
        """ make a simple text dialog. """

def graphics_to_phenix_geo_representation(imol: int, mode: int, g: list) -> None:
        """ phenix GEO bonds representation """

def graphics_to_phenix_geo_representation(imol: int, mode: int, geo_file_name: str) -> None:
        """ phenix GEO bonds representation, read from file """

def set_python_draw_function(command_string: str) -> None:
        """ client/server functions """

def pathology_data(mtz_file_name: str, fp_col: str, sigfp_col: str):
        """ Sphinx-Doc-Placeholder"""
        pass

def encode_ints(i1: int, i2: int) -> int:
        """ encoding of ints """
        return 0

def decode_ints(i: int):
        """ Sphinx-Doc-Placeholder"""
        pass

def store_keyed_user_name(key: str, user_name: str, passwd: str) -> None:
        """ store username and password for the database. """

def py_to_residue_specs(s: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def key_sym_code_py(po: object) -> int:
        """ Sphinx-Doc-Placeholder"""
        return 0

def py_symop_strings_to_space_group(symop_string_list: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def set_use_sounds(state: bool) -> None:
        """ enable or diable sounds (coot needs to have been compiled with sounds of course) """

def curmudgeon_mode() -> None:
        """ no sounds """

def halloween() -> None:
        """ Sphinx-Doc-Placeholder"""

def display_svg_from_file_in_a_dialog(file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def display_svg_from_string_in_a_dialog(string: str, title: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def display_pae_from_file_in_a_dialog(file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def read_interesting_places_json_file(file_name: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def setup_tomo_slider(imol: int) -> int:
        """ return the section index (the middle section currently) """
        return 0

def tomo_section_view(imol: int, axis_id: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_tomo_section_view_section(imol: int, section_index: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def set_tomo_picker_mode_is_active(state: int) -> None:
        """ Sphinx-Doc-Placeholder"""

def tomo_map_analysis(imol_map: int, spot_positions: object) -> None:
        """ Sphinx-Doc-Placeholder"""

def tomo_map_analysis_2(imol_map: int, spot_positions: object) -> None:
        """ Sphinx-Doc-Placeholder"""

def reverse_map(imol_map: int) -> None:
        """ negative becomes positive and positive becomes negative. Apply an offset so that most of the map is above zero. """

def read_positron_metadata(z_data: str, table: str) -> None:
        """ Sphinx-Doc-Placeholder"""

def positron_pathway(map_molecule_list_py: object, pathway_points_py: object):
        """ Sphinx-Doc-Placeholder"""
        pass

def positron_plot_py(fn_z_csv: str, fn_s_csv: str, base_map_index_list: object) -> None:
        """ Sphinx-Doc-Placeholder"""

def positron_plot_internal(fn_z_csv: str, fn_s_csv: str, base_map_index_list: list) -> None:
        """ Sphinx-Doc-Placeholder"""

