Coot Python API
==================

``coot``
-----------------

The Virtual Trackball
---------------------

.. autofunction:: coot.vt_surface

.. autofunction:: coot.vt_surface_status

NCS
---

.. autofunction:: coot.set_draw_ncs_ghosts

.. autofunction:: coot.draw_ncs_ghosts_state

.. autofunction:: coot.set_ncs_ghost_bond_thickness

.. autofunction:: coot.ncs_update_ghosts

.. autofunction:: coot.make_dynamically_transformed_ncs_maps

.. autofunction:: coot.make_ncs_ghosts_maybe

.. autofunction:: coot.add_ncs_matrix

.. autofunction:: coot.clear_ncs_ghost_matrices

.. autofunction:: coot.add_strict_ncs_matrix

.. autofunction:: coot.add_strict_ncs_from_mtrix_from_self_file

.. autofunction:: coot.show_strict_ncs_state

.. autofunction:: coot.set_show_strict_ncs

.. autofunction:: coot.set_ncs_homology_level

.. autofunction:: coot.copy_chain

.. autofunction:: coot.copy_from_ncs_master_to_others

.. autofunction:: coot.copy_residue_range_from_ncs_master_to_others

.. autofunction:: coot.ncs_control_change_ncs_master_to_chain

.. autofunction:: coot.ncs_control_change_ncs_master_to_chain_id

.. autofunction:: coot.ncs_control_display_chain

.. autofunction:: coot.set_ncs_matrix_type

.. autofunction:: coot.get_ncs_matrix_state

Startup Functions
-----------------

.. autofunction:: coot.set_prefer_python

.. autofunction:: coot.prefer_python

File System Functions
---------------------

.. autofunction:: coot.make_directory_maybe

.. autofunction:: coot.set_show_paths_in_display_manager

.. autofunction:: coot.show_paths_in_display_manager_state

.. autofunction:: coot.add_coordinates_glob_extension

.. autofunction:: coot.add_data_glob_extension

.. autofunction:: coot.add_dictionary_glob_extension

.. autofunction:: coot.add_map_glob_extension

.. autofunction:: coot.remove_coordinates_glob_extension

.. autofunction:: coot.remove_data_glob_extension

.. autofunction:: coot.remove_dictionary_glob_extension

.. autofunction:: coot.remove_map_glob_extension

.. autofunction:: coot.set_sticky_sort_by_date

.. autofunction:: coot.unset_sticky_sort_by_date

.. autofunction:: coot.set_filter_fileselection_filenames

.. autofunction:: coot.filter_fileselection_filenames_state

.. autofunction:: coot.file_type_coords

.. autofunction:: coot.open_coords_dialog

.. autofunction:: coot.set_file_chooser_selector

.. autofunction:: coot.file_chooser_selector_state

.. autofunction:: coot.set_file_chooser_overwrite

.. autofunction:: coot.file_chooser_overwrite_state

.. autofunction:: coot.export_map_gui

Widget Utilities
----------------

.. autofunction:: coot.set_main_window_title

MTZ and data handling utilities
-------------------------------

.. autofunction:: coot.manage_column_selector

Molecule Info Functions
-----------------------

.. autofunction:: coot.chain_n_residues

.. autofunction:: coot.molecule_centre_internal

.. autofunction:: coot.seqnum_from_serial_number

.. autofunction:: coot.insertion_code_from_serial_number

.. autofunction:: coot.n_models

.. autofunction:: coot.n_chains

.. autofunction:: coot.is_solvent_chain_p

.. autofunction:: coot.is_protein_chain_p

.. autofunction:: coot.is_nucleotide_chain_p

.. autofunction:: coot.n_residues

.. autofunction:: coot.n_atoms

.. autofunction:: coot.sort_chains

.. autofunction:: coot.sort_residues

.. autofunction:: coot.remarks_dialog

.. autofunction:: coot.print_header_secondary_structure_info

.. autofunction:: coot.add_header_secondary_structure_info

.. autofunction:: coot.write_header_secondary_structure_info

.. autofunction:: coot.copy_molecule

.. autofunction:: coot.add_ligand_delete_residue_copy_molecule

.. autofunction:: coot.exchange_chain_ids_for_seg_ids

.. autofunction:: coot.show_remarks_browswer

Library and Utility Functions
-----------------------------

.. autofunction:: coot.git_revision_count

.. autofunction:: coot.svn_revision

.. autofunction:: coot.molecule_name

.. autofunction:: coot.set_molecule_name

.. autofunction:: coot.coot_checked_exit

.. autofunction:: coot.coot_real_exit

.. autofunction:: coot.coot_no_state_real_exit

.. autofunction:: coot.coot_clear_backup_or_real_exit

.. autofunction:: coot.coot_save_state_and_exit

.. autofunction:: coot.first_coords_imol

.. autofunction:: coot.first_small_coords_imol

.. autofunction:: coot.first_unsaved_coords_imol

.. autofunction:: coot.mmcif_sfs_to_mtz

Graphics Utility Functions
--------------------------

.. autofunction:: coot.set_do_anti_aliasing

.. autofunction:: coot.do_anti_aliasing_state

.. autofunction:: coot.set_do_GL_lighting

.. autofunction:: coot.do_GL_lighting_state

.. autofunction:: coot.use_graphics_interface_state

.. autofunction:: coot.set_use_dark_mode

.. autofunction:: coot.python_at_prompt_at_startup_state

.. autofunction:: coot.reset_view

.. autofunction:: coot.set_view_rotation_scale_factor

.. autofunction:: coot.get_number_of_molecules

.. autofunction:: coot.graphics_n_molecules

.. autofunction:: coot.molecule_has_hydrogens_raw

.. autofunction:: coot.own_molecule_number

.. autofunction:: coot.toggle_idle_spin_function

.. autofunction:: coot.toggle_idle_rock_function

.. autofunction:: coot.get_idle_function_rock_target_angle

.. autofunction:: coot.set_rocking_factors

.. autofunction:: coot.set_idle_function_rotate_angle

.. autofunction:: coot.idle_function_rotate_angle

.. autofunction:: coot.make_updating_model_molecule

.. autofunction:: coot.show_calculate_updating_maps_pythonic_gui

.. autofunction:: coot.allow_duplicate_sequence_numbers

.. autofunction:: coot.set_convert_to_v2_atom_names

.. autofunction:: coot.assign_hetatms

.. autofunction:: coot.hetify_residue

.. autofunction:: coot.residue_has_hetatms

.. autofunction:: coot.het_group_n_atoms

.. autofunction:: coot.replace_fragment

.. autofunction:: coot.copy_residue_range

.. autofunction:: coot.clear_and_update_model_molecule_from_file

.. autofunction:: coot.screendump_image

.. autofunction:: coot.check_for_dark_blue_density

.. autofunction:: coot.set_draw_solid_density_surface

.. autofunction:: coot.set_draw_map_standard_lines

.. autofunction:: coot.set_solid_density_surface_opacity

.. autofunction:: coot.get_solid_density_surface_opacity

.. autofunction:: coot.set_flat_shading_for_solid_density_surface

Interface Preferences
---------------------

.. autofunction:: coot.set_scroll_by_wheel_mouse

.. autofunction:: coot.scroll_by_wheel_mouse_state

.. autofunction:: coot.set_auto_recontour_map

.. autofunction:: coot.get_auto_recontour_map

.. autofunction:: coot.set_default_initial_contour_level_for_map

.. autofunction:: coot.set_default_initial_contour_level_for_difference_map

.. autofunction:: coot.print_view_matrix

.. autofunction:: coot.get_view_matrix_element

.. autofunction:: coot.get_view_quaternion_internal

.. autofunction:: coot.set_view_quaternion

.. autofunction:: coot.apply_ncs_to_view_orientation

.. autofunction:: coot.apply_ncs_to_view_orientation_and_screen_centre

.. autofunction:: coot.set_show_fps

.. autofunction:: coot.set_fps_flag

.. autofunction:: coot.get_fps_flag

.. autofunction:: coot.set_show_origin_marker

.. autofunction:: coot.show_origin_marker_state

.. autofunction:: coot.hide_main_toolbar

.. autofunction:: coot.show_main_toolbar

.. autofunction:: coot.suck_model_fit_dialog

.. autofunction:: coot.suck_model_fit_dialog_bl

.. autofunction:: coot.set_model_fit_refine_dialog_stays_on_top

.. autofunction:: coot.model_fit_refine_dialog_stays_on_top_state

.. autofunction:: coot.set_accept_reject_dialog_docked

.. autofunction:: coot.accept_reject_dialog_docked_state

.. autofunction:: coot.set_accept_reject_dialog_docked_show

.. autofunction:: coot.accept_reject_dialog_docked_show_state

.. autofunction:: coot.set_main_toolbar_style

.. autofunction:: coot.main_toolbar_style_state

Mouse Buttons
-------------

.. autofunction:: coot.quanta_buttons

.. autofunction:: coot.quanta_like_zoom

.. autofunction:: coot.set_control_key_for_rotate

.. autofunction:: coot.control_key_for_rotate_state

.. autofunction:: coot.blob_under_pointer_to_screen_centre

Cursor Function
---------------

.. autofunction:: coot.normal_cursor

.. autofunction:: coot.fleur_cursor

.. autofunction:: coot.pick_cursor_maybe

.. autofunction:: coot.rotate_cursor

.. autofunction:: coot.set_pick_cursor_index

Model/Fit/Refine Functions
--------------------------

.. autofunction:: coot.show_select_map_frame

.. autofunction:: coot.set_model_fit_refine_rotate_translate_zone_label

.. autofunction:: coot.set_model_fit_refine_place_atom_at_pointer_label

.. autofunction:: coot.set_refinement_move_atoms_with_zero_occupancy

.. autofunction:: coot.refinement_move_atoms_with_zero_occupancy_state

Backup Functions
----------------

.. autofunction:: coot.make_backup

.. autofunction:: coot.turn_off_backup

.. autofunction:: coot.turn_on_backup

.. autofunction:: coot.backup_state

.. autofunction:: coot.apply_undo

.. autofunction:: coot.apply_redo

.. autofunction:: coot.set_have_unsaved_changes

.. autofunction:: coot.have_unsaved_changes_p

.. autofunction:: coot.set_undo_molecule

.. autofunction:: coot.show_set_undo_molecule_chooser

.. autofunction:: coot.set_unpathed_backup_file_names

.. autofunction:: coot.unpathed_backup_file_names_state

.. autofunction:: coot.set_decoloned_backup_file_names

.. autofunction:: coot.decoloned_backup_file_names_state

.. autofunction:: coot.backup_compress_files_state

.. autofunction:: coot.set_backup_compress_files

Recover Session Function
------------------------

.. autofunction:: coot.recover_session

Map Functions
-------------

.. autofunction:: coot.calc_phases_generic

.. autofunction:: coot.map_from_mtz_by_refmac_calc_phases

.. autofunction:: coot.map_from_mtz_by_calc_phases

.. autofunction:: coot.calculate_maps_and_stats_py

.. autofunction:: coot.sfcalc_genmap

.. autofunction:: coot.set_auto_updating_sfcalc_genmap

.. autofunction:: coot.set_auto_updating_sfcalc_genmaps

.. autofunction:: coot.set_scroll_wheel_map

.. autofunction:: coot.set_scrollable_map

.. autofunction:: coot.scroll_wheel_map

.. autofunction:: coot.save_previous_map_colour

.. autofunction:: coot.restore_previous_map_colour

.. autofunction:: coot.set_active_map_drag_flag

.. autofunction:: coot.get_active_map_drag_flag

.. autofunction:: coot.set_last_map_colour

.. autofunction:: coot.set_map_colour

.. autofunction:: coot.set_map_hexcolour

.. autofunction:: coot.set_contour_level_absolute

.. autofunction:: coot.set_contour_level_in_sigma

.. autofunction:: coot.get_contour_level_absolute

.. autofunction:: coot.get_contour_level_in_sigma

.. autofunction:: coot.set_last_map_sigma_step

.. autofunction:: coot.set_contour_by_sigma_step_by_mol

.. autofunction:: coot.data_resolution

.. autofunction:: coot.model_resolution

.. autofunction:: coot.export_map

.. autofunction:: coot.export_map_fragment

.. autofunction:: coot.export_map_fragment_with_text_radius

.. autofunction:: coot.export_map_fragment_with_origin_shift

.. autofunction:: coot.export_map_fragment_to_plain_file

.. autofunction:: coot.transform_map_raw

.. autofunction:: coot.difference_map

.. autofunction:: coot.set_map_has_symmetry

.. autofunction:: coot.reinterp_map

.. autofunction:: coot.smooth_map

Density Increment
-----------------

.. autofunction:: coot.get_text_for_iso_level_increment_entry

.. autofunction:: coot.get_text_for_diff_map_iso_level_increment_entry

.. autofunction:: coot.set_iso_level_increment

.. autofunction:: coot.get_iso_level_increment

.. autofunction:: coot.set_iso_level_increment_from_text

.. autofunction:: coot.set_diff_map_iso_level_increment

.. autofunction:: coot.get_diff_map_iso_level_increment

.. autofunction:: coot.set_diff_map_iso_level_increment_from_text

.. autofunction:: coot.set_map_sampling_rate_text

.. autofunction:: coot.set_map_sampling_rate

.. autofunction:: coot.get_text_for_map_sampling_rate_text

.. autofunction:: coot.get_map_sampling_rate

.. autofunction:: coot.change_contour_level

.. autofunction:: coot.set_last_map_contour_level

.. autofunction:: coot.set_last_map_contour_level_by_sigma

.. autofunction:: coot.set_stop_scroll_diff_map

.. autofunction:: coot.set_stop_scroll_iso_map

.. autofunction:: coot.set_stop_scroll_iso_map_level

.. autofunction:: coot.set_stop_scroll_diff_map_level

.. autofunction:: coot.set_residue_density_fit_scale_factor

Density Functions
-----------------

.. autofunction:: coot.set_map_line_width

.. autofunction:: coot.map_line_width_state

.. autofunction:: coot.make_and_draw_map

.. autofunction:: coot.read_mtz

.. autofunction:: coot.make_and_draw_map_with_refmac_params

.. autofunction:: coot.make_and_draw_map_with_reso_with_refmac_params

.. autofunction:: coot.make_updating_map

.. autofunction:: coot.stop_updating_molecule

.. autofunction:: coot.mtz_file_has_phases_p

.. autofunction:: coot.is_mtz_file_p

.. autofunction:: coot.cns_file_has_phases_p

.. autofunction:: coot.wrapped_auto_read_make_and_draw_maps

.. autofunction:: coot.set_auto_read_do_difference_map_too

.. autofunction:: coot.auto_read_do_difference_map_too_state

.. autofunction:: coot.set_auto_read_column_labels

.. autofunction:: coot.get_text_for_density_size_widget

.. autofunction:: coot.set_density_size_from_widget

.. autofunction:: coot.get_text_for_density_size_em_widget

.. autofunction:: coot.set_density_size_em_from_widget

.. autofunction:: coot.set_map_radius

.. autofunction:: coot.set_map_radius_em

.. autofunction:: coot.set_density_size

.. autofunction:: coot.set_map_radius_slider_max

.. autofunction:: coot.set_display_intro_string

.. autofunction:: coot.get_map_radius

.. autofunction:: coot.set_esoteric_depth_cue

.. autofunction:: coot.esoteric_depth_cue_state

.. autofunction:: coot.set_swap_difference_map_colours

.. autofunction:: coot.swap_difference_map_colours_state

.. autofunction:: coot.set_map_is_difference_map

.. autofunction:: coot.map_is_difference_map

.. autofunction:: coot.another_level

.. autofunction:: coot.another_level_from_map_molecule_number

.. autofunction:: coot.residue_density_fit_scale_factor

.. autofunction:: coot.density_at_point

Parameters from map
-------------------

.. autofunction:: coot.mtz_hklin_for_map

.. autofunction:: coot.mtz_fp_for_map

.. autofunction:: coot.mtz_phi_for_map

.. autofunction:: coot.mtz_weight_for_map

.. autofunction:: coot.mtz_use_weight_for_map

PDB Functions
-------------

.. autofunction:: coot.write_pdb_file

.. autofunction:: coot.write_cif_file

.. autofunction:: coot.write_residue_range_to_pdb_file

.. autofunction:: coot.write_chain_to_pdb_file

.. autofunction:: coot.quick_save

.. autofunction:: coot.get_write_conect_record_state

.. autofunction:: coot.set_write_conect_record_state

Info Dialog
-----------

.. autofunction:: coot.info_dialog

.. autofunction:: coot.info_dialog_and_text

.. autofunction:: coot.info_dialog_with_markup

Refmac Functions
----------------

.. autofunction:: coot.set_refmac_counter

.. autofunction:: coot.swap_map_colours

.. autofunction:: coot.set_keep_map_colour_after_refmac

.. autofunction:: coot.keep_map_colour_after_refmac_state

Symmetry Functions
------------------

.. autofunction:: coot.get_text_for_symmetry_size_widget

.. autofunction:: coot.set_symmetry_size_from_widget

.. autofunction:: coot.set_symmetry_size

.. autofunction:: coot.get_symmetry_bonds_colour

.. autofunction:: coot.get_show_symmetry

.. autofunction:: coot.set_show_symmetry_master

.. autofunction:: coot.set_show_symmetry_molecule

.. autofunction:: coot.symmetry_as_calphas

.. autofunction:: coot.get_symmetry_as_calphas_state

.. autofunction:: coot.set_symmetry_molecule_rotate_colour_map

.. autofunction:: coot.symmetry_molecule_rotate_colour_map_state

.. autofunction:: coot.set_symmetry_colour_by_symop

.. autofunction:: coot.set_symmetry_whole_chain

.. autofunction:: coot.set_symmetry_atom_labels_expanded

.. autofunction:: coot.has_unit_cell_state

.. autofunction:: coot.add_symmetry_on_to_preferences_and_apply

.. autofunction:: coot.undo_symmetry_view

.. autofunction:: coot.first_molecule_with_symmetry_displayed

.. autofunction:: coot.save_symmetry_coords

.. autofunction:: coot.new_molecule_by_symmetry

.. autofunction:: coot.new_molecule_by_symmetry_with_atom_selection

.. autofunction:: coot.new_molecule_by_symop

.. autofunction:: coot.n_symops

.. autofunction:: coot.move_reference_chain_to_symm_chain_position

.. autofunction:: coot.setup_save_symmetry_coords

.. autofunction:: coot.set_space_group

.. autofunction:: coot.set_unit_cell_and_space_group

.. autofunction:: coot.set_unit_cell_and_space_group_using_molecule

.. autofunction:: coot.set_symmetry_shift_search_size

History Functions
-----------------

.. autofunction:: coot.print_all_history_in_scheme

.. autofunction:: coot.print_all_history_in_python

.. autofunction:: coot.set_console_display_commands_state

.. autofunction:: coot.set_console_display_commands_hilights

State Functions
---------------

.. autofunction:: coot.save_state

.. autofunction:: coot.save_state_file

.. autofunction:: coot.save_state_file_py

.. autofunction:: coot.set_save_state_file_name

.. autofunction:: coot.save_state_file_name_raw

.. autofunction:: coot.set_run_state_file_status

.. autofunction:: coot.run_state_file

.. autofunction:: coot.run_state_file_py

.. autofunction:: coot.run_state_file_maybe

Clipping Functions
------------------

.. autofunction:: coot.increase_clipping_front

.. autofunction:: coot.increase_clipping_back

.. autofunction:: coot.decrease_clipping_front

.. autofunction:: coot.decrease_clipping_back

.. autofunction:: coot.set_clipping_back

.. autofunction:: coot.set_clipping_front

.. autofunction:: coot.get_clipping_plane_front

.. autofunction:: coot.get_clipping_plane_back

Unit Cell interface
-------------------

.. autofunction:: coot.get_show_unit_cell

.. autofunction:: coot.set_show_unit_cells_all

.. autofunction:: coot.set_show_unit_cell

.. autofunction:: coot.set_unit_cell_colour

Colour
------

.. autofunction:: coot.set_symmetry_colour_merge

.. autofunction:: coot.set_colour_map_rotation_on_read_pdb

.. autofunction:: coot.set_colour_map_rotation_on_read_pdb_flag

.. autofunction:: coot.set_colour_map_rotation_on_read_pdb_c_only_flag

.. autofunction:: coot.set_colour_by_chain

.. autofunction:: coot.set_colour_by_ncs_chain

.. autofunction:: coot.set_colour_by_chain_goodsell_mode

.. autofunction:: coot.set_goodsell_chain_colour_wheel_step

.. autofunction:: coot.set_colour_by_molecule

.. autofunction:: coot.get_colour_map_rotation_on_read_pdb_c_only_flag

.. autofunction:: coot.set_symmetry_colour

Map colour
----------

.. autofunction:: coot.set_colour_map_rotation_for_map

.. autofunction:: coot.set_molecule_bonds_colour_map_rotation

.. autofunction:: coot.get_molecule_bonds_colour_map_rotation

Anisotropic Atoms Interface
---------------------------

.. autofunction:: coot.get_limit_aniso

.. autofunction:: coot.get_show_limit_aniso

.. autofunction:: coot.get_show_aniso

.. autofunction:: coot.set_limit_aniso

.. autofunction:: coot.set_show_aniso

.. autofunction:: coot.set_show_aniso_atoms

.. autofunction:: coot.set_show_aniso_atoms_as_ortep

.. autofunction:: coot.set_aniso_limit_size_from_widget

.. autofunction:: coot.get_text_for_aniso_limit_radius_entry

.. autofunction:: coot.set_aniso_probability

.. autofunction:: coot.get_aniso_probability

Display Functions
-----------------

.. autofunction:: coot.set_graphics_window_size

.. autofunction:: coot.set_graphics_window_size_internal

.. autofunction:: coot.set_graphics_window_position

.. autofunction:: coot.store_graphics_window_position

.. autofunction:: coot.graphics_window_size_and_position_to_preferences

.. autofunction:: coot.graphics_draw

.. autofunction:: coot.zalman_stereo_mode

.. autofunction:: coot.hardware_stereo_mode

.. autofunction:: coot.set_stereo_style

.. autofunction:: coot.stereo_mode_state

.. autofunction:: coot.mono_mode

.. autofunction:: coot.side_by_side_stereo_mode

.. autofunction:: coot.set_dti_stereo_mode

.. autofunction:: coot.set_hardware_stereo_angle_factor

.. autofunction:: coot.hardware_stereo_angle_factor_state

.. autofunction:: coot.set_model_display_radius

.. autofunction:: coot.set_model_fit_refine_dialog_position

.. autofunction:: coot.set_display_control_dialog_position

.. autofunction:: coot.set_go_to_atom_window_position

.. autofunction:: coot.set_delete_dialog_position

.. autofunction:: coot.set_rotate_translate_dialog_position

.. autofunction:: coot.set_accept_reject_dialog_position

.. autofunction:: coot.set_ramachandran_plot_dialog_position

.. autofunction:: coot.set_edit_chi_angles_dialog_position

.. autofunction:: coot.set_rotamer_selection_dialog_position

Smooth Scrolling
----------------

.. autofunction:: coot.set_smooth_scroll_flag

.. autofunction:: coot.get_smooth_scroll

.. autofunction:: coot.set_smooth_scroll_steps_str

.. autofunction:: coot.set_smooth_scroll_steps

.. autofunction:: coot.get_text_for_smooth_scroll_steps

.. autofunction:: coot.set_smooth_scroll_limit_str

.. autofunction:: coot.set_smooth_scroll_limit

.. autofunction:: coot.get_text_for_smooth_scroll_limit

Font Parameters
---------------

.. autofunction:: coot.set_font_size

.. autofunction:: coot.get_font_size

.. autofunction:: coot.set_font_colour

.. autofunction:: coot.set_use_stroke_characters

Rotation Centre
---------------

.. autofunction:: coot.set_rotation_centre_size_from_widget

.. autofunction:: coot.get_text_for_rotation_centre_cube_size

.. autofunction:: coot.set_rotation_centre_size

.. autofunction:: coot.set_rotation_centre_cross_hairs_colour

.. autofunction:: coot.recentre_on_read_pdb

.. autofunction:: coot.set_rotation_centre

.. autofunction:: coot.set_rotation_centre_internal

.. autofunction:: coot.rotation_centre_position

.. autofunction:: coot.go_to_ligand

.. autofunction:: coot.set_go_to_ligand_n_atoms_limit

.. autofunction:: coot.set_reorienting_next_residue_mode

Orthogonal Axes
---------------

.. autofunction:: coot.set_draw_axes

Atom Selection Utilities
------------------------

.. autofunction:: coot.atom_index

.. autofunction:: coot.atom_index_full

.. autofunction:: coot.atom_index_first_atom_in_residue

.. autofunction:: coot.atom_index_first_atom_in_residue_with_altconf

.. autofunction:: coot.min_resno_in_chain

.. autofunction:: coot.max_resno_in_chain

.. autofunction:: coot.median_temperature_factor

.. autofunction:: coot.average_temperature_factor

.. autofunction:: coot.standard_deviation_temperature_factor

.. autofunction:: coot.clear_pending_picks

.. autofunction:: coot.centre_of_mass_string

.. autofunction:: coot.centre_of_mass_string_py

.. autofunction:: coot.set_default_temperature_factor_for_new_atoms

.. autofunction:: coot.default_new_atoms_b_factor

.. autofunction:: coot.set_reset_b_factor_moved_atoms

.. autofunction:: coot.get_reset_b_factor_moved_atoms_state

.. autofunction:: coot.set_atom_attribute

.. autofunction:: coot.set_atom_string_attribute

.. autofunction:: coot.set_residue_name

Skeletonization Interface
-------------------------

.. autofunction:: coot.skel_greer_on

.. autofunction:: coot.skel_greer_off

.. autofunction:: coot.skeletonize_map

.. autofunction:: coot.unskeletonize_map

.. autofunction:: coot.set_initial_map_for_skeletonize

.. autofunction:: coot.set_max_skeleton_search_depth

.. autofunction:: coot.get_text_for_skeletonization_level_entry

.. autofunction:: coot.set_skeletonization_level_from_widget

.. autofunction:: coot.get_text_for_skeleton_box_size_entry

.. autofunction:: coot.set_skeleton_box_size_from_widget

.. autofunction:: coot.set_skeleton_box_size

Save Coordinates
----------------

.. autofunction:: coot.save_coordinates

.. autofunction:: coot.set_save_coordinates_in_original_directory

.. autofunction:: coot.save_molecule_number_from_option_menu

.. autofunction:: coot.set_save_molecule_number

Read Phases File Functions
--------------------------

.. autofunction:: coot.read_phs_and_coords_and_make_map

.. autofunction:: coot.read_phs_and_make_map_using_cell_symm_from_previous_mol

.. autofunction:: coot.read_phs_and_make_map_using_cell_symm_from_mol

.. autofunction:: coot.read_phs_and_make_map_using_cell_symm_from_mol_using_implicit_phs_filename

.. autofunction:: coot.read_phs_and_make_map_using_cell_symm

.. autofunction:: coot.read_phs_and_make_map_with_reso_limits

.. autofunction:: coot.graphics_store_phs_filename

.. autofunction:: coot.possible_cell_symm_for_phs_file

.. autofunction:: coot.get_text_for_phs_cell_chooser

Graphics Move
-------------

.. autofunction:: coot.undo_last_move

.. autofunction:: coot.translate_molecule_by

.. autofunction:: coot.transform_molecule_by

.. autofunction:: coot.transform_zone

Go To Atom Widget Functions
---------------------------

.. autofunction:: coot.post_go_to_atom_window

.. autofunction:: coot.go_to_atom_molecule_number

.. autofunction:: coot.go_to_atom_chain_id

.. autofunction:: coot.go_to_atom_atom_name

.. autofunction:: coot.go_to_atom_residue_number

.. autofunction:: coot.go_to_atom_ins_code

.. autofunction:: coot.go_to_atom_alt_conf

.. autofunction:: coot.set_go_to_atom_chain_residue_atom_name

.. autofunction:: coot.set_go_to_atom_chain_residue_atom_name_full

.. autofunction:: coot.set_go_to_atom_chain_residue_atom_name_no_redraw

.. autofunction:: coot.set_go_to_atom_chain_residue_atom_name_strings

.. autofunction:: coot.update_go_to_atom_from_current_position

.. autofunction:: coot.update_go_to_atom_residue_list

.. autofunction:: coot.atom_spec_to_atom_index

.. autofunction:: coot.full_atom_spec_to_atom_index

.. autofunction:: coot.update_go_to_atom_window_on_changed_mol

.. autofunction:: coot.update_go_to_atom_window_on_new_mol

.. autofunction:: coot.update_go_to_atom_window_on_other_molecule_chosen

.. autofunction:: coot.set_go_to_atom_molecule

.. autofunction:: coot.unset_go_to_atom_widget

.. autofunction:: coot.autobuild_ca_off

.. autofunction:: coot.test_fragment

.. autofunction:: coot.do_skeleton_prune

.. autofunction:: coot.test_function

.. autofunction:: coot.glyco_tree_test

Map and Molecule Control
------------------------

.. autofunction:: coot.post_display_control_window

.. autofunction:: coot.add_map_display_control_widgets

.. autofunction:: coot.add_mol_display_control_widgets

.. autofunction:: coot.add_map_and_mol_display_control_widgets

.. autofunction:: coot.reset_graphics_display_control_window

.. autofunction:: coot.close_graphics_display_control_window

.. autofunction:: coot.set_map_displayed

.. autofunction:: coot.set_mol_displayed

.. autofunction:: coot.set_display_only_model_mol

.. autofunction:: coot.set_mol_active

.. autofunction:: coot.mol_is_displayed

.. autofunction:: coot.mol_is_active

.. autofunction:: coot.map_is_displayed

.. autofunction:: coot.set_all_maps_displayed

.. autofunction:: coot.set_all_models_displayed_and_active

.. autofunction:: coot.set_only_last_model_molecule_displayed

.. autofunction:: coot.display_only_active

.. autofunction:: coot.show_spacegroup

Align and Mutate
----------------

.. autofunction:: coot.align_and_mutate

.. autofunction:: coot.set_alignment_gap_and_space_penalty

Renumber Residue Range
----------------------

.. autofunction:: coot.renumber_residue_range

.. autofunction:: coot.change_residue_number

Change Chain ID
---------------

.. autofunction:: coot.change_chain_id

Scripting Interface
-------------------

.. autofunction:: coot.probe_available_p

.. autofunction:: coot.probe_available_p_py

.. autofunction:: coot.post_scripting_window

.. autofunction:: coot.post_scheme_scripting_window

.. autofunction:: coot.run_command_line_scripts

.. autofunction:: coot.set_guile_gui_loaded_flag

.. autofunction:: coot.set_python_gui_loaded_flag

.. autofunction:: coot.set_found_coot_gui

.. autofunction:: coot.set_found_coot_python_gui

Monomer
-------

.. autofunction:: coot.get_monomer_for_molecule_by_index

.. autofunction:: coot.run_script

.. autofunction:: coot.run_guile_script

.. autofunction:: coot.run_python_script

.. autofunction:: coot.import_python_module

Regularization and Refinement
-----------------------------

.. autofunction:: coot.do_regularize

.. autofunction:: coot.do_refine

.. autofunction:: coot.add_planar_peptide_restraints

.. autofunction:: coot.remove_planar_peptide_restraints

.. autofunction:: coot.make_tight_planar_peptide_restraints

.. autofunction:: coot.planar_peptide_restraints_state

.. autofunction:: coot.set_use_trans_peptide_restraints

.. autofunction:: coot.add_omega_torsion_restriants

.. autofunction:: coot.remove_omega_torsion_restriants

.. autofunction:: coot.set_refine_hydrogen_bonds

.. autofunction:: coot.set_refinement_immediate_replacement

.. autofunction:: coot.refinement_immediate_replacement_state

.. autofunction:: coot.set_refine_use_noughties_physics

.. autofunction:: coot.get_refine_use_noughties_physics_state

.. autofunction:: coot.set_residue_selection_flash_frames_number

.. autofunction:: coot.c_accept_moving_atoms

.. autofunction:: coot.accept_regularizement

.. autofunction:: coot.clear_up_moving_atoms

.. autofunction:: coot.clear_moving_atoms_object

.. autofunction:: coot.stop_refinement_internal

.. autofunction:: coot.set_refinement_use_soft_mode_nbc_restraints

.. autofunction:: coot.shiftfield_b_factor_refinement

.. autofunction:: coot.shiftfield_xyz_factor_refinement

.. autofunction:: coot.set_refine_with_torsion_restraints

.. autofunction:: coot.refine_with_torsion_restraints_state

.. autofunction:: coot.set_matrix

.. autofunction:: coot.matrix_state

.. autofunction:: coot.get_map_weight

.. autofunction:: coot.estimate_map_weight

.. autofunction:: coot.set_refine_auto_range_step

.. autofunction:: coot.set_refine_max_residues

.. autofunction:: coot.refine_zone_atom_index_define

.. autofunction:: coot.refine_zone

.. autofunction:: coot.repeat_refine_zone

.. autofunction:: coot.refine_auto_range

.. autofunction:: coot.regularize_zone

.. autofunction:: coot.set_dragged_refinement_steps_per_frame

.. autofunction:: coot.dragged_refinement_steps_per_frame

.. autofunction:: coot.set_refinement_refine_per_frame

.. autofunction:: coot.refinement_refine_per_frame_state

.. autofunction:: coot.set_refinement_drag_elasticity

.. autofunction:: coot.set_refine_ramachandran_angles

.. autofunction:: coot.set_refine_ramachandran_torsion_angles

.. autofunction:: coot.set_refine_ramachandran_restraints_type

.. autofunction:: coot.set_refine_ramachandran_restraints_weight

.. autofunction:: coot.refine_ramachandran_restraints_weight

.. autofunction:: coot.set_torsion_restraints_weight

.. autofunction:: coot.set_refine_rotamers

.. autofunction:: coot.set_refinement_geman_mcclure_alpha_from_text

.. autofunction:: coot.set_refinement_lennard_jones_epsilon_from_text

.. autofunction:: coot.set_refinement_ramachandran_restraints_weight_from_text

.. autofunction:: coot.set_refinement_overall_weight_from_text

.. autofunction:: coot.set_refinement_torsion_weight_from_text

.. autofunction:: coot.set_refine_params_dialog_more_control_frame_is_active

.. autofunction:: coot.refine_ramachandran_angles_state

.. autofunction:: coot.set_numerical_gradients

.. autofunction:: coot.set_debug_refinement

.. autofunction:: coot.set_fix_chiral_volumes_before_refinement

.. autofunction:: coot.check_chiral_volumes

.. autofunction:: coot.set_show_chiral_volume_errors_dialog

.. autofunction:: coot.set_secondary_structure_restraints_type

.. autofunction:: coot.secondary_structure_restraints_type

.. autofunction:: coot.imol_refinement_map

.. autofunction:: coot.set_imol_refinement_map

.. autofunction:: coot.does_residue_exist_p

.. autofunction:: coot.delete_restraints

.. autofunction:: coot.add_extra_bond_restraint

.. autofunction:: coot.add_extra_geman_mcclure_restraint

.. autofunction:: coot.set_show_extra_distance_restraints

.. autofunction:: coot.add_extra_angle_restraint

.. autofunction:: coot.add_extra_torsion_restraint

.. autofunction:: coot.add_extra_start_pos_restraint

.. autofunction:: coot.add_extra_target_position_restraint

.. autofunction:: coot.delete_all_extra_restraints

.. autofunction:: coot.delete_extra_restraints_for_residue

.. autofunction:: coot.delete_extra_restraints_worse_than

.. autofunction:: coot.add_refmac_extra_restraints

.. autofunction:: coot.set_show_extra_restraints

.. autofunction:: coot.extra_restraints_are_shown

.. autofunction:: coot.set_extra_restraints_prosmart_sigma_limits

.. autofunction:: coot.generate_local_self_restraints

.. autofunction:: coot.generate_self_restraints

.. autofunction:: coot.write_interpolated_extra_restraints

.. autofunction:: coot.write_interpolated_models_and_extra_restraints

.. autofunction:: coot.set_show_parallel_plane_restraints

.. autofunction:: coot.parallel_plane_restraints_are_shown

.. autofunction:: coot.add_parallel_plane_restraint

.. autofunction:: coot.set_extra_restraints_representation_for_bonds_go_to_CA

.. autofunction:: coot.set_use_only_extra_torsion_restraints_for_torsions

.. autofunction:: coot.use_only_extra_torsion_restraints_for_torsions_state

.. autofunction:: coot.clear_all_atom_pull_restraints

.. autofunction:: coot.set_auto_clear_atom_pull_restraint

.. autofunction:: coot.get_auto_clear_atom_pull_restraint_state

.. autofunction:: coot.increase_proportional_editing_radius

.. autofunction:: coot.decrease_proportional_editing_radius

Simplex Refinement Interface
----------------------------

.. autofunction:: coot.fit_residue_range_to_map_by_simplex

.. autofunction:: coot.score_residue_range_fit_to_map

Nomenclature Errors
-------------------

.. autofunction:: coot.fix_nomenclature_errors

.. autofunction:: coot.set_nomenclature_errors_on_read

Atom Info  Interface
--------------------

.. autofunction:: coot.output_atom_info_as_text

Residue Info
------------

.. autofunction:: coot.do_residue_info_dialog

.. autofunction:: coot.output_residue_info_dialog

.. autofunction:: coot.residue_info_dialog

.. autofunction:: coot.residue_info_dialog_is_displayed

.. autofunction:: coot.output_residue_info_as_text

.. autofunction:: coot.do_distance_define

.. autofunction:: coot.do_angle_define

.. autofunction:: coot.do_torsion_define

.. autofunction:: coot.residue_info_apply_all_checkbutton_toggled

.. autofunction:: coot.clear_residue_info_edit_list

.. autofunction:: coot.unset_residue_info_widget

.. autofunction:: coot.clear_measure_distances

.. autofunction:: coot.clear_last_measure_distance

Edit Fuctions
-------------

.. autofunction:: coot.do_edit_copy_molecule

.. autofunction:: coot.do_edit_copy_fragment

.. autofunction:: coot.do_edit_replace_fragment

Residue Environment Functions
-----------------------------

.. autofunction:: coot.set_show_environment_distances

.. autofunction:: coot.set_show_environment_distances_bumps

.. autofunction:: coot.set_show_environment_distances_h_bonds

.. autofunction:: coot.show_environment_distances_state

.. autofunction:: coot.set_environment_distances_distance_limits

.. autofunction:: coot.set_show_environment_distances_as_solid

.. autofunction:: coot.set_environment_distances_label_atom

.. autofunction:: coot.label_neighbours

.. autofunction:: coot.label_atoms_in_residue

.. autofunction:: coot.set_show_local_b_factors

.. autofunction:: coot.add_geometry_distance

Pointer Functions
-----------------

.. autofunction:: coot.set_show_pointer_distances

.. autofunction:: coot.show_pointer_distances_state

Zoom Functions
--------------

.. autofunction:: coot.scale_zoom

.. autofunction:: coot.scale_zoom_internal

.. autofunction:: coot.zoom_factor

.. autofunction:: coot.set_smooth_scroll_do_zoom

.. autofunction:: coot.smooth_scroll_do_zoom

.. autofunction:: coot.smooth_scroll_zoom_limit

.. autofunction:: coot.set_smooth_scroll_zoom_limit

.. autofunction:: coot.set_zoom

CNS Data Functions
------------------

.. autofunction:: coot.handle_cns_data_file

.. autofunction:: coot.handle_cns_data_file_with_cell

mmCIF Functions
---------------

.. autofunction:: coot.auto_read_cif_data_with_phases

.. autofunction:: coot.read_cif_data_with_phases_sigmaa

.. autofunction:: coot.read_cif_data_with_phases_diff_sigmaa

.. autofunction:: coot.read_cif_data

.. autofunction:: coot.read_cif_data_2fofc_map

.. autofunction:: coot.read_cif_data_fofc_map

.. autofunction:: coot.read_cif_data_with_phases_fo_fc

.. autofunction:: coot.read_cif_data_with_phases_2fo_fc

.. autofunction:: coot.read_cif_data_with_phases_nfo_fc

.. autofunction:: coot.read_cif_data_with_phases_fo_alpha_calc

.. autofunction:: coot.write_connectivity

.. autofunction:: coot.open_cif_dictionary_file_selector_dialog

.. autofunction:: coot.import_all_refmac_cifs

.. autofunction:: coot.read_small_molecule_cif

.. autofunction:: coot.read_small_molecule_data_cif

.. autofunction:: coot.read_small_molecule_data_cif_and_make_map_using_coords

Validation Functions
--------------------

.. autofunction:: coot.deviant_geometry

.. autofunction:: coot.is_valid_model_molecule

.. autofunction:: coot.is_valid_map_molecule

.. autofunction:: coot.difference_map_peaks

.. autofunction:: coot.set_difference_map_peaks_max_closeness

.. autofunction:: coot.difference_map_peaks_max_closeness

.. autofunction:: coot.clear_diff_map_peaks

.. autofunction:: coot.gln_asn_b_factor_outliers

.. autofunction:: coot.gln_asn_b_factor_outliers_py

.. autofunction:: coot.clear_multi_residue_torsion_mode

.. autofunction:: coot.set_multi_residue_torsion_reverse_mode

.. autofunction:: coot.show_multi_residue_torsion_dialog

.. autofunction:: coot.setup_multi_residue_torsion

.. autofunction:: coot.atom_overlap_score

.. autofunction:: coot.set_show_chiral_volume_outliers

Ramachandran Plot Functions
---------------------------

.. autofunction:: coot.do_ramachandran_plot

.. autofunction:: coot.set_kleywegt_plot_n_diffs

.. autofunction:: coot.set_ramachandran_plot_contour_levels

.. autofunction:: coot.set_ramachandran_plot_background_block_size

.. autofunction:: coot.set_ramachandran_psi_axis_mode

.. autofunction:: coot.ramachandran_psi_axis_mode

.. autofunction:: coot.set_moving_atoms

.. autofunction:: coot.accept_phi_psi_moving_atoms

.. autofunction:: coot.setup_edit_phi_psi

.. autofunction:: coot.setup_dynamic_distances

.. autofunction:: coot.destroy_edit_backbone_rama_plot

.. autofunction:: coot.ramachandran_plot_differences

.. autofunction:: coot.ramachandran_plot_differences_by_chain

Sequence View Interface
-----------------------

.. autofunction:: coot.sequence_view

.. autofunction:: coot.do_sequence_view

.. autofunction:: coot.update_sequence_view_current_position_highlight_from_active_atom

.. autofunction:: coot.change_peptide_carbonyl_by

.. autofunction:: coot.change_peptide_peptide_by

.. autofunction:: coot.execute_setup_backbone_torsion_edit

.. autofunction:: coot.setup_backbone_torsion_edit

.. autofunction:: coot.set_backbone_torsion_peptide_button_start_pos

.. autofunction:: coot.change_peptide_peptide_by_current_button_pos

.. autofunction:: coot.set_backbone_torsion_carbonyl_button_start_pos

.. autofunction:: coot.change_peptide_carbonyl_by_current_button_pos

Atom Labelling
--------------

.. autofunction:: coot.add_atom_label

.. autofunction:: coot.remove_atom_label

.. autofunction:: coot.remove_all_atom_labels

.. autofunction:: coot.set_label_on_recentre_flag

.. autofunction:: coot.centre_atom_label_status

.. autofunction:: coot.set_brief_atom_labels

.. autofunction:: coot.brief_atom_labels_state

.. autofunction:: coot.set_seg_ids_in_atom_labels

Screen Rotation
---------------

.. autofunction:: coot.rotate_y_scene

.. autofunction:: coot.rotate_x_scene

.. autofunction:: coot.rotate_z_scene

.. autofunction:: coot.spin_zoom_trans

Screen Translation
------------------

.. autofunction:: coot.translate_scene_x

.. autofunction:: coot.translate_scene_y

.. autofunction:: coot.translate_scene_z

Views Interface
---------------

.. autofunction:: coot.add_view_here

.. autofunction:: coot.add_view_raw

.. autofunction:: coot.play_views

.. autofunction:: coot.remove_this_view

.. autofunction:: coot.remove_named_view

.. autofunction:: coot.remove_view

.. autofunction:: coot.go_to_first_view

.. autofunction:: coot.go_to_view_number

.. autofunction:: coot.add_spin_view

.. autofunction:: coot.add_view_description

.. autofunction:: coot.add_action_view

.. autofunction:: coot.insert_action_view_after_view

.. autofunction:: coot.n_views

.. autofunction:: coot.save_views

.. autofunction:: coot.views_play_speed

.. autofunction:: coot.set_views_play_speed

.. autofunction:: coot.clear_all_views

Movies Interface
----------------

.. autofunction:: coot.set_movie_file_name_prefix

.. autofunction:: coot.set_movie_frame_number

.. autofunction:: coot.movie_frame_number

.. autofunction:: coot.set_make_movie_mode

Background Colour
-----------------

.. autofunction:: coot.set_background_colour

.. autofunction:: coot.redraw_background

.. autofunction:: coot.background_is_black_p

Ligand Fitting Functions
------------------------

.. autofunction:: coot.set_ligand_acceptable_fit_fraction

.. autofunction:: coot.set_ligand_cluster_sigma_level

.. autofunction:: coot.set_ligand_flexible_ligand_n_samples

.. autofunction:: coot.set_ligand_verbose_reporting

.. autofunction:: coot.set_find_ligand_n_top_ligands

.. autofunction:: coot.set_find_ligand_do_real_space_refinement

.. autofunction:: coot.set_find_ligand_multi_solutions_per_cluster

.. autofunction:: coot.set_find_ligand_mask_waters

.. autofunction:: coot.set_ligand_search_protein_molecule

.. autofunction:: coot.set_ligand_search_map_molecule

.. autofunction:: coot.add_ligand_search_ligand_molecule

.. autofunction:: coot.add_ligand_search_wiggly_ligand_molecule

.. autofunction:: coot.set_find_ligand_here_cluster

.. autofunction:: coot.execute_ligand_search

.. autofunction:: coot.add_ligand_clear_ligands

.. autofunction:: coot.ligand_expert

.. autofunction:: coot.do_find_ligands_dialog

.. autofunction:: coot.match_ligand_atom_names

.. autofunction:: coot.match_ligand_atom_names_to_comp_id

.. autofunction:: coot.exchange_ligand

.. autofunction:: coot.flip_ligand

.. autofunction:: coot.jed_flip

Water Fitting Functions
-----------------------

.. autofunction:: coot.show_create_find_waters_dialog

.. autofunction:: coot.renumber_waters

.. autofunction:: coot.execute_find_waters_real

.. autofunction:: coot.find_waters

.. autofunction:: coot.move_waters_to_around_protein

.. autofunction:: coot.move_hetgroups_to_around_protein

.. autofunction:: coot.max_water_distance

.. autofunction:: coot.get_text_for_find_waters_sigma_cut_off

.. autofunction:: coot.set_value_for_find_waters_sigma_cut_off

.. autofunction:: coot.set_water_check_spherical_variance_limit

.. autofunction:: coot.set_ligand_water_to_protein_distance_limits

.. autofunction:: coot.set_ligand_water_n_cycles

.. autofunction:: coot.set_write_peaksearched_waters

.. autofunction:: coot.execute_find_blobs

.. autofunction:: coot.split_water

Bond Representation
-------------------

.. autofunction:: coot.set_default_bond_thickness

.. autofunction:: coot.set_bond_thickness

.. autofunction:: coot.set_bond_thickness_intermediate_atoms

.. autofunction:: coot.set_use_variable_bond_thickness

.. autofunction:: coot.set_bond_colour_rotation_for_molecule

.. autofunction:: coot.set_draw_stick_mode_atoms_default

.. autofunction:: coot.get_bond_colour_rotation_for_molecule

.. autofunction:: coot.set_unbonded_atom_star_size

.. autofunction:: coot.set_default_representation_type

.. autofunction:: coot.get_default_bond_thickness

.. autofunction:: coot.set_draw_zero_occ_markers

.. autofunction:: coot.set_draw_cis_peptide_markups

.. autofunction:: coot.set_draw_hydrogens

.. autofunction:: coot.draw_hydrogens_state

.. autofunction:: coot.set_draw_stick_mode_atoms

.. autofunction:: coot.set_draw_missing_residues_loops

.. autofunction:: coot.graphics_to_ca_representation

.. autofunction:: coot.graphics_to_colour_by_chain

.. autofunction:: coot.graphics_to_ca_plus_ligands_representation

.. autofunction:: coot.graphics_to_ca_plus_ligands_and_sidechains_representation

.. autofunction:: coot.graphics_to_bonds_no_waters_representation

.. autofunction:: coot.graphics_to_bonds_representation

.. autofunction:: coot.graphics_to_colour_by_molecule

.. autofunction:: coot.graphics_to_ca_plus_ligands_sec_struct_representation

.. autofunction:: coot.graphics_to_sec_struct_bonds_representation

.. autofunction:: coot.graphics_to_rainbow_representation

.. autofunction:: coot.graphics_to_b_factor_representation

.. autofunction:: coot.graphics_to_b_factor_cas_representation

.. autofunction:: coot.graphics_to_occupancy_representation

.. autofunction:: coot.graphics_to_user_defined_atom_colours_representation

.. autofunction:: coot.graphics_to_user_defined_atom_colours_all_atoms_representation

.. autofunction:: coot.get_graphics_molecule_bond_type

.. autofunction:: coot.set_b_factor_bonds_scale_factor

.. autofunction:: coot.change_model_molecule_representation_mode

.. autofunction:: coot.set_use_grey_carbons_for_molecule

.. autofunction:: coot.set_grey_carbon_colour

.. autofunction:: coot.set_draw_moving_atoms_restraints

.. autofunction:: coot.get_draw_moving_atoms_restraints

.. autofunction:: coot.make_ball_and_stick

.. autofunction:: coot.clear_ball_and_stick

.. autofunction:: coot.set_model_molecule_representation_style

.. autofunction:: coot.set_show_molecular_representation

.. autofunction:: coot.set_show_additional_representation

.. autofunction:: coot.set_show_all_additional_representations

.. autofunction:: coot.all_additional_representations_off_except

.. autofunction:: coot.delete_additional_representation

.. autofunction:: coot.additional_representation_by_string

.. autofunction:: coot.additional_representation_by_attributes

.. autofunction:: coot.set_flev_idle_ligand_interactions

.. autofunction:: coot.toggle_flev_idle_ligand_interactions

.. autofunction:: coot.calculate_hydrogen_bonds

.. autofunction:: coot.set_draw_hydrogen_bonds

Dots Representation
-------------------

.. autofunction:: coot.dots

.. autofunction:: coot.set_dots_colour

.. autofunction:: coot.unset_dots_colour

.. autofunction:: coot.clear_dots

.. autofunction:: coot.clear_dots_by_name

.. autofunction:: coot.n_dots_sets

Pep
---

.. autofunction:: coot.do_pepflip

.. autofunction:: coot.pepflip

.. autofunction:: coot.pepflip_intermediate_atoms

.. autofunction:: coot.pepflip_intermediate_atoms_other_peptide

Rigid Body Refinement Interface
-------------------------------

.. autofunction:: coot.do_rigid_body_refine

.. autofunction:: coot.rigid_body_refine_zone

.. autofunction:: coot.rigid_body_refine_by_atom_selection

.. autofunction:: coot.execute_rigid_body_refine

.. autofunction:: coot.set_rigid_body_fit_acceptable_fit_fraction

Dynamic Map
-----------

.. autofunction:: coot.toggle_dynamic_map_display_size

.. autofunction:: coot.toggle_dynamic_map_sampling

.. autofunction:: coot.set_dynamic_map_size_display_on

.. autofunction:: coot.set_dynamic_map_size_display_off

.. autofunction:: coot.get_dynamic_map_size_display

.. autofunction:: coot.set_dynamic_map_sampling_on

.. autofunction:: coot.set_dynamic_map_sampling_off

.. autofunction:: coot.get_dynamic_map_sampling

.. autofunction:: coot.set_dynamic_map_zoom_offset

Add Terminal Residue Functions
------------------------------

.. autofunction:: coot.do_add_terminal_residue

.. autofunction:: coot.set_add_terminal_residue_n_phi_psi_trials

.. autofunction:: coot.set_add_terminal_residue_add_other_residue_flag

.. autofunction:: coot.set_add_terminal_residue_do_rigid_body_refine

.. autofunction:: coot.set_terminal_residue_do_rigid_body_refine

.. autofunction:: coot.set_add_terminal_residue_debug_trials

.. autofunction:: coot.add_terminal_residue_immediate_addition_state

.. autofunction:: coot.set_add_terminal_residue_immediate_addition

.. autofunction:: coot.add_terminal_residue

.. autofunction:: coot.add_nucleotide

.. autofunction:: coot.add_terminal_residue_using_phi_psi

.. autofunction:: coot.set_add_terminal_residue_default_residue_type

.. autofunction:: coot.set_add_terminal_residue_do_post_refine

.. autofunction:: coot.add_terminal_residue_do_post_refine_state

Delete Residues
---------------

.. autofunction:: coot.delete_atom_by_atom_index

.. autofunction:: coot.delete_residue_by_atom_index

.. autofunction:: coot.delete_residue_hydrogens_by_atom_index

.. autofunction:: coot.delete_residue_range

.. autofunction:: coot.delete_residue

.. autofunction:: coot.delete_residue_with_full_spec

.. autofunction:: coot.delete_residue_hydrogens

.. autofunction:: coot.delete_atom

.. autofunction:: coot.delete_residue_sidechain

.. autofunction:: coot.delete_hydrogen_atoms

.. autofunction:: coot.delete_hydrogens

.. autofunction:: coot.delete_waters

.. autofunction:: coot.post_delete_item_dialog

.. autofunction:: coot.set_delete_atom_mode

.. autofunction:: coot.set_delete_residue_mode

.. autofunction:: coot.set_delete_residue_zone_mode

.. autofunction:: coot.set_delete_residue_hydrogens_mode

.. autofunction:: coot.set_delete_water_mode

.. autofunction:: coot.set_delete_sidechain_mode

.. autofunction:: coot.set_delete_sidechain_range_mode

.. autofunction:: coot.set_delete_chain_mode

.. autofunction:: coot.delete_item_mode_is_atom_p

.. autofunction:: coot.delete_item_mode_is_residue_p

.. autofunction:: coot.delete_item_mode_is_water_p

.. autofunction:: coot.delete_item_mode_is_sidechain_p

.. autofunction:: coot.delete_item_mode_is_sidechain_range_p

.. autofunction:: coot.delete_item_mode_is_chain_p

.. autofunction:: coot.clear_pending_delete_item

.. autofunction:: coot.do_rot_trans_setup

.. autofunction:: coot.rot_trans_reset_previous

.. autofunction:: coot.set_rotate_translate_zone_rotates_about_zone_centre

.. autofunction:: coot.set_rot_trans_object_type

.. autofunction:: coot.get_rot_trans_object_type

.. autofunction:: coot.do_cis_trans_conversion_setup

.. autofunction:: coot.cis_trans_convert

.. autofunction:: coot.cis_trans_convert_intermediate_atoms

Mainchain Building Functions
----------------------------

.. autofunction:: coot.do_db_main

.. autofunction:: coot.db_mainchain

.. autofunction:: coot.db_mainchains_fragment

Close Molecule Functions
------------------------

.. autofunction:: coot.close_molecule

Rotamer Functions
-----------------

.. autofunction:: coot.set_rotamer_search_mode

.. autofunction:: coot.rotamer_search_mode_state

.. autofunction:: coot.setup_rotamers

.. autofunction:: coot.do_rotamers

.. autofunction:: coot.show_rotamers_dialog

.. autofunction:: coot.set_rotamer_lowest_probability

.. autofunction:: coot.set_rotamer_check_clashes

.. autofunction:: coot.auto_fit_best_rotamer

.. autofunction:: coot.auto_fit_rotamer_active_residue

.. autofunction:: coot.set_auto_fit_best_rotamer_clash_flag

.. autofunction:: coot.rotamer_score

.. autofunction:: coot.setup_auto_fit_rotamer

.. autofunction:: coot.n_rotamers

.. autofunction:: coot.set_residue_to_rotamer_number

.. autofunction:: coot.set_residue_to_rotamer_name

.. autofunction:: coot.fill_partial_residues

.. autofunction:: coot.fill_partial_residue

.. autofunction:: coot.simple_fill_partial_residues

180 Flip Side chain
-------------------

.. autofunction:: coot.do_180_degree_side_chain_flip

.. autofunction:: coot.setup_180_degree_flip

.. autofunction:: coot.side_chain_flip_180_intermediate_atoms

Mutate Functions
----------------

.. autofunction:: coot.setup_mutate

.. autofunction:: coot.setup_mutate_auto_fit

.. autofunction:: coot.do_mutation

.. autofunction:: coot.mutate_active_residue

.. autofunction:: coot.progressive_residues_in_chain_check

.. autofunction:: coot.mutate

.. autofunction:: coot.mutate_base

.. autofunction:: coot.nudge_residue_sequence

.. autofunction:: coot.set_mutate_auto_fit_do_post_refine

.. autofunction:: coot.mutate_auto_fit_do_post_refine_state

.. autofunction:: coot.set_rotamer_auto_fit_do_post_refine

.. autofunction:: coot.rotamer_auto_fit_do_post_refine_state

.. autofunction:: coot.mutate_single_residue_by_serial_number

.. autofunction:: coot.mutate_single_residue_by_seqno

.. autofunction:: coot.mutate_and_autofit_residue_range

.. autofunction:: coot.do_base_mutation

.. autofunction:: coot.set_residue_type_chooser_stub_state

.. autofunction:: coot.handle_residue_type_chooser_entry_chose_type

Alternative Conformation
------------------------

.. autofunction:: coot.alt_conf_split_type_number

.. autofunction:: coot.set_add_alt_conf_split_type_number

.. autofunction:: coot.unset_add_alt_conf_dialog

.. autofunction:: coot.unset_add_alt_conf_define

.. autofunction:: coot.altconf

.. autofunction:: coot.set_add_alt_conf_new_atoms_occupancy

.. autofunction:: coot.get_add_alt_conf_new_atoms_occupancy

.. autofunction:: coot.set_show_alt_conf_intermediate_atoms

.. autofunction:: coot.show_alt_conf_intermediate_atoms_state

.. autofunction:: coot.zero_occupancy_residue_range

.. autofunction:: coot.fill_occupancy_residue_range

.. autofunction:: coot.set_occupancy_residue_range

.. autofunction:: coot.set_b_factor_residue_range

.. autofunction:: coot.reset_b_factor_residue_range

Pointer Atom Functions
----------------------

.. autofunction:: coot.place_atom_at_pointer

.. autofunction:: coot.place_atom_at_pointer_by_window

.. autofunction:: coot.place_typed_atom_at_pointer

.. autofunction:: coot.set_pointer_atom_is_dummy

.. autofunction:: coot.display_where_is_pointer

.. autofunction:: coot.create_pointer_atom_molecule_maybe

.. autofunction:: coot.pointer_atom_molecule

.. autofunction:: coot.set_pointer_atom_molecule

Baton Build Interface Functions
-------------------------------

.. autofunction:: coot.set_baton_mode

.. autofunction:: coot.try_set_draw_baton

.. autofunction:: coot.accept_baton_position

.. autofunction:: coot.baton_tip_try_another

.. autofunction:: coot.baton_tip_previous

.. autofunction:: coot.shorten_baton

.. autofunction:: coot.lengthen_baton

.. autofunction:: coot.baton_build_delete_last_residue

.. autofunction:: coot.set_baton_build_params

Terminal OXT Atom
-----------------

.. autofunction:: coot.add_OXT_to_residue

Crosshairs  Interface
---------------------

.. autofunction:: coot.set_draw_crosshairs

.. autofunction:: coot.draw_crosshairs_state

Edit Chi Angles
---------------

.. autofunction:: coot.setup_edit_chi_angles

.. autofunction:: coot.rotate_chi

.. autofunction:: coot.set_find_hydrogen_torsions

.. autofunction:: coot.set_graphics_edit_current_chi

.. autofunction:: coot.unset_moving_atom_move_chis

.. autofunction:: coot.set_moving_atom_move_chis

.. autofunction:: coot.edit_chi_angles

.. autofunction:: coot.set_show_chi_angle_bond

.. autofunction:: coot.set_edit_chi_angles_reverse_fragment_state

.. autofunction:: coot.setup_torsion_general

.. autofunction:: coot.toggle_torsion_general_reverse

.. autofunction:: coot.setup_residue_partial_alt_locs

Backrubbing function
--------------------

.. autofunction:: coot.backrub_rotamer

.. autofunction:: coot.backrub_rotamer_intermediate_atoms

Masks
-----

.. autofunction:: coot.mask_map_by_molecule

.. autofunction:: coot.mask_map_by_atom_selection

.. autofunction:: coot.make_masked_maps_split_by_chain

.. autofunction:: coot.set_map_mask_atom_radius

.. autofunction:: coot.map_mask_atom_radius

check Waters Interface
----------------------

.. autofunction:: coot.set_check_waters_b_factor_limit

.. autofunction:: coot.set_check_waters_map_sigma_limit

.. autofunction:: coot.set_check_waters_min_dist_limit

.. autofunction:: coot.set_check_waters_max_dist_limit

.. autofunction:: coot.delete_checked_waters_baddies

.. autofunction:: coot.check_waters_by_difference_map

.. autofunction:: coot.check_waters_by_difference_map_sigma_level_state

.. autofunction:: coot.set_check_waters_by_difference_map_sigma_level

Least
-----

.. autofunction:: coot.clear_lsq_matches

.. autofunction:: coot.add_lsq_match

.. autofunction:: coot.apply_lsq_matches_simple

.. autofunction:: coot.setup_lsq_deviation

.. autofunction:: coot.setup_lsq_plane_define

.. autofunction:: coot.unset_lsq_plane_dialog

.. autofunction:: coot.remove_last_lsq_plane_atom

Trim
----

.. autofunction:: coot.trim_molecule_by_map

.. autofunction:: coot.trim_molecule_by_b_factor

.. autofunction:: coot.pLDDT_to_b_factor

External Ray
------------

.. autofunction:: coot.raster3d

.. autofunction:: coot.povray

.. autofunction:: coot.renderman

.. autofunction:: coot.make_image_raster3d

.. autofunction:: coot.make_image_povray

.. autofunction:: coot.make_image_raster3d_py

.. autofunction:: coot.make_image_povray_py

.. autofunction:: coot.set_raster3d_bond_thickness

.. autofunction:: coot.set_raster3d_atom_radius

.. autofunction:: coot.set_raster3d_density_thickness

.. autofunction:: coot.set_renderer_show_atoms

.. autofunction:: coot.set_raster3d_bone_thickness

.. autofunction:: coot.set_raster3d_shadows_enabled

.. autofunction:: coot.set_raster3d_water_sphere

.. autofunction:: coot.set_raster3d_font_size

.. autofunction:: coot.raster_screen_shot

.. autofunction:: coot.raster_screen_shot_py

.. autofunction:: coot.citation_notice_off

Superposition (SSM)
-------------------

.. autofunction:: coot.superpose

.. autofunction:: coot.superpose_with_chain_selection

.. autofunction:: coot.superpose_with_atom_selection

Helices and Strands
-------------------

.. autofunction:: coot.place_helix_here

.. autofunction:: coot.place_strand_here

.. autofunction:: coot.set_place_helix_here_fudge_factor

.. autofunction:: coot.place_strand_here_dialog

.. autofunction:: coot.find_helices

.. autofunction:: coot.find_strands

.. autofunction:: coot.find_secondary_structure

.. autofunction:: coot.find_secondary_structure_local

Nucleotides
-----------

.. autofunction:: coot.find_nucleic_acids_local

New Molecule by Section Interface
---------------------------------

.. autofunction:: coot.new_molecule_by_residue_type_selection

.. autofunction:: coot.new_molecule_by_atom_selection

.. autofunction:: coot.new_molecule_by_sphere_selection

RNA/DNA
-------

.. autofunction:: coot.ideal_nucleic_acid

.. autofunction:: coot.watson_crick_pair

.. autofunction:: coot.watson_crick_pair_for_residue_range

.. autofunction:: coot.setup_base_pairing

Sequence File (Assignment/Association)
--------------------------------------

.. autofunction:: coot.print_sequence_chain

.. autofunction:: coot.print_sequence_chain_general

.. autofunction:: coot.assign_fasta_sequence

.. autofunction:: coot.assign_pir_sequence

.. autofunction:: coot.assign_sequence

.. autofunction:: coot.assign_sequence_from_file

.. autofunction:: coot.assign_sequence_from_string

.. autofunction:: coot.delete_all_sequences_from_molecule

.. autofunction:: coot.delete_sequence_by_chain_id

.. autofunction:: coot.associate_sequence_from_file

Surface Interface
-----------------

.. autofunction:: coot.do_surface

.. autofunction:: coot.molecule_is_drawn_as_surface_int

.. autofunction:: coot.make_molecular_surface

.. autofunction:: coot.make_electrostatic_surface

.. autofunction:: coot.set_electrostatic_surface_charge_range

.. autofunction:: coot.get_electrostatic_surface_charge_range

.. autofunction:: coot.set_transparent_electrostatic_surface

.. autofunction:: coot.get_electrostatic_surface_opacity

FFFearing
---------

.. autofunction:: coot.fffear_search

.. autofunction:: coot.set_fffear_angular_resolution

.. autofunction:: coot.fffear_angular_resolution

Remote Control
--------------

.. autofunction:: coot.make_socket_listener_maybe

.. autofunction:: coot.set_coot_listener_socket_state_internal

.. autofunction:: coot.set_socket_string_waiting

.. autofunction:: coot.set_socket_python_string_waiting

.. autofunction:: coot.set_remote_control_port

.. autofunction:: coot.get_remote_control_port_number

.. autofunction:: coot.set_tip_of_the_day_flag

Display Lists for Maps
----------------------

.. autofunction:: coot.set_display_lists_for_maps

.. autofunction:: coot.display_lists_for_maps_state

.. autofunction:: coot.update_maps

Browser Interface
-----------------

.. autofunction:: coot.browser_url

.. autofunction:: coot.set_browser_interface

.. autofunction:: coot.handle_online_coot_search_request

Molprobity Interface
--------------------

.. autofunction:: coot.handle_read_draw_probe_dots

.. autofunction:: coot.handle_read_draw_probe_dots_unformatted

.. autofunction:: coot.set_do_probe_dots_on_rotamers_and_chis

.. autofunction:: coot.do_probe_dots_on_rotamers_and_chis_state

.. autofunction:: coot.set_do_probe_dots_post_refine

.. autofunction:: coot.do_probe_dots_post_refine_state

.. autofunction:: coot.set_do_coot_probe_dots_during_refine

.. autofunction:: coot.get_do_coot_probe_dots_during_refine

.. autofunction:: coot.unmangle_hydrogen_name

.. autofunction:: coot.set_interactive_probe_dots_molprobity_radius

.. autofunction:: coot.interactive_probe_dots_molprobity_radius

Map Sharpening Interface
------------------------

.. autofunction:: coot.sharpen

.. autofunction:: coot.sharpen_with_gompertz_scaling

.. autofunction:: coot.set_map_sharpening_scale_limit

Marking Fixed Atom Interface
----------------------------

.. autofunction:: coot.setup_fixed_atom_pick

.. autofunction:: coot.clear_all_fixed_atoms

.. autofunction:: coot.clear_fixed_atoms_all

.. autofunction:: coot.set_debug_atom_picking

Partial Charges
---------------

.. autofunction:: coot.show_partial_charge_info

EM interface
------------

.. autofunction:: coot.scale_cell

.. autofunction:: coot.segment_map

.. autofunction:: coot.segment_map_multi_scale

.. autofunction:: coot.map_histogram

.. autofunction:: coot.set_ignore_pseudo_zeros_for_map_stats

CCP4mg Interface
----------------

.. autofunction:: coot.set_add_ccp4i_projects_to_file_dialogs

.. autofunction:: coot.write_ccp4mg_picture_description

Dipoles
-------

.. autofunction:: coot.delete_dipole

.. autofunction:: coot.make_and_draw_patterson

.. autofunction:: coot.make_and_draw_patterson_using_intensities

Aux functions
-------------

.. autofunction:: coot.laplacian

SMILES
------

.. autofunction:: coot.do_smiles_gui

.. autofunction:: coot.do_tw

PHENIX Support
--------------

.. autofunction:: coot.set_button_label_for_external_refinement

Graphics Text
-------------

.. autofunction:: coot.place_text

.. autofunction:: coot.remove_text

.. autofunction:: coot.edit_text

.. autofunction:: coot.text_index_near_position

PISA Interaction
----------------

.. autofunction:: coot.pisa_interaction

Jiggle Fit
----------

.. autofunction:: coot.fit_to_map_by_random_jiggle

.. autofunction:: coot.fit_molecule_to_map_by_random_jiggle

.. autofunction:: coot.fit_molecule_to_map_by_random_jiggle_and_blur

.. autofunction:: coot.fit_chain_to_map_by_random_jiggle

.. autofunction:: coot.fit_chain_to_map_by_random_jiggle_and_blur

SBase interface
---------------

.. autofunction:: coot.get_ccp4srs_monomer_and_dictionary

.. autofunction:: coot.get_sbase_monomer

FLE
---

.. autofunction:: coot.fle_view

.. autofunction:: coot.fle_view_internal

.. autofunction:: coot.fle_view_internal_to_png

.. autofunction:: coot.fle_view_with_rdkit

.. autofunction:: coot.fle_view_with_rdkit_to_png

.. autofunction:: coot.fle_view_with_rdkit_to_svg

.. autofunction:: coot.fle_view_with_rdkit_internal

.. autofunction:: coot.fle_view_set_water_dist_max

.. autofunction:: coot.fle_view_set_h_bond_dist_max

.. autofunction:: coot.sprout_hydrogens

LSQ
---

.. autofunction:: coot.lsq_improve

ingle
-----

.. autofunction:: coot.single_model_view_model_number

.. autofunction:: coot.single_model_view_this_model_number

.. autofunction:: coot.single_model_view_next_model_number

.. autofunction:: coot.single_model_view_prev_model_number

graphics 2D ligand view
-----------------------

.. autofunction:: coot.set_show_graphics_ligand_view

Experimental
------------

.. autofunction:: coot.fetch_and_superpose_alphafold_models_using_active_molecule

.. autofunction:: coot.start_ligand_builder_gui

.. autofunction:: coot.globularize

Uncategorized
-------------

.. autofunction:: coot.try_load_python_extras_dir

.. autofunction:: coot.reverse_direction_of_fragment

.. autofunction:: coot.setup_reverse_direction

.. autofunction:: coot.set_axis_orientation_matrix

.. autofunction:: coot.set_axis_orientation_matrix_usage

.. autofunction:: coot.optimal_B_kurtosis

.. autofunction:: coot.add_linked_residue

.. autofunction:: coot.set_add_linked_residue_do_fit_and_refine

.. autofunction:: coot.toolbar_multi_refine_stop

.. autofunction:: coot.toolbar_multi_refine_continue

.. autofunction:: coot.toolbar_multi_refine_cancel

.. autofunction:: coot.set_visible_toolbar_multi_refine_stop_button

.. autofunction:: coot.set_visible_toolbar_multi_refine_continue_button

.. autofunction:: coot.set_visible_toolbar_multi_refine_cancel_button

.. autofunction:: coot.toolbar_multi_refine_button_set_sensitive

.. autofunction:: coot.load_tutorial_model_and_data

.. autofunction:: coot.run_update_self_maybe

.. autofunction:: coot.show_go_to_residue_keyboarding_mode_window

.. autofunction:: coot.handle_go_to_residue_keyboarding_mode

.. autofunction:: coot.git_commit

.. autofunction:: coot.filtered_by_glob

.. autofunction:: coot.string_member

.. autofunction:: coot.compare_strings

.. autofunction:: coot.add_cablam_markup

.. autofunction:: coot.add_cablam_markup_py

.. autofunction:: coot.set_rotation_centre

.. autofunction:: coot.goto_next_atom_maybe_py

.. autofunction:: coot.goto_prev_atom_maybe_py

.. autofunction:: coot.set_go_to_atom_from_spec

.. autofunction:: coot.set_go_to_atom_from_res_spec

.. autofunction:: coot.set_go_to_atom_from_res_spec_py

.. autofunction:: coot.set_go_to_atom_from_atom_spec_py

.. autofunction:: coot.active_atom_spec

.. autofunction:: coot.active_atom_spec_py

.. autofunction:: coot.get_symmetry_py

.. autofunction:: coot.clashes_with_symmetry

.. autofunction:: coot.add_molecular_symmetry

.. autofunction:: coot.add_molecular_symmetry_from_mtrix_from_file

.. autofunction:: coot.add_molecular_symmetry_from_mtrix_from_self_file

.. autofunction:: coot.regen_map_internal

.. autofunction:: coot.make_weighted_map_simple_internal

.. autofunction:: coot.colour_map_by_other_map

.. autofunction:: coot.colour_map_by_other_map_py

.. autofunction:: coot.export_molecule_as_x3d

.. autofunction:: coot.export_molecule_as_obj

.. autofunction:: coot.export_molecule_as_gltf

.. autofunction:: coot.colour_map_by_other_map_turn_off

.. autofunction:: coot.add_density_map_cap

.. autofunction:: coot.recolour_mesh_by_map

.. autofunction:: coot.import_bild

.. autofunction:: coot.servalcat_fofc

.. autofunction:: coot.servalcat_refine

.. autofunction:: coot.run_acedrg_link_generation

.. autofunction:: coot.add_toolbar_subprocess_button

.. autofunction:: coot.compare_mtimes

.. autofunction:: coot.parse_ccp4i_defs

.. autofunction:: coot.ccp4_project_directory

.. autofunction:: coot.add_to_history

.. autofunction:: coot.add_to_history_simple

.. autofunction:: coot.add_to_history_typed

.. autofunction:: coot.single_quote

.. autofunction:: coot.pythonize_command_name

.. autofunction:: coot.schemize_command_name

.. autofunction:: coot.languagize_command

.. autofunction:: coot.add_to_database

.. autofunction:: coot.merge_molecules_by_vector

.. autofunction:: coot.monomer_restraints_py

.. autofunction:: coot.monomer_restraints_for_molecule_py

.. autofunction:: coot.set_monomer_restraints_py

.. autofunction:: coot.show_restraints_editor

.. autofunction:: coot.show_restraints_editor_by_index

.. autofunction:: coot.write_restraints_cif_dictionary

.. autofunction:: coot.list_nomenclature_errors

.. autofunction:: coot.list_nomenclature_errors_py

.. autofunction:: coot.show_fix_nomenclature_errors_gui

.. autofunction:: coot.dipole_to_py

.. autofunction:: coot.coot_has_guile

.. autofunction:: coot.coot_can_do_lidia_p

.. autofunction:: coot.run_scheme_command

.. autofunction:: coot.pyrun_simple_string

.. autofunction:: coot.residue_spec_to_py

.. autofunction:: coot.residue_spec_make_triple_py

.. autofunction:: coot.residue_spec_from_py

.. autofunction:: coot.get_residue_by_type

.. autofunction:: coot.get_residue_specs_in_mol

.. autofunction:: coot.get_residue_specs_in_mol_py

.. autofunction:: coot.resname_from_serial_number

.. autofunction:: coot.residue_name

.. autofunction:: coot.serial_number_from_residue_specs

.. autofunction:: coot.hydrogenate_region

.. autofunction:: coot.add_hydrogens_from_file

.. autofunction:: coot.add_hydrogen_atoms_to_residue

.. autofunction:: coot.add_hydrogen_atoms_to_residue_py

.. autofunction:: coot.atom_info_string_py

.. autofunction:: coot.molecule_to_pdb_string_py

.. autofunction:: coot.residue_info_py

.. autofunction:: coot.residue_name_py

.. autofunction:: coot.residue_centre_from_spec_py

.. autofunction:: coot.chain_fragments_py

.. autofunction:: coot.set_b_factor_residues_py

.. autofunction:: coot.rigid_body_fit_with_residue_ranges

.. autofunction:: coot.morph_fit_all

.. autofunction:: coot.morph_fit_chain

.. autofunction:: coot.morph_fit_residues_py

.. autofunction:: coot.morph_fit_residues

.. autofunction:: coot.morph_fit_by_secondary_structure_elements

.. autofunction:: coot.check_waters_baddies

.. autofunction:: coot.find_blobs

.. autofunction:: coot.find_blobs_py

.. autofunction:: coot.b_factor_distribution_graph

.. autofunction:: coot.monomer_lib_3_letter_codes_matching

.. autofunction:: coot.mutate_residue_range

.. autofunction:: coot.mutate_internal

.. autofunction:: coot.mutate_active_residue_to_single_letter_code

.. autofunction:: coot.show_keyboard_mutate_dialog

.. autofunction:: coot.mutate_by_overlap

.. autofunction:: coot.overlap_ligands_internal

.. autofunction:: coot.do_smiles_to_simple_3d_overlay_frame

.. autofunction:: coot.ligand_search_make_conformers_py

.. autofunction:: coot.ligand_search_make_conformers_internal

.. autofunction:: coot.add_animated_ligand_interaction

.. autofunction:: coot.cootaneer_internal

.. autofunction:: coot.cootaneer_py

.. autofunction:: coot.is_interesting_dots_object_next_p

.. autofunction:: coot.generic_string_vector_to_list_internal_py

.. autofunction:: coot.generic_int_vector_to_list_internal_py

.. autofunction:: coot.generic_list_to_string_vector_internal_py

.. autofunction:: coot.rtop_to_python

.. autofunction:: coot.inverse_rtop_py

.. autofunction:: coot.atom_spec_from_python_expression

.. autofunction:: coot.atom_spec_to_py

.. autofunction:: coot.set_display_control_button_state

.. autofunction:: coot.fullscreen

.. autofunction:: coot.unfullscreen

.. autofunction:: coot.set_use_trackpad

.. autofunction:: coot.set_use_primary_mouse_button_for_view_rotation

.. autofunction:: coot.get_coords_for_accession_code

.. autofunction:: coot.coot_get_url_as_string_internal

.. autofunction:: coot.stop_curl_download

.. autofunction:: coot.get_drug_mdl_via_wikipedia_and_drugbank

.. autofunction:: coot.fetch_and_superpose_alphafold_models

.. autofunction:: coot.fetch_alphafold_model_for_uniprot_id

.. autofunction:: coot.fetch_emdb_map

.. autofunction:: coot.fetch_cod_entry

.. autofunction:: coot.orient_view

.. autofunction:: coot.topological_equivalence_chiral_centres

.. autofunction:: coot.screendump_tga

.. autofunction:: coot.set_framebuffer_scale_factor

.. autofunction:: coot.set_use_perspective_projection

.. autofunction:: coot.use_perspective_projection_state

.. autofunction:: coot.set_perspective_fov

.. autofunction:: coot.set_use_ambient_occlusion

.. autofunction:: coot.use_ambient_occlusion_state

.. autofunction:: coot.set_use_depth_blur

.. autofunction:: coot.use_depth_blur_state

.. autofunction:: coot.set_use_fog

.. autofunction:: coot.use_fog_state

.. autofunction:: coot.set_use_outline

.. autofunction:: coot.use_outline_state

.. autofunction:: coot.set_map_shininess

.. autofunction:: coot.set_map_specular_strength

.. autofunction:: coot.set_draw_normals

.. autofunction:: coot.draw_normals_state

.. autofunction:: coot.set_draw_mesh

.. autofunction:: coot.draw_mesh_state

.. autofunction:: coot.set_map_material_specular

.. autofunction:: coot.set_model_material_specular

.. autofunction:: coot.set_model_material_ambient

.. autofunction:: coot.set_model_material_diffuse

.. autofunction:: coot.set_model_goodselliness

.. autofunction:: coot.set_map_fresnel_settings

.. autofunction:: coot.reload_map_shader

.. autofunction:: coot.reload_model_shader

.. autofunction:: coot.set_atom_radius_scale_factor

.. autofunction:: coot.set_use_fancy_lighting

.. autofunction:: coot.set_use_simple_lines_for_model_molecules

.. autofunction:: coot.set_fresnel_colour

.. autofunction:: coot.set_focus_blur_z_depth

.. autofunction:: coot.set_focus_blur_strength

.. autofunction:: coot.set_shadow_strength

.. autofunction:: coot.set_shadow_resolution

.. autofunction:: coot.set_shadow_box_size

.. autofunction:: coot.set_ssao_kernel_n_samples

.. autofunction:: coot.set_ssao_strength

.. autofunction:: coot.set_ssao_radius

.. autofunction:: coot.set_ssao_bias

.. autofunction:: coot.set_ssao_blur_size

.. autofunction:: coot.set_shadow_softness

.. autofunction:: coot.set_shadow_texture_resolution_multiplier

.. autofunction:: coot.set_effects_shader_output_type

.. autofunction:: coot.set_effects_shader_brightness

.. autofunction:: coot.set_effects_shader_gamma

.. autofunction:: coot.set_bond_smoothness_factor

.. autofunction:: coot.set_draw_gl_ramachandran_plot_during_refinement

.. autofunction:: coot.set_fps_timing_scale_factor

.. autofunction:: coot.set_draw_background_image

.. autofunction:: coot.read_test_gltf_models

.. autofunction:: coot.load_gltf_model

.. autofunction:: coot.scale_model

.. autofunction:: coot.reset_framebuffers

Extra Map Functions
-------------------

.. autofunction:: coot.auto_read_make_and_draw_maps

.. autofunction:: coot.auto_read_make_and_draw_maps_from_mtz

.. autofunction:: coot.auto_read_make_and_draw_maps_from_cns

.. autofunction:: coot.valid_labels

.. autofunction:: coot.sharpen_blur_map

.. autofunction:: coot.sharpen_blur_map_with_resampling

.. autofunction:: coot.sharpen_blur_map_with_resampling_threaded_version

.. autofunction:: coot.multi_sharpen_blur_map_py

.. autofunction:: coot.amplitude_vs_resolution_py

.. autofunction:: coot.flip_hand

.. autofunction:: coot.analyse_map_point_density_change

.. autofunction:: coot.analyse_map_point_density_change_py

.. autofunction:: coot.go_to_map_molecule_centre

.. autofunction:: coot.b_factor_from_map

.. autofunction:: coot.map_colour_components_py

.. autofunction:: coot.read_ccp4_map

.. autofunction:: coot.handle_read_ccp4_map

.. autofunction:: coot.handle_read_emdb_data

.. autofunction:: coot.show_map_partition_by_chain_dialog

.. autofunction:: coot.map_partition_by_chain

.. autofunction:: coot.map_partition_by_chain_threaded

.. autofunction:: coot.set_use_vertex_gradients_for_map_normals

.. autofunction:: coot.use_vertex_gradients_for_map_normals_for_latest_map

.. autofunction:: coot.set_use_vertex_gradients_for_map_normals_for_latest_map

Multi
-----

.. autofunction:: coot.multi_residue_torsion_fit

.. autofunction:: coot.multi_residue_torsion_fit_py

Add an Atom
-----------

.. autofunction:: coot.add_an_atom

Nudge the B
-----------

.. autofunction:: coot.nudge_the_temperature_factors_py

Merge Fragments
---------------

.. autofunction:: coot.merge_fragments

Delete Items
------------

.. autofunction:: coot.delete_chain

.. autofunction:: coot.delete_sidechains_for_chain

Execute Refmac
--------------

.. autofunction:: coot.execute_refmac_real

.. autofunction:: coot.refmac_name

Dictionary Functions
--------------------

.. autofunction:: coot.handle_cif_dictionary

.. autofunction:: coot.read_cif_dictionary

.. autofunction:: coot.handle_cif_dictionary_for_molecule

.. autofunction:: coot.dictionary_entries

.. autofunction:: coot.debug_dictionary

.. autofunction:: coot.SMILES_for_comp_id

.. autofunction:: coot.dictionaries_read_py

.. autofunction:: coot.cif_file_for_comp_id_py

.. autofunction:: coot.dictionary_entries_py

.. autofunction:: coot.SMILES_for_comp_id_py

Using S
-------

.. autofunction:: coot.clear_and_update_molecule_py

.. autofunction:: coot.add_molecule_py

.. autofunction:: coot.active_residue_py

.. autofunction:: coot.closest_atom_simple_py

.. autofunction:: coot.closest_atom_py

.. autofunction:: coot.closest_atom_raw_py

.. autofunction:: coot.residues_near_residue_py

.. autofunction:: coot.residues_near_residues_py

.. autofunction:: coot.residues_near_position_py

.. autofunction:: coot.label_closest_atoms_in_neighbour_residues_py

.. autofunction:: coot.get_bonds_representation

.. autofunction:: coot.set_new_non_drawn_bonds

.. autofunction:: coot.add_to_non_drawn_bonds

.. autofunction:: coot.clear_non_drawn_bonds

.. autofunction:: coot.get_dictionary_radii

.. autofunction:: coot.get_environment_distances_representation_py

.. autofunction:: coot.get_intermediate_atoms_bonds_representation

.. autofunction:: coot.get_continue_updating_refinement_atoms_state

tatus bar string functions
--------------------------

.. autofunction:: coot.atom_info_as_text_for_statusbar

.. autofunction:: coot.atom_info_as_text_for_statusbar

Refinement with specs
---------------------

.. autofunction:: coot.all_residues_with_serial_numbers_py

.. autofunction:: coot.regularize_residues

.. autofunction:: coot.mtz_file_name

.. autofunction:: coot.refine_zone_with_full_residue_spec_py

.. autofunction:: coot.set_draw_moving_atoms_rota_markup

.. autofunction:: coot.set_draw_moving_atoms_rama_markup

.. autofunction:: coot.set_show_intermediate_atoms_rota_markup

.. autofunction:: coot.set_show_intermediate_atoms_rama_markup

.. autofunction:: coot.get_draw_moving_atoms_rota_markup_state

.. autofunction:: coot.get_draw_moving_atoms_rama_markup_state

.. autofunction:: coot.get_show_intermediate_atoms_rota_markup

.. autofunction:: coot.get_show_intermediate_atoms_rama_markup

.. autofunction:: coot.set_cryo_em_refinement

.. autofunction:: coot.get_cryo_em_refinement

.. autofunction:: coot.accept_moving_atoms_py

.. autofunction:: coot.register_post_intermediate_atoms_moved_hook

.. autofunction:: coot.set_regenerate_bonds_needs_make_bonds_type_checked

.. autofunction:: coot.get_regenerate_bonds_needs_make_bonds_type_checked_state

Water Chain Functions
---------------------

.. autofunction:: coot.water_chain_from_shelx_ins_py

.. autofunction:: coot.water_chain_py

Interface Utils
---------------

.. autofunction:: coot.add_status_bar_text

.. autofunction:: coot.set_logging_level

Glyco Tools
-----------

.. autofunction:: coot.print_glyco_tree

Variance Map
------------

.. autofunction:: coot.make_variance_map

.. autofunction:: coot.make_variance_map_py

Spin Search Functions
---------------------

.. autofunction:: coot.spin_search_by_atom_vectors

.. autofunction:: coot.spin_search_py

.. autofunction:: coot.spin_N_py

.. autofunction:: coot.CG_spin_search_py

Sequence from Map
-----------------

.. autofunction:: coot.sequence_from_map

.. autofunction:: coot.apply_sequence_to_fragment

.. autofunction:: coot.assign_sequence_to_active_fragment

Rotamer Scoring
---------------

.. autofunction:: coot.score_rotamers

.. autofunction:: coot.score_rotamers_py

protein
-------

.. autofunction:: coot.protein_db_loops

.. autofunction:: coot.protein_db_loop_specs_to_atom_selection_string

.. autofunction:: coot.protein_db_loops_py

Coot's Hole implementation
--------------------------

.. autofunction:: coot.hole

Coot's Gaussian Surface
-----------------------

.. autofunction:: coot.gaussian_surface

.. autofunction:: coot.set_gaussian_surface_sigma

.. autofunction:: coot.set_gaussian_surface_contour_level

.. autofunction:: coot.set_gaussian_surface_box_radius

.. autofunction:: coot.set_gaussian_surface_grid_scale

.. autofunction:: coot.set_gaussian_surface_fft_b_factor

.. autofunction:: coot.set_gaussian_surface_chain_colour_mode

.. autofunction:: coot.show_gaussian_surface_overlay

.. autofunction:: coot.make_acedrg_dictionary_via_CCD_dictionary

.. autofunction:: coot.make_link

.. autofunction:: coot.make_link_py

.. autofunction:: coot.link_info_py

.. autofunction:: coot.show_acedrg_link_interface_overlay

Drag and Drop Functions
-----------------------

.. autofunction:: coot.handle_drag_and_drop_string

Map Display Control
-------------------

.. autofunction:: coot.undisplay_all_maps_except

Map Contouring Functions
------------------------

.. autofunction:: coot.map_contours

.. autofunction:: coot.map_contours_as_triangles

.. autofunction:: coot.set_radial_map_colouring_enabled

.. autofunction:: coot.set_radial_map_colouring_centre

.. autofunction:: coot.set_radial_map_colouring_min_radius

.. autofunction:: coot.set_radial_map_colouring_max_radius

.. autofunction:: coot.set_radial_map_colouring_invert

.. autofunction:: coot.set_radial_map_colouring_saturation

Map to Model Correlation
------------------------

.. autofunction:: coot.set_map_correlation_atom_radius

.. autofunction:: coot.map_to_model_correlation_py

.. autofunction:: coot.map_to_model_correlation_stats_py

.. autofunction:: coot.map_to_model_correlation_stats_per_residue_range_py

.. autofunction:: coot.map_to_model_correlation

.. autofunction:: coot.map_to_model_correlation_stats

.. autofunction:: coot.map_to_model_correlation_per_residue

.. autofunction:: coot.map_to_model_correlation_stats_per_residue

.. autofunction:: coot.map_to_model_correlation_stats_per_residue_range

.. autofunction:: coot.map_to_model_correlation_per_residue_py

.. autofunction:: coot.qq_plot_map_and_model_py

.. autofunction:: coot.density_score_residue

.. autofunction:: coot.map_mean_py

.. autofunction:: coot.map_sigma_py

.. autofunction:: coot.map_statistics_py

Get Sequence
------------

.. autofunction:: coot.get_sequence_as_fasta_for_chain

.. autofunction:: coot.write_sequence

.. autofunction:: coot.res_tracer

.. autofunction:: coot.register_interesting_positions_list_py

.. autofunction:: coot.molecule_atom_overlaps_py

.. autofunction:: coot.prodrg_import_function

.. autofunction:: coot.sbase_import_function

.. autofunction:: coot.align_to_closest_chain

.. autofunction:: coot.resolve_clashing_sidechains_by_deletion

.. autofunction:: coot.resolve_clashing_sidechains_by_rebuilding

.. autofunction:: coot.simple_text_dialog

.. autofunction:: coot.graphics_to_phenix_geo_representation

.. autofunction:: coot.graphics_to_phenix_geo_representation

.. autofunction:: coot.set_python_draw_function

.. autofunction:: coot.pathology_data

.. autofunction:: coot.encode_ints

.. autofunction:: coot.decode_ints

.. autofunction:: coot.store_keyed_user_name

.. autofunction:: coot.py_to_residue_specs

.. autofunction:: coot.key_sym_code_py

.. autofunction:: coot.py_symop_strings_to_space_group

.. autofunction:: coot.set_use_sounds

.. autofunction:: coot.curmudgeon_mode

.. autofunction:: coot.halloween

.. autofunction:: coot.display_svg_from_file_in_a_dialog

.. autofunction:: coot.display_svg_from_string_in_a_dialog

.. autofunction:: coot.display_pae_from_file_in_a_dialog

.. autofunction:: coot.read_interesting_places_json_file

.. autofunction:: coot.setup_tomo_slider

.. autofunction:: coot.tomo_section_view

.. autofunction:: coot.set_tomo_section_view_section

.. autofunction:: coot.set_tomo_picker_mode_is_active

.. autofunction:: coot.tomo_map_analysis

.. autofunction:: coot.tomo_map_analysis_2

.. autofunction:: coot.reverse_map

.. autofunction:: coot.read_positron_metadata

.. autofunction:: coot.positron_pathway

.. autofunction:: coot.positron_plot_py

.. autofunction:: coot.positron_plot_internal

