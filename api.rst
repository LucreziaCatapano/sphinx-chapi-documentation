Python API
==========================================

.. currentmodule:: chapi

Basics Utilities
----------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: set_make_backups
    .. automethod:: get_make_backups
    .. automethod:: contains_unsaved_models
    .. automethod:: save_unsaved_model_changes
    .. automethod:: set_show_timings
    .. automethod:: set_use_gemmi
    .. automethod:: get_use_gemmi
    .. automethod:: undo
    .. automethod:: redo


Reading and Writing
----------------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: read_coordinates
    .. automethod:: read_pdb
    .. automethod:: read_small_molecule_cif
    .. automethod:: read_mtz
    .. automethod:: auto_read_mtz
    .. automethod:: read_ccp4_map
    .. automethod:: write_coordinates
    .. automethod:: write_map


Molecular Information
---------------------------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: get_header_info
    .. automethod:: get_number_of_molecules
    .. automethod:: print_secondary_structure_info
    .. automethod:: get_molecule_name
    .. automethod:: get_molecule_centre
    .. automethod:: get_molecule_diameter
    .. automethod:: get_chains_in_model
    .. automethod:: get_ncs_related_chains
    .. automethod:: get_single_letter_codes_for_chain
    .. automethod:: get_residue_name
    .. automethod:: get_residue_names_with_no_dictionary
    .. automethod:: get_residues_near_residue
    .. automethod:: get_torsion
    .. automethod:: get_residue_CA_position
    .. automethod:: get_residue_average_position
    .. automethod:: get_residue_sidechain_average_position
    .. automethod:: residues_with_missing_atoms
    .. automethod:: get_missing_residue_ranges
    .. automethod:: get_number_of_atoms
    .. automethod:: get_number_of_hydrogen_atoms
    .. automethod:: get_active_atom
    .. automethod:: get_imol_enc_any
    .. automethod:: get_hb_type
    .. automethod:: get_cell
    .. automethod:: get_symmetry
    .. automethod:: get_mutation_info


Geometry and Dictionaries
----------------------------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: geometry_init_standard
    .. automethod:: non_standard_residue_types_in_model
    .. automethod:: get_non_standard_residues_in_molecule
    .. automethod:: import_cif_dictionary
    .. automethod:: get_cif_file_name
    .. automethod:: get_cif_restraints_as_string
    .. automethod:: copy_dictionary
    .. automethod:: get_monomer
    .. automethod:: get_monomer_from_dictionary
    .. automethod:: get_monomer_and_position_at
    .. automethod:: dictionary_atom_name_map
    .. automethod:: get_group_for_monomer
    .. automethod:: get_groups_for_monomers
    .. automethod:: get_gphl_chem_comp_info
    .. automethod:: get_acedrg_atom_types
    .. automethod:: get_acedrg_atom_types_for_ligand
    .. automethod:: get_dictionary_conformers
    .. automethod:: get_SMILES_for_residue_type


Model Manipulation
----------------------------------------

.. class:: lsq_results_t
   :noindex:

.. class:: molecules_container_t
   :noindex:

    .. automethod:: replace_molecule_by_model_from_file
    .. automethod:: split_multi_model_molecule
    .. automethod:: make_ensemble
    .. automethod:: change_chain_id
    .. automethod:: associate_sequence
    .. automethod:: assign_sequence
    .. automethod:: molecule_to_PDB_string
    .. automethod:: molecule_to_mmCIF_string
    .. automethod:: file_name_to_string
    .. automethod:: create_empty_molecules
    .. automethod:: new_molecule
    .. automethod:: close_molecule
    .. automethod:: end_delete_closed_molecules
    .. automethod:: pop_back
    .. automethod:: clear
    .. automethod:: get_eigenvalues
    .. automethod:: atom_cid_to_atom_spec
    .. automethod:: residue_cid_to_residue_spec
    .. automethod:: set_molecule_name
    .. automethod:: display_molecule_names_table
    .. automethod:: is_valid_model_molecule
    .. automethod:: change_to_next_rotamer
    .. automethod:: change_to_previous_rotamer
    .. automethod:: change_to_first_rotamer
    .. automethod:: add_compound
    .. automethod:: delete_using_cid
    .. automethod:: delete_atom
    .. automethod:: delete_atom_using_cid
    .. automethod:: delete_residue
    .. automethod:: delete_residue_using_cid
    .. automethod:: delete_residue_atoms_using_cid
    .. automethod:: delete_residue_atoms_with_alt_conf
    .. automethod:: delete_side_chain
    .. automethod:: delete_side_chain_using_cid
    .. automethod:: delete_chain_using_cid
    .. automethod:: delete_literal_using_cid
    .. automethod:: add_terminal_residue_directly
    .. automethod:: add_terminal_residue_directly_using_cid
    .. automethod:: set_add_waters_water_to_protein_distance_lim_min
    .. automethod:: set_add_waters_water_to_protein_distance_lim_max
    .. automethod:: set_add_waters_variance_limit
    .. automethod:: set_add_waters_sigma_cutoff
    .. automethod:: add_waters
    .. automethod:: flood
    .. automethod:: add_hydrogen_atoms
    .. automethod:: delete_hydrogen_atoms
    .. automethod:: add_alternative_conformation
    .. automethod:: change_alt_locs
    .. automethod:: split_residue_using_map
    .. automethod:: fill_partial_residue
    .. automethod:: fill_partial_residues
    .. automethod:: fill_partial_residue_using_cid
    .. automethod:: flip_peptide
    .. automethod:: flip_peptide_using_cid
    .. automethod:: eigen_flip_ligand
    .. automethod:: eigen_flip_ligand_using_cid
    .. automethod:: mutate
    .. automethod:: side_chain_180
    .. automethod:: jed_flip
    .. automethod:: move_molecule_to_new_centre
    .. automethod:: copy_molecule
    .. automethod:: copy_fragment_using_cid
    .. automethod:: copy_fragment_for_refinement_using_cid
    .. automethod:: copy_fragment_using_residue_range
    .. automethod:: apply_transformation_to_atom_selection
    .. automethod:: new_positions_for_residue_atoms
    .. automethod:: new_positions_for_atoms_in_residues
    .. automethod:: merge_molecules
    .. automethod:: cis_trans_convert
    .. automethod:: replace_residue
    .. automethod:: replace_fragment
    .. automethod:: SSM_superpose
    .. automethod:: add_lsq_superpose_match
    .. automethod:: add_lsq_superpose_atom_match
    .. automethod:: clear_lsq_matches
    .. automethod:: lsq_superpose
    .. automethod:: get_lsq_matrix
    .. automethod:: transform_map_using_lsq_matrix
    .. automethod:: rotate_around_bond
    

Map Tools
-----------------------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: set_imol_refinement_map
    .. automethod:: is_valid_map_molecule
    .. automethod:: is_a_difference_map
    .. automethod:: get_map_weight
    .. automethod:: set_map_weight
    .. automethod:: get_map_molecule_centre
    .. automethod:: get_map_sampling_rate
    .. automethod:: set_map_sampling_rate
    .. automethod:: get_density_at_position
    .. automethod:: get_diff_diff_map_peaks
    .. automethod:: get_data_set_file_name
    .. automethod:: go_to_blob
    .. automethod:: calculate_new_rail_points
    .. automethod:: rail_points_total
    .. automethod:: replace_map_by_mtz_from_file
    .. automethod:: get_map_mean
    .. automethod:: get_map_rmsd_approx
    .. automethod:: get_map_histogram
    .. automethod:: get_suggested_initial_contour_level
    .. automethod:: is_EM_map
    .. automethod:: sharpen_blur_map
    .. automethod:: sharpen_blur_map_with_resample
    .. automethod:: mask_map_by_atom_selection
    .. automethod:: partition_map_by_chain
    .. automethod:: make_mask
    .. automethod:: flip_hand
    .. automethod:: make_masked_maps_split_by_chain
    .. automethod:: average_map
    .. automethod:: regen_map
    .. automethod:: make_power_scaled_map
    .. automethod:: associate_data_mtz_file_with_map
    .. automethod:: connect_updating_maps
    .. automethod:: fourier_shell_correlation
    .. automethod:: get_q_score
    .. automethod:: get_q_score_for_cid
    .. automethod:: get_map_section_texture
    .. automethod:: get_number_of_map_sections
    .. automethod:: set_map_colour
    .. automethod:: set_map_is_contoured_with_thread_pool
    .. automethod:: get_map_contours_mesh
    .. automethod:: get_map_contours_mesh_using_other_map_for_colours
    .. automethod:: set_map_colour_saturation


Structure Factor
----------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: get_latest_sfcalc_stats
    .. automethod:: get_r_factor_stats
    .. automethod:: r_factor_stats_as_string
    .. automethod:: sfcalc_genmap
    .. automethod:: sfcalc_genmaps_using_bulk_solvent


Real Space Refinement
-----------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: multiply_residue_temperature_factors
    .. automethod:: get_median_temperature_factor
    .. automethod:: set_temperature_factors_using_cid
    .. automethod:: set_occupancy
    .. automethod:: shift_field_b_factor_refinement
    .. automethod:: refine
    .. automethod:: adjust_refinement_residue_name
    .. automethod:: refine_residues_using_atom_cid
    .. automethod:: refine_residues
    .. automethod:: refine_residue_range
    .. automethod:: minimize_energy
    .. automethod:: fix_atom_selection_during_refinement
    .. automethod:: add_target_position_restraint
    .. automethod:: clear_target_position_restraint
    .. automethod:: turn_off_when_close_target_position_restraint
    .. automethod:: get_use_rama_plot_restraints
    .. automethod:: set_use_rama_plot_restraints
    .. automethod:: get_rama_plot_restraints_weight
    .. automethod:: set_rama_plot_restraints_weight
    .. automethod:: get_use_torsion_restraints
    .. automethod:: set_use_torsion_restraints
    .. automethod:: get_torsion_restraints_weight
    .. automethod:: set_torsion_restraints_weight
    .. automethod:: init_refinement_of_molecule_as_fragment_based_on_reference
    .. automethod:: add_target_position_restraint_and_refine
    .. automethod:: clear_target_position_restraints
    .. automethod:: clear_refinement
    .. automethod:: set_refinement_is_verbose
    .. automethod:: get_geman_mcclure_alpha
    .. automethod:: set_refinement_geman_mcclure_alpha
    .. automethod:: generate_self_restraints
    .. automethod:: generate_chain_self_restraints
    .. automethod:: generate_local_self_restraints
    .. automethod:: add_parallel_plane_restraint
    .. automethod:: get_extra_restraints_mesh
    .. automethod:: read_extra_restraints
    .. automethod:: clear_extra_restraints


Servalcat
-----------------------------------

.. class:: molecules_container_t
   :noindex:
   
   .. automethod:: servalcat_refine_xray


Fitting
---------------------------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: auto_fit_rotamer
    .. automethod:: rigid_body_fit
    .. automethod:: fit_ligand_right_here
    .. automethod:: fit_ligand
    .. automethod:: fit_ligand_multi_ligand
    .. automethod:: fit_to_map_by_random_jiggle
    .. automethod:: fit_to_map_by_random_jiggle_using_cid
    .. automethod:: fit_to_map_by_random_jiggle_with_blur_using_cid


Validation
--------------------------------------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: get_rotamer_dodecs
    .. automethod:: get_rotamer_dodecs_instanced
    .. automethod:: fill_rotamer_probability_tables
    .. automethod:: accept_rotamer_probability_tables_compressed_data
    .. automethod:: get_ramachandran_validation_markup_mesh
    .. automethod:: ramachandran_validation
    .. automethod:: contact_dots_for_ligand
    .. automethod:: all_molecule_contact_dots
    .. automethod:: get_simple_molecule
    .. automethod:: make_exportable_environment_bond_box
    .. automethod:: get_h_bonds
    .. automethod:: get_mesh_for_ligand_validation_vs_dictionary
    .. automethod:: get_ligand_validation_vs_dictionary
    .. automethod:: get_validation_vs_dictionary_for_selection
    .. automethod:: match_ligand_torsions
    .. automethod:: get_ligand_distortion
    .. automethod:: match_ligand_torsions_and_position
    .. automethod:: match_ligand_torsions_and_position_using_cid
    .. automethod:: get_overlap_dots
    .. automethod:: get_overlap_dots_for_ligand
    .. automethod:: get_overlaps
    .. automethod:: get_overlaps_for_ligand
    .. automethod:: density_fit_analysis
    .. automethod:: density_correlation_analysis
    .. automethod:: rotamer_analysis
    .. automethod:: ramachandran_analysis
    .. automethod:: ramachandran_analysis_for_chain
    .. automethod:: peptide_omega_analysis
    .. automethod:: get_interesting_places
    .. automethod:: difference_map_peaks
    .. automethod:: pepflips_using_difference_map
    .. automethod:: unmodelled_blobs
    .. automethod:: find_water_baddies
    


Molecular Graphics Representation
----------------------------------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: set_draw_missing_residue_loops
    .. automethod:: get_bonds_mesh
    .. automethod:: get_bonds_mesh_instanced
    .. automethod:: get_bonds_mesh_for_selection_instanced
    .. automethod:: get_goodsell_style_mesh_instanced
    .. automethod:: get_molecular_representation_mesh
    .. automethod:: get_gaussian_surface
    .. automethod:: get_chemical_features_mesh
    .. automethod:: export_map_molecule_as_gltf
    .. automethod:: export_model_molecule_as_gltf
    .. automethod:: export_molecular_representation_as_gltf
    .. automethod:: make_mesh_from_gltf_file
    .. automethod:: get_colour_table
    .. automethod:: set_colour_wheel_rotation_base
    .. automethod:: set_base_colour_for_bonds
    .. automethod:: add_to_non_drawn_bonds
    .. automethod:: clear_non_drawn_bonds
    .. automethod:: print_non_drawn_bonds
    .. automethod:: set_user_defined_bond_colours
    .. automethod:: set_user_defined_atom_colour_by_selection
    .. automethod:: get_colour_rules
    .. automethod:: print_colour_rules
    .. automethod:: add_colour_rule
    .. automethod:: add_colour_rules_multi
    .. automethod:: delete_colour_rules
    .. automethod:: set_use_bespoke_carbon_atom_colour
    .. automethod:: set_bespoke_carbon_atom_colour
    .. automethod:: M2T_updateFloatParameter
    .. automethod:: M2T_updateIntParameter
    .. automethod:: get_octahemisphere
    .. automethod:: get_svg_for_residue_type
    .. automethod:: write_png
    .. automethod:: pae_png


Testing functions
-----------------------------------

.. class:: molecules_container_t
   :noindex:

    .. automethod:: testing_start_long_term_job
    .. automethod:: testing_stop_long_term_job
    .. automethod:: testing_interrogate_long_term_job
    .. automethod:: get_contouring_time
    .. automethod:: set_max_number_of_threads
    .. automethod:: set_max_number_of_threads_in_thread_pool
    .. automethod:: test_the_threading
    .. automethod:: test_launching_threads
    .. automethod:: test_thread_pool_threads
    .. automethod:: mmcif_tests
    .. automethod:: test_origin_cube


Blender functions
--------------------------------------
These functions are for BlendCoot


.. class:: molecules_container_t
   :noindex:

    .. automethod:: make_mesh_for_map_contours_for_blender
    .. automethod:: make_mesh_for_bonds_for_blender
    .. automethod:: make_mesh_for_molecular_representation_for_blender
    .. automethod:: make_mesh_for_gaussian_surface_for_blender
    .. automethod:: make_mesh_for_goodsell_style_for_blender
    .. automethod:: get_colour_table_for_blender
    .. automethod:: get_vertices_for_blender
    .. automethod:: get_triangles_for_blender



