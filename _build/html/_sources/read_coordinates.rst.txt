Reading Files
================================================================

.. _read_coordinates:

Coordinate Files
------------------------------------------


Chapi supports the following coordinate file formats:

    * mmCIF (PDBx/mmCIF)
    * PDB
    * SHELXL
    * MOL file



The :code:`read_coordinates()` is a function provided by the :code:`molecules_container_t` class of the :code:`chapi` module, designed to read files in different format, including mmCIF or PDB format.
Additionally, the function :code:`read_pdb()` can be used as well to read PDB and mmCIF files.

.. code-block:: python

   import chapi

   mc = chapi.molecules_container_t(True)

   # read PDB file
   imol_pdb = mc.read_coordinates("rnase.pdb")

   # read mmCIF file
   imol_mmcif = mc.read_coordinates("rnase.cif")


MTZ and Map Files
------------------------


MTZ file format can be read by using the :code:`read_mtz()` function.

.. code-block:: python

   imol_mtz = mc.read_mtz("rnase-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", False, False)

The latest two arguments are:

* :code:`use_weight` (bool): the flag to use the map weight

* :code:`is_a_difference_map` (bool): the flag to set the map as a difference map

MRC/CCP4 Maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EM maps are typically in MRC/CCP4 map format

.. code-block:: python

   imol_map = mc.read_ccp4_map("emd_16890.map", False)

The latest argument is :code:`is_a_difference_map` (bool): the flag to set the map as a difference map


Writing Files
=========================================

The :code:`write_coordinates()` and :code:`write_map()` functions are used to write coordinates files (PDB and mmCIF) 
and map files (MTZ and map) respectively.


.. code-block:: python

   # read mmCIF file
   imol_mmcif = mc.read_coordinates("rnase.cif")

   # write mmCIF file
   imol_mmcif_new = mc.write_coordinates(imol_mmcif, "rnase-new.cif")

   # read MTZ file
   imol_mtz = mc.read_mtz("rnase-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", False, False)

   # write MTZ file
   imol_mtz_new = mc.write_map(imol_mtz, "rnase-1.8-all_refmac1_new.mtz")

   # read map file
   imol_map = mc.read_ccp4_map("emd_16890.map", False)

   # write map file
   imol_map_new = mc.write_map(imol_map, "emd_16890_new.map")


.. Molecular Models
.. ##################

Molecular Information
========================================


The following functions return information about macromolecular models. For more details see the **Chapi Python API** - `Molecular Information <https://www.mrc-lmb.cam.ac.uk/lucrezia/libcootapi-documentation/chapi_api.html#molecular-information>`_ section.

**General Information**

* :code:`get_molecule_diameter()`
* :code:`get_number_of_atoms()`
* :code:`get_number_of_hydrogen_atoms()`
* :code:`get_cell()`
* :code:`get_symmetry()`
* :code:`get_hb_type()`

**Chains**

* :code:`get_chains_in_model()`
* :code:`get_single_letter_codes_for_chain()`
* :code:`get_ncs_related_chains()`

**Residue**

* :code:`get_residue_name()`
* :code:`get_residue_names_with_no_dictionary()`
* :code:`residues_with_missing_atoms()`
* :code:`get_residues_near_residue()`

Examples
------------------------------------------------------

**1. General Information**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. doctest::
   
   >>> import chapi

   >>> mc = chapi.molecules_container_t(False)
   >>> mc.set_use_gemmi(False)

   >>> imol = mc.read_pdb("rnase.pdb")
   >>> molecule_diameter = mc.get_molecule_diameter(imol)
   >>> print("molecule diameter:", round(molecule_diameter))
   molecule diameter: 56

   >>> number_of_atoms = mc.get_number_of_atoms(imol)
   >>> print("number of atoms:", number_of_atoms)
   number of atoms: 1465

   >>> chains = mc.get_chains_in_model(imol)
   >>> print("chains:", chains)
   chains: ['A', 'B']

   >>> residue45_name = mc.get_residue_name(imol, 'A', 45, '')
   >>> print("residue 45 name:", residue45_name)
   residue 45 name: PRO 


**2. Header Info**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. doctest::

   >>> import chapi

   >>> mc = chapi.molecules_container_t(False)
   >>> mc.set_use_gemmi(False)

   >>> imol = mc.read_pdb("1mcy.pdb")
   >>> header_info = mc.get_header_info(imol)

   INFO:: There are 8 helices and 0 sheets
               Helix info: 
   ------------------------------------------------
   1 1 A 4 A 35 32 
   2 2 A 37 A 42 6 
   3 3 A 52 A 57 6 
   4 4 A 59 A 78 20 
   5 5 A 83 A 95 13 
   6 6 A 101 A 118 18 
   7 7 A 120 A 122 3 
   8 8 A 125 A 149 25 
                  Sheet info: 
   ------------------------------------------------
   ------------------------------------------------

   >>> compound_lines = header_info.compound_lines
   >>> journal_lines  = header_info.journal_lines
   >>> author_lines   = header_info.author_lines
   >>> for l in author_lines:
           print("author_line:", l)
   >>> for l in compound_lines:
           print("compound_line:", l)
   >>> for l in journal_lines:
           print("jounal_line:", l)

   author_line: T.LI,G.N.PHILLIPS JR.
   compound_line: MOL_ID: 1;
   compound_line:  MOLECULE: MYOGLOBIN (CARBONMONOXY);
   compound_line:  CHAIN: A;
   compound_line:  ENGINEERED: YES;
   compound_line:  MUTATION: YES
   jounal_line:   AUTH   X.ZHAO,K.VYAS,B.D.NGUYEN,K.RAJARATHNAM,G.N.LA MAR,T.LI,
   jounal_line:   AUTH 2 G.N.PHILLIPS,R.F.EICH,J.S.OLSON,J.LING
   jounal_line:   TITL   A DOUBLE MUTANT OF SPERM WHALE MYOGLOBIN MIMICS THE
   jounal_line:   TITL 2 STRUCTURE AND FUNCTION OF ELEPHANT MYOGLOBIN.
   jounal_line:   REF    J.BIOL.CHEM.                  V. 270 20763 1995
   jounal_line:   REFN                   ISSN 0021-9258
   jounal_line:   PMID   7657659
   jounal_line:   DOI    10.1074/JBC.270.35.20781


Molecular Editing
========================================
There are many functions in the API that edit molecules, e.g., adding, deleting and moving the atoms.
More functions are documented in the **Chapi Python API** - `Model Manipulation <https://www.mrc-lmb.cam.ac.uk/lucrezia/libcootapi-documentation/chapi_api.html#model-manipulation>`_ section.

Examples
------------------------------------------------------

**1. Adding water molecules**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: chapi-add-waters.py


**2. Deleting water molecules outliers**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: chapi-delete-waters.py


Refinement
=========================

Examples
------------------------------------------------------


**1. Real Space Refinement using Coot**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: chapi-refine.py


**2. Refinement using Servalcat**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Servalcat supports refinement of structures against X-ray data as well. For more details, see:
Yamashita, K., Palmer, C. M., Burnley, T., Murshudov, G. N. (2021) Acta Cryst. D77, 1282-1291
https://servalcat.readthedocs.io/en/latest/overview.html


.. literalinclude:: servalcat-refine.py


The function :code:`servalcat_refine_xray_with_keywords()` can be used to specify the servalcat keywords.


**3. Base pairs Analysis**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, we analyse the base pairs distances in 2PWT after refinement with Coot.
The complete script can be downloaded :download:`here <basepairs.py>`

.. literalinclude:: basepairs.py
  :lines: 1-48


**Plotting the differences after refinement** 
  
  :download:`basepairs-plotting.py <basepairs-plotting.py>`

.. image:: basepairs_distance_differences.png
   :width: 150%
   :align: center



Ligand
=========================

Examples
------------------------------------------------------

**1. Clark and Labute 2D representation**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The Clark and Labute 2D representation is a schematic diagram in which the ligand is displayed in 2D form, and the interactions to and between the residues in its vicinity are summarized in a concise and information-rich manner
-- Clark, A. H., & Labute, P. (2007). J. Chem. Inf. Model., 47(4), 1937-1948.


.. literalinclude:: flev.py
  :language: python

   
.. image:: LZA_2d_ligand_environment.svg
   :width: 80%
   :align: center