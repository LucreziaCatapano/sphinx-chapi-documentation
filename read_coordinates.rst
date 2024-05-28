Reading Files
##################

.. _read_coordinates:

Coordinate Files
================================================================


Chapi supports the following coordinate file formats:

    * mmCIF (PDBx/mmCIF),
    * PDB
    * SHELXL
    * MOL file



The :code:`read_coordinates()` is a function provided by the :code:`molecules_container_t` class of the :code:`chapi` module, designed to read files in different format, including mmCIF or PDB format.
Additionally, the function :code:`read_pdb()` can be used to read PDB files.

.. code-block:: python

   import chapi

   mc = chapi.molecules_container_t(True)

   # read PDB file
   imol_pdb = mc.read_coordinates("/data/tutorial-modern.pdb")  
   
   # read mmCIF file
   imol_mmcif = mc.read_coordinates("/data/tutorial-modern.cif")   


MTZ and Map Files
================================================================


MTZ file format can be read by using the :code:`read_mtz()` function.

.. code-block:: python

   imol_mtz = mc.read_mtz("/data/rnase-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", False, False)

The latest two arguments are: 

* :code:`use_weight` (bool):

* :code:`is_a_difference_map` (bool):

MRC/CCP4 Maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EM maps are typically in MRC/CCP4 map format

.. code-block:: python
   
   imol_map = mc.read_ccp4_map("/data/test-molecule-container-test-data/emd_16890.map", False)

The latest argument is :code:`is_a_difference_map` (bool):

Writing Files
##################

The :code:`write_coordinates()` and :code:`write_map()` functions are used to write coordinates files (PDB and mmCIF) 
and map files (MTZ and map) respectively.


.. code-block:: python

   # read mmCIF file
   imol_mmcif = mc.read_coordinates("/data/tutorial-modern.cif") 

   # write mmCIF file
   imol_mmcif_new = mc.write_coordinates(imol_mmcif, "tutorial-modern-new.cif")

   # read MTZ file
   imol_mtz = mc.read_mtz("/data/rnase-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", False, False)

   # write MTZ file
   imol_mtz_new = mc.write_map(imol_mtz, "rnase-1.8-all_refmac1_new.mtz")

Molecular Models
##################

Molecular Information
================================================================

The following functions return information about macromolecular models.

* :code:`get_molecule_diameter()`: 

* :code:`get_number_of_atoms()`:

* :code:`get_number_of_hydrogen_atoms()`:

* :code:`get_cell()`: 

* :code:`get_symmetry()`: 

* :code:`get_hb_type()`: 

**Chains**

* :code:`get_chains_in_model()`:

* :code:`get_single_letter_codes_for_chain()`:

* :code:`get_ncs_related_chains()`:

**Residue**

* :code:`get_residue_name()`:

* :code:`get_residue_names_with_no_dictionary()`:

* :code:`get_residues_with_missing_atoms()`:

* :code:`get_residues_near_residue()`:


.. doctest::

   >>> imol = mc.read_pdb("tutorial-modern.pdb")
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


**Example #1: header info**

.. doctest::

   >>> import chapi
   >>> mc = chapi.molecules_container_t(False)
   >>> mc.set_use_gemmi(False)

   >>> imol = mc.read_pdb('1mcy.pdb')
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
   INFO:: There are 8 helices and 0 sheets

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
================================================================
There are many functions in the api that edit molecules, e.g., adding, deleting and moving the atoms.


**Example #1: adding water molecules**

.. literalinclude:: chapi-add-waters.py

**Example #2 : deleting water molecules outliers**

.. literalinclude:: chapi-delete-waters.py
