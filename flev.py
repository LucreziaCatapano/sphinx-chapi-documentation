import chapi
import os

# Set the COOT_PREFIX environment variable to the installation path of COOT
os.environ["COOT_PREFIX"] = "/opt/homebrew/Cellar/coot/1.1.15"

# Initialize the COOT headless API with a molecules container
mc = chapi.molecules_container_t(True)

# Disable the use of GEMMI library
mc.set_use_gemmi(False)

# Read the coordinates from the PDB file and get the molecule index (imol)
imol = mc.read_coordinates("2vtq.pdb")

# Retrieve the monomer information for the ligand "LZA"
mc.get_monomer("LZA")

# Add hydrogen atoms to the molecule
mc.add_hydrogen_atoms(imol)

# Generate a 2D SVG representation of the ligand environment for LZA
svg = mc.get_svg_for_2d_ligand_environment_view(imol, "//A/1299", True)

# Write the generated SVG to a file
f = open("LZA_2d_ligand_environment.svg", "w")
f.write(svg)
f.close()