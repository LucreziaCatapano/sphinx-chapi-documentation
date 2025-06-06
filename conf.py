# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Coot API Documentation'
copyright = '2024, Lucrezia Catapano & Paul Emsley'
author = 'Lucrezia Catapano & Paul Emsley'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = ['sphinx.ext.doctest'
# ]
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.inheritance_diagram']

# autodoc_default_options = {
#     'undoc-members': True,
#     'private-members': True,
#     'inherited-members': True,
#     'show-inheritance': True,
# }
#autosummary_generate = True 

autoclass_content = 'class'
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

add_module_names = False
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'



pygments_style = 'sphinx'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = []



html_static_path = ['_static']

html_css_files = [
    'custom.css',
]

# html_css_files = ['custom.css',  # Custom CSS file
#     'https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css', # Example CDN for Font Awesome
# ]
# html_js_files = [
#     'js/clipboard.min.js',
#     'js/copybutton.js',
# ]