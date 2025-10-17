Installation
==========================================

CCP4
------------------------------------------

**chapi** is available as a module in ccp4-python (Python 3.9.18) starting from `CCP4 version 9.0 <https://www.ccp4.ac.uk/download/index.php#os=macos>`_.

To install the latest chapi library, you have several options:

Install via build-it-3-3 script
-------------------------------

1. Download the **build-it-3-3** script from:

   https://github.com/pemsley/coot/blob/main/build-it-3-3

2. Set the environment variable to install only chapi::

     export CHAPI_ONLY=true

3. Run the build script::

     bash build-it-3-3

The script will handle dependencies and attempt to build a Python environment with the latest chapi. Build logs will be saved in `~/public_html/build_logs`.

Install via Homebrew
--------------------

Alternatively, you can install coot (which includes chapi) using **Homebrew**::

   brew install brewsci/bio/coot

Install via Conda
-----------------

To install chapi using conda, follow these steps:

1. Add the conda-forge channel and set strict channel priority::

     conda config --add channels conda-forge

     conda config --set channel_priority strict
   
2. Note: Python bindings for Coot are available for versions 3.9, 3.10, 3.11, 3.12, and 3.13. If no Python version is specified, Python 3.13 will be installed by default.
   Create a new environment and install coot-headless::

     conda create -n coot-env coot-headless -c conda-forge -c bioconda
   

3. Activate the environment::

     conda activate coot-env
   
