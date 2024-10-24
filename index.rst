Libcoot Chapi Documentation
==========================================

Chapi is the alternative name for :code:`coot_headless_api` and is the Pythonic interface to `libcootapi <https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/docs/api/html/>`_.
It is a clear and consistent and easy to use high level interface to the functions of **Coot**.
On creating a new molecule, a *molecule index* will be returned.
Molecules are referred to by this index and using the functions of :code:`molecules_container_t`.
This is unlike many other functions of Python modules, which return a Python representation of the data.


.. toctree::
   :maxdepth: 2

   Introduction<self>
   read_coordinates.rst
   api.rst



