Introduction
============

This code produces a 2D model of a stellar surface that evolves
with time. At each time step in can report the coverage fractions
so that VSPEC can compute a spectrum.

Installation
************

To install using pip:

.. code-block:: shell

    pip install vspec-vsm

or in development mode:

.. code-block:: shell

    git clone https://github.com/VSPEC-collab/vspec-vsm.git
    cd vspec-vsm
    pip install -e .[dev]

The ``[dev]`` flag installs additional dependencies for development.

Additionally, some optional dependencies for plotting can be installed
with:

.. code-block:: shell

    pip install vspec-vsm[plot]

If you have trouble installing cartopy you may need to install from
conda rather than pip:

.. code-block:: shell

    conda install conda-forge cartopy