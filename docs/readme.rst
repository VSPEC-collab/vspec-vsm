VSPEC Template
==============

Template Repository for VSPEC-Associated Packages

Copying
-------

Use this template to create your own repository. Make sure to edit
all necessary files including ``setup.py``, and change the name of
the ``vspec_template`` directory to match the name of your package.

Installation
------------

To install this (or your copy) locally, first clone the repository.
Then navigate to the directory containting ``setup.py`` and 
type

    pip install -e .

This will install the package locally in development mode.

Style
-----

Code written for VSPEC-collab should be documented in the
`numpydoc style <https://numpydoc.readthedocs.io/en/latest/format.html>`_

It is also recommended to use the autopep8 formatting extension
to keep everything readable.

Testing
-------

Tests are located in the ``test`` directory and use the ``pytest`` package.

Documentation
-------------

Documentation is located in the ``docs`` directory and should be
built using ``sphinx``.

Navigate to the ``docs`` directory and issue the command

    sphinx-quickstart

The ``numpydoc`` sphinx extension will be needed, among other extensions.
