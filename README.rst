.. image:: ./img/rhodonite_500_77.png

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis|
        |
..     * - package
..       - | |version| |wheel| |supported-versions| |supported-implementations|
..         | |commits-since|

.. |docs| image:: https://readthedocs.org/projects/rhodonite/badge/?style=flat
    :target: https://readthedocs.org/projects/rhodonite
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/georgerichardson/rhodonite.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/georgerichardson/rhodonite

.. .. |version| image:: https://img.shields.io/pypi/v/rhodonite.svg
..     :alt: PyPI Package latest release
..     :target: https://pypi.python.org/pypi/rhodonite
.. 
.. .. |commits-since| image:: https://img.shields.io/github/commits-since/georgerichardson/rhodonite/v0.1.0.svg
..     :alt: Commits since latest release
..     :target: https://github.com/georgerichardson/rhodonite/compare/v0.1.0...master
.. 
.. .. |wheel| image:: https://img.shields.io/pypi/wheel/rhodonite.svg
..     :alt: PyPI Wheel
..     :target: https://pypi.python.org/pypi/rhodonite
.. 
.. .. |supported-versions| image:: https://img.shields.io/pypi/pyversions/rhodonite.svg
..     :alt: Supported versions
..     :target: https://pypi.python.org/pypi/rhodonite
.. 
.. .. |supported-implementations| image:: https://img.shields.io/pypi/implementation/rhodonite.svg
..     :alt: Supported implementations
..     :target: https://pypi.python.org/pypi/rhodonite


.. end-badges

A Python package for the creation and study of coocurrence networks.


Installation
============

This package is not yet on PyPi. The easiest way to install currently is to
clone this repository and inside the main root directory, run

::

    pip install -e .


Requirements
============

This package is built on ``graph-tool``. Compiling this can take a long time
and 4GB of memory. Pre-compiled versions are available. If you are using Conda
as your package manager then instructions for installation are available

- Linux_
- OSX_

.. _Linux: https://gitlab.com/ostrokach-forge/graph-tool
.. _OSX: https://anaconda.org/ruliana/graph-tool

Some functionality in this package also requires the CFinder tool, which can
be downloaded from http://www.cfinder.org/

.. Documentation
.. =============
.. 
.. https://rhodonite.readthedocs.io/
.. 
.. Development
.. ===========
.. 
.. To run the all tests run::
.. 
..     tox
.. 
.. Note, to combine the coverage data from all the tox environments run:
.. 
.. .. list-table::
..     :widths: 10 90
..     :stub-columns: 1
.. 
..     - - Windows
..       - ::
.. 
..             set PYTEST_ADDOPTS=--cov-append
..             tox
.. 
..     - - Other
..       - ::
.. 
..             PYTEST_ADDOPTS=--cov-append tox
