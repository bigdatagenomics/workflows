bdgenomics.workflows
====================

A package of single-node and distributed workflows for running `Big Data
Genomics`_ tools using `Toil`_.

.. _Big Data Genomics: https://github.com/bigdatagenomics
.. _Toil: https://github.com/BD2KGenomics/toil

For detailed documentation, see our `readthedocs site`_.

.. _readthedocs site: https://bdg-workflows.readthedocs.io

Supported workflows
-------------------

We support Workflows for running:

* Copy number variant calling using `DECA`_.

.. _DECA: https://github.com/bigdatagenomics/deca

Quick Start
-----------

The latest release of `bdgenomics.workflows` can be installed using Pip
into a virtualenv. You must install `Toil`_ first::

  $ virtualenv bdg-workflows
  $
  $ . bdg-workflows/bin/activate
  $
  $ pip install toil==3.10.1
  $
  $ pip install bdgenomics.workflows
  $
  $ bdg-deca

To install from source::

  $ virtualenv bdg-workflows
  $
  $ . bdg-workflows/bin/activate
  $
  $ pip install toil==3.10.1
  $
  $ make develop
  $
  $ bdg-deca
