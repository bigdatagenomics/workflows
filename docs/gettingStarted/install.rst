.. highlight:: console

.. _installation-ref:

Installation
============

This document describes how to install the Big Data Genomics workflows.
First, you must `install Toil`_ into a Python `virtualenv`_.
We support Toil versions 3.8.0 to 3.10.1. If you plan to run your workflow on
`Amazon Web Services`_, you must install Toil with the `AWS extra`_
and `configure your AWS account`_.

.. _virtualenv: https://virtualenv.pypa.io/en/stable/
.. _Install Toil: http://toil.readthedocs.io/en/3.10.1/gettingStarted/install.html
.. _AWS extra: http://toil.readthedocs.io/en/3.10.1/gettingStarted/install.html#installing-extra-features
.. _configure your AWS account: http://toil.readthedocs.io/en/3.10.1/gettingStarted/install.html#preparing-your-aws-environment

Basic Installation
------------------

bdgenomics.workflows can be easily installed using pip::

    $ pip install bdgenomics.workflows

Building from source
--------------------

If extending, modifying, or reusing a workflow, you will need to build from
source. This allows changes you make to the workflows to be reflected
immediately in your runtime environment.

First, clone the source::

   $ git clone https://github.com/bigdatagenomics/workflows
   $ cd workflows

Then, create and activate a virtualenv::

   $ virtualenv venv
   $ . venv/bin/activate

`Install Toil`_ into this virtualenv::

   $ pip install toil==3.10.1
   
From there, you can list all available Make targets by running ``make``.
Now, we can install in development mode (such that changes to the 
source code will immediately affect the virtualenv)::

    $ make develop

To build the docs, run ``make develop`` with all extras followed by

::

    $ make docs
