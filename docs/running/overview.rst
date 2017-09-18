.. highlight:: console

.. _overview-ref:

Overview
========

The workflows package supports running `Big Data Genomics`_ tools using `Toil`_
in both local mode and in distributed mode on AWS using an `autoscaled`_ or
statically provisioned cluster.

.. _Big Data Genomics: https://github.com/bigdatagenomics
.. _Toil: https://github.com/BD2KGenomics/toil
.. _autoscaled: http://toil.readthedocs.io/en/3.10.1/running/amazon.html#running-a-workflow-with-autoscaling

Currently, we support the following workflow:

.. toctree::
   deca

Additional workflows are in development and alpha versions of these workflows
are packaged and made available in this repository. However, since they are
alpha, we have not documented these workflows, and do not currently support
them.

Preparing to Run on AWS
-----------------------

If you plan to run the workflows on your local machine, then you only need to
ensure that `Toil`_ is installed. If you are running on the `Amazon Web
Services`_ cloud, then you need to ensure that Toil is installed with the
AWS extra. You should use `Toil to launch a cluster`_ in AWS. Once you have
launched the cluster, you should `ssh in to the cluster`_ and install the
bdgenomics.workflow package locally in a virtualenv on the cluster leader.
Once you have followed this process, you can run your desired workflow.

.. _Amazon Web Services: https://aws.amazon.com/
.. _Toil to launch a cluster: http://toil.readthedocs.io/en/3.10.1/running/amazon.html
.. _ssh in to the cluster: http://toil.readthedocs.io/en/3.10.1/running/amazon.html#ssh-cluster
