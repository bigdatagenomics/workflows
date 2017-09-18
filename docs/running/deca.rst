.. highlight:: console

.. _deca-ref:

Calling CNVs with DECA
======================

Input Files
-----------

The DECA workflow takes two inputs:

1. A feature file that defines the regions over which to call copy number
   variants. This file can be formatted using any of the BED, GTF/GFF2, GFF3,
   Interval List, or NarrowPeak formats. In the AWS workflow, the ADAM Parquet
   Feature format is also supported.
2. A manifest file that contains paths to a set of sorted BAM files. Each file
   must have a scheme listed. In local mode, the file://, http://, and ftp://
   schemes are supported. On AWS, the s3a://, http://, and ftp:// schemes are
   supported. S3a is an overlay over the AWS Simple Storage System (S3) cloud
   data store which is provided by Apache Hadoop.

Running Locally
---------------

To run locally, we invoke the following command::

  $ bdg-deca \
  $   --targets <regions> \
  $   --samples <manifest> \
  $   --output-dir <path-to-save> \
  $   --memory <memory-in-GB> \
  $   --run-local \
  $   file:<toil-jobstore-path>

This command will run in Toil’s single machine mode, and will save the CNV
calls to `<path-to-save>/cnvs.gff`. `<toil-jobstore-path>` is the path to a
temporary directory where Toil will save intermediate files. The
`<memory-in-GB>` parameter should be specified without units; e.g., to allocate
20GB of memory, pass "--memory 20".

Running on AWS
--------------

Once you have launched a Toil-provisioned cluster on AWS, to run the DECA
workflow, invoke the following command::

  $ bdg-deca \
  $   --targets <regions> \
  $   --samples <manifest> \
  $   --output-dir <path-to-save> \
  $   --memory <memory-in-GB> \
  $   --provisioner aws \
  $   --batchSystem mesos \
  $   --mesosMaster $(hostname -i):5050 \
  $   --nodeType <type> \
  $   --num-nodes <spark-workers + 1> \
  $   --minNodes <spark-workers + 2> \
  $   aws:<region>:<toil-jobstore>
  
Toil will launch a cluster with `spark-workers + 2` worker nodes to run this
workflow. For optimal performance, we recommend choosing a number of Apache
Spark worker nodes such that you have no less than 256MB of data per core.
All file paths used in AWS mode must be files stored in AWS’s S3 storage
system, and must have an s3a:// URI scheme.
