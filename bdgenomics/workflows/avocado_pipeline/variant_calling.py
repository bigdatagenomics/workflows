#!/usr/bin/env python2.7
#
# Licensed to Big Data Genomics (BDG) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The BDG licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import logging
import multiprocessing
import os
import sys
import textwrap
from subprocess import check_call, check_output

from toil_bdg.adam_pipeline.preprocessing import remove_file

import yaml

from toil.job import Job

from toil_lib import require
from toil_lib.files import copy_files, generate_file, move_files
from bdgenomics.workflows.spark import spawn_spark_cluster
from bdgenomics.workflows.tools.spark_tools import call_adam, \
    call_avocado, \
    call_conductor, \
    MasterAddress, \
    HDFS_MASTER_PORT, \
    SPARK_MASTER_PORT

_log = logging.getLogger(__name__)


def download_run_and_upload(job, master_ip, inputs, spark_on_toil):
    """
    Monolithic job that calls data download, conversion, call, upload.
    """
    master_ip = MasterAddress(master_ip)

    bam_name = inputs.sample.split('://')[-1].split('/')[-1]
    sample_name = ".".join(os.path.splitext(bam_name)[:-1])

    hdfs_subdir = sample_name + "-dir"

    if inputs.run_local:
        inputs.local_dir = job.fileStore.getLocalTempDir()
        hdfs_dir = inputs.local_dir
    else:
        inputs.local_dir = None
        hdfs_dir = "hdfs://{0}:{1}/{2}".format(master_ip, HDFS_MASTER_PORT, hdfs_subdir)

    try:
        hdfs_prefix = hdfs_dir + "/" + sample_name
        hdfs_bam = hdfs_dir + "/" + bam_name

        if not inputs.run_local:
            _log.info("Downloading input BAM %s to %s.", bam_name, hdfs_bam)
            call_conductor(job, master_ip, inputs.sample, hdfs_bam,
                           container='fnothaft/conductor',
                           memory=inputs.memory)
        else:
            copy_files([inputs.sample], inputs.local_dir)

        adam_input = hdfs_prefix + ".adam"
        _log.info("Converting input BAM to ADAM.")
        call_adam(job, master_ip,
                  ["transform", hdfs_bam, adam_input],
                  memory=inputs.memory,
                  run_local=inputs.run_local,
                  container='fnothaft/adam')

        avocado_output = hdfs_prefix + ".gt.adam"
        _log.info("Calling variants with avocado.")
        call_avocado(job, master_ip,
                     ["biallelicGenotyper", "-is_not_grc", adam_input, avocado_output],
                     memory=inputs.memory,
                     container='fnothaft/avocado')

        output_vcf = hdfs_prefix + ".vcf"
        _log.info("Converting output ADAM Genotypes to VCF.")
        call_adam(job, master_ip,
                  ["adam2vcf", avocado_output, output_vcf,
                   "-single",
                   "-sort_on_save",
                   "-stringency", "LENIENT"],
                  memory=inputs.memory,
                  run_local=inputs.run_local,
                  container='fnothaft/adam')

        out_file = inputs.output_dir + "/" + sample_name + inputs.suffix + ".bam"

        if not inputs.run_local:
            _log.info("Uploading output VCF %s to %s.", output_vcf, out_file)
            call_conductor(job, master_ip, output_vcf, out_file,
                           memory=inputs.memory,
                           container='fnothaft/conductor')
            remove_file(master_ip, output_vcf, spark_on_toil)
        else:
            local_adam_output = "%s/%s.processed.bam" % (inputs.local_dir, sample_name)
            move_files([local_adam_output], inputs.output_dir)

        remove_file(master_ip, hdfs_subdir, spark_on_toil)
    except:
        remove_file(master_ip, hdfs_subdir, spark_on_toil)
        raise


def static_avocado_dag(job, inputs, sample, output_dir, suffix=''):
    """
    A Toil job function performing Avocado preprocessing on a single sample
    """
    inputs.sample = sample
    inputs.output_dir = output_dir
    inputs.suffix = suffix

    if inputs.master_ip is not None or inputs.run_local:
        # Static, external Spark cluster
        spark_on_toil = False
        spark_work = job.wrapJobFn(download_run_and_upload,
                                   inputs.master_ip, inputs, spark_on_toil)
        job.addChild(spark_work)
    else:
        # Dynamic subclusters, i.e. Spark-on-Toil
        spark_on_toil = True
        cores = multiprocessing.cpu_count()
        master_ip = spawn_spark_cluster(job,
                                        inputs.num_nodes - 1,
                                        cores=cores,
                                        memory=inputs.memory, sparkMasterContainer="fnothaft/apache-spark-master", sparkWorkerContainer = "fnothaft/apache-spark-worker")
        spark_work = job.wrapJobFn(download_run_and_upload,
                                   master_ip, inputs, spark_on_toil)
        job.addChild(spark_work)


def generate_config():
    return textwrap.dedent("""
        # Avocado single sample variant calling pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        num-nodes: 9              # Optional: Number of nodes to use. Do not set if providing master_ip.
        master-ip:                # Optional: IP or hostname of host running for Spark master and HDFS namenode.
                                  # Should be provided instead of num-nodes if pointing at a static (external or
                                  # standalone) Spark cluster. The special value 'auto' indicates the master of
                                  # an externally autoscaled cgcloud spark cluster, i.e. one that is managed by
                                  # the uberscript.
        memory:                   # Required: Amount of memory to allocate for Spark Driver and executor.
                                  # This should be equal to the available memory on each worker node.
        run-local:                # Optional: If true, runs ADAM locally and doesn't connect to a cluster.
    """[1:])


def main():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the avocado variant calling pipeline')
    parser_run.add_argument('--config', default='avocado.config', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--sample', help='The S3 URL or local path to the input SAM or BAM file.'
                            'NOTE: unlike other pipelines, we do not support ftp://, gnos://, etc. schemes.')
    parser_run.add_argument('--output-dir', required=True, default=None,
                            help='full path where final results will be output')
    parser_run.add_argument('-s', '--suffix', default='',
                            help='Additional suffix to add to the names of the output files')

    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()
    cwd = os.getcwd()
    if args.command == 'generate-config':
        generate_file(os.path.join(cwd, 'avocado.config'), generate_config)
    # Pipeline execution
    elif args.command == 'run':
        require(os.path.exists(args.config), '{} not found. Please run '
                                             'generate-config'.format(args.config))
        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        inputs = argparse.Namespace(**parsed_config)

        require(not (inputs.master_ip and (inputs.num_nodes > 0)),
            'Only one of master_ip and num_nodes can be provided.')

        if not hasattr(inputs, 'master_ip'):
            require(inputs.num_nodes > 1,
                'num_nodes allocates one Spark/HDFS master and n-1 workers, and '
                'thus must be greater than 1. %d was passed.' % inputs.num_nodes)

        for arg in [inputs.memory]:
            require(arg, 'Required argument {} missing from config'.format(arg))

        Job.Runner.startToil(Job.wrapJobFn(static_avocado_dag, inputs,
                                           args.sample, args.output_dir), args)

if __name__ == "__main__":
    main()
