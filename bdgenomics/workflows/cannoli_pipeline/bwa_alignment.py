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
    call_cannoli, \
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

    fastq_name = inputs.sample.split('://')[-1].split('/')[-1]
    sample_name = ".".join(os.path.splitext(fastq_name)[:-1])

    hdfs_subdir = sample_name + "-dir"

    if inputs.run_local:
        inputs.local_dir = job.fileStore.getLocalTempDir()
        hdfs_dir = inputs.local_dir
    else:
        inputs.local_dir = None
        hdfs_dir = "hdfs://{0}:{1}/{2}".format(master_ip, HDFS_MASTER_PORT, hdfs_subdir)

    try:
        hdfs_prefix = hdfs_dir + "/" + sample_name
        hdfs_reads = hdfs_dir + "/" + fastq_name

        if not inputs.run_local:
            _log.info("Downloading input reads %s to %s.", fastq_name, hdfs_reads)
            call_conductor(job, master_ip, inputs.sample, hdfs_reads,
                           container='fnothaft/conductor',
                           memory=inputs.memory)
            
            index_exts = ['', '.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']
            hdfs_index = os.path.join(hdfs_prefix, 'reference.fa')
            for ext in index_exts:
                index_path = inputs.index + ext
                hdfs_index_ext = hdfs_index + ext
                _log.info("Downloading index file %s to %s.", index_path, hdfs_index_ext)
                call_conductor(job, master_ip, index_path, hdfs_index_ext,
                               container='fnothaft/conductor',
                               memory=inputs.memory)

            sd_path = inputs.index.replace('.fa', '.dict')
            hdfs_sd = hdfs_index.replace('.fa', '.dict')
            _log.info("Downloading sequence dictionary %s to %s.", sd_path, hdfs_sd)
            call_conductor(job, master_ip, sd_path, hdfs_sd,
                           container='fnothaft/conductor',
                           memory=inputs.memory)
                
        else:
            copy_files([inputs.sample], inputs.local_dir)

        aligned_output = hdfs_prefix + ".bam"
        _log.info("Aligning reads with Cannoli and BWA.")
        call_cannoli(job, master_ip,
                     ["bwa", "-single", hdfs_reads, aligned_output, inputs.sample_id,
                      '-use_docker',
                      '-docker_image', 'fnothaft/bwa:debug-3',
                      '-index', hdfs_index,
                      '-add_indices',
                      '-sequence_dictionary', hdfs_sd],
                     memory=inputs.memory,
                     container='fnothaft/cannoli:1508-1509')
        out_file = inputs.output_dir + "/" + sample_name + inputs.suffix + ".bam"

        if not inputs.run_local:
            _log.info("Uploading output BAM %s to %s.", aligned_output, out_file)
            call_conductor(job, master_ip, aligned_output, out_file,
                           memory=inputs.memory,
                           container='fnothaft/conductor')
            remove_file(master_ip, output_vcf, spark_on_toil)
        else:
            local_adam_output = "%s/%s.bam" % (inputs.local_dir, sample_name)
            move_files([local_adam_output], inputs.output_dir)

        remove_file(master_ip, hdfs_subdir, spark_on_toil)
    except:
        remove_file(master_ip, hdfs_subdir, spark_on_toil)
        raise


def static_cannoli_dag(job, inputs, sample, sample_id, output_dir, suffix=''):
    """
    A Toil job function performing alignment using cannoli on a single sample
    """
    inputs.sample = sample
    inputs.sample_id = sample_id
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
                                        memory=inputs.memory, sparkMasterContainer="fnothaft/apache-spark-master:2.1.0--74e45e9a58550e14db0e1ad48624c839ebd5e8f8", sparkWorkerContainer = "fnothaft/apache-spark-worker:2.1.0--74e45e9a58550e14db0e1ad48624c839ebd5e8f8")
        spark_work = job.wrapJobFn(download_run_and_upload,
                                   master_ip, inputs, spark_on_toil)
        job.addChild(spark_work)


def generate_config():
    return textwrap.dedent("""
        # BWA-on-Cannoli single sample alignment pipeline configuration file
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
        index:                    # Required: Path to the index file.
    """[1:])


def main():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the cannoli alignment pipeline, using BWA')
    parser_run.add_argument('--config', default='cannoli.config', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--sample-id', help='The sample ID.')
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
        generate_file(os.path.join(cwd, 'cannoli.config'), generate_config)
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

        Job.Runner.startToil(Job.wrapJobFn(static_cannoli_dag, inputs,
                                           args.sample, args.sample_id, args.output_dir), args)

if __name__ == "__main__":
    main()
