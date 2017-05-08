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

"""
Toil pipeline for ADAM transform alignments

                        Tree structure of ADAM transform alignments
                                       0
                                       |++(1)
                                       2 --> 4
                                        ++(3)

0 = Start Master
1 = Master Service
2 = Start Workers
3 = Worker Service
4 = Do All The Things (Download from s3, convert to ADAM, upload to s3)

================================================================================
:Dependencies
docker          - apt-get install docker (or 'docker.io' for linux)
toil            - pip install toil
"""

import argparse
import logging
import multiprocessing
import os
import sys
import textwrap
from subprocess import check_call, check_output

import yaml

from toil.job import Job

from toil_lib import require
from toil_lib.files import copy_files, generate_file, move_files
from toil_lib.spark import spawn_spark_cluster
from toil_lib.tools.spark_tools import call_adam, \
    call_conductor, \
    MasterAddress, \
    HDFS_MASTER_PORT, \
    SPARK_MASTER_PORT

_log = logging.getLogger(__name__)


def remove_file(master_ip, filename, spark_on_toil):
    """
    Remove the given file from hdfs with master at the given IP address.

    :type masterIP: MasterAddress
    """
    master_ip = master_ip.actual

    ssh_call = ['ssh', '-o', 'StrictHostKeyChecking=no', master_ip]

    try:
        if spark_on_toil:
            output = check_output(ssh_call + ['docker', 'ps'])
            container_id = next(line.split()[0] for line in output.splitlines() if 'apache-hadoop-master' in line)
            ssh_call += ['docker', 'exec', container_id]

        check_call(ssh_call + ['hdfs', 'dfs', '-rm', '-r', '/' + filename])
    except:
        pass


def truncate_file(master_ip, filename, spark_on_toil):
    """
    Truncate the given hdfs file to 10 bytes with master at the given IP address.

    :type masterIP: MasterAddress
    """
    master_ip = master_ip.actual

    ssh_call = ['ssh', '-o', 'StrictHostKeyChecking=no', master_ip]
    hdfs = ['hdfs']

    if spark_on_toil:
        output = check_output(ssh_call + ['docker', 'ps'])
        container_id = next(line.split()[0] for line in output.splitlines() if 'apache-hadoop-master' in line)
        ssh_call += ['docker', 'exec', container_id]
        hdfs = ['/opt/apache-hadoop/bin/hdfs']

    try:
        check_call(ssh_call + hdfs + ['dfs', '-truncate', '-w', '10', '/' + filename])
    except:
        pass


def download_data(job, master_ip, inputs, bam, hdfs_bam):
    """
    Downloads input data files from S3.

    :type masterIP: MasterAddress
    """

    _log.info("Downloading input BAM %s to %s.", bam, hdfs_bam)
    call_conductor(job, master_ip, bam, hdfs_bam,
                   container='fnothaft/conductor',
                   memory=inputs.memory)


def adam_convert(job, master_ip, inputs, in_file, adam_file, spark_on_toil):
    """
    Convert input sam/bam file into ADAM format.
    """

    _log.info("Converting input BAM to ADAM.")
    call_adam(job, master_ip,
              ["transform", in_file, adam_file],
              memory=inputs.memory,
              run_local=inputs.run_local,
              native_adam_path=inputs.native_adam_path,
              container='fnothaft/adam')

    in_file_name = in_file.split("/")[-1]
    remove_file(master_ip, in_file_name, spark_on_toil)


def upload_data(job, master_ip, inputs, hdfs_name, upload_name, spark_on_toil):
    """
    Upload file hdfs_name from hdfs to upload_name on s3.
    """

    _log.info("Uploading output BAM %s to %s.", hdfs_name, upload_name)
    call_conductor(job, master_ip, hdfs_name, upload_name,
                   memory=inputs.memory,
                   container='fnothaft/conductor')
    remove_file(master_ip, hdfs_name, spark_on_toil)


def download_run_and_upload(job, master_ip, inputs, spark_on_toil):
    """
    Monolithic job that calls data download, conversion, and upload.
    """
    master_ip = MasterAddress(master_ip)

    bam_name = inputs.sample.split('://')[-1].split('/')[-1]
    sample_name = ".".join(os.path.splitext(bam_name)[:-1])

    hdfs_subdir = sample_name + "-dir"

    if inputs.run_local:
        inputs.local_dir = job.fileStore.getLocalTempDir()
        if inputs.native_adam_path is None:
            hdfs_dir = "/data/"
        else:
            hdfs_dir = inputs.local_dir
    else:
        inputs.local_dir = None
        hdfs_dir = "hdfs://{0}:{1}/{2}".format(master_ip, HDFS_MASTER_PORT, hdfs_subdir)

    try:
        hdfs_prefix = hdfs_dir + "/" + sample_name
        hdfs_bam = hdfs_dir + "/" + bam_name

        if not inputs.run_local:
            download_data(job, master_ip, inputs, inputs.sample, hdfs_bam)
        else:
            copy_files([inputs.sample, inputs.dbsnp], inputs.local_dir)

        hdfs_adam = hdfs_prefix + ".alignments.adam"
        adam_convert(job, master_ip, inputs, hdfs_bam, hdfs_adam, spark_on_toil)

        output_adam = inputs.output_dir + "/" + sample_name + inputs.suffix + ".alignments.adam"

        if not inputs.run_local:
            upload_data(job, master_ip, inputs, hdfs_adam, output_adam, spark_on_toil)
        else:
            local_output_adam = "%s/%s.alignments.adam" % (inputs.local_dir, sample_name)
            move_files([local_output_adam], inputs.output_dir)

        remove_file(master_ip, hdfs_subdir, spark_on_toil)
    except:
        remove_file(master_ip, hdfs_subdir, spark_on_toil)
        raise


def static_transform_alignments_dag(job, inputs, sample, output_dir, suffix=''):
    """
    A Toil job function performing ADAM transform alignments on a single sample.
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
        # ADAM transform alignments configuration file
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
        local-dir:                # Required if run-local is true. Sets the local directory to use for input.
        native-adam-path:         # Optional: If set, runs ADAM using the local build of ADAM at this path.
    """[1:])


def main():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the ADAM transform alignments pipeline')
    parser_run.add_argument('--config', default='transform_alignments.config', type=str,
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
        generate_file(os.path.join(cwd, 'transform_alignments.config'), generate_config)
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

        Job.Runner.startToil(Job.wrapJobFn(static_transform_alignments_dag, inputs,
                                           args.sample, args.output_dir), args)

if __name__ == "__main__":
    main()
