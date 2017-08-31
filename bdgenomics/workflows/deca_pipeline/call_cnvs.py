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

from toil.job import Job

from toil_lib import require
from toil_lib.files import generate_file, move_files
from toil_lib.urls import download_url_job

from bdgenomics.workflows.spark import spawn_spark_cluster
from bdgenomics.workflows.tools.spark_tools import call_deca, \
    MasterAddress, \
    HDFS_MASTER_PORT, \
    SPARK_MASTER_PORT

_log = logging.getLogger(__name__)

def setup_deca_state(job,
                     input_files, targets, output,
                     memory,
                     run_local, num_nodes,
                     aws_access_key_id, aws_secret_access_key):

    if run_local:

        # import bams
        loaded_files = []
        for f in input_files:
            
            file_name = os.path.basename(f)
            file_id = job.wrapJobFn(download_url_job, f)
            job.addChild(file_id)
            
            loaded_files.append((file_name, file_id.rv()))
            
        # import target file
        target_id = job.wrapJobFn(download_url_job, targets)
        job.addChild(target_id)
        target = (os.path.basename(targets), target_id.rv())
        
        call_cnvs = job.wrapJobFn(call_deca_cnvs,
                                  loaded_files, target,
                                  output,
                                  memory,
                                  run_local,
                                  None,
                                  aws_access_key_id, aws_secret_access_key)
        job.addFollowOn(call_cnvs)

    else:
        
        # all files must have s3 urls
        def is_s3(f):
            require(f.startswith("s3a"),
                    "url for file %s did not start with s3a scheme" % f)
            
        is_s3(targets)
        is_s3(output)
        for f in input_files:
            is_s3(f)
            
        # launch the spark cluster
        master_ip = spawn_spark_cluster(job,
                                        int(num_nodes) - 1,
                                        cores=multiprocessing.cpu_count(),
                                        memory=memory)
        
        call_cnvs = job.wrapJobFn(call_deca_cnvs,
                                  input_files,
                                  targets,
                                  output,
                                  memory,
                                  False,
                                  master_ip,
                                  aws_access_key_id, aws_secret_access_key)
        job.addChild(call_cnvs)


def call_deca_cnvs(job,
                   input_bams, targets, output,
                   memory,
                   run_local,
                   master_ip,
                   aws_access_key_id, aws_secret_access_key):
    
    if run_local:
        # get work dir
        work_dir = job.fileStore.getLocalTempDir()
        
        # load targets
        job.fileStore.readGlobalFile(targets[1], os.path.join(work_dir, targets[0]))
        
        # load bams
        inputs = []
        for (bam, bam_id) in input_bams:
            
            inputs.append('/data/%s' % bam)
            job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, bam))
            
        call_deca(job, master_ip=None,
                  arguments=['cnv',
                             '-I', ' '.join(inputs),
                             '-L', '/data/%s' % targets[0],
                             '-o', '/data/cnvs.gff'],
                  memory=memory,
                  run_local=True,
                  work_dir=work_dir,
                  aws_access_key_id=aws_access_key_id,
                  aws_secret_access_key=aws_secret_access_key)
        
        # after running deca, move cnvs to output dir
        move_files([os.path.join(work_dir, 'cnvs.gff')], output)

    else:

        call_deca(job, master_ip=master_ip,
                  arguments=['cnv',
                             '-I', ' '.join(input_bams),
                             '-L', targets,
                             '-o', '%s/cnvs.gff' % output],
                  memory=memory,
                  run_local=False,
                  aws_access_key_id=aws_access_key_id,
                  aws_secret_access_key=aws_secret_access_key)


def load_samples(samples):

    fp = open(samples, 'r')
    s_list = []

    for line in fp:
        s_list.append(line.strip().rstrip())

    fp.close()
    return s_list


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', help='Path to a file containing the S3 URL or local paths to the input SAM or BAM file.'
                        'NOTE: unlike other pipelines, we do not support ftp://, gnos://, etc. schemes.',
                        required=True)
    parser.add_argument('--targets', help='Path to a file containing the S3 URL or local path to the target file.'
                        'NOTE: unlike other pipelines, we do not support ftp://, gnos://, etc. schemes.',
                        required=True)
    parser.add_argument('--output-dir', required=True, default=None,
                        help='full path where final results will be output')
    parser.add_argument('--run-local', default=False, action='store_true',
                        help='if specified, runs locally. exclusive of --num-nodes')
    parser.add_argument('--num-nodes', default=None,
                        help='the number of nodes to use for the spark cluster.'
                        'exclusive of --run-local')
    parser.add_argument('--memory', required=True, default=None,
                        help='Amount of memory (in gb) to allocate for DECA')
    parser.add_argument('--aws_access_key', required=False, default=None,
                        help='AWS access key for authenticating with S3')
    parser.add_argument('--aws_secret_key', required=False, default=None,
                        help='AWS secret key for authenticating with S3')

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    cwd = os.getcwd()

    require(not (args.run_local and args.num_nodes),
            'Only one of --run-local and --num-nodes can be provided.')
    require((not args.aws_access_key and not args.aws_secret_key) or
            (args.aws_access_key and args.aws_secret_key),
            'If AWS access key is provided, AWS secret key must also be provided')

    if not args.run_local:
        require(args.num_nodes,
                'neither --run-local or --num-nodes was specified.')
        require(int(args.num_nodes) > 1,
                'num_nodes allocates one Spark/HDFS master and n-1 workers, and '
                'thus must be greater than 1. %s was passed.' % args.num_nodes)

    samples = load_samples(args.samples)    
        
    Job.Runner.startToil(Job.wrapJobFn(setup_deca_state,
                                       samples,
                                       args.targets,
                                       args.output_dir,
                                       args.memory,
                                       args.run_local,
                                       args.num_nodes,
                                       args.aws_access_key,
                                       args.aws_secret_key), args)


if __name__ == "__main__":
    main()
