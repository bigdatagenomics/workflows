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

import yaml

from toil.job import Job

from toil_lib import require
from toil_lib.files import generate_file, move_files
from toil_lib.spark import spawn_spark_cluster
from toil_lib.urls import download_url_job
from toil_lib.tools.spark_tools import call_deca, \
    MasterAddress, \
    HDFS_MASTER_PORT, \
    SPARK_MASTER_PORT

_log = logging.getLogger(__name__)

def setup_deca_state(job, input_files, targets, output, memory, run_local, num_nodes):

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
                                  None)
        job.addFollowOn(call_cnvs)

    else:
        
        # all files must have s3 urls
        def is_s3(f):
            require(f.startswith("s3"),
                    "url for file %s did not start with s3 scheme" % f)                
            
        is_s3(targets)
        is_s3(output)
        for f in input_files:
            is_s3(f)
            
        # launch the spark cluster
        master_ip = spawn_spark_cluster(job,
                                        num_nodes - 1,
                                        cores=cores,
                                        memory=inputs.memory,
                                        sparkMasterContainer="quay.io/ucsc_cgl/apache-spark-master",
                                        sparkWorkerContainer="quay.io/ucsc_cgl/apache-spark-worker")
        
        call_cnvs = job.wrapJobFn(call_deca_cnvs,
                                  input_bams,
                                  targets,
                                  output,
                                  memory,
                                  False,
                                  master_ip)
        job.addChild(call_cnvs)


def call_deca_cnvs(job, input_bams, targets, output, memory, run_local, master_ip):
    
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
                  work_dir=work_dir)
        
        # after running deca, move cnvs to output dir
        move_files([os.path.join(work_dir, 'cnvs.gff')], output)

    else:

        call_deca(job, master_ip=master_ip,
                  arguments=['cnv',
                             '-I', ' '.join(input_bams),
                             '-L', '/data/%s' % targets,
                             '-o', '%s/cnvs.gff' % output],
                  memory=memory,
                  run_local=False)


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

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    cwd = os.getcwd()

    require(not (args.run_local and args.num_nodes),
            'Only one of --run-local and --num-nodes can be provided.')

    if not args.run_local:
        require(args.num_nodes,
                'neither --run-local or --num-nodes was specified.')
        require(args.num_nodes > 1,
                'num_nodes allocates one Spark/HDFS master and n-1 workers, and '
                'thus must be greater than 1. %d was passed.' % inputs.num_nodes)

    samples = load_samples(args.samples)    
        
    Job.Runner.startToil(Job.wrapJobFn(setup_deca_state,
                                       samples,
                                       args.targets,
                                       args.output_dir,
                                       args.memory,
                                       args.run_local,
                                       args.num_nodes), args)


if __name__ == "__main__":
    main()
