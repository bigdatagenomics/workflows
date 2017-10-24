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
from bdgenomics.workflows.tools.functions import is_s3
from bdgenomics.workflows.tools.spark_tools import call_mango_notebook, \
    MasterAddress, \
    HDFS_MASTER_PORT, \
    SPARK_MASTER_PORT

_log = logging.getLogger(__name__)

def setup_mango_state(job,
                     host,
                     port,
                     memory,
                     run_local, run_mac, num_nodes,
                     aws_access_key_id, aws_secret_access_key):

    if run_local:


        run_mango = job.wrapJobFn(run_mango_notebook,
                                  host,
                                  port,
                                  memory,
                                  run_local,
                                  run_mac,
                                  None,
                                  aws_access_key_id, aws_secret_access_key)
        job.addFollowOn(run_mango)

    else:

        # launch the spark cluster
        master_ip = spawn_spark_cluster(job,
                                        int(num_nodes) - 1,
                                        cores=multiprocessing.cpu_count(),
                                        memory=memory)

        run_mango = job.wrapJobFn(run_mango_notebook,
                                  host,
                                  port,
                                  memory,
                                  False,
                                  False,
                                  master_ip,
                                  aws_access_key_id, aws_secret_access_key)
        job.addChild(run_mango)


def run_mango_notebook(job,
                   host,
                   port,
                   memory,
                   run_local,
                   run_mac,
                   master_ip,
                   aws_access_key_id,
                   aws_secret_access_key):

    # get work dir
    work_dir = job.fileStore.getLocalTempDir()

    arguments = []
    arguments.append('--allow-root') # required for npm in docker

    if run_local:

        # TODO: NOT SURE IF WE NEED THIS WHEN NET-HOST IS SET
        arguments.append('--ip=0.0.0.0')

        call_mango_notebook(job, master_ip=None, arguments=arguments,
                  memory=memory,
                  run_local=True,
                  run_mac=run_mac,
                  work_dir=work_dir,
                  aws_access_key_id=aws_access_key_id,
                  aws_secret_access_key=aws_secret_access_key)

    else:

        call_mango_notebook(job, master_ip=master_ip, arguments=arguments,
                  host=host,
                  port=port,
                  memory=memory,
                  run_local=False,
                  run_mac=False,
                  aws_access_key_id=aws_access_key_id,
                  aws_secret_access_key=aws_secret_access_key)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--run-local', default=False, action='store_true',
                        help='if specified, runs locally. exclusive of --num-nodes')
    parser.add_argument('--host', default='localhost', action='store_true',
                            help='host to forward web UI to. Default is localhost.')
    parser.add_argument('--port', default=10000, action='store_true',
                                help='pot to forward web UI to. Default is 10000.')
    parser.add_argument('--run-mac', default=False, action='store_true',
                            help='if specified, runs on mac.')
    parser.add_argument('--num-nodes', default=None,
                        help='the number of nodes to use for the spark cluster.'
                        'exclusive of --run-local')
    parser.add_argument('--memory', required=True, default=None,
                        help='Amount of memory (in gb) to allocate for mango')
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


    _log.info("startToil")

    Job.Runner.startToil(Job.wrapJobFn(setup_mango_state,
                                       args.host,
                                       args.port,
                                       args.memory,
                                       args.run_local,
                                       args.run_mac,
                                       args.num_nodes,
                                       args.aws_access_key,
                                       args.aws_secret_key), args)


if __name__ == "__main__":
    main()
