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

from bdgenomics.workflows.adam_pipeline
from bdgenomics.workflows.spark import spawn_spark_cluster
from bdgenomics.workflows.tools.functions import is_s3a
from bdgenomics.workflows.tools.spark_tools import call_mango_browser, \
    MasterAddress, \
    HDFS_MASTER_PORT, \
    SPARK_MASTER_PORT

_log = logging.getLogger(__name__)

def setup_mango_state(job,
                     reference,
                     genes,
                     reads,
                     variants,
                     features,
                     show_genotypes,
                     host,
                     port,
                     memory,
                     run_local, run_mac, num_nodes,
                     aws_access_key_id, aws_secret_access_key):

    if run_local:
        # import reference
        file_name = os.path.basename(reference)
        file_id = job.wrapJobFn(download_url_job, reference)
        job.addChild(file_id)

        loaded_reference = (file_name, file_id.rv())

        loaded_reads = []

        if reads is not None:
            for f in reads.split(','):
                file_name = os.path.basename(f)
                file_id = job.wrapJobFn(download_url_job, f)
                job.addChild(file_id)

                loaded_reads.append((file_name, file_id.rv()))

                # if file is bam, index is required
                if file_name.endswith('bam'):
                    index_name = file_name + ".bai"
                    index_id = job.wrapJobFn(download_url_job, f + ".bai")
                    job.addChild(index_id)

        loaded_variants = []
        if variants is not None:
            for f in variants.split(','):
                file_name = os.path.basename(f)
                file_id = job.wrapJobFn(download_url_job, f)
                job.addChild(file_id)

                loaded_variants.append((file_name, file_id.rv()))

        loaded_features = []
        if features is not None:
            for f in features.split(','):
                file_name = os.path.basename(f)
                file_id = job.wrapJobFn(download_url_job, f)
                job.addChild(file_id)

                loaded_features.append((file_name, file_id.rv()))

        run_mango = job.wrapJobFn(run_mango_browser,
                                  loaded_reference,
                                  genes,
                                  loaded_reads,
                                  loaded_variants,
                                  loaded_features,
                                  show_genotypes,
                                  host,
                                  port,
                                  memory,
                                  run_local,
                                  run_mac,
                                  None,
                                  aws_access_key_id, aws_secret_access_key)
        job.addFollowOn(run_mango)

    else:

        is_s3a(reference)

        if reads is not None:
            for f in reads.split(','):
                is_s3a(f)
                # browser requires bam files to be indexed
                if f.endswith('bam'):
                    is_s3a(f + '.bai')

        if variants is not None:
            for f in variants.split(','):
                is_s3a(f)

        if features is not None:
            for f in features.split(','):
                is_s3a(f)

        # launch the spark cluster
        master_ip = spawn_spark_cluster(job,
                                        int(num_nodes) - 1,
                                        cores=multiprocessing.cpu_count(),
                                        memory=memory)

        run_mango = job.wrapJobFn(run_mango_browser,
                                  reference,
                                  genes, # usually just url
                                  reads,
                                  variants,
                                  features,
                                  show_genotypes,
                                  host,
                                  port,
                                  memory,
                                  False,
                                  False,
                                  master_ip,
                                  aws_access_key_id, aws_secret_access_key)
        job.addChild(run_mango)


def run_mango_browser(job,
                   reference,
                   genes,
                   reads,
                   variants,
                   features,
                   show_genotypes,
                   host,
                   port,
                   memory,
                   run_local,
                   run_mac,
                   master_ip,
                   aws_access_key_id,
                   aws_secret_access_key):

    if run_local:

        # holds mango arguments
        arguments = []

        # get work dir
        work_dir = job.fileStore.getLocalTempDir()

        # load reference
        job.fileStore.readGlobalFile(reference[1], os.path.join(work_dir, reference[0]))
        arguments.append('/data/%s' % reference[0])

        # load genes
        if genes:
            arguments.extend(['-genes', genes[0]])

        # format reads, variants and features

        # load reads
        formatted_reads = []
        for (f, f_id) in reads:
            formatted_reads.append('/data/%s' % f)
            job.fileStore.readGlobalFile(f_id, os.path.join(work_dir, f))

        if formatted_reads:
            arguments.extend(['-reads', ','.join(formatted_reads)])

        # load variants
        formatted_variants = []
        for (f, f_id) in variants:
            formatted_variants.append('/data/%s' % f)
            job.fileStore.readGlobalFile(f_id, os.path.join(work_dir, f))

        if formatted_variants:
            arguments.extend(['-variants', ','.join(formatted_variants)])

        # load features
        formatted_features = []
        for (f, f_id) in features:
            formatted_features.append('/data/%s' % f)
            job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, f))

        if formatted_features:
            arguments.extend(['-features', ','.join(formatted_features)])

        call_mango_browser(job, master_ip=None,
                  arguments=arguments,
                  memory=memory,
                  run_local=True,
                  run_mac=run_mac,
                  work_dir=work_dir,
                  aws_access_key_id=aws_access_key_id,
                  aws_secret_access_key=aws_secret_access_key)

    else:

        # holds mango arguments
        arguments = [reference]

        if genes:
            arguments.extend(['-genes', genes])

        if reads:
            arguments.extend(['-reads', ','.join(reads)])

        if variants:
            arguments.extend(['-variants', ','.join(variants)])

        if features:
            arguments.extend(['-features', ','.join(features)])

        call_mango_browser(job, master_ip=master_ip,
                  arguments=arguments,
                  host=host,
                  port=port,
                  memory=memory,
                  run_local=False,
                  run_mac=False,
                  aws_access_key_id=aws_access_key_id,
                  aws_secret_access_key=aws_secret_access_key)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', help='Path to a file containing the S3 URL or local paths to the reference .2bit, fasta, or adam file.',
                        required=True)
    parser.add_argument('--genes', help='URL to genes.')
    parser.add_argument('--reads', help='Comma separated (,) list of paths to files containing the S3 URL or local paths to input bam or adam files.')
    parser.add_argument('--variants', help='Comma separated (,) list of paths to files containing the S3 URL or local paths to input vcf or adam files.')
    parser.add_argument('--features', help='Comma separated (,) list of paths to files containing the S3 URL or local paths to input bed, narrowpeak or adam files.')
    parser.add_argument('--show_genotypes', help='If set, shows genotypes from variant files.',default=False)
    parser.add_argument('--run-local', default=False, action='store_true',
                        help='if specified, runs locally. exclusive of --num-nodes')
    parser.add_argument('--host', default='localhost', action='store_true',
                            help='host to forward web UI to. Default is localhost.')
    parser.add_argument('--port', default=8080, action='store_true',
                                help='pot to forward web UI to. Default is 8080.')
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
                                       args.reference,
                                       args.genes,
                                       args.reads,
                                       args.variants,
                                       args.features,
                                       args.show_genotypes,
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
