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
Toil pipeline for ADAM preprocessing

                        Tree structure of ADAM pipeline
                                       0
                                       |++(1)
                                       2 --> 4 --> 5 --> 6 -> 7
                                        ++(3)

0 = Start Master
1 = Master Service
2 = Start Workers
3 = Worker Service
4 = Do All The Things (Download from s3, convert to ADAM, preprocess, upload to s3)

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
from bdgenomics.workflows.tools.preprocessing import run_gatk_preprocessing, \
    run_picard_create_sequence_dictionary, \
    run_samtools_sort, \
    run_samtools_faidx, \
    run_samtools_index
from toil_lib.urls import download_url_job

_log = logging.getLogger(__name__)


def gatk3_transform(job, ref, in_file, snp_file, g1k_indels, mills_indels):

    _log.info("Downloading ref")
    ref_id = download_url_job(job, ref)
    
    _log.info("Indexing reference.")
    faidx = run_samtools_faidx(job, ref_id)

    _log.info("Extracting reference sequence dictionary")
    ref_dict = run_picard_create_sequence_dictionary(job, ref_id)
    
    _log.info("Downloading reads")
    reads_id = download_url_job(job, in_file)

    _log.info("Sorting reads.")
    sorted_bam = run_samtools_sort(job, reads_id)

    _log.info("Indexing reads.")
    bai = run_samtools_index(job, sorted_bam)
    
    _log.info("Downloading resources")
    g1k_id = download_url_job(job, g1k_indels)
    mills_id = download_url_job(job, mills_indels)
    snp_id = download_url_job(job, snp_file)

    _log.info("Running GATK preprocessing")
    return run_gatk_preprocessing(job,
                                  sorted_bam,
                                  bai,
                                  ref_id,
                                  ref_dict,
                                  faidx,
                                  g1k_id,
                                  mills_id,
                                  snp_id,
                                  realign=True)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--sample', help='The S3 URL or local path to the input SAM or BAM file.')
    parser.add_argument('--ref', help='The S3 URL or local path to the input reference FASTA.')
    parser.add_argument('--dbsnp', help='The S3 URL or local path to the input dbsnp VCF.')
    parser.add_argument('--g1k', help='The S3 URL or local path to the input 1kg indel VCF.')
    parser.add_argument('--mills', help='The S3 URL or local path to the input Mills indel VCF.')

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    cwd = os.getcwd()

    Job.Runner.startToil(Job.wrapJobFn(gatk3_transform,
                                       args.ref,
                                       args.sample,
                                       args.dbsnp,
                                       args.g1k,
                                       args.mills), args)

if __name__ == "__main__":
    main()
