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
import os

from toil.job import Job

from bdgenomics.workflows.tools.preprocessing import run_indel_realignment, \
    run_picard_create_sequence_dictionary, \
    run_realigner_target_creator, \
    run_sambamba_sort, \
    run_samtools_faidx, \
    run_samtools_index
from bdgenomics.workflows.tools.spark_tools import call_adam
from toil_lib.urls import download_url_job

_log = logging.getLogger(__name__)


def run_adam_ri(job, sampleId):

    work_dir = job.fileStore.getLocalTempDir()

    # Retrieve file path
    job.fileStore.readGlobalFile(sampleId, os.path.join(work_dir, 'reads.bam'))
    
    add_docker_parameters = ['-v',
                             '{}:/data'.format(work_dir)]
    _log.info("Converting BAM to ADAM format.")
    call_adam(job, None,
              ["transform",
               "/data/reads.bam",
               "/data/reads.adam"],
              memory=str(job.memory),
              run_local=True,
              container='fnothaft/adam',
              add_docker_parameters=add_docker_parameters)

    _log.info("Realiging INDELs using ADAM.")
    call_adam(job, None,
              ["transform",
               "/data/reads.adam",
               "/data/reads.sorted.adam",
               "-realign_indels",
               "-limit_projection"],
              memory=str(job.memory),
              run_local=True,
              container='fnothaft/adam',
              add_docker_parameters=add_docker_parameters)


def run_gatk3_ir(job,
                 reads_id, bam_index, ref_id, faidx, ref_dict,
                 g1k, mills):

    intervals = run_realigner_target_creator(job,
                                             reads_id, bam_index,
                                             ref_id, ref_dict, faidx,
                                             g1k, mills)

    run_indel_realignment(job,
                          intervals,
                          reads_id, bam_index,
                          ref_id, ref_dict, faidx,
                          g1k, mills)


def benchmark_realigners(job, sample, ref, g1k, mills):

    _log.info("Downloading ref")
    ref_id = download_url_job(job, ref)
    
    _log.info("Indexing reference.")
    faidx = run_samtools_faidx(job, ref_id)

    _log.info("Extracting reference sequence dictionary")
    ref_dict = run_picard_create_sequence_dictionary(job, ref_id)
    
    _log.info("Downloading 1000G VCF")
    g1k_id = download_url_job(job, g1k)
    
    _log.info("Downloading Mills VCF")
    mills_id = download_url_job(job, mills)
    
    _log.info("Downloading reads")
    reads_id = download_url_job(job, sample)

    _log.info("Sorting reads by coordinate.")
    coordinate_sorted_bam = run_sambamba_sort(job, reads_id)

    _log.info("Indexing sorted BAM.")
    bam_index = run_samtools_index(job, coordinate_sorted_bam)

    run_adam_ri(job, reads_id)

    run_gatk3_ir(job,
                 reads_id, bam_index,
                 ref_id, faidx, ref_dict,
                 g1k_id, mills_id)

    
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--sample', help='The S3 URL or local path to the input SAM or BAM file.')
    parser.add_argument('--ref', help='The S3 URL or local path to the input reference FASTA.')
    parser.add_argument('--g1k', help='The S3 URL or local path to the 1000G VCF.')
    parser.add_argument('--mills', help='The S3 URL or local path to the Mills VCF.')

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    Job.Runner.startToil(Job.wrapJobFn(benchmark_realigners,
                                       args.sample,
                                       args.ref,
                                       args.g1k,
                                       args.mills), args)
    
    
if __name__ == "__main__":
    main()
