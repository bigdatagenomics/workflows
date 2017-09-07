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

from bdgenomics.workflows.tools.preprocessing import run_picard_sort, \
    run_samtools_sort, \
    run_sambamba_sort
from bdgenomics.workflows.tools.spark_tools import call_adam
from toil_lib.urls import download_url_job

_log = logging.getLogger(__name__)


def run_adam_sort(job, sampleId):

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

    _log.info("Sorting reads using ADAM.")
    call_adam(job, None,
              ["transform",
               "/data/reads.adam",
               "/data/reads.sorted.adam",
               "-sort_reads",
               "-limit_projection"],
              memory=str(job.memory),
              run_local=True,
              container='fnothaft/adam',
              add_docker_parameters=add_docker_parameters)


def benchmark_sorters(job, sample):

    _log.info("Downloading reads")
    reads_id = download_url_job(job, sample)

    _log.info("Sorting reads with picard.")
    picard_sorted_bam = run_picard_sort(job, reads_id)

    _log.info("Sorting reads with samtools.")
    samtools_sorted_bam = run_samtools_sort(job, reads_id)

    _log.info("Sorting reads with sambamba.")
    sambamba_sorted_bam = run_sambamba_sort(job, reads_id)

    run_adam_sort(job, reads_id)
    
    
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--sample', help='The S3 URL or local path to the input SAM or BAM file.')

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    Job.Runner.startToil(Job.wrapJobFn(benchmark_sorters,
                                       args.sample), args)
    
    
if __name__ == "__main__":
    main()
