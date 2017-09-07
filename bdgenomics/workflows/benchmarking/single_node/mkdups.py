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

from bdgenomics.workflows.tools.preprocessing import picard_mark_duplicates, \
    run_sambamba_markdup, \
    run_sambamba_sort, \
    run_samblaster, \
    run_samtools_index, \
    run_samtools_rmdup, \
    run_samtools_view
from bdgenomics.workflows.tools.spark_tools import call_adam
from toil_lib.urls import download_url_job

_log = logging.getLogger(__name__)


def run_adam_markdups(job, sampleId):

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

    _log.info("Marking duplicate reads using ADAM.")
    call_adam(job, None,
              ["transform",
               "/data/reads.adam",
               "/data/reads.sorted.adam",
               "-mark_duplicate_reads",
               "-limit_projection"],
              memory=str(job.memory),
              run_local=True,
              container='fnothaft/adam',
              add_docker_parameters=add_docker_parameters)


def benchmark_duplicate_markers(job, sample):

    _log.info("Downloading reads")
    reads_id = download_url_job(job, sample)

    _log.info("Sorting reads by coordinate.")
    coordinate_sorted_bam = run_sambamba_sort(job, reads_id)

    _log.info("Indexing sorted BAM.")
    bam_index = run_samtools_index(job, coordinate_sorted_bam)

    _log.info("Marking duplicates with picard.")
    picard_bam = picard_mark_duplicates(job,
                                        coordinate_sorted_bam,
                                        bam_index)

    _log.info("Marking duplicates with samtools.")
    samtools_bam = run_samtools_rmdup(job, coordinate_sorted_bam)

    _log.info("Marking duplicates with sambamba.")
    sambamba_bam = run_sambamba_markdup(job, coordinate_sorted_bam)

    run_adam_markdups(job, reads_id)
    
    _log.info("Sorting reads by name.")
    queryname_sorted_bam = run_sambamba_sort(job, reads_id, sort_by_name=True)

    _log.info("Dumping queryname sorted sam to bam.")
    queryname_sorted_sam = run_samtools_view(job, queryname_sorted_bam)

    _log.info("Marking duplicates with SAMBLASTER.")
    samblaster_sam = run_samblaster(job, queryname_sorted_sam)

    
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--sample', help='The S3 URL or local path to the input SAM or BAM file.')

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    Job.Runner.startToil(Job.wrapJobFn(benchmark_duplicate_markers,
                                       args.sample), args)
    
    
if __name__ == "__main__":
    main()
