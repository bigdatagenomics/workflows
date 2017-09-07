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

from bdgenomics.workflows.tools.preprocessing import apply_bqsr_recalibration, \
    run_base_recalibration, \
    run_picard_create_sequence_dictionary, \
    run_sambamba_sort, \
    run_samtools_faidx, \
    run_samtools_index
from bdgenomics.workflows.tools.spark_tools import call_adam
from toil_lib.urls import download_url_job

_log = logging.getLogger(__name__)


def run_adam_bqsr(job, sampleId, dbsnp):

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

    # Retrieve file path
    job.fileStore.readGlobalFile(dbsnp, os.path.join(work_dir, 'dbsnp.vcf'))
    
    add_docker_parameters = ['-v',
                             '{}:/data'.format(work_dir)]
    _log.info("Converting dbSNP VCF to ADAM format.")
    call_adam(job, None,
              ["vcf2adam",
               "/data/dbsnp.vcf",
               "/data/dbsnp.adam",
               "-only_variants"],
              memory=str(job.memory),
              run_local=True,
              container='fnothaft/adam',
              add_docker_parameters=add_docker_parameters)

    _log.info("Recalibrating base qualities using ADAM.")
    call_adam(job, None,
              ["transform",
               "/data/reads.adam",
               "/data/reads.sorted.adam",
               "-recalibrate_base_qualities",
               "-known_snps", "/data/dbsnp.adam",
               "-limit_projection"],
              memory=str(job.memory),
              run_local=True,
              container='fnothaft/adam',
              add_docker_parameters=add_docker_parameters)


def run_gatk3_bqsr(job,
                   reads_id, bam_index, ref_id, faidx, ref_dict,
                   dbsnp, mills):

    table = run_base_recalibration(job,
                                   reads_id, bam_index,
                                   ref_id, ref_dict, faidx,
                                   dbsnp, mills)

    apply_bqsr_recalibration(job,
                             table,
                             reads_id, bam_index,
                             ref_id, ref_dict, faidx)


def benchmark_recalibrators(job, sample, ref, dbsnp, mills):

    _log.info("Downloading ref")
    ref_id = download_url_job(job, ref)
    
    _log.info("Indexing reference.")
    faidx = run_samtools_faidx(job, ref_id)

    _log.info("Extracting reference sequence dictionary")
    ref_dict = run_picard_create_sequence_dictionary(job, ref_id)
    
    _log.info("Downloading dbSNP VCF")
    dbsnp_id = download_url_job(job, dbsnp)
    
    _log.info("Downloading Mills VCF")
    mills_id = download_url_job(job, mills)
    
    _log.info("Downloading reads")
    reads_id = download_url_job(job, sample)

    _log.info("Sorting reads by coordinate.")
    coordinate_sorted_bam = run_sambamba_sort(job, reads_id)

    _log.info("Indexing sorted BAM.")
    bam_index = run_samtools_index(job, coordinate_sorted_bam)

    run_adam_bqsr(job, reads_id, dbsnp_id)

    run_gatk3_bqsr(job,
                   reads_id, bam_index,
                   ref_id, faidx, ref_dict,
                   dbsnp_id, mills_id)

    
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--sample', help='The S3 URL or local path to the input SAM or BAM file.')
    parser.add_argument('--ref', help='The S3 URL or local path to the input reference FASTA.')
    parser.add_argument('--dbsnp', help='The S3 URL or local path to the dbSNP VCF.')
    parser.add_argument('--mills', help='The S3 URL or local path to the Mills VCF.')

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    Job.Runner.startToil(Job.wrapJobFn(benchmark_recalibrators,
                                       args.sample,
                                       args.ref,
                                       args.dbsnp,
                                       args.mills), args)
    
    
if __name__ == "__main__":
    main()
