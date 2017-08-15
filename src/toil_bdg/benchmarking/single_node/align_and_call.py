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

from toil_bdg.benchmarking.single_node.bqsr import run_adam_bqsr
from toil_bdg.benchmarking.single_node.mkdups import run_adam_markdups
from toil_bdg.benchmarking.single_node.realign_indels import run_adam_ri
from toil_bdg.benchmarking.single_node.sort import run_adam_sort

from toil.job import Job

from toil_lib import require
from toil_lib.files import move_files
from toil_lib.urls import download_url_job
from toil_lib.tools.aligners import run_bowtie2, \
    run_bwakit, \
    run_snap
from toil_lib.tools.indexing import run_bowtie2_index, \
    run_bwa_index, \
    run_samtools_faidx, \
    run_snap_index, \
    run_soap3_index
from toil_lib.tools.preprocessing import run_picard_sort, \
    run_picard_create_sequence_dictionary, \
    picard_mark_duplicates, \
    run_base_recalibration, \
    apply_bqsr_recalibration, \
    run_sambamba_index, \
    run_sambamba_markdup, \
    run_sambamba_sort, \
    run_sambamba_view, \
    run_samblaster, \
    run_samtools_view
from toil_lib.tools.spark_tools import call_adam, \
    call_avocado, \
    call_cannoli
from toil_lib.tools.variant_callers import run_freebayes, \
    run_platypus, \
    run_strelka, \
    run_16gt, \
    run_manta, \
    run_samtools_mpileup, \
    run_bcftools_call, \
    run_gatk3_haplotype_caller
from toil_lib.urls import download_url_job

_log = logging.getLogger(__name__)

def align_reads(job,
                sample_name,
                fastq1, fastq2,
                ref,
                ref_fai,
                bowtie2_indices,
                bwa_indices,
                snap_indices):
    '''
    '''

    bwa_config = argparse.Namespace(**bwa_indices)
    bwa_config.ref = ref
    bwa_config.fai = ref_fai
    bwa_config.r1 = fastq1
    bwa_config.r2 = fastq2
    bwa_config.uuid = sample_name
    bwa_config.library = sample_name
    bwa_config.platform = 'ILLUMINA'
    bwa_config.program_unit = 'unknown'
    bwa_bam_id = job.wrapJobFn(run_bwakit,
                               bwa_config,
                               sort=False,
                               trim=False,
                               benchmarking=True)
    job.addChild(bwa_bam_id)

    bowtie_tags = ['LB:%s' % sample_name,
                   'PL:ILLUMINA',
                   'PU:unknown',
                   'SM:%s' % sample_name]
    bowtie2_sam_id = job.wrapJobFn(run_bowtie2,
                                   fastq1,
                                   ref,
                                   bowtie2_indices,
                                   read2=fastq2,
                                   read_group_id=sample_name,
                                   read_group_tags=bowtie_tags,
                                   benchmarking=True)
    job.addChild(bowtie2_sam_id)

    snap_rg_line = 'ID:%s\t%s' % (sample_name, '\t'.join(bowtie_tags))
    snap_bam_id = job.wrapJobFn(run_snap,
                                fastq1,
                                snap_indices,
                                read2=fastq2,
                                read_group_line=snap_rg_line,
                                benchmarking=True)
    job.addChild(snap_bam_id)

    return [('bwa', bwa_bam_id.rv()),
            ('bowtie2', bowtie2_sam_id.rv()),
            ('snap_bam_id', snap_bam_id.rv())]

        
def run_adam_preprocessing(job,
                           ref,
                           dbsnp,
                           bam,
                           aligner,
                           aligner_runtime,
                           is_bam=True,
                           has_md_tags=True):

    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(dbsnp, os.path.join(work_dir, 'dbsnp.vcf'))

    # Retrieve file path
    if is_bam:
        job.fileStore.readGlobalFile(bam, os.path.join(work_dir, 'reads.bam'))
        reads = '/data/reads.bam'
    else:
        job.fileStore.readGlobalFile(bam, os.path.join(work_dir, 'reads.sam'))
        reads = '/data/reads.sam'

    add_docker_parameters = ['-v',
                             '{}:/data'.format(work_dir)]

    if has_md_tags:
        transform_runtime = call_adam(job, None,
                                      ['transformAlignments',
                                       reads,
                                       '/data/reads.adam'],
                                      memory=str((job.memory / (1024 * 1024 * 1024))),
                                      run_local=True,
                                      container='quay.io/ucsc_cgl/adam',
                                      add_docker_parameters=add_docker_parameters,
                                      benchmarking=True)

    else:
        job.fileStore.readGlobalFile(ref, os.path.join(work_dir, 'ref.fa'))
        transform_runtime = call_adam(job, None,
                                      ['transformAlignments',
                                       reads,
                                       '/data/reads.adam',
                                       '-add_md_tags', '/data/ref.fa'],
                                      memory=str((job.memory / (1024 * 1024 * 1024))),
                                      run_local=True,
                                      container='quay.io/ucsc_cgl/adam',
                                      add_docker_parameters=add_docker_parameters,
                                      benchmarking=True)


    markdups_runtime = call_adam(job, None,
                                 ['transformAlignments',
                                  '/data/reads.adam',
                                  '/data/reads.mkdups.adam',
                                  '-limit_projection',
                                  '-mark_duplicate_reads'],
                                 memory=str((job.memory / (1024 * 1024 * 1024))),
                                 run_local=True,
                                 container='quay.io/ucsc_cgl/adam',
                                 add_docker_parameters=add_docker_parameters,
                                 benchmarking=True)
    
    ri_runtime = call_adam(job, None,
                           ['transformAlignments',
                            '/data/reads.mkdups.adam',
                            '/data/reads.ri.adam',
                            '-realign_indels'],
                           memory=str((job.memory / (1024 * 1024 * 1024))),
                           run_local=True,
                           container='quay.io/ucsc_cgl/adam',
                           add_docker_parameters=add_docker_parameters,
                           benchmarking=True)

    bqsr_runtime = call_adam(job, None,
                             ['transformAlignments',
                              '/data/reads.ri.adam',
                              '/data/reads.bqsr.adam',
                              '-recalibrate_base_qualities',
                              '-known_snps', '/data/dbsnp.vcf'],
                             memory=str((job.memory / (1024 * 1024 * 1024))),
                             run_local=True,
                             container='quay.io/ucsc_cgl/adam',
                             add_docker_parameters=add_docker_parameters,
                             benchmarking=True)

    sort_runtime = call_adam(job, None,
                             ['transformAlignments',
                              '/data/reads.bqsr.adam',
                              '/data/reads.sorted.bam',
                              '-sort_reads',
                              '-single'],
                             memory=str((job.memory / (1024 * 1024 * 1024))),
                             run_local=True,
                             container='quay.io/ucsc_cgl/adam',
                             add_docker_parameters=add_docker_parameters,
                             benchmarking=True)
    
    sorted_bam_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'reads.sorted.bam'))
    bai_tuple = job.wrapJobFn(run_sambamba_index,
                              sorted_bam_id,
                              benchmarking=True)
    job.addChild(bai_tuple)

    return (sorted_bam_id, bai_tuple.rv(0), [
        (aligner, aligner_runtime),
        ('ADAM transformAlignments', transform_runtime),
        ('ADAM transformAlignments/mark_duplicate_reads', markdups_runtime),
        ('ADAM transformAlignments/realign_indels', ri_runtime),
        ('ADAM transformAlignments/recalibrate_base_qualities', bqsr_runtime),
        ('ADAM transformAlignments/sort_reads', sort_runtime),
        ('sambamba index', bai_tuple.rv(1))
    ])


def run_gatk3_preprocessing(job,
                            ref,
                            ref_fai,
                            ref_dict,
                            dbsnp,
                            bam,
                            aligner,
                            aligner_runtime,
                            is_bam=True):
    '''
    '''

    if is_bam:
        _log.info('Sorting BAM.')
        sorted_bam = job.wrapJobFn(run_picard_sort,
                                   bam,
                                   benchmarking=True)
        job.addChild(sorted_bam)

    else:
        _log.info('Converting SAM->BAM then sorting')

        sam_as_bam = job.wrapJobFn(run_sambamba_view,
                                   bam,
                                   benchmarking=False)
        job.addChild(sam_as_bam)

        sorted_bam = job.wrapJobFn(run_picard_sort,
                                   sam_as_bam.rv(),
                                   benchmarking=True)
        sam_as_bam.addChild(sorted_bam)
        
    bqsr_tables = sorted_bam.wrapJobFn(run_base_recalibration,
                                       sorted_bam.rv(0),
                                       sorted_bam.rv(1),
                                       ref,
                                       ref_dict,
                                       ref_fai,
                                       dbsnp,
                                       dbsnp,
                                       benchmarking=True)
    sorted_bam.addChild(bqsr_tables)

    bqsred_bam = bqsr_tables.wrapJobFn(apply_bqsr_recalibration,
                                       bqsr_tables.rv(0),
                                       sorted_bam.rv(0),
                                       sorted_bam.rv(1),
                                       ref,
                                       ref_dict,
                                       ref_fai,
                                       benchmarking=True)
    bqsr_tables.addChild(bqsred_bam)

    return (bqsred_bam.rv(0), bqsred_bam.rv(1), [
        (aligner, aligner_runtime),
        ('Picard SortSam', sorted_bam.rv(2)),
        ('GATK3 BQSR', bqsr_tables.rv(1)),
        ('GATK3 BQSR/PrintReads', bqsred_bam.rv(2))
    ])


def run_speedseq(job,
                 sam,
                 aligner,
                 aligner_runtime,
                 is_bam=False):
    '''
    Runs the SpeedSeq preprocessing pipeline. Consists of:
    
    - Duplicate marking using SAMBLASTER
    - Sort using Sambamba
    - Index BAM using Sambamba

    This is equivalent to the post-alignment stages of SpeedSeq align..
    '''
    
    if is_bam:
        _log.info('Converting BAM to SAM before running SpeedSeq.')
        sam_id = job.wrapJobFn(run_samtools_view,
                               sam)
        job.addChild(sam_id)

        deduped_sam_tuple = job.wrapJobFn(run_samblaster,
                                          sam_id.rv(),
                                          benchmarking=True)
        sam_id.addFollowOn(deduped_sam_tuple)
    else:
        _log.info('Input to SpeedSeq is SAM.')
        deduped_sam_tuple = job.wrapJobFn(run_samblaster,
                                          sam,
                                          benchmarking=True)
        job.addChild(deduped_sam_tuple)
        
    deduped_sam_id = deduped_sam_tuple.rv(0)
    samblaster_time = deduped_sam_tuple.rv(1)
    
    deduped_bam_tuple = job.wrapJobFn(run_sambamba_view,
                                      deduped_sam_id,
                                      benchmarking=True)
    deduped_sam_tuple.addFollowOn(deduped_bam_tuple)
    deduped_bam_id = deduped_bam_tuple.rv(0)
    sambamba_view_time = deduped_bam_tuple.rv(1)

    sorted_bam_tuple = job.wrapJobFn(run_sambamba_sort,
                                     deduped_bam_id,
                                     benchmarking=True)
    deduped_bam_tuple.addFollowOn(sorted_bam_tuple)
    sorted_bam_id = sorted_bam_tuple.rv(0)
    sambamba_sort_time = sorted_bam_tuple.rv(1)

    bai_tuple = job.wrapJobFn(run_sambamba_index,
                              sorted_bam_id,
                              benchmarking=True)
    sorted_bam_tuple.addFollowOn(bai_tuple)
    bai_id = bai_tuple.rv(0)
    sambamba_index_time = bai_tuple.rv(1)

    return (sorted_bam_id, bai_id, [
        (aligner, aligner_runtime),
        ('samblaster', samblaster_time),
        ('sambamba view', sambamba_view_time),
        ('sambamba sort', sambamba_sort_time),
        ('sambamba index', sambamba_index_time)
    ])


def generate_plots(job,
                   runs,
                   output_dir):
    '''
    '''
    work_dir = job.fileStore.getLocalTempDir()
    timings_path = os.path.join(work_dir, 'timings.txt')
    fp = open(timings_path, 'w')

    print runs

    for i in range(len(runs)):

        print >> fp, 'Run %d stages:' % i

        (vcf, compressed, job_timings) = runs[i]

        for (job_name, runtime) in job_timings:

            print >> fp, '%s,%s' % (job_name, runtime)

    fp.flush()
    fp.close()

    move_files([timings_path], output_dir)


def generate_indices(job,
                     ref_url):
    
    download_fasta = job.wrapJobFn(download_url_job,
                                   ref_url)
    job.addChild(download_fasta)

    faidx = job.wrapJobFn(run_samtools_faidx,
                          download_fasta.rv())
    ref_dict = job.wrapJobFn(run_picard_create_sequence_dictionary,
                             download_fasta.rv())
    bowtie_indices = job.wrapJobFn(run_bowtie2_index,
                                   download_fasta.rv())
    bwa_indices = job.wrapJobFn(run_bwa_index,
                                download_fasta.rv())
    snap_indices = job.wrapJobFn(run_snap_index,
                                 download_fasta.rv())
    download_fasta.addChild(faidx)
    download_fasta.addChild(ref_dict)
    download_fasta.addChild(bowtie_indices)
    download_fasta.addChild(bwa_indices)
    download_fasta.addChild(snap_indices)

    soap3_indices = job.wrapJobFn(run_soap3_index,
                                  download_fasta.rv(),
                                  bwa_indices.rv())
    bwa_indices.addChild(soap3_indices)

    return (download_fasta.rv(),
            faidx.rv(),
            bowtie_indices.rv(),
            bwa_indices.rv(),
            snap_indices.rv(),
            soap3_indices.rv(),
            ref_dict.rv())
    

def download_fastqs(job,
                    fastq1_url, fastq2_url=None):

    fastq1 = job.wrapJobFn(download_url_job,
                           fastq1_url)
    job.addChild(fastq1)

    if fastq2_url:
        fastq2 = job.wrapJobFn(download_url_job,
                               fastq2_url)
        job.addChild(fastq2)
        return (fastq1.rv(), fastq2.rv())
    else:
        return (fastq1.rv(), None)
    

def run_preprocessing(job,
                      alignment_runs,
                      ref,
                      ref_fai,
                      ref_dict,
                      dbsnp):
    '''
    '''

    pp_runs = []

    for (aligner, (reads_id, runtime)) in alignment_runs:

        is_bam = (aligner != 'bowtie2')
        has_md_tags = False #(aligner != 'bwa')

        speedseq_run = job.wrapJobFn(run_speedseq,
                                     reads_id,
                                     aligner,
                                     runtime,
                                     is_bam=is_bam)
        pp_runs.append(speedseq_run.rv())
        job.addChild(speedseq_run)

        gatk3_run = job.wrapJobFn(run_gatk3_preprocessing,
                                  ref,
                                  ref_fai,
                                  ref_dict,
                                  dbsnp,
                                  reads_id,
                                  aligner,
                                  runtime,
                                  is_bam=is_bam)
        pp_runs.append(gatk3_run.rv())
        job.addChild(gatk3_run)

        adam_run = job.wrapJobFn(run_adam_preprocessing,
                                 ref,
                                 dbsnp,
                                 reads_id,
                                 aligner,
                                 runtime,
                                 is_bam=is_bam,
                                 has_md_tags=has_md_tags)
        pp_runs.append(adam_run.rv())
        job.addChild(adam_run)

    return pp_runs


def call_freebayes(job,
                   ref,
                   ref_fai,
                   bam,
                   bai,
                   run):

    freebayes = job.wrapJobFn(run_freebayes,
                              ref,
                              ref_fai,
                              bam,
                              bai,
                              benchmarking=True)
    job.addChild(freebayes)

    return (freebayes.rv(0), False, run.append(('freebayes', freebayes.rv(1))))


def call_platypus(job,
                  ref,
                  ref_fai,
                  bam,
                  bai,
                  assemble,
                  run):

    platypus = job.wrapJobFn(run_platypus,
                             ref,
                             ref_fai,
                             bam,
                             bai,
                             assemble=assemble,
                             benchmarking=True)
    job.addChild(platypus)
    
    if assemble:
        name = 'platypus-assemble'
    else:
        name = 'platypus'

    return (platypus.rv(0), False, run.append((name, platypus.rv(1))))


def call_strelka(job,
                 ref,
                 ref_fai,
                 bam,
                 bai,
                 manta,
                 run):
    
    if manta:
        
        manta = job.wrapJobFn(run_manta,
                              ref,
                              ref_fai,
                              bam,
                              bai,
                              benchmarking=True)
        job.addChild(manta)

        strelka = job.wrapJobFn(run_strelka,
                                ref,
                                ref_fai,
                                bam,
                                bai,
                                candidate_indels=manta.rv(2),
                                candidate_indels_tbi=manta.rv(3),
                                benchmarking=True)
        manta.addChild(strelka)

        return (strelka.rv(0), True, run.extend([('manta', manta.rv(2)),
                                                 ('strelka', strelka.rv(4))]))

    else:

        strelka = job.wrapJobFn(run_strelka,
                                ref,
                                ref_fai,
                                bam,
                                bai,
                                candidate_indels=None,
                                benchmarking=True)
        job.addChild(strelka)

        return (strelka.rv(0), True, run.append(('strelka', strelka.rv(1))))


def call_avocado_bg(job,
                    bam,
                    run):
    '''
    '''
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(bam, os.path.join(work_dir, 'reads.bam'))

    add_docker_parameters = ['-v',
                             '{}:/data'.format(work_dir)]

    transform_runtime = call_adam(job, None,
                                  ['transformAlignments',
                                   '/data/reads.bam',
                                   '/data/reads.adam'],
                                  memory=str((job.memory / (1024 * 1024 * 1024))),
                                  run_local=True,
                                  container='quay.io/ucsc_cgl/adam',
                                  add_docker_parameters=add_docker_parameters,
                                  benchmarking=True)

    avocado_runtime = call_avocado(job, None,
                                   ['biallelicGenotyper',
                                    '/data/reads.adam',
                                    '/data/gt.adam'],
                                  memory=str((job.memory / (1024 * 1024 * 1024))),
                                  run_local=True,
                                  container='quay.io/ucsc_cgl/avocado',
                                  add_docker_parameters=add_docker_parameters,
                                  benchmarking=True)
    
    vcf_runtime = call_adam(job, None,
                            ['transformGenotypes',
                             '/data/gt.adam',
                             '/data/avocado.vcf',
                             '-sort_on_save',
                             '-single'],
                            memory=str((job.memory / (1024 * 1024 * 1024))),
                            run_local=True,
                            container='quay.io/ucsc_cgl/adam',
                            add_docker_parameters=add_docker_parameters,
                            benchmarking=True)
    vcf_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'avocado.vcf'))

    return (vcf_id, False, run.extend([('ADAM transformAlignments', transform_runtime),
                                           ('avocado', avocado_runtime),
                                           ('ADAM transformGenotypes', vcf_runtime)]))

    
def call_16gt(job,
              ref,
              soap_index,
              bam,
              dbsnp,
              sample_name,
              run):
    '''
    '''
    
    gt16 = job.wrapJobFn(run_16gt,
                         ref,
                         soap_index,
                         bam,
                         dbsnp,
                         sample_name,
                         benchmarking=True)
    job.addChild(gt16)

    return (gt16.rv(0), False, run.extend([('16gt snapshot', gt16.rv(1)),
                                           ('16gt caller', gt16.rv(2)),
                                           ('16gt text2vcf', gt16.rv(3)),
                                           ('16gt filterVCF', gt16.rv(4))]))


def call_gatk3_hc(job,
                  ref,
                  ref_fai,
                  ref_dict,
                  bam,
                  bai,
                  run):
    '''
    '''

    hc = job.wrapJobFn(run_gatk3_haplotype_caller,
                       ref,
                       ref_fai,
                       ref_dict,
                       bam,
                       bai,
                       benchmarking=True)
    job.addChild(hc)

    return (hc.rv(0), True, run.append(('gatk3 hc', hc.rv(1))))


def call_mpileup(job,
                 ref,
                 ref_fai,
                 bam,
                 bai,
                 run):
    '''
    '''

    mpileup = job.wrapJobFn(run_samtools_mpileup,
                            ref,
                            ref_fai,
                            bam,
                            bai,
                            benchmarking=True)
    job.addChild(mpileup)

    bcftools = mpileup.wrapJobFn(run_bcftools_call,
                                 mpileup.rv(0),
                                 benchmarking=True)
    mpileup.addChild(bcftools)

    return (bcftools.rv(0), True, run.extend([('mpileup', mpileup.rv(1)),
                                              ('bcftools', bcftools.rv(1))]))


def call_variants(job,
                  pp_runs,
                  ref,
                  ref_fai,
                  ref_dict,
                  soap_index,
                  dbsnp,
                  sample_name):

    vc_runs = []
    for (bam, bai, timings) in pp_runs:

        freebayes = job.wrapJobFn(call_freebayes,
                                  ref,
                                  ref_fai,
                                  bam,
                                  bai,
                                  timings)
        job.addChild(freebayes)
        vc_runs.append(freebayes.rv())

        platypus = job.wrapJobFn(call_platypus,
                                 ref,
                                 ref_fai,
                                 bam,
                                 bai,
                                 False,
                                 timings)
        job.addChild(platypus)
        vc_runs.append(platypus.rv())

        platypus_assemble = job.wrapJobFn(call_platypus,
                                          ref,
                                          ref_fai,
                                          bam,
                                          bai,
                                          True,
                                          timings)
        job.addChild(platypus_assemble)
        vc_runs.append(platypus_assemble.rv())

        strelka = job.wrapJobFn(call_strelka,
                                ref,
                                ref_fai,
                                bam,
                                bai,
                                False,
                                timings)
        job.addChild(strelka)
        vc_runs.append(strelka.rv())

        manta_strelka = job.wrapJobFn(call_strelka,
                                      ref,
                                      ref_fai,
                                      bam,
                                      bai,
                                      True,
                                      timings)
        job.addChild(manta_strelka)
        vc_runs.append(manta_strelka.rv())

        gt16 = job.wrapJobFn(call_16gt,
                             ref,
                             soap_index,
                             bam,
                             dbsnp,
                             sample_name,
                             timings)
        job.addChild(gt16)
        vc_runs.append(gt16.rv())

        gatk3_hc = job.wrapJobFn(call_gatk3_hc,
                                 ref,
                                 ref_fai,
                                 ref_dict,
                                 bam,
                                 bai,
                                 timings)
        job.addChild(gatk3_hc)
        vc_runs.append(gatk3_hc.rv())

        mpileup = job.wrapJobFn(call_mpileup,
                                ref,
                                ref_fai,
                                bam,
                                bai,
                                timings)
        job.addChild(mpileup)
        vc_runs.append(mpileup.rv())

        avocado = job.wrapJobFn(call_avocado_bg,
                                bam,
                                timings)
        job.addChild(avocado)
        vc_runs.append(avocado.rv())

    return vc_runs


def download_and_index_data(job,
                            fastq1_url,
                            fastq2_url,
                            ref_url,
                            dbsnp_url):

    downloaded_fastqs = job.wrapJobFn(download_fastqs,
                                      fastq1_url,
                                      fastq2_url)
    job.addChild(downloaded_fastqs)

    dbsnp_id = job.wrapJobFn(download_url_job,
                             dbsnp_url)
    job.addChild(dbsnp_id)
    
    ref_indices = job.wrapJobFn(generate_indices,
                                ref_url)
    job.addChild(ref_indices)

    return (downloaded_fastqs.rv(0),
            downloaded_fastqs.rv(1),
            ref_indices.rv(0), # fa
            ref_indices.rv(1), # fai
            ref_indices.rv(2), # bowtie2
            ref_indices.rv(3), # bwa
            ref_indices.rv(4), # snap
            ref_indices.rv(5), # soap3
            ref_indices.rv(6), # ref dict
            dbsnp_id.rv())


def run_benchmarks(job,
                   sample_name,
                   fastq1_url,
                   fastq2_url,
                   ref_url,
                   dbsnp_url,
                   output_dir):


    reads_and_indices = job.wrapJobFn(download_and_index_data,
                                      fastq1_url,
                                      fastq2_url,
                                      ref_url,
                                      dbsnp_url)
    job.addChild(reads_and_indices)
    
    fq1 = reads_and_indices.rv(0)
    fq2 = reads_and_indices.rv(1)
    ref = reads_and_indices.rv(2)
    ref_fai = reads_and_indices.rv(3)
    bowtie2 = reads_and_indices.rv(4)
    bwa = reads_and_indices.rv(5)
    snap = reads_and_indices.rv(6)
    soap3 = reads_and_indices.rv(7)
    ref_dict = reads_and_indices.rv(8)
    dbsnp = reads_and_indices.rv(9)
    aligned_reads = job.wrapJobFn(align_reads,
                                  sample_name,
                                  fq1,
                                  fq2,
                                  ref,
                                  ref_fai,
                                  bowtie2,
                                  bwa,
                                  snap)
    reads_and_indices.addFollowOn(aligned_reads)

    preprocessed_reads = job.wrapJobFn(run_preprocessing,
                                       aligned_reads.rv(),
                                       ref,
                                       ref_fai,
                                       ref_dict,
                                       dbsnp)
    aligned_reads.addFollowOn(preprocessed_reads)

    variant_calls = job.wrapJobFn(call_variants,
                                  preprocessed_reads.rv(),
                                  ref,
                                  ref_fai,
                                  ref_dict,
                                  soap3,
                                  dbsnp,
                                  sample_name)
    preprocessed_reads.addFollowOn(variant_calls)

    plots = job.wrapJobFn(generate_plots,
                          variant_calls.rv(),
                          output_dir)
    variant_calls.addFollowOn(plots)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--reference',
                        help='The S3 URL or local path to the reference FASTA.',
                        required=True)
    parser.add_argument('--sample-name',
                        help='The name of the sample.',
                        required=True)
    parser.add_argument('--fastq1',
                        help='The S3 URL or local path to the first-of-pair'
                        'FASTQ.',
                        required=True)
    parser.add_argument('--fastq2',
                        help='The S3 URL or local path to the second-of-pair'
                        'FASTQ.',
                        default=None)
    parser.add_argument('--dbsnp',
                        help='The S3 URL or local path to a dbSNP VCF.',
                        required=True)
    parser.add_argument('--output-dir', required=True, default=None,
                        help='full path where final results will be output')

    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()

    Job.Runner.startToil(Job.wrapJobFn(run_benchmarks,
                                       args.sample_name,
                                       args.fastq1,
                                       args.fastq2,
                                       args.reference,
                                       args.dbsnp,
                                       args.output_dir), args)

if __name__ == '__main__':
    main()

                                       
