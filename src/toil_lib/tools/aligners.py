import os

import subprocess

from toil_lib.programs import docker_call
from toil_lib.urls import download_url


def run_star(job, r1_id, r2_id, star_index_url, wiggle=False):
    """
    Performs alignment of fastqs to bam via STAR

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str r1_id: FileStoreID of fastq (pair 1)
    :param str r2_id: FileStoreID of fastq (pair 2 if applicable, else pass None)
    :param str star_index_url: STAR index tarball
    :param bool wiggle: If True, will output a wiggle file and return it
    :return: FileStoreID from RSEM
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    download_url(url=star_index_url, name='starIndex.tar.gz', work_dir=work_dir)
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'starIndex.tar.gz'), '-C', work_dir])
    os.remove(os.path.join(work_dir, 'starIndex.tar.gz'))
    # Determine tarball structure - star index contains are either in a subdir or in the tarball itself
    star_index = os.path.join('/data', os.listdir(work_dir)[0]) if len(os.listdir(work_dir)) == 1 else '/data'
    # Parameter handling for paired / single-end data
    parameters = ['--runThreadN', str(job.cores),
                  '--genomeDir', star_index,
                  '--outFileNamePrefix', 'rna',
                  '--outSAMtype', 'BAM', 'SortedByCoordinate',
                  '--outSAMunmapped', 'Within KeepPairs',
                  '--quantMode', 'TranscriptomeSAM',
                  '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                  '--outFilterType', 'BySJout',
                  '--outFilterMultimapNmax', '20',
                  '--outFilterMismatchNmax', '999',
                  '--outFilterMismatchNoverReadLmax', '0.04',
                  '--alignIntronMin', '20',
                  '--alignIntronMax', '1000000',
                  '--alignMatesGapMax', '1000000',
                  '--alignSJoverhangMin', '8',
                  '--alignSJDBoverhangMin', '1',
                  '--sjdbScore', '1']
    if wiggle:
        parameters.extend(['--outWigType', 'bedGraph',
                           '--outWigStrand', 'Unstranded',
                           '--outWigReferencesPrefix', 'chr'])
    if r1_id and r2_id:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
        job.fileStore.readGlobalFile(r2_id, os.path.join(work_dir, 'R2.fastq'))
        parameters.extend(['--readFilesIn', '/data/R1.fastq', '/data/R2.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
        parameters.extend(['--readFilesIn', '/data/R1.fastq'])
    # Call: STAR Mapping
    docker_call(job=job, tool='quay.io/ucsc_cgl/star:2.4.2a--bcbd5122b69ff6ac4ef61958e47bde94001cfe80',
                work_dir=work_dir, parameters=parameters)
    # Write to fileStore
    transcriptome_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.toTranscriptome.out.bam'))
    sorted_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.sortedByCoord.out.bam'))
    log_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaLog.final.out'))
    if wiggle:
        wiggle_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaSignal.UniqueMultiple.str1.out.bg'))
        return transcriptome_id, sorted_id, wiggle_id, log_id
    else:
        return transcriptome_id, sorted_id, log_id


def run_bwakit(job, config, sort=True, trim=False):
    """
    Runs BWA-Kit to align single or paired-end fastq files or realign SAM/BAM files.

    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param Namespace config: A configuration object that holds strings as attributes.
        The attributes must be accessible via the dot operator.
        The config must have:
        config.r1               FileStoreID for FASTQ file, or None if realigning SAM/BAM
        config.r2               FileStoreID for paired FASTQ file, or None if single-ended
        config.bam              FileStoreID for BAM file to be realigned, or None if aligning fastq
        config.sam              FileStoreID for SAM file to be realigned, or None if aligning fastq
        config.ref              FileStoreID for the reference genome
        config.fai              FileStoreID for the reference index file
        config.amb              FileStoreID for the reference amb file
        config.ann              FileStoreID for the reference ann file
        config.bwt              FileStoreID for the reference bwt file
        config.pac              FileStoreID for the reference pac file
        config.sa               FileStoreID for the reference sa file
        config.alt              FileStoreID for the reference alt (or None)
        config.rg_line          The read group value to use (or None -- see below)
        config.library          Read group attribute: library
        config.platform         Read group attribute: platform
        config.program_unit     Read group attribute: program unit
        config.uuid             Read group attribute: sample ID

        If specifying config.rg_line, use the following format:
            BAM read group header line (@RG), as defined on page 3 of the SAM spec.
            Tabs should be escaped, e.g., @RG\\tID:foo\\tLB:bar...
            for the read group "foo" from sequencing library "bar".
            Multiple @RG lines can be defined, but should be split by an escaped newline \\n,
            e.g., @RG\\tID:foo\\t:LB:bar\\n@RG\\tID:santa\\tLB:cruz.

    :param bool sort: If True, sorts the BAM
    :param bool trim: If True, performs adapter trimming
    :return: FileStoreID of BAM
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'ref.fa': config.ref,
              'ref.fa.fai': config.fai,
              'ref.fa.amb': config.amb,
              'ref.fa.ann': config.ann,
              'ref.fa.bwt': config.bwt,
              'ref.fa.pac': config.pac,
              'ref.fa.sa': config.sa}
    samples = []
    realignment = False
    # If a fastq pair was provided
    if getattr(config, 'r1', None):
        inputs['input.1.fq.gz'] = config.r1
        samples.append('input.1.fq.gz')
    if getattr(config, 'r2', None):
        inputs['input.2.fq.gz'] = config.r2
        samples.append('input.2.fq.gz')
    if getattr(config, 'bam', None):
        inputs['input.bam'] = config.bam
        samples.append('input.bam')
        realignment = True
    if getattr(config, 'sam', None):
        inputs['input.sam'] = config.sam
        samples.append('input.sam')
        realignment = True
    # If an alt file was provided
    if getattr(config, 'alt', None):
        inputs['ref.fa.alt'] = config.alt
    for name, fileStoreID in inputs.iteritems():
        job.fileStore.readGlobalFile(fileStoreID, os.path.join(work_dir, name))
    # If a read group line was provided
    if getattr(config, 'rg_line', None):
        rg = config.rg_line
    # Otherwise, generate a read group line to place in the BAM.
    elif all(getattr(config, elem, None) for elem in ['library', 'platform', 'program_unit', 'uuid']):
        rg = "@RG\\tID:{0}".format(config.uuid)  # '\' character is escaped so bwakit gets passed '\t' properly
        rg_attributes = [config.library, config.platform, config.program_unit, config.uuid]
        for tag, info in zip(['LB', 'PL', 'PU', 'SM'], rg_attributes):
            rg += '\\t{0}:{1}'.format(tag, info)
    # If realigning, then bwakit can use pre-existing read group data
    elif realignment:
        rg = None

    # BWA Options
    opt_args = []
    if sort:
        opt_args.append('-s')
    if trim:
        opt_args.append('-a')
    # Call: bwakit
    parameters = ['-t', str(job.cores)] + opt_args + ['-o', '/data/aligned', '/data/ref.fa']
    if rg is not None:
        parameters = ['-R', rg] + parameters
    for sample in samples:
        parameters.append('/data/{}'.format(sample))
    mock_bam = config.uuid + '.bam'
    outputs = {'aligned.aln.bam': mock_bam}

    docker_call(job=job, tool='quay.io/ucsc_cgl/bwakit:0.7.12--528bb9bf73099a31e74a7f5e6e3f2e0a41da486e',
                parameters=parameters, inputs=inputs.keys(), outputs=outputs, work_dir=work_dir)

    # Either write file to local output directory or upload to S3 cloud storage
    job.fileStore.logToMaster('Aligned sample: {}'.format(config.uuid))
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'aligned.aln.bam'))
