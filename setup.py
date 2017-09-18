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

from __future__ import print_function

import sys

from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand
from pkg_resources import parse_version, require, DistributionNotFound

def check_provided(distribution, min_version, max_version=None, optional=False):
    min_version = parse_version(min_version)
    if max_version is not None:
        max_version = parse_version(max_version)

    messages = []

    toil_missing = 'Cannot find a valid installation of Toil.'
    dist_missing = 'Cannot find an installed copy of the %s distribution, typically provided by Toil.' % distribution
    version_too_low = 'The installed copy of %s is out of date. It is typically provided by Toil.' % distribution
    version_too_high = 'The installed copy of %s is too new. It is typically provided by Toil.' % distribution
    required_version = 'Setup requires version %s or higher' % min_version
    required_version += '.' if max_version is None else ', up to but not including %s.' % max_version
    install_toil = 'Installing Toil should fix this problem.'
    upgrade_toil = 'Upgrading Toil should fix this problem.'
    reinstall_dist = 'Uninstalling %s and reinstalling Toil should fix this problem.' % distribution
    reinstall_toil = 'Uninstalling Toil and reinstalling it should fix this problem.'
    footer = ("Setup doesn't install Toil automatically to give you a chance to choose any of the optional extras "
              "that Toil provides. More on installing Toil at http://toil.readthedocs.io/en/latest/installation.html.")
    try:
        # This check will fail if the distribution or any of its dependencies are missing.
        version = require(distribution)[0].version
    except DistributionNotFound:
        version = None
        if not optional:
            messages.extend([toil_missing if distribution == 'toil' else dist_missing, install_toil])
    else:
        if parse_version(version) < min_version:
            messages.extend([version_too_low, required_version,
                             upgrade_toil if distribution == 'toil' else reinstall_dist])
        elif max_version is not None and max_version < parse_version(version):
            messages.extend([version_too_high, required_version,
                             reinstall_toil if distribution == 'toil' else reinstall_dist])
    if messages:
        messages.append(footer)
        raise RuntimeError(' '.join(messages))
    else:
        return version

def importVersion():
    """
    Load and return the module object for bdgenomics/workflows/version.py, generating it from the template if
    required.
    """

    try:
        # Attempt to load the template first. It only exists in a working copy cloned via git.
        import version_template
    except ImportError:
        # If loading the template fails we must be in a unpacked source distribution and
        # src/toil/version.py will already exist.
        pass
    else:
        # Use the template to generate src/toil/version.py
        import os
        import errno
        from tempfile import NamedTemporaryFile

        new = version_template.expand_()

        print(new, sys.stderr)
        
        try:
            with open('bdgenomics/workflows/version.py') as f:
                old = f.read()
        except IOError as e:
            if e.errno == errno.ENOENT:
                old = None
            else:
                raise

        if old != new:
            with NamedTemporaryFile(dir='bdgenomics/workflows', prefix='version.py.', delete=False) as f:
                f.write(new)
            os.rename(f.name, 'bdgenomics/workflows/version.py')

    import bdgenomics.workflows.version
    return bdgenomics.workflows.version

version = importVersion()
print(version, sys.stderr)

toil_version = check_provided('toil', min_version='3.7.0a1.dev392', max_version='3.12.0')

kwargs = dict(
    name='bdgenomics.workflows',
    version=version.distVersion,
    description='A repository of genomic workflows developed by the UC Berkeley AMPLab and UCSC Computational Genomics lab that use Toil to run ADAM/BDG tools',
    author='UC Berkeley AMP Lab',
    author_email='adam-developers@googlegroups.com',
    url="https://github.com/bigdatagenomics/workflows",
    install_requires=[
        'pyyaml==3.11',
        'toil-lib==1.1.8'],
    tests_require=[
        'pytest==2.8.3'],
    entry_points={
        'console_scripts': [
            'bdg-adam = bdgenomics.workflows.adam_pipeline.preprocessing:main',
            'bdg-avocado = bdgenomics.workflows.avocado_pipeline.variant_calling:main',
            'bdg-deca = bdgenomics.workflows.deca_pipeline.call_cnvs:main',
            'bdg-cannoli-bwa = bdgenomics.workflows.cannoli_pipeline.bwa_alignment:main',
            'bdg-gatk3-benchmark = bdgenomics.workflows.benchmarking.gatk3_pipeline.preprocessing:main',
            'bdg-mkdups-benchmark = bdgenomics.workflows.benchmarking.single_node.mkdups:main',
            'bdg-sort-benchmark = bdgenomics.workflows.benchmarking.single_node.sort:main',
            'bdg-ri-benchmark = bdgenomics.workflows.benchmarking.single_node.realign_indels:main',
            'bdg-bqsr-benchmark = bdgenomics.workflows.benchmarking.single_node.bqsr:main']},
    packages=find_packages(),
    classifiers=["License :: OSI Approved :: Apache Software License"])

setup(**kwargs)
