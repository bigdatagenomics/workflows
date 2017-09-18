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

"""This script is a template for bdgenomics/workflows/version.py. Running it without arguments echoes all
globals, i.e. module attributes. Constant assignments will be echoed verbatim while callables
will be invoked and their result echoed as an assignment using the function name as the left-hand
side and the return value of the function as right-hand side. To prevent a module attribute from
being echoed, start or end the attribute name with an underscore. To print the value of a single
symbol, pass the name of that attribute to the script as a command line argument. You can also
import the expand_ function and invoke it directly with either no or exactly one argument."""

# Note to maintainers:
#
#  - don't import at module level unless you want the imported value to be included in the output
#  - only import from the Python standard run-time library (you can't have any dependencies)

baseVersion = '0.0.1a1'

def version():
    """
    A version identifier that includes the full-legth commit SHA1 and an optional suffix to
    indicate that the working copy is dirty.
    """
    return _version()


def shortVersion():
    """
    A version identifier that includes the abbreviated commit SHA1 and an optional suffix to
    indicate that the working copy is dirty.
    """
    return _version(shorten=True)


def _version(shorten=False):
    return '-'.join(filter(None, [distVersion(),
                                  currentCommit()[:7 if shorten else None],
                                  ('dirty' if dirty() else None)]))


def distVersion():
    """
    The distribution version identifying a published release on PyPI.
    """
    from pkg_resources import parse_version
    build_number = buildNumber()
    parsedBaseVersion = parse_version(baseVersion)
    if isinstance(parsedBaseVersion, tuple):
        raise RuntimeError("Setuptools version 8.0 or newer required. Update by running "
                           "'pip install setuptools --upgrade'")

    if build_number is not None and parsedBaseVersion.is_prerelease:
        return baseVersion + '.dev' + build_number
    else:
        return baseVersion


def buildNumber():
    """
    The Jenkins build number, if defined, else None.
    """
    import os
    return os.getenv('BUILD_NUMBER')


def currentCommit():
    from subprocess import check_output
    return check_output('git log --pretty=oneline -n 1 -- $(pwd)', shell=True).split()[0]


def dirty():
    from subprocess import call
    return 0 != call('(git diff --exit-code '
                     '&& git diff --cached --exit-code) > /dev/null', shell=True)


def expand_(name=None):
    variables = {k: v for k, v in globals().items()
                 if not k.startswith('_') and not k.endswith('_')}

    def resolve(k):
        v = variables[k]
        if callable(v):
            v = v()
        return v

    if name is None:
        return ''.join("%s = %s\n" % (k, repr(resolve(k))) for k, v in variables.items())
    else:
        return resolve(name)


def _main():
    import sys
    sys.stdout.write(expand_(*sys.argv[1:]))


if __name__ == '__main__':
    _main()
