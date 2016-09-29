# Copyright (C) 2016 Regents of the University of California
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import tempfile
import unittest
import os

from toil.job import Job


class DockerCallTest(unittest.TestCase):
    """
    This class handles creating a tmpdir and Toil options suitable for a unittest.
    """
    def setUp(self):
        super(DockerCallTest, self).setUp()
        # the test tmpdir needs to be in the home directory so files written onto mounted
        # directories from a Docker container will be visible on the host
        # https://docs.docker.com/docker-for-mac/osxfs/
        home = os.path.expanduser("~") + '/'
        self.tmpdir = tempfile.mkdtemp(prefix=home)
        self.options = Job.Runner.getDefaultOptions(os.path.join(str(self.tmpdir), 'jobstore'))
        self.options.clean = 'always'

    def tearDown(self):
        # delete temp
        super(DockerCallTest, self).tearDown()
        for file in os.listdir(self.tmpdir):
            os.remove(os.path.join(self.tmpdir, file))
        os.removedirs(self.tmpdir)


# this is lifted from toil.test; perhaps refactor into bpl?
try:
    # noinspection PyUnresolvedReferences
    from _pytest.mark import MarkDecorator
except ImportError:
    # noinspection PyUnusedLocal
    def _mark_test(name, test_item):
        return test_item
else:
    def _mark_test(name, test_item):
        return MarkDecorator(name)(test_item)


def needs_spark(test_item):
    """
    Use as a decorator before test classes or methods to only run them if Spark is usable.
    """
    test_item = _mark_test('spark', test_item)

    try:
        # noinspection PyUnresolvedReferences
        import pyspark
    except ImportError:
        return unittest.skip("Skipping test. Install PySpark to include this test.")(test_item)
    except:
        raise
    else:
        return test_item

