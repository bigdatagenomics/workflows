"""
Functions used across pipelines

@author Alyssa Morrow
"""
from toil_lib import require

# all files must have s3 urls
def is_s3(f):
    require(f.startswith("s3a"),
            "url for file %s did not start with s3a scheme" % f)
