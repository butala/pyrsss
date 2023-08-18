import os


try:
    GNSSTK_BUILD_PATH = os.environ['GNSSTK_BUILD']
except KeyError:
    raise RuntimeError('environment variable GNSSTK_BUILD not set')
