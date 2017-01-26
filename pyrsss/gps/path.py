import os


try:
    GPSTK_BUILD_PATH = os.environ['GPSTK_BUILD']
except KeyError:
    raise RuntimeError('environment variable GPSTK_BUILD not set')
