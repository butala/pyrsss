import os


def get_gnsstk_build_path():
    """
    Return the GNSSTk build directory. Raises if the GNSSTK_BUILD
    environment variable is not set (only the modules that need the
    optional gnsstk tools/extension require it).
    """
    try:
        return os.environ['GNSSTK_BUILD']
    except KeyError:
        raise RuntimeError('environment variable GNSSTK_BUILD not set')

