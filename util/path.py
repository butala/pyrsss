import os
import shutil
import tempfile


class TempDirectory(object):
    def __init__(self, delete=True, **kwds):
        """
        Context manager for a temporary directory. If *delete*, then
        cleanup all files in the temporary location upon exiting the
        context manager block.
        """
        self.kwds = kwds
        self.delete = delete

    def __enter__(self):
        self.tempdir = tempfile.mkdtemp(**self.kwds)
        return self.tempdir

    def __exit__(self, type, value, traceback):
        if self.delete:
            shutil.rmtree(self.tempdir)


class SmartTempDir(TempDirectory):
    def __init__(self, dirname=None, **kwds):
        """
        Context manager that creates a temporary directory if *dirname* is
        `None` and, otherwise, returns the given path.
        """
        super(SmartTempDir, self).__init__(**kwds)
        self.dirname = dirname

    def __enter__(self):
        if self.dirname is None:
            return super(SmartTempDir, self).__enter__()
        else:
            return self.dirname

    def __exit__(self, type, value, traceback):
        if self.dirname is None:
            super(SmartTempDir, self).__exit__(type, value, traceback)


def touch_path(path):
    """
    If *path* does not exist, create it. Returns *path*.
    """
    if not os.path.isdir(path):
        os.makedirs(path)
    return path


def replace_path(path, fname):
    """
    Return the replacement the path of *fname* with *path*.
    """
    return os.path.join(path, os.path.basename(fname))
