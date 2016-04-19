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


def tail(f, window=20):
    """
    Returns the last *window* lines of file object *f* as a list.

    From http://stackoverflow.com/a/7047765/6226265.
    """
    if window == 0:
        return []
    BUFSIZ = 1024
    f.seek(0, 2)
    bytes_i = f.tell()
    size = window + 1
    block = -1
    data = []
    while size > 0 and bytes_i > 0:
        if bytes_i - BUFSIZ > 0:
            # Seek back one whole BUFSIZ
            f.seek(block * BUFSIZ, 2)
            # read BUFFER
            data.insert(0, f.read(BUFSIZ))
        else:
            # file too small, start from beginning
            f.seek(0,0)
            # only read what was not read
            data.insert(0, f.read(bytes_i))
        linesFound = data[0].count('\n')
        size -= linesFound
        bytes_i -= BUFSIZ
        block -= 1
    return ''.join(data).splitlines()[-window:]
