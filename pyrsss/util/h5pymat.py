from pathlib import Path
import datetime
import sys

import h5py
import numpy as np
import scipy as sp


def get_scalar(dset):
    """
    """
    assert dset.shape == (1, 1)
    return dset[0][0]


def get_vector(dset):
    """
    """
    assert dset.shape[0] == 1 or dset.shape[1] == 1
    v = np.empty(dset.shape, dtype=dset.dtype)
    dset.read_direct(v)
    return v.squeeze()


def get_string(dset):
    """
    """
    return ''.join(map(chr, get_vector(dset)))


def loadmat(mat_fname, tosparse=True, toint=True):
    """
    """
    d = {}
    with h5py.File(mat_fname, 'r') as h5_mat:
        for k, v in h5_mat.items():
            if v.dtype == np.dtype('uint16') and (v.attrs['MATLAB_class'] == b'char' or v.attrs['MATLAB_class'] == 'char') and v.attrs['MATLAB_int_decode'] == 2:
                v = get_string(v)
            elif v.shape == (1, 1):
                v = get_scalar(v)
                if toint and 'uint' in str(v.dtype):
                    v = int(v)
            elif v.shape[0] == 1 or v.shape[1] == 1:
                v = get_vector(v)
                if toint and 'uint' in str(v.dtype):
                    v = v.astype(int)
            d[k] = v
        if tosparse:
            sparse_keys = [k for k in h5_mat.keys() if k.endswith('_csr_data')]
            for k in sparse_keys:
                mat_name = k[:-9]
                assert mat_name not in d
                csr_data = d.pop(f'{mat_name}_csr_data')
                csr_indices = d.pop(f'{mat_name}_csr_indices')
                csr_indptr = d.pop(f'{mat_name}_csr_indptr')
                csr_shape = d.pop(f'{mat_name}_csr_shape')
                d[mat_name] = sp.sparse.csr_array((csr_data,
                                                   csr_indices,
                                                   csr_indptr),
                                                  shape=csr_shape)
    return d


################################################################################


def write_matlab_7_3_header(mat_fname):
    """
    Copied from hdf5storage __init__.py
    """
    # Get the time.
    now = datetime.datetime.now()

    # Construct the leading string. The MATLAB one looks like
    #
    # s = 'MATLAB 7.3 MAT-file, Platform: GLNXA64, Created on: ' \
    #     + now.strftime('%a %b %d %H:%M:%S %Y') \
    #     + ' HDF5 schema 1.00 .'
    #
    # Platform is going to be changed to CPython version. The
    # version is just gotten from sys.version_info, which is a class
    # for Python >= 2.7, but a tuple before that.

    v = sys.version_info
    if sys.hexversion >= 0x02070000:
        v = {'major': v.major, 'minor': v.minor, 'micro': v.micro}
    else:
        v = {'major': v[0], 'minor': v[1], 'micro': v[1]}

    s = 'MATLAB 7.3 MAT-file, Platform: CPython ' \
        + '{0}.{1}.{2}'.format(v['major'], v['minor'], v['micro']) \
        + ', Created on: {0} {1}'.format(
        ('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')[ \
        now.weekday()], \
        ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', \
         'Sep', 'Oct', 'Nov', 'Dec')[now.month - 1]) \
        + now.strftime(' %d %H:%M:%S %Y') \
        + ' HDF5 schema 1.00 .'

    # Make the bytearray while padding with spaces up to 128-12
    # (the minus 12 is there since the last 12 bytes are special.

    b = bytearray(s + (128-12-len(s))*' ', encoding='utf-8')

    # Add 8 nulls (0) and the magic number (or something) that
    # MATLAB uses. Lengths must be gone to to make sure the argument
    # to fromhex is unicode because Python 2.6 requires it.

    b.extend(bytearray.fromhex(
             b'00000000 00000000 0002494D'.decode()))


    with h5py.File(mat_fname, 'w', userblock_size=512):
        pass

    with open(mat_fname, 'r+b') as fid:
        fid.write(b)


class File:
    """
    """
    def __init__(self, name, mode='r', **kwds):
        if mode in ['w-', 'x']:
            assert not Path(name).exists()
        if mode in ['w', 'w-', 'x']:
            write_matlab_7_3_header(name)
            self.h5_file = h5py.File(name, mode='a', **kwds)
        else:
            self.h5_file = h5pt.File(name, mode=mode, **kwds)

    def __enter__(self):
        return self.h5_file

    def __exit__(self, type, value, traceback):
        self.h5_file.close()


def savemat(mat_fname, m):
    """
    """
    with File(mat_fname, 'w') as h5_mat:
        for k in [x for x in m.keys() if sp.sparse.issparse(m[x])]:
            A_csr = m.pop(k).tocsr()
            m[f'{k}_csr_data'] = A_csr.data
            m[f'{k}_csr_indices'] = A_csr.indices
            m[f'{k}_csr_indptr'] = A_csr.indptr
            m[f'{k}_csr_shape'] = A_csr.shape
        for k, v in m.items():
            if isinstance(v, str):
                shape = len(v), 1
                v = np.array(list(v.encode('utf-8')), dtype='uint16')
                v.shape = v.size, 1
                dset = h5_mat.create_dataset(k, shape, dtype='uint16')
                dset.attrs.create('MATLAB_class', b'char')
                dset.attrs.create('MATLAB_int_decode', 2)
            else:
                if not isinstance(v, np.ndarray):
                    v = np.array(v)
                match v.ndim:
                    case 0:
                        dset = h5_mat.create_dataset(k, (1, 1), dtype=v.dtype)
                    case 1:
                        dset = h5_mat.create_dataset(k, (1, len(v)), dtype=v.dtype)
                    case 2:
                        dset = h5_mat.create_dataset(k, v.shape, dtype=v.dtype)
                    case _:
                        assert False
            dset.write_direct(v)
