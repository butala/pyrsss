from dataclasses import dataclass
from typing import List

import numpy as np


@dataclass
class ZRecord:
    T: List[float]
    Zxx: List[complex]
    Zxy: List[complex]
    Zyx: List[complex]
    Zyy: List[complex]

    @property
    def f(self):
        return 1/np.array(self.T)


def parse_zfile(fname):
    """
    """
    with open(fname) as fid:
        # reading through header
        for line in fid:
            if line.startswith('number of channels'):
                toks = line.split()
                assert len(toks) == 8
                n_chan = int(toks[3])
                assert n_chan == 5
                n_freq = int(toks[7])
            if len(line.strip()) == 0:
                # end of header
                break
        T = []
        Zxx = []
        Zxy = []
        Zyx = []
        Zyy = []
        # read data blocks
        for _ in range(n_freq):
            # period line
            line = next(fid)
            assert line.startswith('period :')
            toks = line.split()
            assert len(toks) == 12
            T.append(float(toks[2]))
            # data points line
            line = next(fid)
            # start of transfer function block
            line = next(fid)
            assert line.startswith(' Transfer Functions')
            # tipper line
            line = next(fid)
            # Zxx / Zxy line
            line = next(fid)
            toks = line.split()
            toks = list(map(float, line.split()))
            assert len(toks) == 4
            Zxx.append(complex(toks[0], toks[1]))
            Zxy.append(complex(toks[2], toks[3]))
            # Zyx / Zyy line
            line = next(fid)
            toks = list(map(float, line.split()))
            assert len(toks) == 4
            Zyx.append(complex(toks[0], toks[1]))
            Zyy.append(complex(toks[2], toks[3]))
            # signal power block
            line = next(fid)
            assert line.startswith(' Inverse Coherent Signal Power Matrix')
            line = next(fid)
            line = next(fid)
            # residual covariance block
            line = next(fid)
            assert line.startswith(' Residual Covariance')
            line = next(fid)
            line = next(fid)
            line = next(fid)
    return ZRecord(T, Zxx, Zxy, Zyx, Zyy)
