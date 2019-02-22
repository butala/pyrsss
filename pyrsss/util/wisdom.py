import logging
import os
from cPickle import load, dump

import pyfftw

logger = logging.getLogger('pyrsss.util.wisdom')



DIRECTION_MAP = {'forward':  'FFTW_FORWARD',
                 'backward': 'FFTW_BACKWARD'}


RIGOR_MAP = {'measure': 'FFTW_MEASURE',
             'patient': 'FFTW_PATIENT'}



def wisdom_fname(path,
                 plan_type,
                 pow2,
                 threads,
                 direction,
                 rigor):
    """
    """
    return os.path.join(path,
                        'wisdom_{plan_type}_{pow2:02d}_{threads:02d}_'
                        '{direction}_{rigor}.pkl'.format(plan_type=plan_type,
                                                         pow2=pow2,
                                                         threads=threads,
                                                         direction=direction,
                                                         rigor=rigor))


def compute_pow2_real_wisdom(path,
                             pow2=range(20),
                             rigor='measure',
                             threads=16):
    """
    ???

    If you plan with FFTW_PATIENT, it will automatically disable
    threads for sizes that don't benefit from parallelization.
    """
    flags = [RIGOR_MAP[rigor],
             'FFTW_DESTROY_INPUT']
    wisdom_fnames = []
    for pow2_i in pow2:
        N = 2**pow2_i
        for direction in ['forward', 'backward']:
            logger.info('building wisdom for real 2**{} {}'.format(pow2_i,
                                                                   direction))
            if direction == 'forward':
                x_input  = pyfftw.empty_aligned(N, dtype='float64')
                x_output = pyfftw.empty_aligned(int(N // 2) + 1, dtype='complex128')
            else:
                x_output = pyfftw.empty_aligned(N, dtype='float64')
                x_input  = pyfftw.empty_aligned(int(N // 2) + 1, dtype='complex128')
            plan = pyfftw.FFTW(x_input,
                               x_output,
                               direction=DIRECTION_MAP[direction],
                               flags=flags,
                               threads=threads)
            wisdom_fnames.append(wisdom_fname(path,
                                              'real',
                                              pow2_i,
                                              threads,
                                              direction,
                                              rigor))
            logger.info('writing to {}'.format(wisdom_fnames[-1]))
            with open(wisdom_fnames[-1], 'w') as fid:
                dump(pyfftw.export_wisdom(), fid, -1)
            pyfftw.forget_wisdom()
    return wisdom_fnames


def compute_pow2_complex_wisdom(path,
                                pow2=range(20),
                                rigor='measure',
                                threads=16):
    """
    ???

    If you plan with FFTW_PATIENT, it will automatically disable
    threads for sizes that don't benefit from parallelization.
    """
    flags = [RIGOR_MAP[rigor],
             'FFTW_DESTROY_INPUT']
    wisdom_fnames = []
    for pow2_i in pow2:
        N = 2**pow2_i
        x_input = pyfftw.empty_aligned(N, dtype='complex128')
        x_output = pyfftw.empty_aligned(N, dtype='complex128')
        for direction in ['forward', 'backward']:
            logger.info('building wisdom for complex 2**{} {}'.format(pow2_i,
                                                                      direction))
            plan = pyfftw.FFTW(x_input,
                               x_output,
                               direction=DIRECTION_MAP[direction],
                               flags=flags,
                               threads=threads)
            wisdom_fnames.append(wisdom_fname(path,
                                              'complex',
                                              pow2_i,
                                              threads,
                                              direction,
                                              rigor))
            logger.info('writing to {}'.format(wisdom_fnames[-1]))
            with open(wisdom_fnames[-1], 'w') as fid:
                dump(pyfftw.export_wisdom(), fid, -1)
            pyfftw.forget_wisdom()
    return wisdom_fnames


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    path = '/rdata/airglow/butala/data/fftw_wisdom'

    # compute_pow2_real_wisdom(path,
    #                          pow2=[20])

    compute_pow2_complex_wisdom('/rdata/airglow/butala/data/fftw_wisdom',
                                pow2=[20, 21, 22, 23, 24, 25, 26])
