#!/usr/bin/env python3
import argparse
from itertools import islice
from multiprocessing import Pool
import warnings

import numpy as np


class Positive(object):
    """Create an argparse argument type that accepts only positive values

    :param mytype: Type function for type to accept, e.g. `int` or `float`
    """

    def __init__(self, mytype):
        self.mytype = mytype

    def __repr__(self):
        return "positive {}".format(self.mytype)

    def __call__(self, y):
        yt = self.mytype(y)
        if yt <= 0:
            raise argparse.ArgumentTypeError('Argument must be {}'.format(self))
        return yt


parser = argparse.ArgumentParser()
parser.add_argument('--limit', default=None, type=Positive(int),
                    help='Limit number of reads processed')
parser.add_argument('--rlc', default=False, action='store_true',
                    help='Call run-length compressed sequence')
parser.add_argument('--no-rlc', dest='rlc', action='store_false',
                    help="Don't call run-length compressed sequence")
parser.add_argument('--run_max', default=50, type=Positive(int),
                    help='Maximum run for mean approximation')
parser.add_argument('--scale', default=(1.00, 1.00, 1.00, 1.00), nargs=4,
                    metavar=('scaleA', 'scaleC', 'scaleG', 'scaleT'), type=Positive(float),
                    help='Factors for per-base scale parameter')
parser.add_argument('--shape', default=(1.00, 1.00, 1.00, 1.00), nargs=4,
                    metavar=('shapeA', 'shapeC', 'shapeG', 'shapeT'), type=Positive(float),
                    help='Factors for per-base shape parameter')
parser.add_argument('-t', '--threads', default=1, type=Positive(int),
                    help='Number of threads to use')
parser.add_argument('--width', default=60, type=Positive(int),
                    help='Line width for Fasta output')
parser.add_argument('file', default='/dev/stdin', nargs='?')


def weibull_pmf(x, shape, scale):
    """ weibull pmf for x - 1 """
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=RuntimeWarning)
        a = np.power((x - 1) / scale, shape)
        b = np.power(x / scale, shape)
        pmf = -np.exp(-a) * np.expm1(a - b)
    return pmf


def pow1p(x, y):
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=RuntimeWarning)
        res = np.exp(y * np.log1p(x))
    return res


def run_estimate_modes(shape, scale, imax=50, discrete_correction=True):
    #  Estimate run-length using mode of continuous distribution
    inv_shape = np.reciprocal(shape)

    pdf_mode = scale * np.where(shape > 1.0, pow1p(-inv_shape, inv_shape), 0.0)
    run_mode = np.int32(1 + np.floor(pdf_mode))

    if discrete_correction:
        #  Extend run.
        pmf_mode = weibull_pmf(run_mode, shape, scale)
        pmf_mode_p1 = weibull_pmf(run_mode + 1, shape, scale)
        run_mode = np.where(pmf_mode_p1 > pmf_mode, run_mode + 1, run_mode)
        #  Contract run, if possible.
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=RuntimeWarning)
            pmf_mode_m1 = weibull_pmf(run_mode - 1, shape, scale)
            run_mode = np.where((run_mode > 1) & (pmf_mode_m1 > pmf_mode),
                                run_mode - 1,
                                run_mode)

    return run_mode.astype(int)


ALPHABET = 'ACGT'
def runlength_basecall(read_data, shapef, scalef, imax=50):
    base_vec = np.array([elt[0] for elt in read_data])
    for i, b in enumerate(ALPHABET):
        #  Convert bases from characters to numerical index
        base_vec[base_vec == b] = i
    base_vec = base_vec.astype('i4')
    shape_vec = np.array([float.fromhex(elt[1]) for elt in read_data])
    scale_vec = np.array([float.fromhex(elt[2]) for elt in read_data])
    runlen_est = run_estimate_modes(shape_vec * shapef[base_vec],
                                    scale_vec * scalef[base_vec], imax=imax)
    return ''.join([ALPHABET[b] * r for b, r in zip(base_vec, runlen_est)])


def read_generator(fh):
    first_read = True
    for line in fh:
        if line.startswith('#'):
            if not first_read:
                yield read_name, read_data
            #  New record
            first_read = False
            read_name = line[2:-1]
            read_data = []
        else:
            read_data.append(line)
    yield read_name, read_data


gbl_args = None
def init_basecall_worker(*args):
    global gbl_args
    if len(args) > 0:
        gbl_args = {'shape' : np.array(args[0]),
                    'scale' : np.array(args[1]),
                    'run_max' : args[2]}


def basecall_worker(indata):
    read_name, read_lines = indata
    read_data = [line.split('\t') for line in read_lines]
    if gbl_args is None:
        basecall = ''.join([elt[0] for elt in read_data])
    else:
        basecall = runlength_basecall(read_data, gbl_args['shape'],
                                      gbl_args['scale'], gbl_args['run_max'])

    return read_name, basecall


if __name__ == '__main__':
    args = parser.parse_args()

    if args.rlc:
        init_params = []
    else:
        init_params = [args.shape, args.scale, args.run_max]

    with open(args.file, 'r') as fh:
        with Pool(processes=args.threads, initializer=init_basecall_worker, initargs=init_params) as pool:
            for res in pool.imap(basecall_worker, islice(read_generator(fh), args.limit)):
                read_name, basecall = res
                print('>{}'.format(read_name))
                print('\n'.join([basecall[st : st + args.width]
                                 for st in range(0, len(basecall), args.width)]))
